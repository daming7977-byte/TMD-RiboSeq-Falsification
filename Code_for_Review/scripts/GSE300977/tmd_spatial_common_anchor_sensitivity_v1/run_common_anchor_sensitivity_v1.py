from pathlib import Path
import numpy as np
import pandas as pd

PAIR = Path(
    "crossmap_window_preregister/output/"
    "GSE300977_PAIRED_REAL_PSEUDO_CROSSMAP_MASK_AUDIT_v1.tsv"
)

CATALOG = Path(
    "tmd_spatial_secondary_scan_v1/output/"
    "SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv"
)

COUNTDIR = Path(
    "tmd_spatial_secondary_scan_v1/geometry_counts"
)

OUTDIR = Path(
    "tmd_spatial_common_anchor_sensitivity_v1/output"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

CENTERS = [-30,-15,0,15,30,45,60,75,90,105,120]

LIBRARIES = [
    "B01","B02","B03","B04","B05","B06",
    "B07","B08","B09","B10","B11","B12",
]

CONTRASTS = [
    ("ER",1,"B01","B02"),
    ("ER",2,"B11","B12"),
    ("mitochondria",1,"B09","B10"),
    ("mitochondria",2,"B07","B08"),
    ("cytosol",1,"B05","B06"),
    ("cytosol",2,"B03","B04"),
]

BOOT_N = 10000
BOOT_SEED = 30097775
BOOT_BATCH = 200

# --------------------------------------------------
# Eligible paired rows
# --------------------------------------------------

d = pd.read_csv(PAIR, sep="\t", dtype=str)

d["center"] = pd.to_numeric(
    d["center"]
).astype(int)

d["real_crossmap_unmasked_nt"] = pd.to_numeric(
    d["real_crossmap_unmasked_nt"]
)

d["pseudo_min_unmasked_nt"] = pd.to_numeric(
    d["pseudo_min_unmasked_nt"]
)

d["real_anchor_genomic_start0"] = pd.to_numeric(
    d["real_anchor_genomic_start0"]
).astype(int)

d = d[
    (d["real_crossmap_unmasked_nt"] >= 1)
    &
    (d["pseudo_min_unmasked_nt"] >= 1)
].copy()

KEY = [
    "chrom",
    "strand",
    "real_anchor_genomic_start0",
]

# --------------------------------------------------
# Same genomic anchor must occur exactly once
# at every one of the 11 centers
# --------------------------------------------------

unique_keys_by_center = {}

for c in CENTERS:

    z = d[d["center"].eq(c)]

    sizes = z.groupby(KEY).size()

    keys = set(
        sizes[
            sizes.eq(1)
        ].index.tolist()
    )

    unique_keys_by_center[c] = keys

common = set.intersection(
    *[
        unique_keys_by_center[c]
        for c in CENTERS
    ]
)

common = sorted(common)

N = len(common)

print("COMMON_ANCHOR_UNIVERSE =", N)

if N == 0:
    raise SystemExit(
        "COMMON_ANCHOR_GATE = FAIL"
    )

print("COMMON_ANCHOR_GATE = PASS")

common_set = set(common)

# --------------------------------------------------
# Geometry catalog -> geometry ID
# --------------------------------------------------

cat = pd.read_csv(
    CATALOG,
    sep="\t",
    dtype=str
)

key_to_gid = {}

for _,r in cat.iterrows():

    key = (
        str(r["strand"]),
        str(r["geometry"]),
        str(r["mask_bits"]),
    )

    if key in key_to_gid:
        raise SystemExit(
            "CATALOG_KEY_DUPLICATION"
        )

    key_to_gid[key] = str(
        r["geometry_id"]
    )

print(
    "GEOMETRY_CATALOG_ROWS =",
    len(cat)
)

print(
    "GEOMETRY_CATALOG_GATE = PASS"
)

# --------------------------------------------------
# Reuse completed secondary density checkpoints
# --------------------------------------------------

density = {}

for lib in LIBRARIES:

    p = (
        COUNTDIR /
        f"{lib}.SECONDARY_GEOMETRY_DENSITIES_v1.tsv"
    )

    if not p.exists():
        raise SystemExit(
            f"MISSING_COUNT_CHECKPOINT:{lib}:{p}"
        )

    x = pd.read_csv(
        p,
        sep="\t"
    )

    density[lib] = dict(
        zip(
            x["geometry_id"].astype(str),
            x["density"].astype(float)
        )
    )

print(
    "SECONDARY_COUNT_CHECKPOINTS = REUSED"
)

# --------------------------------------------------
# Helpers
# --------------------------------------------------

def splitsemi(x):
    return [
        z
        for z in str(x).split(";")
        if z != ""
    ]

# --------------------------------------------------
# Fixed common-anchor matrix at each center
# --------------------------------------------------

matrices = {}

anchor_order = common

for c in CENTERS:

    z = d[
        d["center"].eq(c)
    ].copy()

    z["anchor_key"] = list(
        zip(
            z["chrom"],
            z["strand"],
            z["real_anchor_genomic_start0"],
        )
    )

    z = z[
        z["anchor_key"].isin(
            common_set
        )
    ].copy()

    if len(z) != N:
        raise SystemExit(
            f"CENTER_COMMON_ROW_GATE_FAIL:"
            f"{c}:{len(z)}:{N}"
        )

    z = z.set_index(
        "anchor_key",
        verify_integrity=True
    ).loc[anchor_order]

    rows = []

    for _,r in z.iterrows():

        real_key = (
            str(r["strand"]),
            str(r["real_window_geometry"]),
            str(r["real_crossmap_mask_bits"]),
        )

        if real_key not in key_to_gid:
            raise SystemExit(
                f"REAL_GEOMETRY_NOT_IN_CATALOG:{c}"
            )

        real_gid = key_to_gid[
            real_key
        ]

        pgeoms = splitsemi(
            r["pseudo_window_geometries"]
        )

        pbits = splitsemi(
            r["pseudo_crossmap_mask_bits"]
        )

        if len(pgeoms) != len(pbits):
            raise SystemExit(
                "PSEUDO_CARDINALITY_FAIL"
            )

        pseudo_gids = []

        for geom,bits in zip(
            pgeoms,
            pbits
        ):

            pkey = (
                str(r["strand"]),
                geom,
                bits,
            )

            if pkey not in key_to_gid:
                raise SystemExit(
                    f"PSEUDO_GEOMETRY_NOT_IN_CATALOG:{c}"
                )

            pseudo_gids.append(
                key_to_gid[pkey]
            )

        row = []

        for lib in LIBRARIES:

            real_density = density[
                lib
            ][real_gid]

            pseudo_density = float(
                np.mean(
                    [
                        density[lib][gid]
                        for gid in pseudo_gids
                    ]
                )
            )

            row.extend(
                [
                    real_density,
                    pseudo_density,
                ]
            )

        rows.append(row)

    X = np.asarray(
        rows,
        dtype=float
    )

    if X.shape != (N,24):
        raise SystemExit(
            f"MATRIX_SHAPE_FAIL:"
            f"{c}:{X.shape}"
        )

    if (
        not np.isfinite(X).all()
        or
        (X < 0).any()
    ):
        raise SystemExit(
            f"MATRIX_VALUE_GATE_FAIL:{c}"
        )

    matrices[c] = X

print(
    "COMMON_ANCHOR_MATRIX_GATE = PASS"
)

# --------------------------------------------------
# Column map
# --------------------------------------------------

cols = []

for lib in LIBRARIES:
    cols.extend(
        [
            f"{lib}_real",
            f"{lib}_pseudo",
        ]
    )

ci = {
    name:i
    for i,name in enumerate(cols)
}

# --------------------------------------------------
# D estimator
# --------------------------------------------------

def dvec(M,pdlib,inlib):

    Rpd = M[:,ci[f"{pdlib}_real"]]
    Rin = M[:,ci[f"{inlib}_real"]]

    Ppd = M[:,ci[f"{pdlib}_pseudo"]]
    Pin = M[:,ci[f"{inlib}_pseudo"]]

    if (
        (Rpd <= 0).any()
        or
        (Rin <= 0).any()
        or
        (Ppd <= 0).any()
        or
        (Pin <= 0).any()
    ):
        raise SystemExit(
            "ZERO_AGGREGATE_GATE_FAIL"
        )

    return (
        np.log2(Rpd/Rin)
        -
        np.log2(Ppd/Pin)
    )

# --------------------------------------------------
# Observed effects
# --------------------------------------------------

observed_rows = []

for c in CENTERS:

    mean = matrices[c].mean(
        axis=0,
        keepdims=True
    )

    for loc,rep,pdlib,inlib in CONTRASTS:

        D = float(
            dvec(
                mean,
                pdlib,
                inlib
            )[0]
        )

        observed_rows.append({
            "center":c,
            "location":loc,
            "replicate":rep,
            "D":D,
        })

obs = pd.DataFrame(
    observed_rows
)

# --------------------------------------------------
# Joint common-anchor bootstrap
# --------------------------------------------------

print(
    "BOOTSTRAP_START =",
    BOOT_N
)

rng = np.random.default_rng(
    BOOT_SEED
)

prob = np.full(
    N,
    1.0/N
)

boot = {
    (c,loc):[]
    for c in CENTERS
    for loc in [
        "ER",
        "mitochondria",
        "cytosol",
    ]
}

done = 0

while done < BOOT_N:

    b = min(
        BOOT_BATCH,
        BOOT_N-done
    )

    W = rng.multinomial(
        N,
        prob,
        size=b
    )

    for c in CENTERS:

        M = (
            W @ matrices[c]
        ) / N

        vals = {}

        for loc,rep,pdlib,inlib in CONTRASTS:

            vals[(loc,rep)] = dvec(
                M,
                pdlib,
                inlib
            )

        for loc in [
            "ER",
            "mitochondria",
            "cytosol",
        ]:

            boot[(c,loc)].append(
                (
                    vals[(loc,1)]
                    +
                    vals[(loc,2)]
                ) / 2
            )

    done += b

    if (
        done % 1000 == 0
        or
        done == BOOT_N
    ):
        print(
            "BOOTSTRAP_COMPLETED =",
            done,
            "/",
            BOOT_N
        )

# --------------------------------------------------
# Summary table
# --------------------------------------------------

rows = []

for c in CENTERS:

    for loc in [
        "ER",
        "mitochondria",
        "cytosol",
    ]:

        z = obs[
            (obs["center"].eq(c))
            &
            (obs["location"].eq(loc))
        ]

        d1 = float(
            z[
                z["replicate"].eq(1)
            ]["D"].iloc[0]
        )

        d2 = float(
            z[
                z["replicate"].eq(2)
            ]["D"].iloc[0]
        )

        vals = np.concatenate(
            boot[(c,loc)]
        )

        lo = float(
            np.percentile(
                vals,
                2.5
            )
        )

        hi = float(
            np.percentile(
                vals,
                97.5
            )
        )

        rows.append({
            "center":c,
            "N_common_anchors":N,
            "location":loc,
            "D_rep1":d1,
            "D_rep2":d2,
            "D_mean":(d1+d2)/2,
            "CI_2.5":lo,
            "CI_97.5":hi,
            "replicate_concordance":
                (
                    "SAME_SIGN"
                    if np.sign(d1)==np.sign(d2)
                    else "OPPOSITE_SIGN"
                ),
            "CI_excludes_zero":
                (
                    "YES"
                    if lo > 0 or hi < 0
                    else "NO"
                ),
        })

out = pd.DataFrame(
    rows
)

out.to_csv(
    OUTDIR /
    "COMMON_ANCHOR_11CENTER_PROFILE_v1.tsv",
    sep="\t",
    index=False
)

pd.DataFrame(
    common,
    columns=[
        "chrom",
        "strand",
        "real_anchor_genomic_start0",
    ]
).to_csv(
    OUTDIR /
    "COMMON_ANCHOR_UNIVERSE_v1.tsv",
    sep="\t",
    index=False
)

# --------------------------------------------------
# Reveal
# --------------------------------------------------

print()
print(
    "==========================================="
)
print(
    "COMMON-ANCHOR POSITIONAL SENSITIVITY"
)
print(
    "==========================================="
)

for loc in [
    "ER",
    "mitochondria",
    "cytosol",
]:

    print()
    print(
        "LOCATION =",
        loc
    )

    z = out[
        out["location"].eq(loc)
    ]

    for _,r in z.iterrows():

        print(
            f"CENTER={int(r['center']):+d}",
            f"N={N}",
            f"D1={r['D_rep1']:.12f}",
            f"D2={r['D_rep2']:.12f}",
            f"Dmean={r['D_mean']:.12f}",
            f"CI95=[{r['CI_2.5']:.12f},"
            f"{r['CI_97.5']:.12f}]",
            f"CONCORDANCE={r['replicate_concordance']}",
            f"CI_EXCLUDES_0={r['CI_excludes_zero']}",
        )

print()
print(
    "COMMON_ANCHOR_SENSITIVITY_GATE = PASS"
)
