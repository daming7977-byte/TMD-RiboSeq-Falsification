from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

FORMAL = Path(
    "tmd_anchor_preregister/bridge/"
    "FORMAL_TMD_ANCHOR_CANDIDATES_v1p1.tsv"
)

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
    "tmd_spatial_multipass_first_internal_v1/output"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

CENTERS = [
    -30,-15,0,15,30,45,60,75,90,105,120
]

CLASSES = [
    "FIRST",
    "INTERNAL",
]

LOCATIONS = [
    "ER",
    "mitochondria",
    "cytosol",
]

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
BOOT_SEED_BASE = 30097775
BOOT_BATCH = 200


# ============================================================
# 1. AUTHORITATIVE target transcript -> pass_class map
# ============================================================

f = pd.read_csv(
    FORMAL,
    sep="\t",
    dtype=str
)

required_formal = {
    "target_transcript",
    "pass_class",
}

missing = required_formal - set(f.columns)

if missing:
    raise SystemExit(
        f"FORMAL_COLUMN_GATE_FAIL:{sorted(missing)}"
    )

ft = f[
    [
        "target_transcript",
        "pass_class",
    ]
].dropna().copy()

allowed_pass = {
    "SINGLE_PASS",
    "MULTI_PASS",
}

observed_pass = set(
    ft["pass_class"]
    .astype(str)
    .unique()
)

print(
    "FORMAL_PASS_CLASSES =",
    sorted(observed_pass)
)

if not observed_pass.issubset(
    allowed_pass
):
    raise SystemExit(
        "PASS_CLASS_VALUE_GATE_FAIL"
    )

class_nunique = (
    ft.groupby(
        "target_transcript"
    )["pass_class"]
    .nunique()
)

bad_tx = int(
    (class_nunique > 1).sum()
)

print(
    "TARGET_TRANSCRIPT_PASS_CLASS_CONFLICTS =",
    bad_tx
)

if bad_tx != 0:
    raise SystemExit(
        "TRANSCRIPT_PASS_CLASS_CONSISTENCY_GATE_FAIL"
    )

tx_class = (
    ft.drop_duplicates(
        "target_transcript"
    )
    .set_index(
        "target_transcript"
    )["pass_class"]
    .to_dict()
)

print(
    "TARGET_TRANSCRIPTS_WITH_PASS_CLASS =",
    len(tx_class)
)

print(
    "PASS_CLASS_MAPPING_GATE = PASS"
)


# ============================================================
# 2. LOAD FROZEN PAIRED/CROSSMAP UNITS
# ============================================================

d = pd.read_csv(
    PAIR,
    sep="\t",
    dtype=str
)

d["center"] = pd.to_numeric(
    d["center"]
).astype(int)

d["real_crossmap_unmasked_nt"] = pd.to_numeric(
    d["real_crossmap_unmasked_nt"]
)

d["pseudo_min_unmasked_nt"] = pd.to_numeric(
    d["pseudo_min_unmasked_nt"]
)

# Frozen secondary eligibility
d = d[
    (d["real_crossmap_unmasked_nt"] >= 1)
    &
    (d["pseudo_min_unmasked_nt"] >= 1)
].copy()


def split_csv(x):
    return [
        y
        for y in str(x).split(",")
        if y != ""
    ]


def topology_class(row):

    txs = split_csv(
        row["target_transcripts"]
    )

    if len(txs) == 0:
        raise ValueError(
            "EMPTY_TARGET_TRANSCRIPT_LIST"
        )

    classes = []

    for tx in txs:

        if tx not in tx_class:
            raise KeyError(
                f"NO_PASS_CLASS_FOR_TARGET:{tx}"
            )

        classes.append(
            tx_class[tx]
        )

    classes = sorted(
        set(classes)
    )

    if classes == ["MULTI_PASS"]:
        return "MULTI_PASS"

    if classes == ["SINGLE_PASS"]:
        return "SINGLE_PASS"

    return "MIXED_TOPOLOGY"


d["topology_class"] = d.apply(
    topology_class,
    axis=1
)

print()
print(
    "=== TOPOLOGY ACCOUNTING BY CENTER ==="
)

account_rows = []

for c in CENTERS:

    z = d[
        d["center"].eq(c)
    ]

    nmulti = int(
        z["topology_class"]
        .eq("MULTI_PASS")
        .sum()
    )

    nsingle = int(
        z["topology_class"]
        .eq("SINGLE_PASS")
        .sum()
    )

    nmixedtop = int(
        z["topology_class"]
        .eq("MIXED_TOPOLOGY")
        .sum()
    )

    zm = z[
        z["topology_class"]
        .eq("MULTI_PASS")
    ]

    nfirst = int(
        zm[
            "position_class_consensus"
        ]
        .eq("FIRST")
        .sum()
    )

    ninternal = int(
        zm[
            "position_class_consensus"
        ]
        .eq("INTERNAL")
        .sum()
    )

    nmixedpos = int(
        zm[
            "position_class_consensus"
        ]
        .eq("MIXED")
        .sum()
    )

    if (
        nmulti
        + nsingle
        + nmixedtop
        != len(z)
    ):
        raise SystemExit(
            f"TOPOLOGY_ACCOUNTING_FAIL:{c}"
        )

    if (
        nfirst
        + ninternal
        + nmixedpos
        != nmulti
    ):
        raise SystemExit(
            f"MULTIPASS_POSITION_ACCOUNTING_FAIL:{c}"
        )

    print(
        f"CENTER={c:+d}",
        f"TOTAL={len(z)}",
        f"MULTI_PASS={nmulti}",
        f"SINGLE_PASS={nsingle}",
        f"MIXED_TOPOLOGY={nmixedtop}",
        f"MULTI_FIRST={nfirst}",
        f"MULTI_INTERNAL={ninternal}",
        f"MULTI_MIXED_POSITION={nmixedpos}",
    )

    account_rows.append({
        "center":c,
        "total_eligible":len(z),
        "MULTI_PASS":nmulti,
        "SINGLE_PASS":nsingle,
        "MIXED_TOPOLOGY":nmixedtop,
        "MULTI_FIRST":nfirst,
        "MULTI_INTERNAL":ninternal,
        "MULTI_MIXED_POSITION":
            nmixedpos,
    })


account = pd.DataFrame(
    account_rows
)

account.to_csv(
    OUTDIR /
    "MULTIPASS_TOPOLOGY_ACCOUNTING_v1.tsv",
    sep="\t",
    index=False
)

print(
    "TOPOLOGY_ACCOUNTING_GATE = PASS"
)


# ============================================================
# 3. RETAIN MULTI-PASS FIRST / INTERNAL ONLY
# ============================================================

m = d[
    d["topology_class"]
    .eq("MULTI_PASS")
].copy()

m = m[
    m["position_class_consensus"]
    .isin(
        ["FIRST","INTERNAL"]
    )
].copy()

if len(m) == 0:
    raise SystemExit(
        "EMPTY_MULTIPASS_UNIVERSE"
    )

print()
print(
    "MULTIPASS_ANALYTIC_ROWS_ALL_CENTERS =",
    len(m)
)


# ============================================================
# 4. GEOMETRY CATALOG
# ============================================================

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


# ============================================================
# 5. REUSE SECONDARY DENSITY CHECKPOINTS
# ============================================================

density = {}

for lib in LIBRARIES:

    p = (
        COUNTDIR /
        f"{lib}.SECONDARY_GEOMETRY_DENSITIES_v1.tsv"
    )

    if not p.exists():
        raise SystemExit(
            f"MISSING_COUNT_CHECKPOINT:{lib}"
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


# ============================================================
# 6. HELPERS
# ============================================================

def splitsemi(x):

    return [
        z
        for z in str(x).split(";")
        if z != ""
    ]


columns = []

for lib in LIBRARIES:
    columns.extend([
        f"{lib}_real",
        f"{lib}_pseudo",
    ])

col_idx = {
    name:i
    for i,name in enumerate(
        columns
    )
}


def build_matrix(z):

    rows = []

    for _,r in z.iterrows():

        real_key = (
            str(r["strand"]),
            str(r["real_window_geometry"]),
            str(r["real_crossmap_mask_bits"]),
        )

        if real_key not in key_to_gid:
            raise SystemExit(
                "REAL_GEOMETRY_NOT_IN_CATALOG"
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
                "PSEUDO_GEOMETRY_MASK_CARDINALITY_FAIL"
            )

        pseudo_gids = []

        for geom,bits in zip(
            pgeoms,
            pbits
        ):

            key = (
                str(r["strand"]),
                geom,
                bits,
            )

            if key not in key_to_gid:
                raise SystemExit(
                    "PSEUDO_GEOMETRY_NOT_IN_CATALOG"
                )

            pseudo_gids.append(
                key_to_gid[key]
            )

        row = []

        for lib in LIBRARIES:

            real_density = float(
                density[lib][real_gid]
            )

            pseudo_density = float(
                np.mean(
                    [
                        density[lib][gid]
                        for gid in pseudo_gids
                    ]
                )
            )

            row.extend([
                real_density,
                pseudo_density,
            ])

        rows.append(row)

    X = np.asarray(
        rows,
        dtype=float
    )

    if (
        not np.isfinite(X).all()
        or
        (X < 0).any()
    ):
        raise SystemExit(
            "MATRIX_VALUE_GATE_FAIL"
        )

    return X


def D_vector(
    M,
    pd_lib,
    in_lib,
):

    Rpd = M[
        :,
        col_idx[
            f"{pd_lib}_real"
        ]
    ]

    Rin = M[
        :,
        col_idx[
            f"{in_lib}_real"
        ]
    ]

    Ppd = M[
        :,
        col_idx[
            f"{pd_lib}_pseudo"
        ]
    ]

    Pin = M[
        :,
        col_idx[
            f"{in_lib}_pseudo"
        ]
    ]

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
        np.log2(
            Rpd/Rin
        )
        -
        np.log2(
            Ppd/Pin
        )
    )


# ============================================================
# 7. BUILD MULTI-PASS FIRST / INTERNAL MATRICES
# ============================================================

matrices = {}

for c in CENTERS:

    for cls in CLASSES:

        z = m[
            (
                m["center"].eq(c)
            )
            &
            (
                m[
                    "position_class_consensus"
                ].eq(cls)
            )
        ].copy()

        N = len(z)

        if N == 0:
            raise SystemExit(
                f"EMPTY_MULTIPASS_STRATUM:"
                f"{c}:{cls}"
            )

        X = build_matrix(
            z
        )

        if X.shape != (
            N,
            24
        ):
            raise SystemExit(
                f"MATRIX_SHAPE_FAIL:"
                f"{c}:{cls}:{X.shape}"
            )

        matrices[
            (c,cls)
        ] = X


print(
    "MULTIPASS_MATRIX_GATE = PASS"
)


# ============================================================
# 8. OBSERVED D
# ============================================================

observed = {}

for c in CENTERS:

    for cls in CLASSES:

        X = matrices[
            (c,cls)
        ]

        mean = X.mean(
            axis=0,
            keepdims=True
        )

        for loc,rep,pdlib,inlib in CONTRASTS:

            D = float(
                D_vector(
                    mean,
                    pdlib,
                    inlib
                )[0]
            )

            observed[
                (c,cls,loc,rep)
            ] = D


print(
    "MULTIPASS_OBSERVED_EFFECT_GATE = PASS"
)


# ============================================================
# 9. BOOTSTRAP
# ============================================================

boot_location = {}

for ci_center,c in enumerate(
    CENTERS
):

    for ci_class,cls in enumerate(
        CLASSES
    ):

        X = matrices[
            (c,cls)
        ]

        N = X.shape[0]

        seed = (
            BOOT_SEED_BASE
            +
            1000*ci_center
            +
            100*ci_class
        )

        rng = np.random.default_rng(
            seed
        )

        prob = np.full(
            N,
            1.0/N
        )

        store = {
            loc:[]
            for loc in LOCATIONS
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

            M = (
                W @ X
            ) / N

            vals = {}

            for loc,rep,pdlib,inlib in CONTRASTS:

                vals[
                    (loc,rep)
                ] = D_vector(
                    M,
                    pdlib,
                    inlib
                )

            for loc in LOCATIONS:

                store[
                    loc
                ].append(
                    (
                        vals[(loc,1)]
                        +
                        vals[(loc,2)]
                    ) / 2
                )

            done += b

        for loc in LOCATIONS:

            boot_location[
                (c,cls,loc)
            ] = np.concatenate(
                store[loc]
            )

    print(
        "BOOTSTRAP_CENTER_COMPLETED =",
        f"{c:+d}"
    )


print(
    "MULTIPASS_BOOTSTRAP_GATE = PASS"
)


# ============================================================
# 10. CLASS-SPECIFIC SUMMARY
# ============================================================

summary_rows = []

for c in CENTERS:

    for cls in CLASSES:

        N = matrices[
            (c,cls)
        ].shape[0]

        for loc in LOCATIONS:

            d1 = observed[
                (c,cls,loc,1)
            ]

            d2 = observed[
                (c,cls,loc,2)
            ]

            dm = (
                d1+d2
            ) / 2

            vals = boot_location[
                (c,cls,loc)
            ]

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

            summary_rows.append({
                "center":c,
                "TMD_class":cls,
                "N_units":N,
                "location":loc,
                "D_rep1":d1,
                "D_rep2":d2,
                "D_mean":dm,
                "CI_2.5":lo,
                "CI_97.5":hi,
                "replicate_concordance":
                    (
                        "SAME_SIGN"
                        if np.sign(d1)
                        ==
                        np.sign(d2)
                        else
                        "OPPOSITE_SIGN"
                    ),
                "CI_excludes_zero":
                    (
                        "YES"
                        if lo > 0
                        or hi < 0
                        else
                        "NO"
                    ),
            })


summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    OUTDIR /
    "MULTIPASS_FIRST_INTERNAL_11CENTER_PROFILE_v1.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 11. FIRST - INTERNAL DIRECT CONTRAST
# ============================================================

delta_rows = []

for c in CENTERS:

    for loc in LOCATIONS:

        f1 = observed[
            (c,"FIRST",loc,1)
        ]

        f2 = observed[
            (c,"FIRST",loc,2)
        ]

        i1 = observed[
            (c,"INTERNAL",loc,1)
        ]

        i2 = observed[
            (c,"INTERNAL",loc,2)
        ]

        df1 = f1-i1
        df2 = f2-i2

        dfm = (
            df1+df2
        ) / 2

        bf = boot_location[
            (c,"FIRST",loc)
        ]

        bi = boot_location[
            (c,"INTERNAL",loc)
        ]

        delta = bf-bi

        lo = float(
            np.percentile(
                delta,
                2.5
            )
        )

        hi = float(
            np.percentile(
                delta,
                97.5
            )
        )

        delta_rows.append({
            "center":c,
            "location":loc,
            "FIRST_D_rep1":f1,
            "FIRST_D_rep2":f2,
            "INTERNAL_D_rep1":i1,
            "INTERNAL_D_rep2":i2,
            "FIRST_minus_INTERNAL_rep1":
                df1,
            "FIRST_minus_INTERNAL_rep2":
                df2,
            "FIRST_minus_INTERNAL_mean":
                dfm,
            "CI_2.5":lo,
            "CI_97.5":hi,
            "CI_excludes_zero":
                (
                    "YES"
                    if lo > 0
                    or hi < 0
                    else
                    "NO"
                ),
        })


delta = pd.DataFrame(
    delta_rows
)

delta.to_csv(
    OUTDIR /
    "MULTIPASS_FIRST_MINUS_INTERNAL_11CENTER_v1.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 12. REVEAL
# ============================================================

print()
print(
    "================================================"
)
print(
    "MULTI-PASS-ONLY FIRST VS INTERNAL SENSITIVITY"
)
print(
    "================================================"
)

for loc in LOCATIONS:

    print()
    print(
        "LOCATION =",
        loc
    )

    for c in CENTERS:

        fs = summary[
            (
                summary["center"].eq(c)
            )
            &
            (
                summary["TMD_class"].eq(
                    "FIRST"
                )
            )
            &
            (
                summary["location"].eq(
                    loc
                )
            )
        ].iloc[0]

        ins = summary[
            (
                summary["center"].eq(c)
            )
            &
            (
                summary["TMD_class"].eq(
                    "INTERNAL"
                )
            )
            &
            (
                summary["location"].eq(
                    loc
                )
            )
        ].iloc[0]

        de = delta[
            (
                delta["center"].eq(c)
            )
            &
            (
                delta["location"].eq(
                    loc
                )
            )
        ].iloc[0]

        print(
            f"CENTER={c:+d}",
            f"NFIRST_MULTI={int(fs['N_units'])}",
            f"FIRST_MULTI={fs['D_mean']:.12f}",
            f"FCI95=[{fs['CI_2.5']:.12f},"
            f"{fs['CI_97.5']:.12f}]",
            f"NINTERNAL_MULTI={int(ins['N_units'])}",
            f"INTERNAL_MULTI={ins['D_mean']:.12f}",
            f"ICI95=[{ins['CI_2.5']:.12f},"
            f"{ins['CI_97.5']:.12f}]",
            f"FIRST_MINUS_INTERNAL="
            f"{de['FIRST_minus_INTERNAL_mean']:.12f}",
            f"DELTA95=[{de['CI_2.5']:.12f},"
            f"{de['CI_97.5']:.12f}]",
            f"DELTA_CI_EXCLUDES_0="
            f"{de['CI_excludes_zero']}",
        )


print()
print(
    "MULTIPASS_FIRST_INTERNAL_SENSITIVITY_GATE = PASS"
)


# ============================================================
# 13. OUTPUT SHA
# ============================================================

for path in [
    OUTDIR /
    "MULTIPASS_TOPOLOGY_ACCOUNTING_v1.tsv",

    OUTDIR /
    "MULTIPASS_FIRST_INTERNAL_11CENTER_PROFILE_v1.tsv",

    OUTDIR /
    "MULTIPASS_FIRST_MINUS_INTERNAL_11CENTER_v1.tsv",
]:

    sha = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()

    print(
        "SHA256",
        sha,
        path
    )
