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
    "tmd_spatial_first_internal_v1/output"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

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
# Load frozen paired/Crossmap audit
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

# Same eligibility rule as secondary scan
d = d[
    (d["real_crossmap_unmasked_nt"] >= 1)
    &
    (d["pseudo_min_unmasked_nt"] >= 1)
].copy()

allowed_classes = {
    "FIRST",
    "INTERNAL",
    "MIXED",
}

observed_classes = set(
    d["position_class_consensus"]
    .dropna()
    .astype(str)
    .unique()
)

if not observed_classes.issubset(
    allowed_classes
):
    raise SystemExit(
        "POSITION_CLASS_GATE_FAIL:"
        + str(sorted(observed_classes))
    )

print(
    "POSITION_CLASS_GATE = PASS"
)

# ============================================================
# Annotation-only class accounting
# ============================================================

account_rows = []

print()
print(
    "=== FIRST / INTERNAL ELIGIBILITY ==="
)

for c in CENTERS:

    z = d[
        d["center"].eq(c)
    ]

    nf = int(
        z[
            "position_class_consensus"
        ].eq("FIRST").sum()
    )

    ni = int(
        z[
            "position_class_consensus"
        ].eq("INTERNAL").sum()
    )

    nm = int(
        z[
            "position_class_consensus"
        ].eq("MIXED").sum()
    )

    total = len(z)

    if nf + ni + nm != total:
        raise SystemExit(
            f"CLASS_ACCOUNTING_FAIL:{c}"
        )

    print(
        f"CENTER={c:+d}",
        f"TOTAL={total}",
        f"FIRST={nf}",
        f"INTERNAL={ni}",
        f"MIXED={nm}",
    )

    account_rows.append({
        "center":c,
        "total_eligible":total,
        "FIRST":nf,
        "INTERNAL":ni,
        "MIXED":nm,
    })


account = pd.DataFrame(
    account_rows
)

account.to_csv(
    OUTDIR /
    "FIRST_INTERNAL_ELIGIBILITY_ACCOUNTING_v1.tsv",
    sep="\t",
    index=False
)

print(
    "CLASS_ACCOUNTING_GATE = PASS"
)

# ============================================================
# Frozen geometry catalog
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
# Reuse frozen secondary density checkpoints
# ============================================================

density = {}

for lib in LIBRARIES:

    p = (
        COUNTDIR /
        f"{lib}.SECONDARY_GEOMETRY_DENSITIES_v1.tsv"
    )

    if not p.exists():
        raise SystemExit(
            f"MISSING_CHECKPOINT:{lib}:{p}"
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
# Helpers
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
    for i,name in enumerate(columns)
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

            pkey = (
                str(r["strand"]),
                geom,
                bits,
            )

            if pkey not in key_to_gid:
                raise SystemExit(
                    "PSEUDO_GEOMETRY_NOT_IN_CATALOG"
                )

            pseudo_gids.append(
                key_to_gid[pkey]
            )

        row = []

        for lib in LIBRARIES:

            rd = float(
                density[lib][real_gid]
            )

            pdens = float(
                np.mean(
                    [
                        density[lib][gid]
                        for gid in pseudo_gids
                    ]
                )
            )

            row.extend(
                [
                    rd,
                    pdens,
                ]
            )

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


def d_vector(
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
            Rpd / Rin
        )
        -
        np.log2(
            Ppd / Pin
        )
    )


# ============================================================
# Build class/center matrices
# ============================================================

matrices = {}

for c in CENTERS:

    for cls in CLASSES:

        z = d[
            (
                d["center"].eq(c)
            )
            &
            (
                d[
                    "position_class_consensus"
                ].eq(cls)
            )
        ].copy()

        N = len(z)

        if N == 0:
            raise SystemExit(
                f"EMPTY_STRATUM:{c}:{cls}"
            )

        X = build_matrix(z)

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
    "FIRST_INTERNAL_MATRIX_GATE = PASS"
)

# ============================================================
# Observed estimates
# ============================================================

observed = {}

for c in CENTERS:

    for cls in CLASSES:

        X = matrices[
            (c,cls)
        ]

        M = X.mean(
            axis=0,
            keepdims=True
        )

        for loc,rep,pdlib,inlib in CONTRASTS:

            D = float(
                d_vector(
                    M,
                    pdlib,
                    inlib
                )[0]
            )

            observed[
                (c,cls,loc,rep)
            ] = D


print(
    "OBSERVED_STRATIFIED_EFFECT_GATE = PASS"
)

# ============================================================
# Bootstrap each center and class
# ============================================================

boot_location = {}

for center_index,c in enumerate(
    CENTERS
):

    for class_index,cls in enumerate(
        CLASSES
    ):

        X = matrices[
            (c,cls)
        ]

        N = X.shape[0]

        # Deterministic independent stream
        seed = (
            BOOT_SEED_BASE
            +
            1000 * center_index
            +
            100 * class_index
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
                ] = d_vector(
                    M,
                    pdlib,
                    inlib
                )

            for loc in LOCATIONS:

                store[loc].append(
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
    "FIRST_INTERNAL_BOOTSTRAP_GATE = PASS"
)

# ============================================================
# Class-specific summary
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
    "FIRST_INTERNAL_11CENTER_PROFILE_v1.tsv",
    sep="\t",
    index=False
)

# ============================================================
# FIRST minus INTERNAL contrast
# ============================================================

fi_rows = []

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

        delta1 = f1-i1
        delta2 = f2-i2

        delta_mean = (
            delta1+delta2
        ) / 2

        bf = boot_location[
            (c,"FIRST",loc)
        ]

        bi = boot_location[
            (c,"INTERNAL",loc)
        ]

        if len(bf) != len(bi):
            raise SystemExit(
                "BOOTSTRAP_LENGTH_GATE_FAIL"
            )

        bdelta = bf-bi

        lo = float(
            np.percentile(
                bdelta,
                2.5
            )
        )

        hi = float(
            np.percentile(
                bdelta,
                97.5
            )
        )

        fi_rows.append({
            "center":c,
            "location":loc,

            "FIRST_D_rep1":
                f1,

            "FIRST_D_rep2":
                f2,

            "INTERNAL_D_rep1":
                i1,

            "INTERNAL_D_rep2":
                i2,

            "FIRST_minus_INTERNAL_rep1":
                delta1,

            "FIRST_minus_INTERNAL_rep2":
                delta2,

            "FIRST_minus_INTERNAL_mean":
                delta_mean,

            "CI_2.5":
                lo,

            "CI_97.5":
                hi,

            "CI_excludes_zero":
                (
                    "YES"
                    if lo > 0
                    or hi < 0
                    else
                    "NO"
                ),
        })


fi = pd.DataFrame(
    fi_rows
)

fi.to_csv(
    OUTDIR /
    "FIRST_MINUS_INTERNAL_11CENTER_CONTRAST_v1.tsv",
    sep="\t",
    index=False
)

# ============================================================
# Formal reveal
# ============================================================

print()
print(
    "============================================"
)
print(
    "FIRST VS INTERNAL TMD STRATIFIED REVEAL"
)
print(
    "============================================"
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

        delta = fi[
            (
                fi["center"].eq(c)
            )
            &
            (
                fi["location"].eq(
                    loc
                )
            )
        ].iloc[0]

        print(
            f"CENTER={c:+d}",
            f"NFIRST={int(fs['N_units'])}",
            f"FIRST={fs['D_mean']:.12f}",
            f"FI95=[{fs['CI_2.5']:.12f},"
            f"{fs['CI_97.5']:.12f}]",
            f"NINTERNAL={int(ins['N_units'])}",
            f"INTERNAL={ins['D_mean']:.12f}",
            f"II95=[{ins['CI_2.5']:.12f},"
            f"{ins['CI_97.5']:.12f}]",
            f"FIRST_MINUS_INTERNAL="
            f"{delta['FIRST_minus_INTERNAL_mean']:.12f}",
            f"DELTA95=[{delta['CI_2.5']:.12f},"
            f"{delta['CI_97.5']:.12f}]",
            f"DELTA_CI_EXCLUDES_0="
            f"{delta['CI_excludes_zero']}",
        )

print()
print(
    "FIRST_INTERNAL_STRATIFICATION_GATE = PASS"
)

# ============================================================
# Output SHA
# ============================================================

import hashlib

for path in [
    OUTDIR /
    "FIRST_INTERNAL_ELIGIBILITY_ACCOUNTING_v1.tsv",

    OUTDIR /
    "FIRST_INTERNAL_11CENTER_PROFILE_v1.tsv",

    OUTDIR /
    "FIRST_MINUS_INTERNAL_11CENTER_CONTRAST_v1.tsv",
]:

    sha = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()

    print(
        "SHA256",
        sha,
        path
    )
