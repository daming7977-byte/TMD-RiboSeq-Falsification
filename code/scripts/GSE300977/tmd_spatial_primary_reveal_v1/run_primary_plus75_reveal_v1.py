from pathlib import Path
import hashlib
import gc
import math

import numpy as np
import pandas as pd

from plastid import (
    BAMGenomeArray,
    CenterMapFactory,
    SizeFilterFactory,
    GenomicSegment,
    SegmentChain,
)

# ============================================================
# FROZEN INPUTS
# ============================================================

ELIG = Path(
    "crossmap_window_preregister/output/"
    "GSE300977_PRIMARY_PLUS75_PAIRED_ELIGIBLE_v1.tsv"
)

ESTIMAND = Path(
    "tmd_spatial_analysis_preregister/"
    "GSE300977_PRIMARY_SPATIAL_ESTIMAND_CONTRACT_v1.txt"
)

CLOSURE = Path(
    "tmd_spatial_analysis_preregister/audit/"
    "GSE300977_PRIMARY_PREREGISTRATION_CLOSURE_v1.txt"
)

IMPLEMENTATION = Path(
    "tmd_spatial_analysis_preregister/"
    "GSE300977_PRIMARY_COUNTING_IMPLEMENTATION_CONTRACT_v1.txt"
)

OUTDIR = Path("tmd_spatial_primary_reveal_v1")
COUNTDIR = OUTDIR / "geometry_counts"
RESULTDIR = OUTDIR / "output"
AUDITDIR = OUTDIR / "audit"

for p in [COUNTDIR, RESULTDIR, AUDITDIR]:
    p.mkdir(parents=True, exist_ok=True)

EXPECTED_ELIG_SHA = (
    "ba47f5a42f9923743723cf7d6a0b1483a27acaa091b0e32e8954f3f0005f7ab5"
)

EXPECTED_ESTIMAND_SHA = (
    "6faa60107a5034f38aff5d1b586784dbb315335e84897567d16620a0d051f7eb"
)

EXPECTED_CLOSURE_SHA = (
    "48f952b8719c0e21ee168b6aa74c7f33c4c8264abc8e06a281cd099b1ae43bb1"
)

N_UNITS = 9098
BOOTSTRAP_N = 10000
BOOTSTRAP_SEED = 30097775
BOOTSTRAP_BATCH = 200

LIBRARIES = [
    "B01","B02","B03","B04","B05","B06",
    "B07","B08","B09","B10","B11","B12",
]

BAMS = {
    b: Path(
        f"unique_mapping_blinded/bam/{b}.NH1.unique.bam"
    )
    for b in LIBRARIES
}

CONTRASTS = [
    ("ER",           1, "B01", "B02"),
    ("ER",           2, "B11", "B12"),
    ("mitochondria", 1, "B09", "B10"),
    ("mitochondria", 2, "B07", "B08"),
    ("cytosol",      1, "B05", "B06"),
    ("cytosol",      2, "B03", "B04"),
]


def sha256_file(path):
    h = hashlib.sha256()

    with path.open("rb") as fh:
        for block in iter(
            lambda: fh.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


# ============================================================
# HARD PRE-REVEAL SHA GATES
# ============================================================

print("=== FROZEN INPUT SHA GATES ===")

checks = [
    ("ELIGIBLE", ELIG, EXPECTED_ELIG_SHA),
    ("ESTIMAND", ESTIMAND, EXPECTED_ESTIMAND_SHA),
    ("CLOSURE", CLOSURE, EXPECTED_CLOSURE_SHA),
]

for name, path, expected in checks:

    actual = sha256_file(path)

    print(
        name,
        "SHA256 =",
        actual
    )

    if actual != expected:
        raise SystemExit(
            f"{name}_SHA_GATE = FAIL"
        )

    print(
        f"{name}_SHA_GATE = PASS"
    )

implementation_sha = sha256_file(
    IMPLEMENTATION
)

print(
    "IMPLEMENTATION_CONTRACT_SHA256 =",
    implementation_sha
)

# ============================================================
# LOAD FROZEN 9098-UNIT UNIVERSE
# ============================================================

u = pd.read_csv(
    ELIG,
    sep="\t",
    dtype=str,
)

if len(u) != N_UNITS:
    raise SystemExit(
        f"UNIT_COUNT_GATE_FAIL:{len(u)}"
    )

if u["paired_real_unit_id"].duplicated().any():
    raise SystemExit(
        "UNIT_ID_DUPLICATION_GATE_FAIL"
    )

required = {
    "paired_real_unit_id",
    "strand",
    "real_window_geometry",
    "pseudo_window_geometries",
    "real_crossmap_mask_bits",
    "pseudo_crossmap_mask_bits",
    "real_crossmap_unmasked_nt",
    "pseudo_min_unmasked_nt",
}

missing = required - set(u.columns)

if missing:
    raise SystemExit(
        f"ELIGIBLE_COLUMN_GATE_FAIL:{sorted(missing)}"
    )

print()
print("PRIMARY_UNIT_ROWS =", len(u))
print("PRIMARY_UNIT_ID_GATE = PASS")

# ============================================================
# GEOMETRY CATALOG
# ============================================================

def split_semicolon(x):
    return [
        z
        for z in str(x).split(";")
        if z != ""
    ]


geometry_keys = set()
unit_specs = []

for _, r in u.iterrows():

    strand = str(r["strand"])

    real_geom = str(
        r["real_window_geometry"]
    )

    real_bits = str(
        r["real_crossmap_mask_bits"]
    )

    if (
        len(real_bits) != 31
        or set(real_bits) - {"0","1"}
    ):
        raise SystemExit(
            "REAL_MASK_VECTOR_GATE_FAIL"
        )

    real_key = (
        strand,
        real_geom,
        real_bits,
    )

    geometry_keys.add(real_key)

    pseudo_geoms = split_semicolon(
        r["pseudo_window_geometries"]
    )

    pseudo_bits = split_semicolon(
        r["pseudo_crossmap_mask_bits"]
    )

    if len(pseudo_geoms) != len(pseudo_bits):
        raise SystemExit(
            "PSEUDO_GEOMETRY_MASK_CARDINALITY_FAIL"
        )

    if len(pseudo_geoms) not in (1,2):
        raise SystemExit(
            "PSEUDO_CARDINALITY_GATE_FAIL:"
            f"{len(pseudo_geoms)}"
        )

    pseudo_keys = []

    for geom, bits in zip(
        pseudo_geoms,
        pseudo_bits,
    ):

        if (
            len(bits) != 31
            or set(bits) - {"0","1"}
        ):
            raise SystemExit(
                "PSEUDO_MASK_VECTOR_GATE_FAIL"
            )

        key = (
            strand,
            geom,
            bits,
        )

        geometry_keys.add(key)
        pseudo_keys.append(key)

    unit_specs.append({
        "paired_real_unit_id":
            r["paired_real_unit_id"],

        "real_key":
            real_key,

        "pseudo_keys":
            pseudo_keys,
    })


geometry_keys = sorted(
    geometry_keys
)

key_to_gid = {
    key: f"GEO{i:06d}"
    for i,key in enumerate(
        geometry_keys,
        start=1
    )
}

catalog_rows = []

for key in geometry_keys:

    strand, geom, bits = key

    unmasked_nt = bits.count("0")

    if unmasked_nt < 1:
        raise SystemExit(
            "FULLY_MASKED_GEOMETRY_IN_ELIGIBLE_UNIVERSE"
        )

    catalog_rows.append({
        "geometry_id":
            key_to_gid[key],

        "strand":
            strand,

        "geometry":
            geom,

        "crossmap_mask_bits":
            bits,

        "unmasked_nt":
            unmasked_nt,
    })

catalog = pd.DataFrame(
    catalog_rows
)

CATALOG_OUT = (
    RESULTDIR /
    "PRIMARY_PLUS75_GEOMETRY_CATALOG_v1.tsv"
)

catalog.to_csv(
    CATALOG_OUT,
    sep="\t",
    index=False
)

unit_map_rows = []

for spec in unit_specs:

    unit_map_rows.append({
        "paired_real_unit_id":
            spec["paired_real_unit_id"],

        "real_geometry_id":
            key_to_gid[
                spec["real_key"]
            ],

        "pseudo_geometry_ids":
            ";".join(
                key_to_gid[k]
                for k in spec[
                    "pseudo_keys"
                ]
            ),
    })

unit_map = pd.DataFrame(
    unit_map_rows
)

UNIT_MAP_OUT = (
    RESULTDIR /
    "PRIMARY_PLUS75_UNIT_GEOMETRY_MAP_v1.tsv"
)

unit_map.to_csv(
    UNIT_MAP_OUT,
    sep="\t",
    index=False
)

print(
    "UNIQUE_MASKED_GEOMETRIES_TO_COUNT =",
    len(catalog)
)

print(
    "PRECOUT_GEOMETRY_CATALOG_GATE = PASS"
)

# ============================================================
# GEOMETRY -> SegmentChain
# ============================================================

def geometry_to_chain(
    geometry,
    strand,
):

    segments = []

    for token in str(
        geometry
    ).split("^"):

        chrom, coords = token.split(
            ":",
            1
        )

        start, end = map(
            int,
            coords.split("-")
        )

        segments.append(
            GenomicSegment(
                chrom,
                start,
                end,
                strand
            )
        )

    chain = SegmentChain(
        *segments
    )

    if chain.length != 31:
        raise ValueError(
            f"CHAIN_LENGTH_FAIL:"
            f"{geometry}:{chain.length}"
        )

    return chain


chains = {}

for _, r in catalog.iterrows():

    chains[
        r["geometry_id"]
    ] = geometry_to_chain(
        r["geometry"],
        r["strand"],
    )

print(
    "SEGMENTCHAIN_CONSTRUCTION_GATE = PASS"
)

# ============================================================
# COUNT ONE LIBRARY
# ============================================================

def count_library(
    library,
):

    bam = BAMS[library]

    if not bam.exists():
        raise SystemExit(
            f"BAM_MISSING:{library}:{bam}"
        )

    bai = Path(str(bam) + ".bai")

    if not bai.exists():
        raise SystemExit(
            f"BAI_MISSING:{library}:{bai}"
        )

    outfile = (
        COUNTDIR /
        f"{library}.PRIMARY_PLUS75_GEOMETRY_COUNTS_v1.tsv"
    )

    # Safe restart:
    # only accept a prior complete table with exact geometry IDs.
    if outfile.exists():

        old = pd.read_csv(
            outfile,
            sep="\t",
        )

        if (
            len(old) != len(catalog)
            or
            old[
                "geometry_id"
            ].tolist()
            !=
            catalog[
                "geometry_id"
            ].tolist()
        ):
            raise SystemExit(
                f"STALE_OR_INCOMPLETE_CHECKPOINT:"
                f"{library}"
            )

        if old[
            ["unmasked_mass","density"]
        ].isna().any().any():
            raise SystemExit(
                f"BAD_CHECKPOINT_VALUES:{library}"
            )

        print(
            f"{library} CHECKPOINT = REUSED"
        )

        return old

    print()
    print(
        f"=== COUNTING {library} ==="
    )

    ga = BAMGenomeArray(
        str(bam)
    )

    ga.set_mapping(
        CenterMapFactory(
            nibble=13
        )
    )

    ga.add_filter(
        "length_27_36",
        SizeFilterFactory(
            min=27,
            max=36
        )
    )

    rows = []

    for j, r in catalog.iterrows():

        gid = r["geometry_id"]
        bits = str(
            r["crossmap_mask_bits"]
        )

        chain = chains[gid]

        counts = np.asarray(
            chain.get_counts(
                ga,
                stranded=True
            ),
            dtype=float
        )

        if len(counts) != 31:
            raise SystemExit(
                f"COUNT_VECTOR_LENGTH_FAIL:"
                f"{library}:{gid}:{len(counts)}"
            )

        if (
            not np.isfinite(counts).all()
            or
            (counts < 0).any()
        ):
            raise SystemExit(
                f"COUNT_VECTOR_VALUE_FAIL:"
                f"{library}:{gid}"
            )

        keep = np.fromiter(
            (
                b == "0"
                for b in bits
            ),
            dtype=bool,
            count=31
        )

        unmasked_nt = int(
            keep.sum()
        )

        if unmasked_nt != int(
            r["unmasked_nt"]
        ):
            raise SystemExit(
                f"UNMASKED_NT_GATE_FAIL:"
                f"{library}:{gid}"
            )

        mass = float(
            counts[keep].sum()
        )

        density = (
            mass /
            unmasked_nt
        )

        rows.append({
            "geometry_id":
                gid,

            "unmasked_nt":
                unmasked_nt,

            "unmasked_mass":
                mass,

            "density":
                density,
        })

        if (
            (j + 1) % 1000 == 0
            or
            (j + 1) == len(catalog)
        ):
            print(
                library,
                "GEOMETRIES_COUNTED =",
                j + 1,
                "/",
                len(catalog)
            )

    out = pd.DataFrame(rows)

    if len(out) != len(catalog):
        raise SystemExit(
            f"LIBRARY_OUTPUT_ROW_FAIL:{library}"
        )

    if out[
        ["unmasked_mass","density"]
    ].isna().any().any():
        raise SystemExit(
            f"LIBRARY_NAN_FAIL:{library}"
        )

    out.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    print(
        f"{library} COUNTING_GATE = PASS"
    )

    del ga
    gc.collect()

    return out


# ============================================================
# FORMAL 12-LIBRARY COUNTING
# ============================================================

geometry_density = {}

for lib in LIBRARIES:

    libdf = count_library(
        lib
    )

    geometry_density[lib] = dict(
        zip(
            libdf[
                "geometry_id"
            ],
            libdf[
                "density"
            ].astype(float)
        )
    )


print()
print(
    "FORMAL_12_LIBRARY_GEOMETRY_COUNTING = PASS"
)

# ============================================================
# BUILD 9098 x 24 UNIT DENSITY MATRIX
# ============================================================

unit_rows = []

for spec in unit_specs:

    uid = spec[
        "paired_real_unit_id"
    ]

    real_gid = key_to_gid[
        spec["real_key"]
    ]

    pseudo_gids = [
        key_to_gid[k]
        for k in spec[
            "pseudo_keys"
        ]
    ]

    row = {
        "paired_real_unit_id":
            uid,

        "real_geometry_id":
            real_gid,

        "pseudo_geometry_ids":
            ";".join(
                pseudo_gids
            ),

        "n_distinct_pseudos":
            len(pseudo_gids),
    }

    for lib in LIBRARIES:

        rd = float(
            geometry_density[
                lib
            ][real_gid]
        )

        pds = [
            float(
                geometry_density[
                    lib
                ][gid]
            )
            for gid in pseudo_gids
        ]

        pd_mean = float(
            np.mean(pds)
        )

        row[
            f"{lib}_real_density"
        ] = rd

        row[
            f"{lib}_pseudo_density"
        ] = pd_mean

        row[
            f"{lib}_pseudo_component_densities"
        ] = ";".join(
            f"{x:.17g}"
            for x in pds
        )

    unit_rows.append(row)


unit_density = pd.DataFrame(
    unit_rows
)

if len(unit_density) != N_UNITS:
    raise SystemExit(
        "UNIT_DENSITY_ROW_GATE_FAIL"
    )

UNIT_DENSITY_OUT = (
    RESULTDIR /
    "PRIMARY_PLUS75_UNIT_DENSITIES_v1.tsv"
)

unit_density.to_csv(
    UNIT_DENSITY_OUT,
    sep="\t",
    index=False
)

print(
    "UNIT_DENSITY_MATRIX_GATE = PASS"
)

# ============================================================
# NUMERIC MATRIX
# ============================================================

matrix_cols = []

for lib in LIBRARIES:
    matrix_cols.extend([
        f"{lib}_real_density",
        f"{lib}_pseudo_density",
    ])

X = unit_density[
    matrix_cols
].to_numpy(
    dtype=float
)

if X.shape != (
    N_UNITS,
    24
):
    raise SystemExit(
        f"NUMERIC_MATRIX_SHAPE_FAIL:{X.shape}"
    )

if (
    not np.isfinite(X).all()
    or
    (X < 0).any()
):
    raise SystemExit(
        "NUMERIC_MATRIX_VALUE_GATE_FAIL"
    )

col_index = {
    col: i
    for i,col in enumerate(
        matrix_cols
    )
}

observed_means = (
    X.mean(axis=0)
)

if (
    observed_means <= 0
).any():
    zero_cols = [
        matrix_cols[i]
        for i,v in enumerate(
            observed_means
        )
        if v <= 0
    ]

    raise SystemExit(
        "PRIMARY_AGGREGATE_ZERO_GATE_FAIL:"
        + ",".join(zero_cols)
    )

print(
    "PRIMARY_AGGREGATE_ZERO_GATE = PASS"
)

# ============================================================
# EFFECT FUNCTIONS
# ============================================================

def get_mean(
    means,
    lib,
    kind,
):
    return means[
        col_index[
            f"{lib}_{kind}_density"
        ]
    ]


def one_contrast(
    means,
    pd_lib,
    in_lib,
):

    R_pd = get_mean(
        means,
        pd_lib,
        "real"
    )

    R_in = get_mean(
        means,
        in_lib,
        "real"
    )

    P_pd = get_mean(
        means,
        pd_lib,
        "pseudo"
    )

    P_in = get_mean(
        means,
        in_lib,
        "pseudo"
    )

    vals = [
        R_pd,R_in,P_pd,P_in
    ]

    if any(v <= 0 for v in vals):
        raise ValueError(
            "ZERO_AGGREGATE_IN_CONTRAST"
        )

    E_real = math.log2(
        R_pd / R_in
    )

    E_pseudo = math.log2(
        P_pd / P_in
    )

    D = E_real - E_pseudo

    return {
        "R_pulldown":
            R_pd,

        "R_input":
            R_in,

        "P_pulldown":
            P_pd,

        "P_input":
            P_in,

        "E_real":
            E_real,

        "E_pseudo":
            E_pseudo,

        "D":
            D,
    }


contrast_rows = []

D_obs = {}

for loc,rep,pd_lib,in_lib in CONTRASTS:

    z = one_contrast(
        observed_means,
        pd_lib,
        in_lib,
    )

    key = (
        loc,
        rep
    )

    D_obs[key] = z["D"]

    contrast_rows.append({
        "location":
            loc,

        "replicate":
            rep,

        "pulldown_library":
            pd_lib,

        "input_library":
            in_lib,

        **z,
    })


contrasts = pd.DataFrame(
    contrast_rows
)

CONTRAST_OUT = (
    RESULTDIR /
    "PRIMARY_PLUS75_REPLICATE_CONTRASTS_v1.tsv"
)

contrasts.to_csv(
    CONTRAST_OUT,
    sep="\t",
    index=False
)

# ============================================================
# OBSERVED LOCATION SUMMARIES
# ============================================================

location_rows = []

for loc in [
    "ER",
    "mitochondria",
    "cytosol",
]:

    d1 = D_obs[
        (loc,1)
    ]

    d2 = D_obs[
        (loc,2)
    ]

    same_sign = (
        np.sign(d1)
        ==
        np.sign(d2)
    )

    location_rows.append({
        "location":
            loc,

        "D_rep1":
            d1,

        "D_rep2":
            d2,

        "D_location_mean":
            (d1+d2)/2,

        "replicate_concordance":
            (
                "SAME_SIGN"
                if same_sign
                else "OPPOSITE_SIGN"
            ),
    })


location_summary = pd.DataFrame(
    location_rows
)

# ============================================================
# OBSERVED CROSS-LOCATION CONTRASTS
# ============================================================

cross_defs = [
    ("ER_minus_mito",
     "ER","mitochondria"),

    ("ER_minus_cytosol",
     "ER","cytosol"),

    ("mito_minus_cytosol",
     "mitochondria","cytosol"),
]

cross_rows = []

for name,a,b in cross_defs:

    c1 = (
        D_obs[(a,1)]
        -
        D_obs[(b,1)]
    )

    c2 = (
        D_obs[(a,2)]
        -
        D_obs[(b,2)]
    )

    cross_rows.append({
        "contrast":
            name,

        "rep1":
            c1,

        "rep2":
            c2,

        "mean":
            (c1+c2)/2,
    })


cross_summary = pd.DataFrame(
    cross_rows
)

# ============================================================
# UNIT BOOTSTRAP
# ============================================================

print()
print(
    "BOOTSTRAP_START =",
    BOOTSTRAP_N
)

rng = np.random.default_rng(
    BOOTSTRAP_SEED
)

prob = np.full(
    N_UNITS,
    1.0 / N_UNITS,
    dtype=float
)

metric_names = [
    "ER_rep1",
    "ER_rep2",
    "mitochondria_rep1",
    "mitochondria_rep2",
    "cytosol_rep1",
    "cytosol_rep2",

    "ER_location_mean",
    "mitochondria_location_mean",
    "cytosol_location_mean",

    "ER_minus_mito_rep1",
    "ER_minus_mito_rep2",
    "ER_minus_mito_mean",

    "ER_minus_cytosol_rep1",
    "ER_minus_cytosol_rep2",
    "ER_minus_cytosol_mean",

    "mito_minus_cytosol_rep1",
    "mito_minus_cytosol_rep2",
    "mito_minus_cytosol_mean",
]

boot = np.empty(
    (
        BOOTSTRAP_N,
        len(metric_names)
    ),
    dtype=float
)


def vector_D(
    means2d,
    pd_lib,
    in_lib,
):

    R_pd = means2d[
        :,
        col_index[
            f"{pd_lib}_real_density"
        ]
    ]

    R_in = means2d[
        :,
        col_index[
            f"{in_lib}_real_density"
        ]
    ]

    P_pd = means2d[
        :,
        col_index[
            f"{pd_lib}_pseudo_density"
        ]
    ]

    P_in = means2d[
        :,
        col_index[
            f"{in_lib}_pseudo_density"
        ]
    ]

    if (
        (R_pd <= 0).any()
        or
        (R_in <= 0).any()
        or
        (P_pd <= 0).any()
        or
        (P_in <= 0).any()
    ):
        raise SystemExit(
            "BOOTSTRAP_ZERO_AGGREGATE_GATE_FAIL"
        )

    return (
        np.log2(
            R_pd / R_in
        )
        -
        np.log2(
            P_pd / P_in
        )
    )


done = 0

while done < BOOTSTRAP_N:

    b = min(
        BOOTSTRAP_BATCH,
        BOOTSTRAP_N-done
    )

    weights = rng.multinomial(
        N_UNITS,
        prob,
        size=b
    )

    means = (
        weights @ X
    ) / N_UNITS

    er1 = vector_D(
        means,
        "B01","B02"
    )

    er2 = vector_D(
        means,
        "B11","B12"
    )

    mi1 = vector_D(
        means,
        "B09","B10"
    )

    mi2 = vector_D(
        means,
        "B07","B08"
    )

    cy1 = vector_D(
        means,
        "B05","B06"
    )

    cy2 = vector_D(
        means,
        "B03","B04"
    )

    er = (
        er1+er2
    ) / 2

    mi = (
        mi1+mi2
    ) / 2

    cy = (
        cy1+cy2
    ) / 2

    em1 = er1-mi1
    em2 = er2-mi2

    ec1 = er1-cy1
    ec2 = er2-cy2

    mc1 = mi1-cy1
    mc2 = mi2-cy2

    batch_metrics = np.column_stack([
        er1,
        er2,
        mi1,
        mi2,
        cy1,
        cy2,

        er,
        mi,
        cy,

        em1,
        em2,
        (em1+em2)/2,

        ec1,
        ec2,
        (ec1+ec2)/2,

        mc1,
        mc2,
        (mc1+mc2)/2,
    ])

    boot[
        done:done+b,
        :
    ] = batch_metrics

    done += b

    if (
        done % 1000 == 0
        or
        done == BOOTSTRAP_N
    ):
        print(
            "BOOTSTRAP_COMPLETED =",
            done,
            "/",
            BOOTSTRAP_N
        )


if not np.isfinite(
    boot
).all():
    raise SystemExit(
        "BOOTSTRAP_FINITE_GATE_FAIL"
    )


bootstrap_df = pd.DataFrame(
    boot,
    columns=metric_names
)

BOOT_OUT = (
    RESULTDIR /
    "PRIMARY_PLUS75_BOOTSTRAP_10000_v1.tsv"
)

bootstrap_df.to_csv(
    BOOT_OUT,
    sep="\t",
    index=False
)

# ============================================================
# BOOTSTRAP CI SUMMARY
# ============================================================

ci_rows = []

for j,name in enumerate(
    metric_names
):

    vals = boot[:,j]

    ci_rows.append({
        "metric":
            name,

        "bootstrap_mean":
            float(
                np.mean(vals)
            ),

        "ci_2.5":
            float(
                np.percentile(
                    vals,
                    2.5
                )
            ),

        "ci_97.5":
            float(
                np.percentile(
                    vals,
                    97.5
                )
            ),
    })


ci = pd.DataFrame(
    ci_rows
)

CI_OUT = (
    RESULTDIR /
    "PRIMARY_PLUS75_BOOTSTRAP_CI_v1.tsv"
)

ci.to_csv(
    CI_OUT,
    sep="\t",
    index=False
)

ci_lookup = (
    ci.set_index("metric")
)

# ============================================================
# ADD PRIMARY SUPPORT LABELS
# ============================================================

final_location_rows = []

metric_for_location = {
    "ER":
        "ER_location_mean",

    "mitochondria":
        "mitochondria_location_mean",

    "cytosol":
        "cytosol_location_mean",
}

for _,r in location_summary.iterrows():

    loc = r["location"]

    metric = metric_for_location[
        loc
    ]

    lo = float(
        ci_lookup.loc[
            metric,
            "ci_2.5"
        ]
    )

    hi = float(
        ci_lookup.loc[
            metric,
            "ci_97.5"
        ]
    )

    boot_supported = (
        lo > 0
        or
        hi < 0
    )

    replicate_concordant = (
        r["replicate_concordance"]
        ==
        "SAME_SIGN"
    )

    strong = (
        replicate_concordant
        and
        boot_supported
    )

    final_location_rows.append({
        **r.to_dict(),

        "bootstrap_ci_2.5":
            lo,

        "bootstrap_ci_97.5":
            hi,

        "unit_bootstrap_support":
            (
                "SUPPORTED"
                if boot_supported
                else "NOT_SUPPORTED"
            ),

        "strong_support":
            (
                "YES"
                if strong
                else "NO"
            ),
    })


final_location = pd.DataFrame(
    final_location_rows
)

LOCATION_OUT = (
    RESULTDIR /
    "PRIMARY_PLUS75_LOCATION_SUMMARY_v1.tsv"
)

final_location.to_csv(
    LOCATION_OUT,
    sep="\t",
    index=False
)

CROSS_OUT = (
    RESULTDIR /
    "PRIMARY_PLUS75_CROSS_LOCATION_CONTRASTS_v1.tsv"
)

cross_summary.to_csv(
    CROSS_OUT,
    sep="\t",
    index=False
)

# ============================================================
# LIBRARY AGGREGATES
# ============================================================

agg_rows = []

for lib in LIBRARIES:

    R = float(
        observed_means[
            col_index[
                f"{lib}_real_density"
            ]
        ]
    )

    P = float(
        observed_means[
            col_index[
                f"{lib}_pseudo_density"
            ]
        ]
    )

    agg_rows.append({
        "library":
            lib,

        "mean_real_density":
            R,

        "mean_pseudo_density":
            P,
    })


aggregates = pd.DataFrame(
    agg_rows
)

AGG_OUT = (
    RESULTDIR /
    "PRIMARY_PLUS75_LIBRARY_AGGREGATES_v1.tsv"
)

aggregates.to_csv(
    AGG_OUT,
    sep="\t",
    index=False
)

# ============================================================
# REVEAL SUMMARY
# ============================================================

print()
print("============================================")
print("FIRST FORMAL TMD-SPATIAL BIOLOGICAL REVEAL")
print("============================================")

for _,r in contrasts.iterrows():

    print(
        f"{r['location']} rep{int(r['replicate'])}",
        f"D = {r['D']:.12f}",
        f"E_real = {r['E_real']:.12f}",
        f"E_pseudo = {r['E_pseudo']:.12f}",
    )

print()
print("=== LOCATION SUMMARY ===")

for _,r in final_location.iterrows():

    print(
        r["location"],
        f"D_rep1={r['D_rep1']:.12f}",
        f"D_rep2={r['D_rep2']:.12f}",
        f"D_mean={r['D_location_mean']:.12f}",
        f"CI95=[{r['bootstrap_ci_2.5']:.12f},"
        f"{r['bootstrap_ci_97.5']:.12f}]",
        f"CONCORDANCE={r['replicate_concordance']}",
        f"BOOTSTRAP={r['unit_bootstrap_support']}",
        f"STRONG_SUPPORT={r['strong_support']}",
    )

print()
print("PRIMARY_PLUS75_REVEAL_GATE = PASS")
print("BIOLOGICAL_TMD_SIGNAL = VIEWED")

# ============================================================
# OUTPUT SHA MANIFEST
# ============================================================

manifest_paths = [
    CATALOG_OUT,
    UNIT_MAP_OUT,
    UNIT_DENSITY_OUT,
    AGG_OUT,
    CONTRAST_OUT,
    LOCATION_OUT,
    CROSS_OUT,
    BOOT_OUT,
    CI_OUT,
]

with (
    AUDITDIR /
    "PRIMARY_PLUS75_REVEAL_OUTPUT_SHA256_v1.txt"
).open("w") as fh:

    for path in manifest_paths:

        digest = sha256_file(
            path
        )

        fh.write(
            f"{digest}  {path}\n"
        )

print(
    "OUTPUT_SHA_MANIFEST =",
    AUDITDIR /
    "PRIMARY_PLUS75_REVEAL_OUTPUT_SHA256_v1.txt"
)
