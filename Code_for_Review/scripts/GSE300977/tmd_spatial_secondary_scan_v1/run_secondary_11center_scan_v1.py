from pathlib import Path
from collections import defaultdict
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
# Frozen inputs
# ============================================================

PAIR_AUDIT = Path(
    "crossmap_window_preregister/output/"
    "GSE300977_PAIRED_REAL_PSEUDO_CROSSMAP_MASK_AUDIT_v1.tsv"
)

CONTRACT = Path(
    "tmd_spatial_analysis_preregister/"
    "GSE300977_SECONDARY_POSITIONAL_SCAN_CONTRACT_v1.txt"
)

OUTDIR = Path("tmd_spatial_secondary_scan_v1")
COUNTDIR = OUTDIR / "geometry_counts"
RESULTDIR = OUTDIR / "output"
AUDITDIR = OUTDIR / "audit"

for p in [COUNTDIR, RESULTDIR, AUDITDIR]:
    p.mkdir(parents=True, exist_ok=True)

CENTERS = [-30,-15,0,15,30,45,60,75,90,105,120]

EXPECTED_ELIGIBLE_N = {
    -30:6503,
    -15:7526,
    0:8478,
    15:9082,
    30:9081,
    45:9087,
    60:9083,
    75:9098,
    90:7625,
    105:6388,
    120:5570,
}

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
    ("ER",1,"B01","B02"),
    ("ER",2,"B11","B12"),
    ("mitochondria",1,"B09","B10"),
    ("mitochondria",2,"B07","B08"),
    ("cytosol",1,"B05","B06"),
    ("cytosol",2,"B03","B04"),
]

PRIMARY_REFERENCE = {
    ("ER",1): 0.082889583159,
    ("ER",2): 0.095748374758,
    ("mitochondria",1): 0.124859849072,
    ("mitochondria",2): -0.039749358711,
    ("cytosol",1): 0.100608092444,
    ("cytosol",2): 0.067930823332,
}

BOOT_N = 10000
BOOT_SEED = 30097775
BOOT_BATCH = 200

# ============================================================
# Load and freeze center-specific eligible units
# ============================================================

d = pd.read_csv(PAIR_AUDIT, sep="\t", dtype=str)

d["center"] = pd.to_numeric(d["center"]).astype(int)
d["real_crossmap_unmasked_nt"] = pd.to_numeric(
    d["real_crossmap_unmasked_nt"]
)
d["pseudo_min_unmasked_nt"] = pd.to_numeric(
    d["pseudo_min_unmasked_nt"]
)

d = d[
    (d["real_crossmap_unmasked_nt"] >= 1)
    &
    (d["pseudo_min_unmasked_nt"] >= 1)
].copy()

print("=== SECONDARY ELIGIBILITY ===")

for c in CENTERS:
    n = int((d["center"] == c).sum())
    print(f"CENTER={c:+d} ELIGIBLE={n}")

    if n != EXPECTED_ELIGIBLE_N[c]:
        raise SystemExit(
            f"ELIGIBILITY_GATE_FAIL:{c}:{n}"
        )

print("SECONDARY_ELIGIBILITY_GATE = PASS")

ELIG_OUT = (
    RESULTDIR /
    "SECONDARY_11CENTER_ELIGIBLE_UNITS_v1.tsv"
)

d.to_csv(
    ELIG_OUT,
    sep="\t",
    index=False,
)

# ============================================================
# Geometry helpers
# ============================================================

def split_semicolon(x):
    return [z for z in str(x).split(";") if z != ""]


def parse_geometry(s):
    blocks = []

    for token in str(s).split("^"):
        chrom, coords = token.split(":",1)
        a,b = map(int, coords.split("-"))
        blocks.append((chrom,a,b))

    total = sum(b-a for _,a,b in blocks)

    if total != 31:
        raise ValueError(
            f"GEOMETRY_NOT_31:{s}:{total}"
        )

    return blocks


def geometry_positions(geometry, strand, bits):

    blocks = parse_geometry(geometry)

    positions = []

    for chrom,a,b in blocks:

        if strand == "+":
            it = range(a,b)
        elif strand == "-":
            it = range(b-1,a-1,-1)
        else:
            raise ValueError(
                f"BAD_STRAND:{strand}"
            )

        for pos in it:
            positions.append((chrom,pos))

    if len(positions) != 31:
        raise ValueError("POSITION_LENGTH_FAIL")

    bits = str(bits)

    if (
        len(bits) != 31
        or set(bits) - {"0","1"}
    ):
        raise ValueError("MASK_BITS_FAIL")

    keep = [
        p
        for p,bit in zip(positions,bits)
        if bit == "0"
    ]

    if len(keep) < 1:
        raise ValueError(
            "FULLY_MASKED_GEOMETRY_IN_ELIGIBLE_SET"
        )

    return keep


# ============================================================
# Build unique masked geometry catalog
# ============================================================

geometry_keys = set()
unit_specs_by_center = defaultdict(list)

for _,r in d.iterrows():

    c = int(r["center"])
    strand = str(r["strand"])

    real_geom = str(r["real_window_geometry"])
    real_bits = str(r["real_crossmap_mask_bits"])

    real_key = (
        strand,
        real_geom,
        real_bits,
    )

    geometry_keys.add(real_key)

    pgeoms = split_semicolon(
        r["pseudo_window_geometries"]
    )

    pbits = split_semicolon(
        r["pseudo_crossmap_mask_bits"]
    )

    if len(pgeoms) != len(pbits):
        raise SystemExit(
            "PSEUDO_GEOMETRY_MASK_CARDINALITY_FAIL"
        )

    pseudo_keys = []

    for geom,bits in zip(pgeoms,pbits):

        key = (
            strand,
            geom,
            bits,
        )

        geometry_keys.add(key)
        pseudo_keys.append(key)

    unit_specs_by_center[c].append({
        "unit_id":
            r["paired_real_unit_id"],

        "real_key":
            real_key,

        "pseudo_keys":
            pseudo_keys,
    })


geometry_keys = sorted(geometry_keys)

key_to_gid = {
    key:f"GEO{i:06d}"
    for i,key in enumerate(
        geometry_keys,
        start=1
    )
}

geom_positions = {}
wanted = defaultdict(set)

catalog_rows = []

for key in geometry_keys:

    strand,geom,bits = key

    pos = geometry_positions(
        geom,
        strand,
        bits,
    )

    gid = key_to_gid[key]

    geom_positions[gid] = pos

    for chrom,p in pos:
        wanted[(chrom,strand)].add(p)

    catalog_rows.append({
        "geometry_id":gid,
        "strand":strand,
        "geometry":geom,
        "mask_bits":bits,
        "unmasked_nt":len(pos),
    })


catalog = pd.DataFrame(catalog_rows)

catalog.to_csv(
    RESULTDIR /
    "SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv",
    sep="\t",
    index=False,
)

print(
    "UNIQUE_ELIGIBLE_MASKED_GEOMETRIES =",
    len(catalog)
)

# ============================================================
# Merge required genomic positions
# ============================================================

merged_intervals = []

for (chrom,strand),positions in sorted(
    wanted.items()
):

    xs = sorted(positions)

    if not xs:
        continue

    start = xs[0]
    prev = xs[0]

    for pos in xs[1:]:

        if pos == prev + 1:
            prev = pos
            continue

        merged_intervals.append(
            (chrom,start,prev+1,strand)
        )

        start = pos
        prev = pos

    merged_intervals.append(
        (chrom,start,prev+1,strand)
    )


print(
    "MERGED_GENOMIC_INTERVALS =",
    len(merged_intervals)
)

print(
    "UNIQUE_UNMASKED_GENOMIC_POSITIONS =",
    sum(
        len(x)
        for x in wanted.values()
    )
)

# ============================================================
# Count one BAM using merged intervals
# ============================================================

def count_library(lib):

    outfile = (
        COUNTDIR /
        f"{lib}.SECONDARY_GEOMETRY_DENSITIES_v1.tsv"
    )

    if outfile.exists():

        old = pd.read_csv(
            outfile,
            sep="\t"
        )

        if (
            len(old) != len(catalog)
            or old["geometry_id"].tolist()
            != catalog["geometry_id"].tolist()
        ):
            raise SystemExit(
                f"BAD_CHECKPOINT:{lib}"
            )

        print(
            f"{lib} CHECKPOINT = REUSED"
        )

        return dict(
            zip(
                old["geometry_id"],
                old["density"].astype(float)
            )
        )

    print()
    print(f"=== COUNTING {lib} ===")

    ga = BAMGenomeArray(
        str(BAMS[lib])
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

    position_counts = {}

    for i,(chrom,a,b,strand) in enumerate(
        merged_intervals,
        start=1
    ):

        seg = GenomicSegment(
            chrom,a,b,strand
        )

        chain = SegmentChain(seg)

        counts = np.asarray(
            chain.get_counts(
                ga,
                stranded=True
            ),
            dtype=float
        )

        if len(counts) != b-a:
            raise SystemExit(
                f"MERGED_COUNT_LENGTH_FAIL:"
                f"{lib}:{chrom}:{a}-{b}:{strand}"
            )

        if strand == "+":
            positions = range(a,b)
        else:
            positions = range(b-1,a-1,-1)

        for pos,val in zip(
            positions,
            counts
        ):
            position_counts[
                (chrom,strand,pos)
            ] = float(val)

        if (
            i % 2000 == 0
            or i == len(merged_intervals)
        ):
            print(
                lib,
                "MERGED_INTERVALS_COUNTED =",
                i,
                "/",
                len(merged_intervals)
            )

    rows = []
    density_map = {}

    for _,r in catalog.iterrows():

        gid = r["geometry_id"]

        vals = [
            position_counts[
                (chrom,r["strand"],pos)
            ]
            for chrom,pos
            in geom_positions[gid]
        ]

        density = float(
            np.sum(vals)
            /
            len(vals)
        )

        density_map[gid] = density

        rows.append({
            "geometry_id":gid,
            "unmasked_nt":
                len(vals),
            "density":
                density,
        })

    out = pd.DataFrame(rows)

    out.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    print(
        f"{lib} COUNTING_GATE = PASS"
    )

    del ga
    del position_counts
    gc.collect()

    return density_map


geometry_density = {}

for lib in LIBRARIES:
    geometry_density[lib] = count_library(lib)

print()
print(
    "SECONDARY_12_LIBRARY_COUNTING = PASS"
)

# ============================================================
# Center analysis
# ============================================================

def calculate_D(
    means,
    pd_lib,
    in_lib,
    col_idx
):

    Rpd = means[
        col_idx[
            f"{pd_lib}_real"
        ]
    ]

    Rin = means[
        col_idx[
            f"{in_lib}_real"
        ]
    ]

    Ppd = means[
        col_idx[
            f"{pd_lib}_pseudo"
        ]
    ]

    Pin = means[
        col_idx[
            f"{in_lib}_pseudo"
        ]
    ]

    if min(Rpd,Rin,Ppd,Pin) <= 0:
        raise ValueError(
            "ZERO_AGGREGATE"
        )

    Er = math.log2(Rpd/Rin)
    Ep = math.log2(Ppd/Pin)

    return Er,Ep,Er-Ep


replicate_rows = []
location_rows = []
bootstrap_ci_rows = []

# do not print non-primary results until +75 validation passes
stored_results = {}

for center in CENTERS:

    specs = unit_specs_by_center[center]

    N = len(specs)

    matrix_cols = []

    for lib in LIBRARIES:
        matrix_cols.extend([
            f"{lib}_real",
            f"{lib}_pseudo",
        ])

    rows = []

    for spec in specs:

        row = []

        real_gid = key_to_gid[
            spec["real_key"]
        ]

        pseudo_gids = [
            key_to_gid[x]
            for x in spec["pseudo_keys"]
        ]

        for lib in LIBRARIES:

            rd = geometry_density[
                lib
            ][real_gid]

            pdens = float(
                np.mean(
                    [
                        geometry_density[
                            lib
                        ][gid]
                        for gid in pseudo_gids
                    ]
                )
            )

            row.extend([rd,pdens])

        rows.append(row)

    X = np.asarray(
        rows,
        dtype=float
    )

    if X.shape != (N,24):
        raise SystemExit(
            f"CENTER_MATRIX_SHAPE_FAIL:"
            f"{center}:{X.shape}"
        )

    if (
        not np.isfinite(X).all()
        or
        (X < 0).any()
    ):
        raise SystemExit(
            f"CENTER_MATRIX_VALUE_FAIL:{center}"
        )

    col_idx = {
        name:i
        for i,name in enumerate(
            matrix_cols
        )
    }

    means = X.mean(axis=0)

    observed = {}

    for loc,rep,pd_lib,in_lib in CONTRASTS:

        Er,Ep,D = calculate_D(
            means,
            pd_lib,
            in_lib,
            col_idx
        )

        observed[(loc,rep)] = (
            Er,Ep,D
        )

    # bootstrap
    rng = np.random.default_rng(
        BOOT_SEED
    )

    prob = np.full(
        N,
        1.0/N
    )

    boot_location = {
        "ER":[],
        "mitochondria":[],
        "cytosol":[],
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

        Dmat = {}

        for loc,rep,pd_lib,in_lib in CONTRASTS:

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
                    f"BOOT_ZERO_GATE_FAIL:{center}"
                )

            Dmat[(loc,rep)] = (
                np.log2(Rpd/Rin)
                -
                np.log2(Ppd/Pin)
            )

        for loc in [
            "ER",
            "mitochondria",
            "cytosol",
        ]:

            boot_location[loc].append(
                (
                    Dmat[(loc,1)]
                    +
                    Dmat[(loc,2)]
                ) / 2
            )

        done += b

    for loc in boot_location:
        boot_location[loc] = np.concatenate(
            boot_location[loc]
        )

    center_rep_rows = []
    center_loc_rows = []

    for loc,rep,pd_lib,in_lib in CONTRASTS:

        Er,Ep,D = observed[(loc,rep)]

        center_rep_rows.append({
            "center":center,
            "N_units":N,
            "location":loc,
            "replicate":rep,
            "pulldown":pd_lib,
            "input":in_lib,
            "E_real":Er,
            "E_pseudo":Ep,
            "D":D,
        })

    for loc in [
        "ER",
        "mitochondria",
        "cytosol",
    ]:

        d1 = observed[(loc,1)][2]
        d2 = observed[(loc,2)][2]
        dm = (d1+d2)/2

        vals = boot_location[loc]

        lo = float(
            np.percentile(vals,2.5)
        )

        hi = float(
            np.percentile(vals,97.5)
        )

        center_loc_rows.append({
            "center":center,
            "N_units":N,
            "location":loc,
            "D_rep1":d1,
            "D_rep2":d2,
            "D_mean":dm,
            "replicate_concordance":
                (
                    "SAME_SIGN"
                    if np.sign(d1)==np.sign(d2)
                    else "OPPOSITE_SIGN"
                ),
            "pointwise_bootstrap_ci_2.5":
                lo,
            "pointwise_bootstrap_ci_97.5":
                hi,
            "pointwise_ci_excludes_zero":
                (
                    "YES"
                    if lo > 0 or hi < 0
                    else "NO"
                ),
        })

    stored_results[center] = (
        center_rep_rows,
        center_loc_rows
    )

# ============================================================
# Mandatory +75 reproduction gate
# ============================================================

for row in stored_results[75][0]:

    key = (
        row["location"],
        int(row["replicate"])
    )

    diff = abs(
        row["D"]
        -
        PRIMARY_REFERENCE[key]
    )

    print(
        "PRIMARY_REPRO_CHECK",
        key[0],
        "rep"+str(key[1]),
        "D=",
        f"{row['D']:.12f}",
        "ABS_DIFF=",
        f"{diff:.3e}"
    )

    if diff > 1e-10:
        raise SystemExit(
            "PRIMARY_PLUS75_REPRODUCTION_GATE = FAIL"
        )

print(
    "PRIMARY_PLUS75_REPRODUCTION_GATE = PASS"
)

# Only now make non-primary signal visible
for c in CENTERS:

    replicate_rows.extend(
        stored_results[c][0]
    )

    location_rows.extend(
        stored_results[c][1]
    )


repdf = pd.DataFrame(
    replicate_rows
)

locdf = pd.DataFrame(
    location_rows
)

repdf.to_csv(
    RESULTDIR /
    "SECONDARY_11CENTER_REPLICATE_D_v1.tsv",
    sep="\t",
    index=False
)

locdf.to_csv(
    RESULTDIR /
    "SECONDARY_11CENTER_LOCATION_PROFILE_v1.tsv",
    sep="\t",
    index=False
)

# ============================================================
# Formal reveal
# ============================================================

print()
print("============================================")
print("FORMAL SECONDARY 11-CENTER POSITIONAL REVEAL")
print("============================================")

for loc in [
    "ER",
    "mitochondria",
    "cytosol",
]:

    print()
    print("LOCATION =", loc)

    z = locdf[
        locdf["location"] == loc
    ]

    for _,r in z.iterrows():

        print(
            f"CENTER={int(r['center']):+d}",
            f"N={int(r['N_units'])}",
            f"D1={r['D_rep1']:.12f}",
            f"D2={r['D_rep2']:.12f}",
            f"Dmean={r['D_mean']:.12f}",
            f"CI95=[{r['pointwise_bootstrap_ci_2.5']:.12f},"
            f"{r['pointwise_bootstrap_ci_97.5']:.12f}]",
            f"CONCORDANCE={r['replicate_concordance']}",
            f"POINTWISE_CI_EXCLUDES_0={r['pointwise_ci_excludes_zero']}",
        )

print()
print("SECONDARY_11CENTER_REVEAL_GATE = PASS")
print("SECONDARY_POSITIONAL_SIGNAL = VIEWED")

# ============================================================
# SHA
# ============================================================

for path in [
    ELIG_OUT,
    RESULTDIR /
    "SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv",
    RESULTDIR /
    "SECONDARY_11CENTER_REPLICATE_D_v1.tsv",
    RESULTDIR /
    "SECONDARY_11CENTER_LOCATION_PROFILE_v1.tsv",
]:

    digest = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()

    print(
        "SHA256",
        digest,
        path
    )
