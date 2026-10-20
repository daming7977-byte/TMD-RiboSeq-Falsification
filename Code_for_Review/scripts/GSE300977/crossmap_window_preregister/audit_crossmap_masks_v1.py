from pathlib import Path
from collections import defaultdict
import hashlib
import pandas as pd

from plastid import BigBedReader, GenomicSegment

BB = Path(
    "periodicity_blinded/crossmap/"
    "GRCh38_release104_29_1_crossmap.sorted.bb"
)

ALL_REAL = Path(
    "tmd_anchor_preregister/bridge/"
    "GSE300977_TMD_GEOMETRY_RESOLVED_UNITS_ALL_CENTERS_v2.tsv"
)

PAIRED = Path(
    "pseudo_transfer_preregister/bridge/"
    "GSE300977_REAL_PSEUDO_PAIRED_GEOMETRY_UNITS_ALL_CENTERS_v1.tsv"
)

OUTDIR = Path(
    "crossmap_window_preregister/output"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

EXPECTED_BB_SHA = (
    "b6a0f58dcbc0c2503d5549234db33f89683e5a1aeb8b62ef8646257e83fe3fd8"
)

CENTERS = [-30,-15,0,15,30,45,60,75,90,105,120]


def sha256_file(path):
    h = hashlib.sha256()

    with path.open("rb") as fh:
        for block in iter(
            lambda: fh.read(1024*1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


actual_sha = sha256_file(BB)

print("BIGBED_SHA256 =", actual_sha)

if actual_sha != EXPECTED_BB_SHA:
    raise SystemExit(
        "BIGBED_SHA_GATE = FAIL"
    )

print("BIGBED_SHA_GATE = PASS")

bb = BigBedReader(str(BB))

# --------------------------------------------------
# Geometry helpers
# --------------------------------------------------

def parse_geometry(s):

    blocks = []

    for token in str(s).split("^"):

        chrom, coords = token.split(":",1)

        start,end = map(
            int,
            coords.split("-")
        )

        blocks.append(
            (chrom,start,end)
        )

    total = sum(
        end-start
        for chrom,start,end in blocks
    )

    if total != 31:
        raise ValueError(
            f"WINDOW_LENGTH_NOT_31:{s}:{total}"
        )

    return blocks


def transcript_positions(
    blocks,
    strand,
):

    positions = []

    # Geometry strings were frozen in transcript order.
    for chrom,start,end in blocks:

        if strand == "+":
            iterator = range(start,end)
        elif strand == "-":
            iterator = range(
                end-1,
                start-1,
                -1
            )
        else:
            raise ValueError(
                f"BAD_STRAND:{strand}"
            )

        for pos in iterator:
            positions.append(
                (chrom,pos)
            )

    if len(positions) != 31:
        raise ValueError(
            f"POSITION_LENGTH_FAIL:{len(positions)}"
        )

    return positions


cache = {}


def audit_geometry(
    geometry,
    strand,
):

    key = (
        str(geometry),
        str(strand)
    )

    if key in cache:
        return cache[key]

    blocks = parse_geometry(
        geometry
    )

    tx_positions = transcript_positions(
        blocks,
        strand
    )

    masked_intervals = defaultdict(list)

    for chrom,start,end in blocks:

        query = GenomicSegment(
            chrom,
            start,
            end,
            strand
        )

        hits = bb.get(
            query,
            stranded=True
        )

        for feature in hits:

            for seg in feature.segments:

                a = max(
                    start,
                    int(seg.start)
                )

                b = min(
                    end,
                    int(seg.end)
                )

                if a < b:
                    masked_intervals[chrom].append(
                        (a,b)
                    )

    bits = []

    for chrom,pos in tx_positions:

        is_masked = any(
            a <= pos < b
            for a,b
            in masked_intervals.get(
                chrom,
                []
            )
        )

        bits.append(
            "1" if is_masked else "0"
        )

    bitstr = "".join(bits)

    if len(bitstr) != 31:
        raise ValueError(
            "MASK_VECTOR_LENGTH_FAIL"
        )

    masked_nt = bitstr.count("1")
    unmasked_nt = 31-masked_nt

    result = {
        "crossmap_mask_bits":
            bitstr,

        "crossmap_masked_nt":
            masked_nt,

        "crossmap_unmasked_nt":
            unmasked_nt,

        "crossmap_fully_masked":
            "YES"
            if masked_nt == 31
            else "NO",

        "crossmap_any_masked":
            "YES"
            if masked_nt > 0
            else "NO",
    }

    cache[key] = result

    return result


# --------------------------------------------------
# ALL-TMD real windows
# --------------------------------------------------

a = pd.read_csv(
    ALL_REAL,
    sep="\t",
    dtype=str
)

print()
print(
    "ALL_TMD_INPUT_ROWS =",
    len(a)
)

all_rows = []

for _,r in a.iterrows():

    z = audit_geometry(
        r[
            "window_genomic_blocks_transcript_order"
        ],
        r["strand"]
    )

    row = r.to_dict()
    row.update(z)

    all_rows.append(row)

ao = pd.DataFrame(all_rows)

ao.to_csv(
    OUTDIR /
    "GSE300977_ALL_TMD_CROSSMAP_MASK_AUDIT_v1.tsv",
    sep="\t",
    index=False
)

# --------------------------------------------------
# Paired real-pseudo windows
# --------------------------------------------------

p = pd.read_csv(
    PAIRED,
    sep="\t",
    dtype=str
)

print(
    "PAIRED_INPUT_ROWS =",
    len(p)
)

paired_rows = []

for _,r in p.iterrows():

    real = audit_geometry(
        r["real_window_geometry"],
        r["strand"]
    )

    pseudo_geoms = [
        x
        for x in str(
            r["pseudo_window_geometries"]
        ).split(";")
        if x != ""
    ]

    if len(pseudo_geoms) == 0:
        raise SystemExit(
            "NO_PSEUDO_GEOMETRY"
        )

    pseudo_audits = [
        audit_geometry(
            g,
            r["strand"]
        )
        for g in pseudo_geoms
    ]

    row = r.to_dict()

    row.update({
        "real_crossmap_mask_bits":
            real[
                "crossmap_mask_bits"
            ],

        "real_crossmap_masked_nt":
            real[
                "crossmap_masked_nt"
            ],

        "real_crossmap_unmasked_nt":
            real[
                "crossmap_unmasked_nt"
            ],

        "real_crossmap_any_masked":
            real[
                "crossmap_any_masked"
            ],

        "real_crossmap_fully_masked":
            real[
                "crossmap_fully_masked"
            ],

        "pseudo_crossmap_mask_bits":
            ";".join(
                z["crossmap_mask_bits"]
                for z in pseudo_audits
            ),

        "pseudo_crossmap_masked_nt":
            ";".join(
                str(
                    z["crossmap_masked_nt"]
                )
                for z in pseudo_audits
            ),

        "pseudo_crossmap_unmasked_nt":
            ";".join(
                str(
                    z["crossmap_unmasked_nt"]
                )
                for z in pseudo_audits
            ),

        "pseudo_min_unmasked_nt":
            min(
                z["crossmap_unmasked_nt"]
                for z in pseudo_audits
            ),

        "pseudo_max_unmasked_nt":
            max(
                z["crossmap_unmasked_nt"]
                for z in pseudo_audits
            ),

        "pseudo_any_fully_masked":
            (
                "YES"
                if any(
                    z["crossmap_fully_masked"]
                    == "YES"
                    for z in pseudo_audits
                )
                else "NO"
            ),

        "pseudo_any_masked":
            (
                "YES"
                if any(
                    z["crossmap_any_masked"]
                    == "YES"
                    for z in pseudo_audits
                )
                else "NO"
            ),
    })

    paired_rows.append(row)

po = pd.DataFrame(
    paired_rows
)

po.to_csv(
    OUTDIR /
    "GSE300977_PAIRED_REAL_PSEUDO_CROSSMAP_MASK_AUDIT_v1.tsv",
    sep="\t",
    index=False
)

# --------------------------------------------------
# Counts-blind summaries
# --------------------------------------------------

print()
print("=== ALL-TMD CROSSMAP SUMMARY BY CENTER ===")

for c in CENTERS:

    z = ao[
        pd.to_numeric(
            ao["center"]
        ).eq(c)
    ]

    masked = pd.to_numeric(
        z["crossmap_masked_nt"]
    )

    print(
        f"CENTER={c:+d}",
        f"N={len(z)}",
        f"ANY_MASK={int((masked>0).sum())}",
        f"FULLY_MASKED={int((masked==31).sum())}",
        f"MIN_UNMASKED={int((31-masked).min())}",
        f"MEDIAN_UNMASKED={float((31-masked).median()):.1f}",
    )

print()
print("=== PAIRED CROSSMAP SUMMARY BY CENTER ===")

for c in CENTERS:

    z = po[
        pd.to_numeric(
            po["center"]
        ).eq(c)
    ]

    real_unmasked = pd.to_numeric(
        z["real_crossmap_unmasked_nt"]
    )

    pseudo_min = pd.to_numeric(
        z["pseudo_min_unmasked_nt"]
    )

    print(
        f"CENTER={c:+d}",
        f"N_REAL={len(z)}",
        f"REAL_ANY_MASK={int((real_unmasked<31).sum())}",
        f"REAL_FULLY_MASKED={int((real_unmasked==0).sum())}",
        f"PSEUDO_ANY_FULLY_MASKED={int(z['pseudo_any_fully_masked'].eq('YES').sum())}",
        f"JOINT_ALL_HAVE_AT_LEAST_1_UNMASKED={int(((real_unmasked>0)&(pseudo_min>0)).sum())}",
        f"REAL_MIN_UNMASKED={int(real_unmasked.min())}",
        f"PSEUDO_MIN_UNMASKED={int(pseudo_min.min())}",
    )

# --------------------------------------------------
# Primary +75 detailed summary
# --------------------------------------------------

p75 = po[
    pd.to_numeric(
        po["center"]
    ).eq(75)
].copy()

ru = pd.to_numeric(
    p75["real_crossmap_unmasked_nt"]
)

pu = pd.to_numeric(
    p75["pseudo_min_unmasked_nt"]
)

print()
print("=== PRIMARY +75 CROSSMAP SUMMARY ===")

print(
    "PRIMARY_PLUS75_REAL_UNITS =",
    len(p75)
)

print(
    "PRIMARY_REAL_ANY_MASKED =",
    int((ru < 31).sum())
)

print(
    "PRIMARY_REAL_FULLY_UNMASKED =",
    int((ru == 31).sum())
)

print(
    "PRIMARY_REAL_FULLY_MASKED =",
    int((ru == 0).sum())
)

print(
    "PRIMARY_PSEUDO_ANY_FULLY_MASKED =",
    int(
        p75[
            "pseudo_any_fully_masked"
        ].eq("YES").sum()
    )
)

print(
    "PRIMARY_JOINT_ALL_HAVE_AT_LEAST_1_UNMASKED =",
    int(
        (
            (ru > 0)
            &
            (pu > 0)
        ).sum()
    )
)

print(
    "PRIMARY_REAL_MIN_UNMASKED_NT =",
    int(ru.min())
)

print(
    "PRIMARY_REAL_MEDIAN_UNMASKED_NT =",
    float(ru.median())
)

print(
    "PRIMARY_PSEUDO_MIN_UNMASKED_NT =",
    int(pu.min())
)

print(
    "PRIMARY_PSEUDO_MEDIAN_MIN_UNMASKED_NT =",
    float(pu.median())
)

print()
print(
    "UNIQUE_GEOMETRIES_QUERIED =",
    len(cache)
)

print(
    "CROSSMAP_WINDOW_MASK_AUDIT = PASS"
)

print(
    "BIOLOGICAL_TMD_SIGNAL = NOT_VIEWED"
)
