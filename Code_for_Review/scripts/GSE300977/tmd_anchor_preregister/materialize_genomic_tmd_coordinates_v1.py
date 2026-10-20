from pathlib import Path
from collections import defaultdict
import re
import pandas as pd

GTF = Path(
    "mapping_reference_blinded/"
    "Homo_sapiens.GRCh38.104.gtf"
)

INP = Path(
    "tmd_anchor_preregister/bridge/"
    "FORMAL_TMD_ANCHOR_CANDIDATES_v1p1.tsv"
)

OUT = Path(
    "tmd_anchor_preregister/bridge/"
    "GSE300977_FORMAL_GENOMIC_TMD_ANCHORS_v1.tsv"
)

# --------------------------------------------------
# Parse CDS segments from Ensembl104
# --------------------------------------------------

attr_re = re.compile(r'(\S+) "([^"]*)"')

cds = defaultdict(list)

with GTF.open() as fh:
    for line in fh:

        if line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")

        if len(f) != 9 or f[2] != "CDS":
            continue

        chrom = f[0]
        start0 = int(f[3]) - 1
        end0 = int(f[4])
        strand = f[6]
        phase = f[7]

        attrs = dict(attr_re.findall(f[8]))

        tx = attrs.get("transcript_id")

        if tx is None:
            continue

        cds[tx].append({
            "chrom": chrom,
            "start0": start0,
            "end0": end0,
            "strand": strand,
            "phase": phase,
        })

# --------------------------------------------------
# Helpers
# --------------------------------------------------

def ordered_segments(tx):

    segs = cds.get(tx, [])

    if not segs:
        raise ValueError(f"NO_CDS:{tx}")

    chroms = {x["chrom"] for x in segs}
    strands = {x["strand"] for x in segs}

    if len(chroms) != 1:
        raise ValueError(f"MULTI_CHROM:{tx}:{chroms}")

    if len(strands) != 1:
        raise ValueError(f"MULTI_STRAND:{tx}:{strands}")

    strand = next(iter(strands))

    if strand == "+":
        segs = sorted(
            segs,
            key=lambda x: (x["start0"], x["end0"])
        )
    elif strand == "-":
        segs = sorted(
            segs,
            key=lambda x: (x["start0"], x["end0"]),
            reverse=True
        )
    else:
        raise ValueError(f"BAD_STRAND:{tx}:{strand}")

    return segs


def cds_offset_to_genome(segs, offset0):

    cumulative = 0

    for seg in segs:

        L = seg["end0"] - seg["start0"]

        if cumulative <= offset0 < cumulative + L:

            local = offset0 - cumulative

            if seg["strand"] == "+":
                return seg["start0"] + local
            else:
                return seg["end0"] - 1 - local

        cumulative += L

    raise ValueError(
        f"OFFSET_OUTSIDE_CDS:{offset0}/{cumulative}"
    )


def cds_interval_to_blocks(segs, start0, end0):

    if not (0 <= start0 < end0):
        raise ValueError(
            f"BAD_CDS_INTERVAL:{start0}-{end0}"
        )

    blocks = []
    cumulative = 0
    mapped_nt = 0

    for seg in segs:

        L = seg["end0"] - seg["start0"]

        seg_cds_start = cumulative
        seg_cds_end = cumulative + L

        ov0 = max(start0, seg_cds_start)
        ov1 = min(end0, seg_cds_end)

        if ov0 < ov1:

            local0 = ov0 - seg_cds_start
            local1 = ov1 - seg_cds_start

            if seg["strand"] == "+":
                g0 = seg["start0"] + local0
                g1 = seg["start0"] + local1
            else:
                g0 = seg["end0"] - local1
                g1 = seg["end0"] - local0

            blocks.append((g0, g1))
            mapped_nt += (g1 - g0)

        cumulative += L

    if mapped_nt != end0 - start0:
        raise ValueError(
            f"INTERVAL_INCOMPLETE:"
            f"{start0}-{end0}:"
            f"expected={end0-start0}:mapped={mapped_nt}"
        )

    return blocks


# --------------------------------------------------
# Input
# --------------------------------------------------

d = pd.read_csv(INP, sep="\t")

required = {
    "target_transcript",
    "gene_id",
    "gene_name",
    "TMD_index",
    "TMD_position_class",
    "aa_start0",
    "aa_end0_exclusive",
    "target_tmd_cds_nt_start0",
    "target_tmd_cds_nt_end0_exclusive",
    "target_cds_length",
    "target_protein_length",
}

missing = required - set(d.columns)

if missing:
    raise SystemExit(
        f"REQUIRED_COLUMNS_FAIL:{sorted(missing)}"
    )

if len(d) != 19113:
    raise SystemExit(
        f"INPUT_ROW_COUNT_FAIL:{len(d)}"
    )

# --------------------------------------------------
# Materialize
# --------------------------------------------------

rows = []

for _, r in d.iterrows():

    tx = str(r["target_transcript"])

    segs = ordered_segments(tx)

    chrom = segs[0]["chrom"]
    strand = segs[0]["strand"]

    cds_len_gtf = sum(
        x["end0"] - x["start0"]
        for x in segs
    )

    expected_cds_len = int(r["target_cds_length"])

    if cds_len_gtf != expected_cds_len:
        raise SystemExit(
            f"CDS_LENGTH_GATE_FAIL:"
            f"{tx}:gtf={cds_len_gtf}:"
            f"expected={expected_cds_len}"
        )

    tmd_start = int(r["target_tmd_cds_nt_start0"])
    tmd_end = int(r["target_tmd_cds_nt_end0_exclusive"])

    anchor_g0 = cds_offset_to_genome(
        segs,
        tmd_start
    )

    blocks = cds_interval_to_blocks(
        segs,
        tmd_start,
        tmd_end
    )

    block_nt = sum(b-a for a,b in blocks)

    if block_nt != tmd_end - tmd_start:
        raise SystemExit(
            f"TMD_BLOCK_LENGTH_FAIL:{tx}"
        )

    # blocks are stored in transcript 5'->3' order
    block_string = "^".join(
        f"{chrom}:{a}-{b}"
        for a,b in blocks
    )

    span_start0 = min(a for a,b in blocks)
    span_end0 = max(b for a,b in blocks)

    rows.append({
        "target_transcript": tx,
        "gene_id": r["gene_id"],
        "gene_name": r["gene_name"],
        "chrom": chrom,
        "strand": strand,

        "TMD_index": int(r["TMD_index"]),
        "TMD_position_class": r["TMD_position_class"],

        "aa_start0": int(r["aa_start0"]),
        "aa_end0_exclusive": int(r["aa_end0_exclusive"]),

        "target_protein_length":
            int(r["target_protein_length"]),

        "target_cds_length":
            expected_cds_len,

        "tmd_cds_nt_start0": tmd_start,
        "tmd_cds_nt_end0_exclusive": tmd_end,

        "anchor_genomic_start0": anchor_g0,
        "anchor_genomic_end0": anchor_g0 + 1,
        "anchor_genomic_1based": anchor_g0 + 1,

        "tmd_genomic_block_count": len(blocks),
        "tmd_genomic_blocks_transcript_order":
            block_string,

        "tmd_genomic_span_start0": span_start0,
        "tmd_genomic_span_end0": span_end0,

        "mapped_tmd_nt": block_nt,
    })

out = pd.DataFrame(rows)

# --------------------------------------------------
# Formal gates
# --------------------------------------------------

if len(out) != 19113:
    raise SystemExit(
        f"OUTPUT_ROW_COUNT_FAIL:{len(out)}"
    )

if out["target_transcript"].nunique() != 4502:
    raise SystemExit(
        "OUTPUT_TRANSCRIPT_COUNT_FAIL:"
        f"{out['target_transcript'].nunique()}"
    )

first_n = int(
    (out["TMD_position_class"] == "FIRST").sum()
)

internal_n = int(
    (out["TMD_position_class"] == "INTERNAL").sum()
)

if first_n != 4502:
    raise SystemExit(
        f"FIRST_COUNT_FAIL:{first_n}"
    )

if internal_n != 14611:
    raise SystemExit(
        f"INTERNAL_COUNT_FAIL:{internal_n}"
    )

# each transcript/TMD_index must be unique
dups = int(
    out.duplicated(
        subset=["target_transcript","TMD_index"]
    ).sum()
)

if dups != 0:
    raise SystemExit(
        f"TRANSCRIPT_TMD_INDEX_DUPLICATES:{dups}"
    )

# anchor must lie inside at least one TMD genomic block
bad_anchor = 0

for _, r in out.iterrows():

    anchor = int(r["anchor_genomic_start0"])

    good = False

    for token in str(
        r["tmd_genomic_blocks_transcript_order"]
    ).split("^"):

        coords = token.split(":",1)[1]
        a,b = map(int, coords.split("-"))

        if a <= anchor < b:
            good = True
            break

    if not good:
        bad_anchor += 1

if bad_anchor != 0:
    raise SystemExit(
        f"ANCHOR_INSIDE_TMD_GATE_FAIL:{bad_anchor}"
    )

out.to_csv(
    OUT,
    sep="\t",
    index=False
)

print("INPUT_TMD_ROWS =", len(d))
print("OUTPUT_TMD_ROWS =", len(out))
print("OUTPUT_TMD_TRANSCRIPTS =", out["target_transcript"].nunique())

print("FIRST_TMD_ROWS =", first_n)
print("INTERNAL_TMD_ROWS =", internal_n)

print(
    "SINGLE_BLOCK_TMD_ROWS =",
    int((out["tmd_genomic_block_count"] == 1).sum())
)

print(
    "MULTIBLOCK_TMD_ROWS =",
    int((out["tmd_genomic_block_count"] > 1).sum())
)

print(
    "MAX_TMD_BLOCK_COUNT =",
    int(out["tmd_genomic_block_count"].max())
)

print(
    "TRANSCRIPT_TMD_INDEX_DUPLICATES =",
    dups
)

print(
    "ANCHOR_INSIDE_TMD_FAILURES =",
    bad_anchor
)

# annotation-only duplicate diagnostics
dup_anchor_rows = int(
    out.duplicated(
        subset=[
            "chrom",
            "strand",
            "anchor_genomic_start0",
        ],
        keep=False
    ).sum()
)

dup_anchor_groups = (
    out.groupby(
        ["chrom","strand","anchor_genomic_start0"]
    )
    .size()
)

dup_anchor_groups = int(
    (dup_anchor_groups > 1).sum()
)

multi_tmd_genes = int(
    (
        out.groupby("gene_id").size() > 1
    ).sum()
)

print(
    "ROWS_IN_DUPLICATE_GENOMIC_ANCHOR_GROUPS =",
    dup_anchor_rows
)

print(
    "DUPLICATE_GENOMIC_ANCHOR_GROUPS =",
    dup_anchor_groups
)

print(
    "GENES_WITH_MULTIPLE_TMD_ROWS =",
    multi_tmd_genes
)

print()
print("GENOMIC_TMD_COORDINATE_GATE = PASS")
print("BIOLOGICAL_TMD_SIGNAL = NOT_VIEWED")
