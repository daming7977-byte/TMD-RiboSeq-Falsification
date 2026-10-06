from pathlib import Path
from collections import defaultdict
import hashlib
import re
import pandas as pd

PSEUDO = Path(
    "REPRO_WORKSPACE/GSE133111_stage0C/reference/tmd/"
    "primary_tmd_pseudo_anchor_pairs.tsv"
)

FORMAL = Path(
    "tmd_anchor_preregister/bridge/"
    "FORMAL_TMD_ANCHOR_CANDIDATES_v1p1.tsv"
)

GENOMIC = Path(
    "tmd_anchor_preregister/bridge/"
    "GSE300977_FORMAL_GENOMIC_TMD_ANCHORS_v1.tsv"
)

GTF = Path(
    "mapping_reference_blinded/"
    "Homo_sapiens.GRCh38.104.gtf"
)

OUTDIR = Path("pseudo_transfer_preregister/bridge")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT = (
    OUTDIR /
    "GSE300977_TRANSFERRED_PSEUDO_PAIRS_v2.tsv"
)

EXCLUDED_OUT = (
    OUTDIR /
    "GSE300977_PSEUDO_MATCHED_EXCLUDED_NO_FORMAL_REAL_v2.tsv"
)

CENTERS = [-30,-15,0,15,30,45,60,75,90,105,120]

EXPECTED_PSEUDO_SHA = (
    "1f6150d9985c82460d9f12ccc565edb416a1539ef70c761149f35b8fefba60a5"
)

EXPECTED_GTF_SHA = (
    "f3de2bccc02877709272a0ce66af59a2922d407a4076541bc545311ff77d2fa"
)

# NOTE: replace the above only if your already-frozen Ensembl104 GTF SHA
# printed previously differs. This must equal the frozen project value.

def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()

if sha256_file(PSEUDO) != EXPECTED_PSEUDO_SHA:
    raise SystemExit("PSEUDO_SHA_GATE_FAIL")

actual_gtf_sha = sha256_file(GTF)

# authoritative frozen project GTF SHA
FROZEN_GTF_SHA = (
    "f3de2bccc02877709272a0ce66af59a2922d407a4076541bc5a65c9d3518581e"
)

if actual_gtf_sha != FROZEN_GTF_SHA:
    raise SystemExit(
        f"GTF_SHA_GATE_FAIL:{actual_gtf_sha}"
    )

print("SOURCE_SHA_GATES = PASS")

# --------------------------------------------------
# Load
# --------------------------------------------------

p = pd.read_csv(PSEUDO, sep="\t", dtype=str)
f = pd.read_csv(FORMAL, sep="\t")
gg = pd.read_csv(GENOMIC, sep="\t")

p["TMD_index"] = p["TMD_index"].astype(int)
f["TMD_index"] = f["TMD_index"].astype(int)
gg["TMD_index"] = gg["TMD_index"].astype(int)

matched = p[
    p["match_status"].eq("MATCHED")
].copy()

if len(matched) != 9829:
    raise SystemExit(
        f"MATCHED_SOURCE_GATE_FAIL:{len(matched)}"
    )

formal_cols = [
    "transcript",
    "TMD_index",
    "match_class",
    "target_transcript",
    "gene_id",
    "gene_name",
    "TMD_position_class",
    "target_cds_length",
    "target_protein_length",
    "target_tmd_cds_nt_start0",
    "target_tmd_cds_nt_end0_exclusive",
]

fm = f[formal_cols].copy().rename(
    columns={
        "TMD_position_class":
            "target_TMD_position_class"
    }
)

x = matched.merge(
    fm,
    on=["transcript","TMD_index"],
    how="left",
    validate="one_to_one"
)

x["formal_real_present"] = (
    x["target_transcript"].notna()
)

excluded = x[
    ~x["formal_real_present"]
].copy()

x = x[
    x["formal_real_present"]
].copy()

print("SOURCE_MATCHED_PAIRS =", len(matched))
print(
    "MATCHED_WITH_FORMAL_REAL_TMD =",
    len(x)
)
print(
    "MATCHED_EXCLUDED_NO_FORMAL_REAL =",
    len(excluded)
)

# position class
pos_mismatch = int(
    (
        x["TMD_position_class"]
        != x["target_TMD_position_class"]
    ).sum()
)

print(
    "SOURCE_TARGET_POSITION_CLASS_MISMATCHES =",
    pos_mismatch
)

if pos_mismatch != 0:
    raise SystemExit(
        "SOURCE_TARGET_POSITION_CLASS_GATE = FAIL"
    )

print(
    "SOURCE_TARGET_POSITION_CLASS_GATE = PASS"
)

# --------------------------------------------------
# Attach frozen real genomic coordinates
# --------------------------------------------------

gx = gg[
    [
        "target_transcript",
        "TMD_index",
        "chrom",
        "strand",
        "anchor_genomic_start0",
    ]
].copy()

x = x.merge(
    gx,
    on=["target_transcript","TMD_index"],
    how="left",
    validate="one_to_one"
)

if x[
    ["chrom","strand","anchor_genomic_start0"]
].isna().any().any():
    raise SystemExit(
        "REAL_GENOMIC_JOIN_GATE_FAIL"
    )

# numeric source fields
for c in [
    "CDS_length_nt",
    "real_TMD_length_nt",
    "pseudo_anchor_cds0",
    "pseudo_segment_cds_start0",
    "pseudo_segment_cds_end0_exclusive",
    "pseudo_segment_length_nt",
]:
    x[c] = pd.to_numeric(x[c])

x["target_cds_length"] = pd.to_numeric(
    x["target_cds_length"]
)

x["source_target_span_delta"] = (
    x["CDS_length_nt"]
    - x["target_cds_length"]
)

observed_delta = set(
    x["source_target_span_delta"]
    .astype(int)
    .unique()
)

if not observed_delta.issubset({0,3}):
    raise SystemExit(
        f"SPAN_DELTA_GATE_FAIL:{sorted(observed_delta)}"
    )

print()
print("=== SOURCE-TARGET SPAN DELTA ===")
print(
    x["source_target_span_delta"]
    .value_counts()
    .sort_index()
    .to_string()
)

# --------------------------------------------------
# Parse Ensembl104 CDS + stop_codon
# --------------------------------------------------

attr_re = re.compile(r'(\S+) "([^"]*)"')

features = {
    "CDS": defaultdict(list),
    "stop_codon": defaultdict(list),
}

with GTF.open() as fh:
    for line in fh:

        if line.startswith("#"):
            continue

        z = line.rstrip("\n").split("\t")

        if len(z) != 9:
            continue

        ft = z[2]

        if ft not in features:
            continue

        attrs = dict(attr_re.findall(z[8]))
        tx = attrs.get("transcript_id")

        if tx is None:
            continue

        features[ft][tx].append({
            "chrom": z[0],
            "start0": int(z[3])-1,
            "end0": int(z[4]),
            "strand": z[6],
            "feature": ft,
        })


def ordered(segs):

    if not segs:
        return []

    strands = {z["strand"] for z in segs}
    chroms = {z["chrom"] for z in segs}

    if len(strands) != 1 or len(chroms) != 1:
        raise ValueError("BAD_TRANSCRIPT_GEOMETRY")

    strand = next(iter(strands))

    return sorted(
        segs,
        key=lambda z: (z["start0"],z["end0"]),
        reverse=(strand=="-")
    )


def make_path(tx, delta):

    cds = ordered(
        features["CDS"].get(tx, [])
    )

    if not cds:
        raise ValueError(f"NO_CDS:{tx}")

    translated_len = sum(
        z["end0"]-z["start0"]
        for z in cds
    )

    if delta == 0:
        return cds, translated_len, 0

    stop = ordered(
        features["stop_codon"].get(tx, [])
    )

    stop_len = sum(
        z["end0"]-z["start0"]
        for z in stop
    )

    if stop_len != 3:
        raise ValueError(
            f"STOP_LENGTH_FAIL:{tx}:{stop_len}"
        )

    allseg = ordered(cds + stop)
    tags = [z["feature"] for z in allseg]

    first_stop = tags.index("stop_codon")

    if any(t != "CDS" for t in tags[:first_stop]):
        raise ValueError(f"CDS_STOP_ORDER_FAIL:{tx}")

    if any(t != "stop_codon" for t in tags[first_stop:]):
        raise ValueError(f"CDS_STOP_ORDER_FAIL:{tx}")

    return allseg, translated_len, stop_len


def offset_to_genome(path, offset0):

    cum = 0

    for seg in path:

        L = seg["end0"]-seg["start0"]

        if cum <= offset0 < cum+L:

            local = offset0-cum

            if seg["strand"] == "+":
                return seg["start0"] + local

            return seg["end0"] - 1 - local

        cum += L

    raise ValueError(
        f"OFFSET_OUTSIDE_SPAN:{offset0}/{cum}"
    )


def interval_to_blocks(path,a,b):

    if a < 0 or a >= b:
        return None

    total = sum(
        z["end0"]-z["start0"]
        for z in path
    )

    if b > total:
        return None

    out = []
    cum = 0
    mapped = 0

    for seg in path:

        L = seg["end0"]-seg["start0"]

        x0 = max(a,cum)
        x1 = min(b,cum+L)

        if x0 < x1:

            l0 = x0-cum
            l1 = x1-cum

            if seg["strand"] == "+":
                g0 = seg["start0"]+l0
                g1 = seg["start0"]+l1
            else:
                g0 = seg["end0"]-l1
                g1 = seg["end0"]-l0

            out.append(
                (seg["chrom"],g0,g1)
            )

            mapped += g1-g0

        cum += L

    if mapped != b-a:
        raise ValueError(
            f"INCOMPLETE_WINDOW:{a}-{b}:{mapped}"
        )

    return out


def block_sig(blocks):
    return "^".join(
        f"{chrom}:{a}-{b}"
        for chrom,a,b in blocks
    )

# --------------------------------------------------
# Formal target TMD intervals
# --------------------------------------------------

tmds = defaultdict(list)

for _,r in f.iterrows():

    tx = str(r["target_transcript"])

    tmds[tx].append({
        "index": int(r["TMD_index"]),
        "start":
            int(r["target_tmd_cds_nt_start0"]),
        "end":
            int(r["target_tmd_cds_nt_end0_exclusive"]),
    })


def overlaps_tmd(tx,a,b,exclude_index=None):

    for z in tmds[tx]:

        if (
            exclude_index is not None
            and z["index"] == exclude_index
        ):
            continue

        if a < z["end"] and b > z["start"]:
            return True

    return False


def window(
    tx,
    path,
    anchor,
    center,
    exclude_real_index=None,
):

    a = anchor + center - 15
    b = anchor + center + 16

    blocks = interval_to_blocks(
        path,a,b
    )

    if blocks is None:
        return "0","OUTSIDE_FROZEN_HUMAN_SPAN"

    if overlaps_tmd(
        tx,a,b,
        exclude_index=exclude_real_index
    ):
        return "0","OVERLAP_OTHER_TMD"

    return "1","VALID|" + block_sig(blocks)

# --------------------------------------------------
# Transfer
# --------------------------------------------------

rows = []

real_anchor_fail = 0
pseudo_segment_overlap_fail = 0
source_real_mask_mismatch = 0
source_pseudo_mask_mismatch = 0
target_pair_mask_mismatch = 0
stop_feature_fail = 0

for _,r in x.iterrows():

    tx = str(r["target_transcript"])
    idx = int(r["TMD_index"])
    delta = int(r["source_target_span_delta"])

    try:
        path, translated_len, stop_len = (
            make_path(tx,delta)
        )
    except Exception:
        if delta == 3:
            stop_feature_fail += 1
        raise

    if translated_len != int(
        r["target_cds_length"]
    ):
        raise SystemExit(
            f"TRANSLATED_CDS_LENGTH_FAIL:{tx}"
        )

    if translated_len + stop_len != int(
        r["CDS_length_nt"]
    ):
        raise SystemExit(
            f"FROZEN_SPAN_LENGTH_FAIL:{tx}"
        )

    real_anchor = int(
        r["target_tmd_cds_nt_start0"]
    )

    # real anchor is always inside translated CDS
    real_anchor_g = offset_to_genome(
        path,real_anchor
    )

    if real_anchor_g != int(
        r["anchor_genomic_start0"]
    ):
        real_anchor_fail += 1

    pseudo_anchor = int(
        r["pseudo_anchor_cds0"]
    )

    ps0 = int(
        r["pseudo_segment_cds_start0"]
    )

    ps1 = int(
        r["pseudo_segment_cds_end0_exclusive"]
    )

    plen = int(
        r["pseudo_segment_length_nt"]
    )

    rlen = int(
        r["real_TMD_length_nt"]
    )

    total_span = (
        translated_len + stop_len
    )

    if not (
        0 <= pseudo_anchor < total_span
        and 0 <= ps0 < ps1 <= total_span
    ):
        raise SystemExit(
            f"PSEUDO_SPAN_BOUNDS_FAIL:{tx}|{idx}"
        )

    if ps1-ps0 != plen or plen != rlen:
        raise SystemExit(
            f"PSEUDO_LENGTH_GATE_FAIL:{tx}|{idx}"
        )

    if ps0 % 3 != 0:
        raise SystemExit(
            f"PSEUDO_CODON_ALIGNMENT_FAIL:{tx}|{idx}"
        )

    if pseudo_anchor != ps0:
        raise SystemExit(
            f"PSEUDO_ANCHOR_GATE_FAIL:{tx}|{idx}"
        )

    if overlaps_tmd(
        tx,ps0,ps1,
        exclude_index=None
    ):
        pseudo_segment_overlap_fail += 1

    pseudo_anchor_g = offset_to_genome(
        path,pseudo_anchor
    )

    real_bits = []
    pseudo_bits = []

    real_sigs = []
    pseudo_sigs = []

    for c in CENTERS:

        br,sr = window(
            tx,path,real_anchor,c,
            exclude_real_index=idx
        )

        bp,sp = window(
            tx,path,pseudo_anchor,c,
            exclude_real_index=None
        )

        real_bits.append(br)
        pseudo_bits.append(bp)

        real_sigs.append(
            f"{c}:{sr}"
        )

        pseudo_sigs.append(
            f"{c}:{sp}"
        )

    real_mask = "".join(real_bits)
    pseudo_mask = "".join(pseudo_bits)

    source_real = str(
        r["real_scan_validity_mask"]
    )

    source_pseudo = str(
        r["pseudo_scan_validity_mask"]
    )

    if real_mask != source_real:
        source_real_mask_mismatch += 1

    if pseudo_mask != source_pseudo:
        source_pseudo_mask_mismatch += 1

    if real_mask != pseudo_mask:
        target_pair_mask_mismatch += 1

    i75 = CENTERS.index(75)

    rows.append({
        "source_transcript":
            r["transcript"],

        "target_transcript":
            tx,

        "match_class":
            r["match_class"],

        "gene_id":
            r["gene_id"],

        "gene_name":
            r["gene_name"],

        "chrom":
            r["chrom"],

        "strand":
            r["strand"],

        "TMD_index":
            idx,

        "TMD_position_class":
            r["TMD_position_class"],

        "source_CDS_length_nt":
            int(r["CDS_length_nt"]),

        "target_translated_cds_length":
            translated_len,

        "target_stop_codon_nt_used":
            stop_len,

        "real_anchor_cds0":
            real_anchor,

        "real_anchor_genomic_start0":
            real_anchor_g,

        "pseudo_anchor_cds0":
            pseudo_anchor,

        "pseudo_anchor_genomic_start0":
            pseudo_anchor_g,

        "pseudo_segment_cds_start0":
            ps0,

        "pseudo_segment_cds_end0_exclusive":
            ps1,

        "real_scan_validity_mask_source":
            source_real,

        "pseudo_scan_validity_mask_source":
            source_pseudo,

        "real_scan_validity_mask_target":
            real_mask,

        "pseudo_scan_validity_mask_target":
            pseudo_mask,

        "real_primary_signature":
            real_sigs[i75].split(":",1)[1],

        "pseudo_primary_signature":
            pseudo_sigs[i75].split(":",1)[1],

        "real_full_11center_geometry":
            ";".join(real_sigs),

        "pseudo_full_11center_geometry":
            ";".join(pseudo_sigs),

        "positional_mismatch_nt":
            int(r["positional_mismatch_nt"]),
    })

out = pd.DataFrame(rows)

print()
print("TRANSFERRED_PAIRS =", len(out))
print(
    "TRANSFERRED_TRANSCRIPTS =",
    out["target_transcript"].nunique()
)

print(
    "TRANSFERRED_FIRST_PAIRS =",
    int(
        out["TMD_position_class"]
        .eq("FIRST")
        .sum()
    )
)

print(
    "TRANSFERRED_INTERNAL_PAIRS =",
    int(
        out["TMD_position_class"]
        .eq("INTERNAL")
        .sum()
    )
)

print()
print(
    "REAL_ANCHOR_REPRODUCTION_FAILURES =",
    real_anchor_fail
)

print(
    "PSEUDO_SEGMENT_TMD_OVERLAP_FAILURES =",
    pseudo_segment_overlap_fail
)

print(
    "DELTA3_STOP_FEATURE_FAILURES =",
    stop_feature_fail
)

print(
    "SOURCE_REAL_MASK_MISMATCHES =",
    source_real_mask_mismatch
)

print(
    "SOURCE_PSEUDO_MASK_MISMATCHES =",
    source_pseudo_mask_mismatch
)

print(
    "TARGET_REAL_PSEUDO_MASK_MISMATCHES =",
    target_pair_mask_mismatch
)

hard = [
    real_anchor_fail,
    pseudo_segment_overlap_fail,
    stop_feature_fail,
    source_real_mask_mismatch,
    source_pseudo_mask_mismatch,
    target_pair_mask_mismatch,
]

if any(v != 0 for v in hard):
    raise SystemExit(
        "PSEUDO_TRANSFER_V2_HARD_GATE = FAIL"
    )

print()
print("PSEUDO_TRANSFER_V2_HARD_GATE = PASS")

# --------------------------------------------------
# Primary pair geometry accounting
# --------------------------------------------------

REALKEY = [
    "chrom",
    "strand",
    "real_anchor_genomic_start0",
    "real_primary_signature",
]

real_groups = (
    out.groupby(REALKEY)
)

sizes = real_groups.size()

dup_real_groups = int(
    (sizes > 1).sum()
)

rows_in_dup_real = int(
    sizes[sizes > 1].sum()
)

multiple_pseudo_geometry = 0

for _,grp in real_groups:
    if len(grp) > 1:
        if grp[
            "pseudo_primary_signature"
        ].nunique() > 1:
            multiple_pseudo_geometry += 1

pair_geom = (
    out[
        [
            "chrom",
            "strand",
            "real_anchor_genomic_start0",
            "real_primary_signature",
            "pseudo_primary_signature",
        ]
    ]
    .drop_duplicates()
)

print()
print(
    "DUPLICATE_REAL_PRIMARY_GEOMETRY_GROUPS =",
    dup_real_groups
)

print(
    "ROWS_IN_DUPLICATE_REAL_PRIMARY_GEOMETRY_GROUPS =",
    rows_in_dup_real
)

print(
    "REAL_GEOMETRY_GROUPS_WITH_MULTIPLE_PSEUDO_GEOMETRIES =",
    multiple_pseudo_geometry
)

print(
    "UNIQUE_PRIMARY_REAL_PSEUDO_GEOMETRY_COMBINATIONS =",
    len(pair_geom)
)

# deterministic final output
out = out.sort_values(
    [
        "target_transcript",
        "TMD_index",
    ]
).reset_index(drop=True)

out.insert(
    0,
    "pseudo_pair_id",
    [
        f"TMDPS{i:05d}"
        for i in range(1,len(out)+1)
    ]
)

out.to_csv(
    OUT,
    sep="\t",
    index=False
)

excluded[
    [
        "transcript",
        "TMD_index",
        "match_status",
    ]
].to_csv(
    EXCLUDED_OUT,
    sep="\t",
    index=False
)

print()
print("PSEUDO_TRANSFER_V2_GATE = PASS")
print("BIOLOGICAL_TMD_SIGNAL = NOT_VIEWED")
