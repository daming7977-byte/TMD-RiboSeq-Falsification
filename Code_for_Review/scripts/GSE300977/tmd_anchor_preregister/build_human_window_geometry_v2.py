from pathlib import Path
from collections import defaultdict
import re
import pandas as pd

GTF = Path(
    "mapping_reference_blinded/"
    "Homo_sapiens.GRCh38.104.gtf"
)

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

OUTDIR = Path("tmd_anchor_preregister/bridge")

GEOM_OUT = (
    OUTDIR /
    "GSE300977_TRANSCRIPT_TMD_WINDOW_GEOMETRY_v2.tsv"
)

DUP_OUT = (
    OUTDIR /
    "DUPLICATE_ANCHOR_WINDOW_GEOMETRY_AUDIT_v2.tsv"
)

CENTERS = [-30,-15,0,15,30,45,60,75,90,105,120]

# --------------------------------------------------
# Load formal/source metadata
# --------------------------------------------------

f = pd.read_csv(FORMAL, sep="\t")
gg = pd.read_csv(GENOMIC, sep="\t")
p = pd.read_csv(PSEUDO, sep="\t", dtype=str)

f["TMD_index"] = f["TMD_index"].astype(int)
gg["TMD_index"] = gg["TMD_index"].astype(int)
p["TMD_index"] = p["TMD_index"].astype(int)

src = p[
    [
        "transcript",
        "TMD_index",
        "CDS_length_nt",
        "real_primary_window_valid",
        "real_scan_validity_mask",
    ]
].copy()

src["source_CDS_length_nt"] = pd.to_numeric(
    src["CDS_length_nt"]
)

src = src.drop(columns=["CDS_length_nt"])

x = f.merge(
    src,
    on=["transcript","TMD_index"],
    how="left",
    validate="one_to_one"
)

gcols = [
    "target_transcript",
    "TMD_index",
    "chrom",
    "strand",
    "anchor_genomic_start0",
]

x = x.merge(
    gg[gcols],
    on=["target_transcript","TMD_index"],
    how="left",
    validate="one_to_one"
)

if len(x) != 19113:
    raise SystemExit(
        f"FORMAL_ROW_GATE_FAIL:{len(x)}"
    )

if x[
    [
        "source_CDS_length_nt",
        "real_scan_validity_mask",
        "anchor_genomic_start0",
    ]
].isna().any().any():
    raise SystemExit("SOURCE_METADATA_GATE_FAIL")

x["target_cds_length"] = pd.to_numeric(
    x["target_cds_length"]
)

x["source_target_span_delta"] = (
    x["source_CDS_length_nt"]
    - x["target_cds_length"]
)

allowed = {0,3}

observed = set(
    x["source_target_span_delta"]
    .astype(int)
    .unique()
)

if not observed.issubset(allowed):
    raise SystemExit(
        f"SPAN_DELTA_GATE_FAIL:{sorted(observed)}"
    )

print("FORMAL_TMD_ROWS =", len(x))

print()
print("=== SOURCE-TARGET WINDOW SPAN DELTA ===")

print(
    x["source_target_span_delta"]
    .value_counts()
    .sort_index()
    .to_string()
)

# --------------------------------------------------
# Parse CDS + stop codon
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

        feature = z[2]

        if feature not in features:
            continue

        attrs = dict(attr_re.findall(z[8]))
        tx = attrs.get("transcript_id")

        if tx is None:
            continue

        features[feature][tx].append({
            "chrom": z[0],
            "start0": int(z[3])-1,
            "end0": int(z[4]),
            "strand": z[6],
            "feature": feature,
        })


def order_segments(segs):

    if not segs:
        return []

    chroms = {a["chrom"] for a in segs}
    strands = {a["strand"] for a in segs}

    if len(chroms) != 1:
        raise ValueError(
            f"MULTI_CHROM:{chroms}"
        )

    if len(strands) != 1:
        raise ValueError(
            f"MULTI_STRAND:{strands}"
        )

    strand = next(iter(strands))

    if strand == "+":
        return sorted(
            segs,
            key=lambda a: (a["start0"],a["end0"])
        )

    return sorted(
        segs,
        key=lambda a: (a["start0"],a["end0"]),
        reverse=True
    )


def transcript_path(tx, delta):

    cds = order_segments(
        features["CDS"].get(tx, [])
    )

    if not cds:
        raise ValueError(f"NO_CDS:{tx}")

    cds_len = sum(
        a["end0"]-a["start0"]
        for a in cds
    )

    if delta == 0:
        return cds, cds_len, 0

    if delta != 3:
        raise ValueError(
            f"BAD_DELTA:{tx}:{delta}"
        )

    stop = order_segments(
        features["stop_codon"].get(tx, [])
    )

    stop_len = sum(
        a["end0"]-a["start0"]
        for a in stop
    )

    if stop_len != 3:
        raise ValueError(
            f"STOP_LENGTH_FAIL:{tx}:{stop_len}"
        )

    # Combined transcript order must contain all CDS before stop codon.
    allseg = order_segments(cds + stop)

    tags = [a["feature"] for a in allseg]

    first_stop = tags.index("stop_codon")

    if any(
        t != "CDS"
        for t in tags[:first_stop]
    ):
        raise ValueError(
            f"CDS_STOP_ORDER_FAIL:{tx}"
        )

    if any(
        t != "stop_codon"
        for t in tags[first_stop:]
    ):
        raise ValueError(
            f"CDS_STOP_ORDER_FAIL:{tx}"
        )

    # No genomic overlap between translated CDS and stop feature.
    for c in cds:
        for s in stop:
            if (
                c["chrom"] == s["chrom"]
                and c["start0"] < s["end0"]
                and c["end0"] > s["start0"]
            ):
                raise ValueError(
                    f"CDS_STOP_OVERLAP:{tx}"
                )

    return cds + stop, cds_len, stop_len


def interval_to_blocks(path,start0,end0):

    if start0 < 0 or start0 >= end0:
        return None

    total = sum(
        a["end0"]-a["start0"]
        for a in path
    )

    if end0 > total:
        return None

    blocks = []
    cum = 0
    mapped = 0

    for seg in path:

        L = seg["end0"]-seg["start0"]

        a = max(start0,cum)
        b = min(end0,cum+L)

        if a < b:

            l0 = a-cum
            l1 = b-cum

            if seg["strand"] == "+":
                g0 = seg["start0"] + l0
                g1 = seg["start0"] + l1
            else:
                g0 = seg["end0"] - l1
                g1 = seg["end0"] - l0

            blocks.append(
                (seg["chrom"],g0,g1)
            )

            mapped += g1-g0

        cum += L

    if mapped != end0-start0:
        raise ValueError(
            f"INCOMPLETE_MAPPING:"
            f"{start0}-{end0}:{mapped}"
        )

    return blocks


def block_signature(blocks):

    return "^".join(
        f"{chrom}:{a}-{b}"
        for chrom,a,b in blocks
    )

# --------------------------------------------------
# Formal TMD intervals
# --------------------------------------------------

tmd_intervals = defaultdict(list)

for _,r in x.iterrows():

    tx = str(r["target_transcript"])

    tmd_intervals[tx].append({
        "index": int(r["TMD_index"]),
        "start":
            int(r["target_tmd_cds_nt_start0"]),
        "end":
            int(r["target_tmd_cds_nt_end0_exclusive"]),
    })


def overlaps_other_tmd(
    tx,
    focal_index,
    start0,
    end0,
):

    for t in tmd_intervals[tx]:

        if t["index"] == focal_index:
            continue

        if (
            start0 < t["end"]
            and end0 > t["start"]
        ):
            return True

    return False

# --------------------------------------------------
# Build geometry
# --------------------------------------------------

rows = []

source_mask_mismatch = 0
source_primary_mismatch = 0

delta3_stop_fail = 0

for _,r in x.iterrows():

    tx = str(r["target_transcript"])
    idx = int(r["TMD_index"])

    delta = int(
        r["source_target_span_delta"]
    )

    try:
        path, translated_len, stop_len = (
            transcript_path(tx,delta)
        )
    except Exception:
        if delta == 3:
            delta3_stop_fail += 1
        raise

    expected_translated = int(
        r["target_cds_length"]
    )

    if translated_len != expected_translated:
        raise SystemExit(
            f"TARGET_TRANSLATED_CDS_FAIL:"
            f"{tx}:{translated_len}:"
            f"{expected_translated}"
        )

    total_span = translated_len + stop_len

    expected_span = int(
        r["source_CDS_length_nt"]
    )

    if total_span != expected_span:
        raise SystemExit(
            f"WINDOW_SPAN_LENGTH_FAIL:"
            f"{tx}:{total_span}:{expected_span}"
        )

    anchor = int(
        r["target_tmd_cds_nt_start0"]
    )

    record = {
        "source_transcript":
            r["transcript"],
        "target_transcript":
            tx,
        "TMD_index":
            idx,
        "gene_id":
            r["gene_id"],
        "gene_name":
            r["gene_name"],
        "TMD_position_class":
            r["TMD_position_class"],
        "chrom":
            r["chrom"],
        "strand":
            r["strand"],
        "anchor_genomic_start0":
            int(r["anchor_genomic_start0"]),
        "anchor_cds0":
            anchor,
        "target_translated_cds_length":
            translated_len,
        "target_stop_codon_nt_used":
            stop_len,
        "frozen_source_window_span_nt":
            expected_span,
        "source_target_span_delta":
            delta,
    }

    bits = []
    full_sig = []

    for c in CENTERS:

        w0 = anchor + c - 15
        w1 = anchor + c + 16

        reason = "VALID"
        blocks = None

        if (
            w0 < 0
            or w1 > total_span
        ):
            reason = (
                "OUTSIDE_FROZEN_HUMAN_SPAN"
            )
        else:
            blocks = interval_to_blocks(
                path,w0,w1
            )

            if blocks is None:
                reason = (
                    "OUTSIDE_FROZEN_HUMAN_SPAN"
                )

        if reason == "VALID":

            if overlaps_other_tmd(
                tx,idx,w0,w1
            ):
                reason = (
                    "OVERLAP_OTHER_TMD"
                )

        if reason == "VALID":
            sig = (
                "VALID|"
                + block_signature(blocks)
            )
            bit = "1"
        else:
            sig = reason
            bit = "0"

        record[
            f"center_{c:+d}_signature"
        ] = sig

        record[
            f"center_{c:+d}_valid"
        ] = (
            "YES" if bit=="1" else "NO"
        )

        bits.append(bit)
        full_sig.append(
            f"{c}:{sig}"
        )

    mask = "".join(bits)

    record[
        "target_scan_validity_mask"
    ] = mask

    record[
        "source_scan_validity_mask"
    ] = str(
        r["real_scan_validity_mask"]
    )

    record[
        "full_11center_signature"
    ] = ";".join(full_sig)

    record[
        "primary_signature"
    ] = record[
        "center_+75_signature"
    ]

    record[
        "primary_valid"
    ] = record[
        "center_+75_valid"
    ]

    if (
        mask
        !=
        str(r["real_scan_validity_mask"])
    ):
        source_mask_mismatch += 1

    source_primary = str(
        r["real_primary_window_valid"]
    )

    if (
        record["primary_valid"]
        != source_primary
    ):
        source_primary_mismatch += 1

    rows.append(record)

geom = pd.DataFrame(rows)

print()
print(
    "SOURCE_REAL_MASK_MISMATCHES =",
    source_mask_mismatch
)

print(
    "SOURCE_PRIMARY_VALIDITY_MISMATCHES =",
    source_primary_mismatch
)

print(
    "DELTA3_STOP_FEATURE_FAILURES =",
    delta3_stop_fail
)

if source_mask_mismatch != 0:
    raise SystemExit(
        "SOURCE_MASK_REPRODUCTION_GATE = FAIL"
    )

if source_primary_mismatch != 0:
    raise SystemExit(
        "SOURCE_PRIMARY_REPRODUCTION_GATE = FAIL"
    )

print("SOURCE_MASK_REPRODUCTION_GATE = PASS")
print("SOURCE_PRIMARY_REPRODUCTION_GATE = PASS")

geom.to_csv(
    GEOM_OUT,
    sep="\t",
    index=False
)

# --------------------------------------------------
# Structural summary
# --------------------------------------------------

primary_valid_n = int(
    (geom["primary_valid"]=="YES").sum()
)

primary_invalid_n = (
    len(geom)-primary_valid_n
)

print()
print(
    "PRIMARY_STRUCTURALLY_VALID_ROWS =",
    primary_valid_n
)

print(
    "PRIMARY_STRUCTURALLY_INVALID_ROWS =",
    primary_invalid_n
)

# --------------------------------------------------
# Duplicate-anchor geometry audit
# --------------------------------------------------

KEY = [
    "chrom",
    "strand",
    "anchor_genomic_start0",
]

audit_rows = []

for key,grp in geom.groupby(KEY):

    if len(grp) <= 1:
        continue

    full_n = grp[
        "full_11center_signature"
    ].nunique()

    primary_n = grp[
        "primary_signature"
    ].nunique()

    valid_cols = [
        f"center_{c:+d}_valid"
        for c in CENTERS
    ]

    masks = (
        grp[valid_cols]
        .astype(str)
        .agg("|".join,axis=1)
    )

    mask_n = masks.nunique()

    audit_rows.append({
        "chrom": key[0],
        "strand": key[1],
        "anchor_genomic_start0":
            int(key[2]),
        "n_source_rows":
            len(grp),
        "n_primary_geometry_signatures":
            primary_n,
        "n_full_11center_geometry_signatures":
            full_n,
        "n_validity_masks":
            mask_n,
        "primary_geometry_consensus":
            "YES" if primary_n==1 else "NO",
        "full_geometry_consensus":
            "YES" if full_n==1 else "NO",
        "validity_mask_consensus":
            "YES" if mask_n==1 else "NO",
    })

audit = pd.DataFrame(audit_rows)

audit.to_csv(
    DUP_OUT,
    sep="\t",
    index=False
)

dup_groups = len(audit)

primary_discord = int(
    (
        audit["primary_geometry_consensus"]
        =="NO"
    ).sum()
)

full_discord = int(
    (
        audit["full_geometry_consensus"]
        =="NO"
    ).sum()
)

mask_discord = int(
    (
        audit["validity_mask_consensus"]
        =="NO"
    ).sum()
)

print()
print("DUPLICATE_ANCHOR_GROUPS =", dup_groups)

print(
    "DUP_GROUPS_PRIMARY_GEOMETRY_DISCORDANT =",
    primary_discord
)

print(
    "DUP_GROUPS_FULL_11CENTER_GEOMETRY_DISCORDANT =",
    full_discord
)

print(
    "DUP_GROUPS_VALIDITY_MASK_DISCORDANT =",
    mask_discord
)

print()
print("HUMAN_WINDOW_GEOMETRY_V2_GATE = PASS")
print("BIOLOGICAL_TMD_SIGNAL = NOT_VIEWED")
