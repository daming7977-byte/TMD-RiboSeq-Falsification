from pathlib import Path
import pandas as pd

INP = Path(
    "tmd_anchor_preregister/bridge/"
    "GSE300977_FORMAL_GENOMIC_TMD_ANCHORS_v1.tsv"
)

OUTDIR = Path("tmd_anchor_preregister/bridge")

GROUP_OUT = OUTDIR / "DUPLICATE_GENOMIC_ANCHOR_GROUPS_v1.tsv"
UNIQUE_OUT = OUTDIR / "GSE300977_FORMAL_UNIQUE_GENOMIC_TMD_ANCHORS_v1.tsv"

d = pd.read_csv(INP, sep="\t")

KEY = [
    "chrom",
    "strand",
    "anchor_genomic_start0",
]

if len(d) != 19113:
    raise SystemExit(f"INPUT_ROW_COUNT_FAIL:{len(d)}")

gsize = d.groupby(KEY).size()

dup_keys = gsize[gsize > 1]

if len(dup_keys) != 105:
    raise SystemExit(
        f"DUPLICATE_GROUP_COUNT_FAIL:{len(dup_keys)}"
    )

if int(dup_keys.sum()) != 217:
    raise SystemExit(
        f"DUPLICATE_ROW_COUNT_FAIL:{int(dup_keys.sum())}"
    )

def join_unique(x):
    return ",".join(
        sorted({
            str(v)
            for v in x
            if pd.notna(v) and str(v) != ""
        })
    )

rows = []
duplicate_rows = []

for key, grp in d.groupby(KEY, sort=False):

    chrom, strand, anchor = key

    txs = sorted(set(grp["target_transcript"].astype(str)))
    genes = sorted(set(grp["gene_id"].dropna().astype(str)))
    gene_names = sorted(set(grp["gene_name"].dropna().astype(str)))
    classes = sorted(set(grp["TMD_position_class"].astype(str)))
    indices = sorted(set(grp["TMD_index"].astype(int)))

    blocks = sorted(
        set(
            grp[
                "tmd_genomic_blocks_transcript_order"
            ].astype(str)
        )
    )

    if len(classes) == 1:
        pos_consensus = classes[0]
        stratifiable = "YES"
    else:
        pos_consensus = "MIXED"
        stratifiable = "NO"

    if len(genes) == 1:
        gene_assignment = genes[0]
        gene_unique = "YES"
    else:
        gene_assignment = "AMBIGUOUS_MULTI_GENE"
        gene_unique = "NO"

    if len(grp) == 1:
        category = "SINGLETON"
    elif len(genes) == 1 and len(classes) == 1:
        category = "DUP_SAME_GENE_SAME_CLASS"
    elif len(genes) == 1 and len(classes) > 1:
        category = "DUP_SAME_GENE_MIXED_CLASS"
    elif len(genes) > 1 and len(classes) == 1:
        category = "DUP_MULTI_GENE_SAME_CLASS"
    else:
        category = "DUP_MULTI_GENE_MIXED_CLASS"

    row = {
        "chrom": chrom,
        "strand": strand,
        "anchor_genomic_start0": int(anchor),
        "anchor_genomic_end0": int(anchor) + 1,
        "anchor_genomic_1based": int(anchor) + 1,

        "n_source_rows": len(grp),

        "n_target_transcripts": len(txs),
        "target_transcripts": ",".join(txs),

        "n_gene_ids": len(genes),
        "gene_ids": ",".join(genes),
        "gene_names": ",".join(gene_names),
        "gene_assignment": gene_assignment,
        "unique_gene_assignment": gene_unique,

        "n_position_classes": len(classes),
        "position_classes": ",".join(classes),
        "position_class_consensus": pos_consensus,
        "position_stratifiable": stratifiable,

        "n_tmd_indices": len(indices),
        "tmd_indices": ",".join(map(str, indices)),

        "n_distinct_tmd_blocks": len(blocks),
        "tmd_blocks": ";".join(blocks),

        "collapse_category": category,
    }

    rows.append(row)

    if len(grp) > 1:
        duplicate_rows.append(row)

unique = pd.DataFrame(rows)
dups = pd.DataFrame(duplicate_rows)

# deterministic order
unique = unique.sort_values(
    ["chrom","anchor_genomic_start0","strand"]
).reset_index(drop=True)

dups = dups.sort_values(
    ["chrom","anchor_genomic_start0","strand"]
).reset_index(drop=True)

unique.insert(
    0,
    "anchor_id",
    [f"TMDGA{i:05d}" for i in range(1, len(unique)+1)]
)

# --------------------------------------------------
# Formal checks
# --------------------------------------------------

if len(unique) != 19001:
    raise SystemExit(
        f"UNIQUE_ANCHOR_COUNT_FAIL:{len(unique)}"
    )

if len(dups) != 105:
    raise SystemExit(
        f"DUP_GROUP_OUTPUT_FAIL:{len(dups)}"
    )

if unique.duplicated(KEY).any():
    raise SystemExit("UNIQUE_KEY_DUPLICATION_FAIL")

if int(unique["n_source_rows"].sum()) != 19113:
    raise SystemExit(
        "SOURCE_ROW_ACCOUNTING_FAIL:"
        f"{int(unique['n_source_rows'].sum())}"
    )

unique.to_csv(
    UNIQUE_OUT,
    sep="\t",
    index=False
)

dups.to_csv(
    GROUP_OUT,
    sep="\t",
    index=False
)

# --------------------------------------------------
# Annotation-only summary
# --------------------------------------------------

print("INPUT_TMD_ROWS =", len(d))
print("UNIQUE_GENOMIC_ANCHORS =", len(unique))
print("COLLAPSED_EXCESS_ROWS =", len(d) - len(unique))

print()
print("DUPLICATE_GROUPS =", len(dups))
print(
    "ROWS_IN_DUPLICATE_GROUPS =",
    int(dups["n_source_rows"].sum())
)

print()
print("=== DUPLICATE CATEGORY COUNTS ===")
print(
    dups["collapse_category"]
    .value_counts()
    .to_string()
)

print()
print(
    "MIXED_FIRST_INTERNAL_ANCHORS =",
    int(
        (
            unique["position_class_consensus"]
            == "MIXED"
        ).sum()
    )
)

print(
    "AMBIGUOUS_MULTI_GENE_ANCHORS =",
    int(
        (
            unique["unique_gene_assignment"]
            == "NO"
        ).sum()
    )
)

print()
print(
    "CONSENSUS_FIRST_ANCHORS =",
    int(
        (
            unique["position_class_consensus"]
            == "FIRST"
        ).sum()
    )
)

print(
    "CONSENSUS_INTERNAL_ANCHORS =",
    int(
        (
            unique["position_class_consensus"]
            == "INTERNAL"
        ).sum()
    )
)

print()
print(
    "SOURCE_ROW_ACCOUNTING =",
    int(unique["n_source_rows"].sum())
)

print("UNIQUE_GENOMIC_ANCHOR_GATE = PASS")
print("BIOLOGICAL_TMD_SIGNAL = NOT_VIEWED")
