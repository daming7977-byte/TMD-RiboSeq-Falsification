from pathlib import Path
import pandas as pd

GEOM = Path(
    "tmd_anchor_preregister/bridge/"
    "GSE300977_TRANSCRIPT_TMD_WINDOW_GEOMETRY_v2.tsv"
)

OUTDIR = Path("tmd_anchor_preregister/bridge")

CENTERS = [-30,-15,0,15,30,45,60,75,90,105,120]

g = pd.read_csv(GEOM, sep="\t")

if len(g) != 19113:
    raise SystemExit(
        f"GEOMETRY_ROW_GATE_FAIL:{len(g)}"
    )

required = {
    "target_transcript",
    "TMD_index",
    "gene_id",
    "gene_name",
    "TMD_position_class",
    "chrom",
    "strand",
    "anchor_genomic_start0",
}

missing = required - set(g.columns)

if missing:
    raise SystemExit(
        f"REQUIRED_COLUMN_GATE_FAIL:{sorted(missing)}"
    )

all_units = []
summaries = []

for center in CENTERS:

    sigcol = f"center_{center:+d}_signature"
    validcol = f"center_{center:+d}_valid"

    valid = g[
        g[validcol].eq("YES")
    ].copy()

    invalid = g[
        ~g[validcol].eq("YES")
    ].copy()

    outside = int(
        invalid[sigcol]
        .eq("OUTSIDE_FROZEN_HUMAN_SPAN")
        .sum()
    )

    overlap = int(
        invalid[sigcol]
        .eq("OVERLAP_OTHER_TMD")
        .sum()
    )

    if outside + overlap != len(invalid):
        raise SystemExit(
            f"INVALID_REASON_ACCOUNTING_FAIL:"
            f"{center}:"
            f"{outside}+{overlap}!={len(invalid)}"
        )

    key = [
        "chrom",
        "strand",
        "anchor_genomic_start0",
        sigcol,
    ]

    rows = []

    for k, grp in valid.groupby(
        key,
        sort=False
    ):

        chrom, strand, anchor, signature = k

        if not str(signature).startswith("VALID|"):
            raise SystemExit(
                f"VALID_SIGNATURE_FAIL:{center}"
            )

        blocks = str(signature)[len("VALID|"):]

        txs = sorted(
            set(
                grp["target_transcript"]
                .astype(str)
            )
        )

        genes = sorted(
            set(
                grp["gene_id"]
                .dropna()
                .astype(str)
            )
        )

        gene_names = sorted(
            set(
                grp["gene_name"]
                .dropna()
                .astype(str)
            )
        )

        classes = sorted(
            set(
                grp["TMD_position_class"]
                .astype(str)
            )
        )

        indices = sorted(
            set(
                grp["TMD_index"]
                .astype(int)
            )
        )

        if len(classes) == 1:
            pos_class = classes[0]
            pos_stratifiable = "YES"
        else:
            pos_class = "MIXED"
            pos_stratifiable = "NO"

        if len(genes) == 1:
            gene_assignment = genes[0]
            unique_gene = "YES"
        else:
            gene_assignment = "AMBIGUOUS_MULTI_GENE"
            unique_gene = "NO"

        rows.append({
            "center": center,
            "window_rel_start":
                center - 15,
            "window_rel_end_exclusive":
                center + 16,

            "chrom": chrom,
            "strand": strand,
            "anchor_genomic_start0":
                int(anchor),
            "anchor_genomic_end0":
                int(anchor)+1,

            "window_genomic_blocks_transcript_order":
                blocks,

            "n_source_rows":
                len(grp),

            "n_target_transcripts":
                len(txs),
            "target_transcripts":
                ",".join(txs),

            "n_gene_ids":
                len(genes),
            "gene_ids":
                ",".join(genes),
            "gene_names":
                ",".join(gene_names),

            "gene_assignment":
                gene_assignment,
            "unique_gene_assignment":
                unique_gene,

            "position_class_consensus":
                pos_class,
            "position_stratifiable":
                pos_stratifiable,

            "tmd_indices":
                ",".join(
                    map(str,indices)
                ),
        })

    c = pd.DataFrame(rows)

    if len(c) == 0:
        raise SystemExit(
            f"NO_VALID_UNIT_FAIL:{center}"
        )

    c = c.sort_values(
        [
            "chrom",
            "anchor_genomic_start0",
            "strand",
            "window_genomic_blocks_transcript_order",
        ]
    ).reset_index(drop=True)

    c.insert(
        0,
        "window_unit_id",
        [
            f"TMDW_{center:+04d}_{i:05d}"
            for i in range(1,len(c)+1)
        ]
    )

    per_anchor = (
        c.groupby(
            [
                "chrom",
                "strand",
                "anchor_genomic_start0",
            ]
        )
        .size()
    )

    unique_anchor_n = len(per_anchor)

    split_anchor_n = int(
        (per_anchor > 1).sum()
    )

    max_geometries = int(
        per_anchor.max()
    )

    mixed_n = int(
        c["position_class_consensus"]
        .eq("MIXED")
        .sum()
    )

    ambiguous_gene_n = int(
        c["unique_gene_assignment"]
        .eq("NO")
        .sum()
    )

    summaries.append({
        "center": center,

        "source_valid_rows":
            len(valid),

        "source_invalid_rows":
            len(invalid),

        "invalid_outside_frozen_human_span":
            outside,

        "invalid_overlap_other_tmd":
            overlap,

        "unique_valid_genomic_anchors":
            unique_anchor_n,

        "geometry_resolved_units":
            len(c),

        "collapsed_excess_source_rows":
            len(valid)-len(c),

        "anchors_with_multiple_valid_geometries":
            split_anchor_n,

        "max_valid_geometries_per_anchor":
            max_geometries,

        "mixed_position_class_units":
            mixed_n,

        "ambiguous_multi_gene_units":
            ambiguous_gene_n,
    })

    all_units.append(c)

long = pd.concat(
    all_units,
    ignore_index=True
)

summary = pd.DataFrame(summaries)

primary = long[
    long["center"].eq(75)
].copy()

# ----------------------------
# Hard accounting
# ----------------------------

if (
    summary["source_valid_rows"]
    + summary["source_invalid_rows"]
).ne(19113).any():

    raise SystemExit(
        "CENTER_SOURCE_ROW_ACCOUNTING_FAIL"
    )

s75 = summary[
    summary["center"].eq(75)
].iloc[0]

if int(
    s75["source_valid_rows"]
) != 15010:
    raise SystemExit(
        "PRIMARY_VALID_ROW_GATE_FAIL:"
        f"{int(s75['source_valid_rows'])}"
    )

if int(
    s75["source_invalid_rows"]
) != 4103:
    raise SystemExit(
        "PRIMARY_INVALID_ROW_GATE_FAIL:"
        f"{int(s75['source_invalid_rows'])}"
    )

if len(primary) != int(
    s75["geometry_resolved_units"]
):
    raise SystemExit(
        "PRIMARY_UNIT_ACCOUNTING_FAIL"
    )

long.to_csv(
    OUTDIR /
    "GSE300977_TMD_GEOMETRY_RESOLVED_UNITS_ALL_CENTERS_v2.tsv",
    sep="\t",
    index=False
)

summary.to_csv(
    OUTDIR /
    "GSE300977_TMD_GEOMETRY_RESOLVED_UNIT_SUMMARY_v2.tsv",
    sep="\t",
    index=False
)

primary.to_csv(
    OUTDIR /
    "GSE300977_TMD_PRIMARY_PLUS75_GEOMETRY_UNITS_v2.tsv",
    sep="\t",
    index=False
)

print("TMD_SOURCE_ROWS = 19113")

print()
print("=== CENTER SUMMARY V2 ===")
print(
    summary.to_string(index=False)
)

print()
print(
    "PRIMARY_PLUS75_SOURCE_VALID_ROWS =",
    int(s75["source_valid_rows"])
)

print(
    "PRIMARY_PLUS75_SOURCE_INVALID_ROWS =",
    int(s75["source_invalid_rows"])
)

print(
    "PRIMARY_PLUS75_UNIQUE_VALID_ANCHORS =",
    int(
        s75[
            "unique_valid_genomic_anchors"
        ]
    )
)

print(
    "PRIMARY_PLUS75_GEOMETRY_UNITS =",
    int(
        s75[
            "geometry_resolved_units"
        ]
    )
)

print(
    "PRIMARY_PLUS75_SPLIT_ANCHORS =",
    int(
        s75[
            "anchors_with_multiple_valid_geometries"
        ]
    )
)

print(
    "PRIMARY_PLUS75_MIXED_UNITS =",
    int(
        s75[
            "mixed_position_class_units"
        ]
    )
)

print(
    "PRIMARY_PLUS75_AMBIGUOUS_GENE_UNITS =",
    int(
        s75[
            "ambiguous_multi_gene_units"
        ]
    )
)

print()
print("GEOMETRY_RESOLVED_UNIT_V2_GATE = PASS")
print("BIOLOGICAL_TMD_SIGNAL = NOT_VIEWED")
