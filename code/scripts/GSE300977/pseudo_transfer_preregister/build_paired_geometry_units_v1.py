from pathlib import Path
import pandas as pd

INP = Path(
    "pseudo_transfer_preregister/bridge/"
    "GSE300977_TRANSFERRED_PSEUDO_PAIRS_v2.tsv"
)

OUTDIR = Path("pseudo_transfer_preregister/bridge")

LONG_OUT = (
    OUTDIR /
    "GSE300977_REAL_PSEUDO_PAIRED_GEOMETRY_UNITS_ALL_CENTERS_v1.tsv"
)

SUMMARY_OUT = (
    OUTDIR /
    "GSE300977_REAL_PSEUDO_PAIRED_GEOMETRY_SUMMARY_v1.tsv"
)

PRIMARY_OUT = (
    OUTDIR /
    "GSE300977_REAL_PSEUDO_PRIMARY_PLUS75_UNITS_v1.tsv"
)

CENTERS = [-30,-15,0,15,30,45,60,75,90,105,120]

d = pd.read_csv(INP, sep="\t", dtype=str)

if len(d) != 9624:
    raise SystemExit(
        f"SOURCE_PAIR_ROW_GATE_FAIL:{len(d)}"
    )

def parse_geometry(s):
    out = {}

    for token in str(s).split(";"):
        center, sig = token.split(":",1)
        out[int(center)] = sig

    return out


def consensus(values, mixed_label="MIXED"):

    vals = sorted(
        {
            str(v)
            for v in values
            if pd.notna(v) and str(v) != ""
        }
    )

    if len(vals) == 1:
        return vals[0]

    return mixed_label


expanded = []

for _,r in d.iterrows():

    rm = str(r["real_scan_validity_mask_target"])
    pm = str(r["pseudo_scan_validity_mask_target"])

    if (
        len(rm) != 11
        or len(pm) != 11
        or rm != pm
    ):
        raise SystemExit(
            "PAIR_VALIDITY_MASK_GATE_FAIL"
        )

    rg = parse_geometry(
        r["real_full_11center_geometry"]
    )

    pg = parse_geometry(
        r["pseudo_full_11center_geometry"]
    )

    for i,c in enumerate(CENTERS):

        if rm[i] == "0":
            continue

        rsig = rg[c]
        psig = pg[c]

        if not rsig.startswith("VALID|"):
            raise SystemExit(
                f"REAL_VALID_SIGNATURE_FAIL:{c}"
            )

        if not psig.startswith("VALID|"):
            raise SystemExit(
                f"PSEUDO_VALID_SIGNATURE_FAIL:{c}"
            )

        expanded.append({
            "center": c,

            "chrom": r["chrom"],
            "strand": r["strand"],

            "real_anchor_genomic_start0":
                int(r["real_anchor_genomic_start0"]),

            "real_window_geometry":
                rsig[len("VALID|"):],

            "pseudo_window_geometry":
                psig[len("VALID|"):],

            "source_transcript":
                r["source_transcript"],

            "target_transcript":
                r["target_transcript"],

            "gene_id":
                r["gene_id"],

            "gene_name":
                r["gene_name"],

            "TMD_index":
                int(r["TMD_index"]),

            "TMD_position_class":
                r["TMD_position_class"],
        })

e = pd.DataFrame(expanded)

if len(e) == 0:
    raise SystemExit("NO_VALID_PAIR_WINDOWS")

all_units = []
summary_rows = []

for c in CENTERS:

    z = e[e["center"].eq(c)].copy()

    REALKEY = [
        "chrom",
        "strand",
        "real_anchor_genomic_start0",
        "real_window_geometry",
    ]

    rows = []

    for key,grp in z.groupby(
        REALKEY,
        sort=False
    ):

        chrom,strand,anchor,real_geom = key

        pseudo_geoms = sorted(
            set(
                grp["pseudo_window_geometry"]
                .astype(str)
            )
        )

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
            stratifiable = "YES"
        else:
            pos_class = "MIXED"
            stratifiable = "NO"

        if len(genes) == 1:
            gene_assignment = genes[0]
            unique_gene = "YES"
        else:
            gene_assignment = "AMBIGUOUS_MULTI_GENE"
            unique_gene = "NO"

        # Number of distinct real-pseudo geometry combinations
        pair_combo_n = (
            grp[
                [
                    "real_window_geometry",
                    "pseudo_window_geometry",
                ]
            ]
            .drop_duplicates()
            .shape[0]
        )

        rows.append({
            "center": c,

            "chrom": chrom,
            "strand": strand,

            "real_anchor_genomic_start0":
                int(anchor),

            "real_window_geometry":
                real_geom,

            "n_source_pair_rows":
                len(grp),

            "n_distinct_pair_geometries":
                pair_combo_n,

            "n_distinct_pseudo_geometries":
                len(pseudo_geoms),

            "pseudo_window_geometries":
                ";".join(pseudo_geoms),

            "pseudo_control_aggregation":
                (
                    "SINGLE"
                    if len(pseudo_geoms) == 1
                    else "UNWEIGHTED_MEAN_OF_DISTINCT_PSEUDO_EFFECTS"
                ),

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
                stratifiable,

            "tmd_indices":
                ",".join(map(str,indices)),
        })

    u = pd.DataFrame(rows)

    u = u.sort_values(
        [
            "chrom",
            "real_anchor_genomic_start0",
            "strand",
            "real_window_geometry",
        ]
    ).reset_index(drop=True)

    u.insert(
        0,
        "paired_real_unit_id",
        [
            f"TMDPAIR_{c:+04d}_{i:05d}"
            for i in range(1,len(u)+1)
        ]
    )

    pair_combo_n = (
        z[
            REALKEY
            + ["pseudo_window_geometry"]
        ]
        .drop_duplicates()
        .shape[0]
    )

    multi_pseudo_n = int(
        (
            u["n_distinct_pseudo_geometries"] > 1
        ).sum()
    )

    max_pseudo_n = int(
        u["n_distinct_pseudo_geometries"].max()
    )

    exact_duplicate_excess = (
        len(z) - pair_combo_n
    )

    real_collapse_excess = (
        len(z) - len(u)
    )

    summary_rows.append({
        "center": c,

        "source_valid_pair_rows":
            len(z),

        "unique_real_geometry_units":
            len(u),

        "unique_real_pseudo_geometry_combinations":
            pair_combo_n,

        "exact_duplicate_pair_row_excess":
            exact_duplicate_excess,

        "source_rows_collapsed_to_real_units":
            real_collapse_excess,

        "real_units_with_multiple_pseudo_geometries":
            multi_pseudo_n,

        "max_distinct_pseudo_geometries_per_real":
            max_pseudo_n,

        "mixed_position_class_units":
            int(
                u["position_class_consensus"]
                .eq("MIXED")
                .sum()
            ),

        "ambiguous_multi_gene_units":
            int(
                u["unique_gene_assignment"]
                .eq("NO")
                .sum()
            ),
    })

    all_units.append(u)

long = pd.concat(
    all_units,
    ignore_index=True
)

summary = pd.DataFrame(summary_rows)

primary = long[
    long["center"].eq(75)
].copy()

# --------------------------------------------------
# Primary hard accounting
# --------------------------------------------------

s75 = summary[
    summary["center"].eq(75)
].iloc[0]

if int(
    s75["source_valid_pair_rows"]
) != 9624:
    raise SystemExit(
        "PRIMARY_SOURCE_PAIR_GATE_FAIL:"
        f"{int(s75['source_valid_pair_rows'])}"
    )

if int(
    s75[
        "unique_real_pseudo_geometry_combinations"
    ]
) != 9574:
    raise SystemExit(
        "PRIMARY_PAIR_COMBO_GATE_FAIL:"
        f"{int(s75['unique_real_pseudo_geometry_combinations'])}"
    )

long.to_csv(
    LONG_OUT,
    sep="\t",
    index=False
)

summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False
)

primary.to_csv(
    PRIMARY_OUT,
    sep="\t",
    index=False
)

print("SOURCE_TRANSFERRED_PAIRS = 9624")

print()
print("=== PAIRED CENTER SUMMARY ===")
print(summary.to_string(index=False))

print()
print(
    "PRIMARY_PLUS75_SOURCE_PAIR_ROWS =",
    int(s75["source_valid_pair_rows"])
)

print(
    "PRIMARY_PLUS75_UNIQUE_REAL_UNITS =",
    int(s75["unique_real_geometry_units"])
)

print(
    "PRIMARY_PLUS75_UNIQUE_REAL_PSEUDO_COMBINATIONS =",
    int(
        s75[
            "unique_real_pseudo_geometry_combinations"
        ]
    )
)

print(
    "PRIMARY_PLUS75_EXACT_DUPLICATE_PAIR_EXCESS =",
    int(
        s75[
            "exact_duplicate_pair_row_excess"
        ]
    )
)

print(
    "PRIMARY_PLUS75_REAL_UNITS_WITH_MULTIPLE_PSEUDOS =",
    int(
        s75[
            "real_units_with_multiple_pseudo_geometries"
        ]
    )
)

print(
    "PRIMARY_PLUS75_MAX_PSEUDOS_PER_REAL =",
    int(
        s75[
            "max_distinct_pseudo_geometries_per_real"
        ]
    )
)

print()
print("PAIRED_GEOMETRY_UNIT_GATE = PASS")
print("BIOLOGICAL_TMD_SIGNAL = NOT_VIEWED")
