from pathlib import Path
import pandas as pd

ANCHOR = Path(
    "REPRO_WORKSPACE/GSE133111_stage0C/reference/tmd/"
    "primary_tmd_anchors.tsv"
)

BRIDGE = Path(
    "tmd_anchor_preregister/bridge/"
    "APPRIS19192_TO_ENSEMBL104_EXACT_PROTEIN_BRIDGE_v1.tsv"
)

TARGET = Path(
    "tmd_anchor_preregister/bridge/"
    "ENSEMBL104_TRANSLATED_PROTEIN_SHA_v1.tsv"
)

OUTDIR = Path("tmd_anchor_preregister/bridge")
OUTDIR.mkdir(parents=True, exist_ok=True)

a = pd.read_csv(ANCHOR, sep="\t")
b = pd.read_csv(BRIDGE, sep="\t")
t = pd.read_csv(TARGET, sep="\t")

admissible_classes = {
    "SAME_ENST_EXACT",
    "UNIQUE_EXACT_OTHER_ENST",
}

# --------------------------------------------------
# Transcript-level bridge universe
# --------------------------------------------------

bt = b[b["match_class"].isin(admissible_classes)].copy()

target_source_n = (
    bt.groupby("target_transcript")["source_transcript"]
      .nunique()
)

colliding_targets = set(
    target_source_n[target_source_n > 1].index
)

bt["target_collision"] = (
    bt["target_transcript"].isin(colliding_targets)
)

tm = t[
    [
        "target_transcript",
        "gene_id",
        "gene_name",
        "transcript_biotype",
        "cds_length",
        "protein_length",
        "protein_sha256",
        "internal_stop",
    ]
].copy()

tm = tm.rename(columns={
    "cds_length": "target_cds_length",
    "protein_length": "target_protein_length",
    "protein_sha256": "target_protein_sha256",
    "internal_stop": "target_internal_stop",
})

bt = bt.merge(
    tm,
    on="target_transcript",
    how="left",
    validate="many_to_one"
)

bt["target_metadata_present"] = (
    bt["target_protein_sha256"].notna()
)

bt["protein_length_match"] = (
    bt["source_protein_length"].astype(int)
    ==
    bt["target_protein_length"].fillna(-1).astype(int)
)

bt["protein_sha_match"] = (
    bt["source_protein_sha256"]
    ==
    bt["target_protein_sha256"]
)

bt["complete_codon_cds"] = (
    bt["target_cds_length"].fillna(-1).astype(int)
    ==
    3 * bt["target_protein_length"].fillna(-1).astype(int)
)

bt["no_internal_stop"] = (
    bt["target_internal_stop"].eq("NO")
)

bt["transcript_structural_pass"] = (
    (~bt["target_collision"])
    & bt["target_metadata_present"]
    & bt["protein_length_match"]
    & bt["protein_sha_match"]
    & bt["complete_codon_cds"]
    & bt["no_internal_stop"]
)

bt.to_csv(
    OUTDIR / "TMD_BRIDGE_TRANSCRIPT_STRUCTURAL_QC_v1.tsv",
    sep="\t",
    index=False
)

# --------------------------------------------------
# Materialize TMD rows
# --------------------------------------------------

x = a.merge(
    bt[
        [
            "source_transcript",
            "match_class",
            "target_transcript",
            "gene_id",
            "gene_name",
            "transcript_biotype",
            "target_cds_length",
            "target_protein_length",
            "target_collision",
            "transcript_structural_pass",
        ]
    ],
    left_on="transcript",
    right_on="source_transcript",
    how="left",
    validate="many_to_one"
)

x["bridge_class_admissible"] = (
    x["match_class"].isin(admissible_classes)
)

x["tmd_bounds_pass"] = (
    (x["aa_start0"].astype(int) >= 0)
    &
    (
        x["aa_end0_exclusive"].astype(int)
        >
        x["aa_start0"].astype(int)
    )
    &
    (
        x["aa_end0_exclusive"].astype(int)
        <=
        x["target_protein_length"].fillna(-1).astype(int)
    )
)

x["tmhelix_label_pass"] = (
    x["raw_deeptmhmm_label"].eq("TMhelix")
)

x["roundtrip_pass"] = (
    x["roundtrip_match"].eq("YES")
)

x["formal_anchor_candidate"] = (
    x["bridge_class_admissible"]
    & x["transcript_structural_pass"].fillna(False)
    & x["tmd_bounds_pass"]
    & x["tmhelix_label_pass"]
    & x["roundtrip_pass"]
)

formal = x[x["formal_anchor_candidate"]].copy()

formal["target_tmd_cds_nt_start0"] = (
    formal["aa_start0"].astype(int) * 3
)

formal["target_tmd_cds_nt_end0_exclusive"] = (
    formal["aa_end0_exclusive"].astype(int) * 3
)

formal.to_csv(
    OUTDIR / "FORMAL_TMD_ANCHOR_CANDIDATES_v1.tsv",
    sep="\t",
    index=False
)

# --------------------------------------------------
# Audit summaries
# --------------------------------------------------

print("SOURCE_TMD_TRANSCRIPTS =", a["transcript"].nunique())
print("SOURCE_TMD_ROWS =", len(a))

print()
print("PRESTRUCTURAL_ADMISSIBLE_TRANSCRIPTS =", len(bt))
print("TARGET_TRANSCRIPT_COLLISIONS =", len(colliding_targets))

print(
    "STRUCTURAL_PASS_TRANSCRIPTS =",
    int(bt["transcript_structural_pass"].sum())
)

print(
    "STRUCTURAL_FAIL_TRANSCRIPTS =",
    int((~bt["transcript_structural_pass"]).sum())
)

print()
print("FORMAL_TMD_ROWS =", len(formal))
print(
    "FORMAL_TMD_TRANSCRIPTS =",
    formal["target_transcript"].nunique()
)

print(
    "FORMAL_FIRST_TMD_ROWS =",
    int((formal["TMD_position_class"] == "FIRST").sum())
)

print(
    "FORMAL_INTERNAL_TMD_ROWS =",
    int((formal["TMD_position_class"] == "INTERNAL").sum())
)

print()
print("=== STRUCTURAL FAILURE COMPONENTS ===")

print(
    "TARGET_COLLISION_TRANSCRIPTS =",
    int(bt["target_collision"].sum())
)

print(
    "PROTEIN_LENGTH_MISMATCH_TRANSCRIPTS =",
    int((~bt["protein_length_match"]).sum())
)

print(
    "PROTEIN_SHA_MISMATCH_TRANSCRIPTS =",
    int((~bt["protein_sha_match"]).sum())
)

print(
    "INCOMPLETE_CODON_CDS_TRANSCRIPTS =",
    int((~bt["complete_codon_cds"]).sum())
)

print(
    "INTERNAL_STOP_TRANSCRIPTS =",
    int((~bt["no_internal_stop"]).sum())
)

print()
print(
    "FORMAL_TMD_ANCHOR_CANDIDATE_GATE =",
    "PASS"
)
