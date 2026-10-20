#!/usr/bin/env python3
"""Build and validate the frozen annotation-only primary TMD anchor table.

Audit-only constructor for Stage 0C. This is not production Module 05 and does
not read Ribo-seq data or create pseudo-anchors/null coordinates.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path("REPRO_WORKSPACE/GSE133111_stage0C")
SEGMENTS = ROOT / "reference/deeptmhmm_raw/deeptmhmm_full_normalized_topology_segments.tsv"
FEATURES = ROOT / "reference/author_reference/human_features.tsv"
INPUT_MANIFEST = ROOT / "reference/deeptmhmm_raw/deeptmhmm_full_input_manifest.tsv"
PREDICTION_CLASSES = ROOT / "reference/deeptmhmm_raw/deeptmhmm_full_raw_prediction_classes.tsv"
OUTPUT = ROOT / "reference/tmd/primary_tmd_anchors.tsv"
QC_DIR = ROOT / "qc/tmd_anchor_table"
SUMMARY = QC_DIR / "primary_tmd_anchor_table_validation_summary.json"

EXPECTED_HASHES = {
    SEGMENTS: "ca60babf576a8547077d68be1f52d29835ce8829a5d2b671c4ef3fd175c94040",
    FEATURES: "c54965e9bca3cca852f422402180b8f0811a9ed58e8c9fd7a18b92b0e3099f36",
    INPUT_MANIFEST: "e08a6caa640c57474199319f8edba10bd797cb5c006f5ee03a2b67dc73d72dc8",
    PREDICTION_CLASSES: "7c621e7081d1d5a91ba990cd8dd268b40a579eda12ff7181c8d8b169f17ec49b",
}

DEEPTMHMM_SOURCE_FREEZE = "stage0c_deeptmhmm_v1_2026-08-24"
RAW_ANNOTATION_FREEZE = "stage0c_deeptmhmm_raw_v1_2026-08-24"
TMD_LABEL_POLICY_FREEZE = "stage0c_tmd_label_policy_v1_2026-08-24"
TMD_ANCHOR_POLICY_FREEZE = "stage0c_tmd_anchor_policy_v1_2026-08-24"
ANCHOR_TABLE_FREEZE = "stage0c_primary_tmd_anchors_v1_2026-08-24"
APP_VERSION = "DTU/DeepTMHMM:1.0.24"
ANCHOR_DEFINITION = "TMD_N_TERMINAL_CODON_FIRST_NT_TRANSCRIPT_0BASED"

EXPECTED_ROWS = 19_784
EXPECTED_ZERO = 14_225
EXPECTED_SINGLE = 1_993
EXPECTED_MULTI = 2_627
EXPECTED_INELIGIBLE = 347

OUTPUT_FIELDS = [
    "transcript",
    "protein_length",
    "protein_sha256",
    "segment_index_raw_order",
    "TMD_index",
    "n_primary_TMDs",
    "pass_class",
    "TMD_position_class",
    "raw_deeptmhmm_label",
    "source_residue_start",
    "source_residue_end",
    "aa_start0",
    "aa_end0_exclusive",
    "cds_start",
    "cds_nt_start0",
    "cds_nt_end0_exclusive",
    "transcript_nt_start0",
    "transcript_nt_end0_exclusive",
    "primary_anchor_coordinate",
    "primary_anchor_cds_coordinate",
    "anchor_coordinate_system",
    "anchor_definition",
    "source_job_id",
    "source_batch_id",
    "deeptmhmm_app_version",
    "deeptmhmm_source_freeze_version",
    "raw_annotation_freeze_version",
    "tmd_label_policy_freeze_version",
    "tmd_anchor_policy_freeze_version",
    "anchor_table_freeze_version",
    "frozen_reference_fasta_sha256",
    "human_features_sha256",
    "normalized_raw_annotation_sha256",
    "roundtrip_match",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    for path, expected in EXPECTED_HASHES.items():
        observed = sha256(path)
        require(observed == expected, f"Frozen input hash mismatch: {path}: {observed}")

    feature_rows = read_tsv(FEATURES)
    require(len(feature_rows) == 19_192, "human_features.tsv row count is not 19,192")
    features = {row["Transcript"]: row for row in feature_rows}
    require(len(features) == 19_192, "Duplicate transcript IDs in human_features.tsv")

    manifest_rows = read_tsv(INPUT_MANIFEST)
    require(len(manifest_rows) == 19_192, "Full input manifest row count is not 19,192")
    eligible_ids = {r["transcript"] for r in manifest_rows if r["deeptmhmm_input_eligible"] == "YES"}
    ineligible_ids = {r["transcript"] for r in manifest_rows if r["deeptmhmm_input_eligible"] == "NO"}
    require(len(eligible_ids) == 18_845, "Eligible transcript count changed")
    require(len(ineligible_ids) == EXPECTED_INELIGIBLE, "Ineligible transcript count changed")
    require(not eligible_ids.intersection(ineligible_ids), "Eligibility sets overlap")

    class_rows = read_tsv(PREDICTION_CLASSES)
    class_ids = [r["transcript"] for r in class_rows]
    require(len(class_ids) == 18_845 and len(set(class_ids)) == 18_845, "Prediction-class accounting failed")
    require(set(class_ids) == eligible_ids, "Prediction-class IDs do not equal eligible IDs")

    segment_rows = read_tsv(SEGMENTS)
    accepted = [r for r in segment_rows if r["raw_deeptmhmm_label"] == "TMhelix"]
    require(len(accepted) == EXPECTED_ROWS, "Accepted TMhelix row count changed")

    by_transcript: dict[str, list[dict[str, str]]] = defaultdict(list)
    coordinate_mismatches = 0
    invalid_intervals = 0
    for row in accepted:
        tx = row["transcript"]
        require(tx in eligible_ids, f"Accepted TMD belongs to noneligible transcript: {tx}")
        require(tx in features, f"Missing feature row: {tx}")
        require(row["deeptmhmm_app_version"] == APP_VERSION, f"App version mismatch: {tx}")
        require(row["freeze_version"] == DEEPTMHMM_SOURCE_FREEZE, f"Source freeze mismatch: {tx}")
        require(row["roundtrip_match"] == "YES", f"Frozen round trip failed: {tx}")

        cds_start = int(features[tx]["UTR5_Length"])
        aa_start = int(row["normalized_aa_start0"])
        aa_end = int(row["normalized_aa_end0_exclusive"])
        cds_start_rel = int(row["cds_nt_start0"])
        cds_end_rel = int(row["cds_nt_end0_exclusive"])
        tx_start = int(row["transcript_nt_start0"])
        tx_end = int(row["transcript_nt_end0_exclusive"])
        protein_length = int(row["protein_length"])
        source_start = int(row["source_residue_start"])
        source_end = int(row["source_residue_end"])

        coordinate_mismatches += sum(
            [
                tx_start != cds_start + cds_start_rel,
                cds_start_rel != 3 * aa_start,
                cds_end_rel != 3 * aa_end,
                tx_end != cds_start + cds_end_rel,
                aa_start != source_start - 1,
                aa_end != source_end,
            ]
        )
        invalid_intervals += int(
            not (
                0 <= aa_start < aa_end <= protein_length
                and 0 <= cds_start_rel < cds_end_rel
                and 0 <= tx_start < tx_end
            )
        )
        by_transcript[tx].append(row)

    require(coordinate_mismatches == 0, f"Coordinate invariant mismatches: {coordinate_mismatches}")
    require(invalid_intervals == 0, f"Invalid intervals: {invalid_intervals}")

    duplicate_intervals = 0
    overlapping_intervals = 0
    adjacent_intervals = 0
    duplicate_anchors = 0
    output_rows: list[dict[str, str | int]] = []

    for tx in sorted(by_transcript):
        rows = sorted(
            by_transcript[tx],
            key=lambda r: (
                int(r["transcript_nt_start0"]),
                int(r["transcript_nt_end0_exclusive"]),
                int(r["segment_index_raw_order"]),
            ),
        )
        intervals = [
            (int(r["transcript_nt_start0"]), int(r["transcript_nt_end0_exclusive"])) for r in rows
        ]
        duplicate_intervals += len(intervals) - len(set(intervals))
        anchors = [start for start, _ in intervals]
        duplicate_anchors += len(anchors) - len(set(anchors))
        for left, right in zip(intervals, intervals[1:]):
            overlapping_intervals += int(right[0] < left[1])
            adjacent_intervals += int(right[0] == left[1])

        n_tmds = len(rows)
        pass_class = "SINGLE_PASS" if n_tmds == 1 else "MULTI_PASS"
        for index, row in enumerate(rows, 1):
            cds_start = int(features[tx]["UTR5_Length"])
            tx_start = int(row["transcript_nt_start0"])
            cds_start_rel = int(row["cds_nt_start0"])
            out = {
                "transcript": tx,
                "protein_length": int(row["protein_length"]),
                "protein_sha256": row["protein_sha256"],
                "segment_index_raw_order": int(row["segment_index_raw_order"]),
                "TMD_index": index,
                "n_primary_TMDs": n_tmds,
                "pass_class": pass_class,
                "TMD_position_class": "FIRST" if index == 1 else "INTERNAL",
                "raw_deeptmhmm_label": row["raw_deeptmhmm_label"],
                "source_residue_start": int(row["source_residue_start"]),
                "source_residue_end": int(row["source_residue_end"]),
                "aa_start0": int(row["normalized_aa_start0"]),
                "aa_end0_exclusive": int(row["normalized_aa_end0_exclusive"]),
                "cds_start": cds_start,
                "cds_nt_start0": cds_start_rel,
                "cds_nt_end0_exclusive": int(row["cds_nt_end0_exclusive"]),
                "transcript_nt_start0": tx_start,
                "transcript_nt_end0_exclusive": int(row["transcript_nt_end0_exclusive"]),
                "primary_anchor_coordinate": tx_start,
                "primary_anchor_cds_coordinate": cds_start_rel,
                "anchor_coordinate_system": "transcript_nt_0based",
                "anchor_definition": ANCHOR_DEFINITION,
                "source_job_id": row["source_job_id"],
                "source_batch_id": row["source_batch_id"],
                "deeptmhmm_app_version": row["deeptmhmm_app_version"],
                "deeptmhmm_source_freeze_version": DEEPTMHMM_SOURCE_FREEZE,
                "raw_annotation_freeze_version": RAW_ANNOTATION_FREEZE,
                "tmd_label_policy_freeze_version": TMD_LABEL_POLICY_FREEZE,
                "tmd_anchor_policy_freeze_version": TMD_ANCHOR_POLICY_FREEZE,
                "anchor_table_freeze_version": ANCHOR_TABLE_FREEZE,
                "frozen_reference_fasta_sha256": "ecda9babdf126cffa6b91e5731c355db9e40db96e2c5c5f737125342961dbf56",
                "human_features_sha256": EXPECTED_HASHES[FEATURES],
                "normalized_raw_annotation_sha256": EXPECTED_HASHES[SEGMENTS],
                "roundtrip_match": row["roundtrip_match"],
            }
            output_rows.append(out)

    require(duplicate_intervals == 0, f"Duplicate accepted intervals: {duplicate_intervals}")
    require(overlapping_intervals == 0, f"Overlapping accepted intervals: {overlapping_intervals}")
    require(adjacent_intervals == 0, f"Adjacent accepted intervals: {adjacent_intervals}")
    require(duplicate_anchors == 0, f"Duplicate anchors within transcript: {duplicate_anchors}")

    output_rows.sort(
        key=lambda r: (
            str(r["transcript"]),
            int(r["transcript_nt_start0"]),
            int(r["transcript_nt_end0_exclusive"]),
        )
    )
    require(len(output_rows) == EXPECTED_ROWS, "Final output row count mismatch")

    keys = [(str(r["transcript"]), int(r["TMD_index"])) for r in output_rows]
    duplicate_tmd_index = len(keys) - len(set(keys))
    classification_mismatches = 0
    for tx, rows in by_transcript.items():
        out = [r for r in output_rows if r["transcript"] == tx]
        n = len(rows)
        classification_mismatches += int([int(r["TMD_index"]) for r in out] != list(range(1, n + 1)))
        classification_mismatches += sum(int(r["n_primary_TMDs"]) != n for r in out)
        classification_mismatches += sum((r["TMD_position_class"] == "FIRST") != (int(r["TMD_index"]) == 1) for r in out)
        classification_mismatches += sum((r["pass_class"] == "SINGLE_PASS") != (n == 1) for r in out)

    first_rows = sum(r["TMD_position_class"] == "FIRST" for r in output_rows)
    internal_rows = sum(r["TMD_position_class"] == "INTERNAL" for r in output_rows)
    single_transcripts = sum(len(rows) == 1 for rows in by_transcript.values())
    multi_transcripts = sum(len(rows) >= 2 for rows in by_transcript.values())
    zero_transcripts = len(eligible_ids - set(by_transcript))

    require(duplicate_tmd_index == 0, "Duplicate transcript/TMD_index keys")
    require(classification_mismatches == 0, f"Classification mismatches: {classification_mismatches}")
    require(first_rows == len(by_transcript), "FIRST count mismatch")
    require(internal_rows == EXPECTED_ROWS - first_rows, "INTERNAL count mismatch")
    require(single_transcripts == EXPECTED_SINGLE, "Single-pass transcript count changed")
    require(multi_transcripts == EXPECTED_MULTI, "Multi-pass transcript count changed")
    require(zero_transcripts == EXPECTED_ZERO, "Zero-accepted-TMD transcript count changed")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    require(not OUTPUT.exists(), f"Refusing to overwrite existing authoritative output: {OUTPUT}")
    temporary = OUTPUT.with_suffix(".tsv.tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    os.replace(temporary, OUTPUT)
    output_hash = sha256(OUTPUT)

    summary = {
        "accepted_primary_tmd_segments": EXPECTED_ROWS,
        "adjacent_accepted_intervals": adjacent_intervals,
        "anchor_definition": ANCHOR_DEFINITION,
        "anchor_table": str(OUTPUT),
        "anchor_table_freeze_version": ANCHOR_TABLE_FREEZE,
        "anchor_table_sha256": output_hash,
        "biological_comparison": "NOT PERFORMED",
        "classification_invariant_mismatches": classification_mismatches,
        "coordinate_invariant_mismatches": coordinate_mismatches,
        "duplicate_accepted_intervals": duplicate_intervals,
        "duplicate_anchors_within_transcript": duplicate_anchors,
        "duplicate_transcript_tmd_index": duplicate_tmd_index,
        "first_tmd_rows": first_rows,
        "internal_tmd_rows": internal_rows,
        "invalid_intervals": invalid_intervals,
        "module_05": "NOT CREATED",
        "multi_pass_transcripts": multi_transcripts,
        "overlapping_accepted_intervals": overlapping_intervals,
        "pseudo_anchor": "NOT YET FROZEN",
        "ribo_seq_join": "NOT PERFORMED",
        "rows": len(output_rows),
        "single_pass_transcripts": single_transcripts,
        "technically_ineligible_unclassified": len(ineligible_ids),
        "tmd_blind": "PRESERVED",
        "unique_tmd_transcripts": len(by_transcript),
        "zero_accepted_tmd_successful_transcripts": zero_transcripts,
    }
    summary_tmp = SUMMARY.with_suffix(".json.tmp")
    summary_tmp.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    os.replace(summary_tmp, SUMMARY)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
