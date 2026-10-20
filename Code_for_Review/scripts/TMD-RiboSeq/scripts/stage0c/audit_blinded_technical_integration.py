#!/usr/bin/env python3
"""Annotation-only GSE133111 Stage 0C blinded integration audit.

This script reads frozen reference/topology/anchor/control manifests plus the
pre-reveal design and analysis-contract records. It does not read BAM,
calibrated P-site tables, profile matrices, benchmark labels, enrichment
outputs, or biological-result files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


CENTERS = (-30, -15, 0, 15, 30, 45, 60, 75, 90, 105, 120)
PRIMARY_CENTER = 75
PRIMARY_INDEX = CENTERS.index(PRIMARY_CENTER)
EXPECTED_OFFSETS = {27: 9, 28: 9, 29: 9, 30: 12, 32: 13, 33: 13, 34: 13, 35: 13}
EXPECTED_UNRESOLVED = {25, 26, 31}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--offsets", type=Path, required=True)
    parser.add_argument("--raw-input", type=Path, required=True)
    parser.add_argument("--raw-topology", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--pseudo-pairs", type=Path, required=True)
    parser.add_argument("--null-universe", type=Path, required=True)
    parser.add_argument("--null-columns", type=Path, required=True)
    parser.add_argument("--null-validation", type=Path, required=True)
    parser.add_argument("--null-addendum", type=Path, required=True)
    parser.add_argument("--decision-table", type=Path, required=True)
    parser.add_argument("--contract-freeze", type=Path, required=True)
    parser.add_argument("--design-matrix", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--universe-tsv", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_key_value_tsv(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    if not rows or set(rows[0]) != {"field", "value"}:
        raise ValueError(f"Expected field/value TSV: {path}")
    return {row["field"]: row["value"] for row in rows}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def overlaps(a0: int, a1: int, b0: int, b1: int) -> bool:
    return a0 < b1 and b0 < a1


def fail_if(counter: Counter[str]) -> None:
    nonzero = {key: value for key, value in counter.items() if value}
    if nonzero:
        raise RuntimeError(f"Integration invariant failure: {nonzero}")


def main() -> None:
    args = parse_args()
    for path in (
        args.fasta,
        args.features,
        args.offsets,
        args.raw_input,
        args.raw_topology,
        args.anchors,
        args.pseudo_pairs,
        args.null_universe,
        args.null_columns,
        args.null_validation,
        args.null_addendum,
        args.decision_table,
        args.contract_freeze,
        args.design_matrix,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)

    violations: Counter[str] = Counter()

    # Final pre-reveal analysis contract and outcome-independent design.
    decision_rows = read_tsv(args.decision_table)
    allowed_statuses = {
        "INHERITED / FROZEN",
        "DETERMINISTIC ADAPTATION / FROZEN",
        "NEW PAPER3 DECISION / HUMAN APPROVED AND FROZEN",
    }
    decision_by_id = {row["decision_id"]: row for row in decision_rows}
    violations["decision_row_count"] = int(len(decision_rows) != 25)
    violations["decision_duplicate_ids"] = len(decision_rows) - len(decision_by_id)
    violations["decision_id_set"] = int(set(decision_by_id) != {f"P3-{i:02d}" for i in range(1, 26)})
    violations["unresolved_decision_rows"] = sum(
        row["status"] not in allowed_statuses for row in decision_rows
    )
    p3_08 = decision_by_id.get("P3-08", {})
    p3_14 = decision_by_id.get("P3-14", {})
    violations["p3_08_contract"] = int(
        "TMhelix, signal, and Beta sheet" not in p3_08.get("recommended_option", "")
        or p3_08.get("status") != "NEW PAPER3 DECISION / HUMAN APPROVED AND FROZEN"
    )
    violations["p3_14_contract"] = int(
        "[A-30,max(E,A+136))" not in p3_14.get("recommended_option", "")
        or p3_14.get("status") != "NEW PAPER3 DECISION / HUMAN APPROVED AND FROZEN"
    )

    contract = read_key_value_tsv(args.contract_freeze)
    expected_contract = {
        "freeze_version": "stage0c_paper3_analysis_contract_v1_2026-08-25",
        "status": "PASS",
        "decision_table_total": "25",
        "decision_table_resolved": "25",
        "decision_table_unresolved": "0",
        "coordinate_null_authoritative_version": "V2",
        "coordinate_null_policy": "stage0c_coordinate_null_policy_v2_2026-08-25",
        "coordinate_null_production": "stage0c_coordinate_null_realizations_v2_2026-08-25",
        "coordinate_null_anchors": "13060",
        "coordinate_null_transcripts": "3997",
        "coordinate_null_realizations": "10000",
        "coordinate_null_reproduction_mismatches": "0",
        "coordinate_null_v1": "SUPERSEDED_BEFORE_UNBLINDING_PROHIBITED_FOR_INFERENCE",
    }
    for field, expected in expected_contract.items():
        violations[f"contract_{field}"] = int(contract.get(field) != expected)
    violations["contract_decision_table_hash"] = int(
        contract.get("decision_table_sha256") != sha256(args.decision_table)
    )
    violations["contract_design_matrix_hash"] = int(
        contract.get("design_matrix_sha256") != sha256(args.design_matrix)
    )
    violations["contract_null_addendum_hash"] = int(
        contract.get("coordinate_null_addendum_sha256") != sha256(args.null_addendum)
    )

    design_rows = read_tsv(args.design_matrix)
    expected_gsms = {f"GSM3900{i}" for i in range(208, 214)}
    observed_gsms = {row["GSM"] for row in design_rows}
    condition_counts = Counter(row["condition"] for row in design_rows)
    violations["design_row_count"] = int(len(design_rows) != 6)
    violations["design_gsm_set"] = int(observed_gsms != expected_gsms)
    violations["design_condition_balance"] = int(condition_counts != {"Control": 3, "Knockdown": 3})
    violations["design_contract_fields"] = sum(
        row["assay"] != "Ribo-seq"
        or row["primary_model_included"] != "YES"
        or row["model_formula"] != "~ condition"
        or row["reference_level"] != "Control"
        or row["coefficient"] != "conditionKnockdown"
        or row["contrast"] != "Knockdown-Control"
        or row["freeze_version"] != "stage0c_paper3_analysis_contract_v1_2026-08-25"
        for row in design_rows
    )

    # Frozen reference identity.
    feature_rows = read_tsv(args.features)
    feature_by_tx = {row["Transcript"]: row for row in feature_rows}
    fasta = fasta_ids(args.fasta)
    violations["feature_duplicate_ids"] = len(feature_rows) - len(feature_by_tx)
    violations["fasta_duplicate_ids"] = len(fasta) - len(set(fasta))
    violations["feature_fasta_symmetric_difference"] = len(set(feature_by_tx) ^ set(fasta))
    if len(feature_rows) != 19_192 or len(fasta) != 19_192:
        violations["frozen_reference_count"] += 1

    # P-site freeze identity without reading any calibrated P-site output.
    offset_rows = read_tsv(args.offsets)
    observed_offsets: dict[int, int] = {}
    observed_unresolved: set[int] = set()
    for row in offset_rows:
        length = int(row["read_length"])
        value = row["frozen_offset_distance"]
        if value == "NA":
            observed_unresolved.add(length)
        else:
            observed_offsets[length] = int(value)
        if row["freeze_version"] != "stage0c_psite_v1_2026-08-24":
            violations["psite_freeze_version"] += 1
    violations["psite_offset_contract"] = int(observed_offsets != EXPECTED_OFFSETS)
    violations["psite_unresolved_contract"] = int(observed_unresolved != EXPECTED_UNRESOLVED)

    # DeepTMHMM input eligibility and normalized topology.
    input_rows = read_tsv(args.raw_input)
    input_by_tx = {row["transcript"]: row for row in input_rows}
    eligible_tx = {tx for tx, row in input_by_tx.items() if row["deeptmhmm_input_eligible"] == "YES"}
    ineligible_tx = set(input_by_tx) - eligible_tx
    violations["raw_input_duplicate_transcripts"] = len(input_rows) - len(input_by_tx)
    violations["raw_input_total"] = int(len(input_rows) != 19_192)
    violations["raw_input_eligible"] = int(len(eligible_tx) != 18_845)
    violations["raw_input_ineligible"] = int(len(ineligible_tx) != 347)

    raw_rows = read_tsv(args.raw_topology)
    raw_by_identity: dict[tuple[str, int], dict[str, str]] = {}
    topology_intervals: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    primary_intervals: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for row in raw_rows:
        identity = (row["transcript"], int(row["segment_index_raw_order"]))
        if identity in raw_by_identity:
            violations["raw_duplicate_segment_identity"] += 1
        raw_by_identity[identity] = row
        label = row["raw_deeptmhmm_label"]
        start = int(row["cds_nt_start0"])
        end = int(row["cds_nt_end0_exclusive"])
        if label in {"TMhelix", "signal", "Beta sheet"}:
            topology_intervals[row["transcript"]].append((start, end, label))
        if row["roundtrip_match"] != "YES":
            violations["raw_roundtrip"] += 1
    violations["raw_topology_rows"] = int(len(raw_rows) != 61_861)

    # Primary anchor contract and exact raw-source mapping.
    anchors = read_tsv(args.anchors)
    anchors.sort(key=lambda row: (row["transcript"], int(row["TMD_index"])))
    anchor_by_key = {(row["transcript"], int(row["TMD_index"])): row for row in anchors}
    violations["anchor_rows"] = int(len(anchors) != 19_784)
    violations["anchor_duplicate_keys"] = len(anchors) - len(anchor_by_key)
    for row in anchors:
        tx = row["transcript"]
        tmd_index = int(row["TMD_index"])
        if row["raw_deeptmhmm_label"] != "TMhelix":
            violations["anchor_label_policy"] += 1
        expected_versions = {
            "deeptmhmm_source_freeze_version": "stage0c_deeptmhmm_v1_2026-08-24",
            "raw_annotation_freeze_version": "stage0c_deeptmhmm_raw_v1_2026-08-24",
            "tmd_label_policy_freeze_version": "stage0c_tmd_label_policy_v1_2026-08-24",
            "tmd_anchor_policy_freeze_version": "stage0c_tmd_anchor_policy_v1_2026-08-24",
            "anchor_table_freeze_version": "stage0c_primary_tmd_anchors_v1_2026-08-24",
        }
        for field, value in expected_versions.items():
            if row[field] != value:
                violations[f"anchor_version_{field}"] += 1
        aa0 = int(row["aa_start0"])
        aa1 = int(row["aa_end0_exclusive"])
        cds0 = int(row["cds_nt_start0"])
        cds1 = int(row["cds_nt_end0_exclusive"])
        tx0 = int(row["transcript_nt_start0"])
        tx1 = int(row["transcript_nt_end0_exclusive"])
        cds_start = int(row["cds_start"])
        if not (
            cds0 == 3 * aa0
            and cds1 == 3 * aa1
            and tx0 == cds_start + cds0
            and tx1 == cds_start + cds1
            and int(row["primary_anchor_coordinate"]) == tx0
            and int(row["primary_anchor_cds_coordinate"]) == cds0
        ):
            violations["anchor_coordinate_contract"] += 1
        raw = raw_by_identity.get((tx, int(row["segment_index_raw_order"])))
        if raw is None:
            violations["anchor_missing_raw_segment"] += 1
        else:
            comparable = (
                raw["raw_deeptmhmm_label"] == "TMhelix"
                and int(raw["normalized_aa_start0"]) == aa0
                and int(raw["normalized_aa_end0_exclusive"]) == aa1
                and int(raw["cds_nt_start0"]) == cds0
                and int(raw["cds_nt_end0_exclusive"]) == cds1
                and int(raw["transcript_nt_start0"]) == tx0
                and int(raw["transcript_nt_end0_exclusive"]) == tx1
                and raw["protein_sha256"] == row["protein_sha256"]
            )
            if not comparable:
                violations["anchor_raw_coordinate_mismatch"] += 1
        if tx not in eligible_tx:
            violations["ineligible_transcript_in_primary_universe"] += 1
        primary_intervals[tx].append((cds0, cds1, tmd_index))
    raw_primary_count = sum(row["raw_deeptmhmm_label"] == "TMhelix" for row in raw_rows)
    violations["raw_primary_segment_count"] = int(raw_primary_count != 19_784)

    # Rebuild all 11 structural-validity masks and the primary-valid universe.
    masks_by_key: dict[tuple[str, int], str] = {}
    primary_valid_keys: set[tuple[str, int]] = set()
    for row in anchors:
        tx = row["transcript"]
        key = (tx, int(row["TMD_index"]))
        anchor_cds = int(row["primary_anchor_cds_coordinate"])
        cds_length = int(feature_by_tx[tx]["CDS_Length"])
        bits: list[str] = []
        for center in CENTERS:
            start = anchor_cds + center - 15
            end = anchor_cds + center + 16
            valid = start >= 0 and end <= cds_length
            if valid:
                valid = not any(
                    other_index != key[1] and overlaps(start, end, other_start, other_end)
                    for other_start, other_end, other_index in primary_intervals[tx]
                )
            bits.append("1" if valid else "0")
        mask = "".join(bits)
        masks_by_key[key] = mask
        if mask[PRIMARY_INDEX] == "1":
            primary_valid_keys.add(key)
    violations["primary_valid_count"] = int(len(primary_valid_keys) != 15_597)

    # Pseudo production table and frozen-policy invariants.
    pairs = read_tsv(args.pseudo_pairs)
    pair_by_key = {(row["transcript"], int(row["TMD_index"])): row for row in pairs}
    violations["pseudo_rows"] = int(len(pairs) != 19_784)
    violations["pseudo_duplicate_keys"] = len(pairs) - len(pair_by_key)
    violations["pseudo_real_anchor_key_set"] = len(set(pair_by_key) ^ set(anchor_by_key))
    matched_keys: set[tuple[str, int]] = set()
    pseudo_reuse: set[tuple[str, int]] = set()
    status_counts = Counter(row["match_status"] for row in pairs)
    expected_status = {
        "MATCHED": 9_829,
        "PRIMARY_WINDOW_INVALID": 4_187,
        "NO_PSEUDO_CANDIDATE": 5_736,
        "ASSIGNMENT_CONFLICT": 32,
    }
    violations["pseudo_status_counts"] = int(dict(status_counts) != expected_status)
    for key, row in pair_by_key.items():
        anchor = anchor_by_key[key]
        if (
            int(row["real_anchor_coordinate"]) != int(anchor["primary_anchor_coordinate"])
            or int(row["real_anchor_cds0"]) != int(anchor["primary_anchor_cds_coordinate"])
            or int(row["real_TMD_cds_start0"]) != int(anchor["cds_nt_start0"])
            or int(row["real_TMD_cds_end0_exclusive"]) != int(anchor["cds_nt_end0_exclusive"])
        ):
            violations["pseudo_real_anchor_reference"] += 1
        recomputed_mask = masks_by_key[key]
        if row["real_scan_validity_mask"] != recomputed_mask:
            violations["pseudo_real_mask_rebuild"] += 1
        recomputed_primary = "YES" if key in primary_valid_keys else "NO"
        if row["real_primary_window_valid"] != recomputed_primary:
            violations["pseudo_real_primary_rebuild"] += 1
        if row["match_status"] != "MATCHED":
            continue
        matched_keys.add(key)
        if key not in primary_valid_keys or row["pseudo_primary_window_valid"] != "YES":
            violations["matched_primary_validity"] += 1
        if row["pseudo_transcript"] != key[0]:
            violations["matched_same_transcript"] += 1
        if row["pseudo_scan_validity_mask"] != recomputed_mask:
            violations["matched_validity_mask"] += 1
        pseudo_start = int(row["pseudo_segment_cds_start0"])
        pseudo_end = int(row["pseudo_segment_cds_end0_exclusive"])
        real_length = int(row["real_TMD_length_nt"])
        if pseudo_end - pseudo_start != real_length or int(row["pseudo_segment_length_nt"]) != real_length:
            violations["matched_tmd_length"] += 1
        if pseudo_start % 3 or pseudo_end % 3:
            violations["matched_codon_alignment"] += 1
        if any(overlaps(pseudo_start, pseudo_end, start, end) for start, end, _ in topology_intervals[key[0]]):
            violations["matched_topology_overlap"] += 1
        real_starts = {start for start, _end, _idx in primary_intervals[key[0]]}
        if pseudo_start in real_starts:
            violations["matched_real_start_collision"] += 1
        pseudo_key = (key[0], int(row["pseudo_anchor_coordinate"]))
        if pseudo_key in pseudo_reuse:
            violations["matched_pseudo_reuse"] += 1
        pseudo_reuse.add(pseudo_key)
        if int(row["positional_mismatch_nt"]) != abs(pseudo_start - int(row["real_anchor_cds0"])):
            violations["matched_positional_metric"] += 1
    violations["matched_count"] = int(len(matched_keys) != 9_829)

    # Coordinate-null v2 universe/subset and production reconciliation.
    null_rows = read_tsv(args.null_universe)
    null_by_key = {(row["transcript"], int(row["TMD_index"])): row for row in null_rows}
    violations["null_v2_audit_rows"] = int(len(null_rows) != 13_480)
    violations["null_v2_audit_duplicate_keys"] = len(null_rows) - len(null_by_key)
    violations["null_v2_audit_not_primary_valid"] = len(set(null_by_key) - primary_valid_keys)
    null_eligible_keys = {
        key for key, row in null_by_key.items()
        if row["v2_coordinate_null_status"] == "COORDINATE_NULL_V2_ELIGIBLE"
    }
    null_v2_ineligible_keys = set(null_by_key) - null_eligible_keys
    null_all_ineligible_keys = primary_valid_keys - null_eligible_keys
    violations["null_v2_eligible_count"] = int(len(null_eligible_keys) != 13_060)
    violations["null_v2_newly_ineligible_count"] = int(len(null_v2_ineligible_keys) != 420)
    violations["null_v2_total_ineligible_count"] = int(len(null_all_ineligible_keys) != 2_537)
    violations["null_v2_partition"] = len((null_eligible_keys | null_all_ineligible_keys) ^ primary_valid_keys)
    violations["null_v2_policy_versions"] = sum(
        row["coordinate_null_policy_v2_freeze_version"]
        != "stage0c_coordinate_null_policy_v2_2026-08-25"
        for row in null_rows
    )
    null_columns = read_tsv(args.null_columns)
    null_column_keys = {(row["transcript"], int(row["TMD_index"])) for row in null_columns}
    violations["null_column_rows"] = int(len(null_columns) != 13_060)
    violations["null_column_duplicate_keys"] = len(null_columns) - len(null_column_keys)
    violations["null_column_key_set"] = len(null_column_keys ^ null_eligible_keys)
    violations["null_column_policy_versions"] = sum(
        row["coordinate_null_policy_freeze_version"]
        != "stage0c_coordinate_null_policy_v2_2026-08-25"
        for row in null_columns
    )
    for row in null_columns:
        key = (row["transcript"], int(row["TMD_index"]))
        anchor = anchor_by_key[key]
        if (
            int(row["real_anchor_coordinate"]) != int(anchor["primary_anchor_coordinate"])
            or int(row["real_anchor_cds0"]) != int(anchor["primary_anchor_cds_coordinate"])
        ):
            violations["null_column_anchor_reference"] += 1

    null_validation = json.loads(args.null_validation.read_text(encoding="utf-8"))
    null_counts = null_validation["counts"]
    null_storage = null_validation["storage"]
    violations["null_validation_realizations"] = int(null_counts["realizations"] != 10_000)
    violations["null_validation_anchors"] = int(null_counts["anchors_per_realization"] != 13_060)
    violations["null_validation_transcripts"] = int(null_counts["transcripts"] != 3_997)
    violations["null_validation_assignments"] = int(null_counts["total_null_assignments"] != 130_600_000)
    violations["null_validation_reproduction"] = int(
        null_validation["independent_replay"]["coordinate_assignment_mismatches"] != 0
        or null_validation["independent_replay"]["status"] != "PASS"
    )
    violations["null_validation_matrix_hash"] = int(
        null_storage["matrix_file_sha256"]
        != "7983d6dc1c382e57521e676669891a03271a4f1ffda0dda8d9ecf101e01071b3"
    )
    violations["null_validation_content_hash"] = int(
        null_storage["canonical_logical_C_order_int32_content_sha256"]
        != "c98f4b899a3f01398000b246b04c22002560b4b6f787f87801673b61e287aaf4"
    )
    violations["null_validation_v2_versions"] = int(
        null_validation["policy_freeze_version"]
        != "stage0c_coordinate_null_policy_v2_2026-08-25"
        or null_validation["production_freeze_version"]
        != "stage0c_coordinate_null_realizations_v2_2026-08-25"
    )
    violations["null_validation_any_violation"] = sum(null_validation["violations"].values())

    # Synthetic/annotation-only P-site-to-anchor algebra.
    for row in anchors:
        anchor_tx = int(row["primary_anchor_coordinate"])
        anchor_cds = int(row["primary_anchor_cds_coordinate"])
        cds_start = int(row["cds_start"])
        if anchor_tx != cds_start + anchor_cds:
            violations["psite_anchor_base_equivalence"] += 1
        for d_nt in (-45, 59, 60, 75, 90, 91, 135):
            psite = anchor_tx + d_nt
            direct = psite - anchor_tx
            relative = (psite - cds_start) - anchor_cds
            if direct != d_nt or relative != d_nt:
                violations["psite_anchor_synthetic_algebra"] += 1
            primary_member = 60 <= direct < 91
            if primary_member != (d_nt in {60, 75, 90}):
                violations["primary_membership_synthetic"] += 1
        for center in CENTERS:
            for d_nt, expected in ((center - 15, True), (center + 15, True), (center + 16, False)):
                observed = center - 15 <= d_nt < center + 16
                if observed != expected:
                    violations["scan_membership_synthetic"] += 1

    fail_if(violations)

    universe_rows: list[dict[str, Any]] = []
    for row in anchors:
        key = (row["transcript"], int(row["TMD_index"]))
        pair = pair_by_key[key]
        primary_valid = key in primary_valid_keys
        in_pseudo = key in matched_keys
        in_null = key in null_eligible_keys
        if not primary_valid:
            null_reason = "PRIMARY_WINDOW_INVALID"
        elif in_null:
            null_reason = "IN_COORDINATE_NULL_V2_UNIVERSE"
        elif key in null_v2_ineligible_keys:
            null_reason = "COORDINATE_NULL_V2_INELIGIBLE"
        else:
            null_reason = "V1_CAPACITY_INELIGIBLE_NOT_IN_V2_INPUT"
        universe_rows.append({
            "transcript": key[0],
            "TMD_index": key[1],
            "n_primary_TMDs": int(row["n_primary_TMDs"]),
            "pass_class": row["pass_class"],
            "TMD_position_class": row["TMD_position_class"],
            "primary_anchor_coordinate": int(row["primary_anchor_coordinate"]),
            "primary_anchor_cds0": int(row["primary_anchor_cds_coordinate"]),
            "in_full_primary_tmd_universe": "YES",
            "primary_window_valid": "YES" if primary_valid else "NO",
            "primary_window_status": "PRIMARY_WINDOW_VALID" if primary_valid else "PRIMARY_WINDOW_INVALID",
            "in_strict_pseudo_matched_universe": "YES" if in_pseudo else "NO",
            "pseudo_match_status": pair["match_status"],
            "pseudo_exclusion_reason": "IN_STRICT_MATCHED_UNIVERSE" if in_pseudo else pair["match_status"],
            "in_coordinate_null_universe": "YES" if in_null else "NO",
            "coordinate_null_status": null_reason,
            "primary_anchor_table_freeze_version": "stage0c_primary_tmd_anchors_v1_2026-08-24",
            "metagene_window_policy_freeze_version": "stage0c_metagene_window_policy_v1_2026-08-24",
            "pseudo_anchor_table_freeze_version": "stage0c_pseudo_anchor_table_v1_2026-08-25",
            "coordinate_null_policy_freeze_version": "stage0c_coordinate_null_policy_v2_2026-08-25",
            "coordinate_null_production_freeze_version": "stage0c_coordinate_null_realizations_v2_2026-08-25",
            "analysis_contract_freeze_version": "stage0c_paper3_analysis_contract_v1_2026-08-25",
            "scope": "ANNOTATION_ONLY_NO_RIBOSEQ",
        })

    args.universe_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.universe_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(universe_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(universe_rows)

    report = {
        "scope": "annotation_only_plus_outcome_independent_design_no_riboseq_no_profiles_no_benchmarks",
        "analysis_contract": {
            "freeze_version": "stage0c_paper3_analysis_contract_v1_2026-08-25",
            "status": "PASS",
            "decision_table": "25/25 RESOLVED",
            "coordinate_null_authoritative_version": "V2",
            "coordinate_null_v1": "SUPERSEDED_BEFORE_UNBLINDING_NOT_USED_FOR_INFERENCE",
        },
        "counts": {
            "frozen_reference_transcripts": len(fasta),
            "deeptmhmm_input_eligible": len(eligible_tx),
            "deeptmhmm_technically_ineligible": len(ineligible_tx),
            "raw_topology_segments": len(raw_rows),
            "full_primary_tmd_universe": len(anchors),
            "primary_window_valid_universe": len(primary_valid_keys),
            "strict_pseudo_matched_universe": len(matched_keys),
            "coordinate_null_universe": len(null_eligible_keys),
            "coordinate_null_v2_newly_ineligible": len(null_v2_ineligible_keys),
            "coordinate_null_total_ineligible": len(null_all_ineligible_keys),
            "decision_rows_resolved": len(decision_rows),
            "design_libraries": len(design_rows),
        },
        "cross_freeze_checks": {
            "all_anchor_labels_follow_policy": True,
            "all_anchor_coordinates_follow_raw_contract": True,
            "all_anchors_exist_in_raw_source": True,
            "all_pseudo_rows_reference_authoritative_anchors": True,
            "all_matched_pseudo_rows_obey_frozen_annotation_geometry": True,
            "primary_window_valid_universe_rebuilt_exactly": True,
            "coordinate_null_is_exact_subset": True,
            "coordinate_null_v2_partition_reconciles": True,
            "technically_ineligible_transcripts_absent_from_primary_universe": True,
            "module04_anchor_coordinate_algebra": True,
            "decision_table_25_of_25_resolved": True,
            "design_matrix_frozen": True,
            "p3_08_topology_mask_frozen": True,
            "p3_14_pseudo_focal_exclusion_frozen": True,
            "coordinate_null_v2_only_authoritative": True,
            "coordinate_null_v1_not_used_for_inference": True,
        },
        "violations": dict(violations),
        "source_sha256": {
            "fasta": sha256(args.fasta),
            "features": sha256(args.features),
            "offsets": sha256(args.offsets),
            "raw_input": sha256(args.raw_input),
            "raw_topology": sha256(args.raw_topology),
            "anchors": sha256(args.anchors),
            "pseudo_pairs": sha256(args.pseudo_pairs),
            "null_universe": sha256(args.null_universe),
            "null_columns": sha256(args.null_columns),
            "null_validation": sha256(args.null_validation),
            "null_addendum": sha256(args.null_addendum),
            "decision_table": sha256(args.decision_table),
            "contract_freeze": sha256(args.contract_freeze),
            "design_matrix": sha256(args.design_matrix),
        },
        "psite_anchor_interface": {
            "d_nt": "psite_coordinate-primary_anchor_coordinate",
            "equivalent_cds_formula": "relative_psite-primary_anchor_cds0",
            "primary_membership": "60<=d_nt<91",
            "scan_membership": "center-15<=d_nt<center+16",
            "synthetic_and_annotation_checks": "PASS",
        },
        "riboseq_used": False,
        "biological_results_inspected": False,
        "benchmark_outcomes_inspected": False,
        "stage1_analysis_run": False,
        "tmd_blind": "PRESERVED",
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["counts"], sort_keys=True))
    print("INTEGRATION_INVARIANTS=PASS")


if __name__ == "__main__":
    main()
