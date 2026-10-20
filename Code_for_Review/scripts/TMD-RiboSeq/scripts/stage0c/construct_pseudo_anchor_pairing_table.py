#!/usr/bin/env python3
"""Construct the frozen annotation-only primary-TMD/pseudo-anchor pairing table.

This is a reference-table constructor, not Module 05. It reads no Ribo-seq,
sample, condition, benchmark, or biological-result data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_only_pseudo_anchor_feasibility import (  # noqa: E402
    CENTERS,
    MEMBRANE_FEATURE_LABELS,
    PRIMARY_CENTER,
    RealAnchor,
    build_mask,
    deterministic_min_cost_assignment,
    overlaps,
    read_tsv,
    sha256,
)


SCRIPT_VERSION = "construct_pseudo_anchor_pairing_table_v1_2026-08-25"
TABLE_FREEZE_VERSION = "stage0c_pseudo_anchor_table_v1_2026-08-25"
PSEUDO_POLICY_FREEZE = "stage0c_pseudo_anchor_policy_v1_2026-08-24"
WINDOW_POLICY_FREEZE = "stage0c_metagene_window_policy_v1_2026-08-24"
EXPECTED_ANCHOR_FREEZE = "stage0c_primary_tmd_anchors_v1_2026-08-24"
EXPECTED_LABEL_FREEZE = "stage0c_tmd_label_policy_v1_2026-08-24"
EXPECTED_ANCHORS = 19_784
EXPECTED_COUNTS = {
    "MATCHED": 9_829,
    "PRIMARY_WINDOW_INVALID": 4_187,
    "NO_PSEUDO_CANDIDATE": 5_736,
    "ASSIGNMENT_CONFLICT": 32,
}
NA = "NA"


OUTPUT_FIELDS = [
    "transcript",
    "protein_length",
    "protein_sha256",
    "TMD_index",
    "n_primary_TMDs",
    "pass_class",
    "TMD_position_class",
    "real_TMD_length_nt",
    "real_TMD_aa_start0",
    "real_TMD_aa_end0_exclusive",
    "real_TMD_cds_start0",
    "real_TMD_cds_end0_exclusive",
    "real_anchor_coordinate",
    "real_anchor_cds0",
    "CDS_length_nt",
    "real_normalized_CDS_position",
    "real_primary_window_valid",
    "real_scan_validity_mask",
    "pseudo_transcript",
    "pseudo_anchor_coordinate",
    "pseudo_anchor_cds0",
    "pseudo_segment_start0",
    "pseudo_segment_end0_exclusive",
    "pseudo_segment_cds_start0",
    "pseudo_segment_cds_end0_exclusive",
    "pseudo_segment_aa_start0",
    "pseudo_segment_aa_end0_exclusive",
    "pseudo_segment_length_nt",
    "pseudo_primary_window_valid",
    "pseudo_scan_validity_mask",
    "positional_mismatch_nt",
    "pseudo_candidate_count",
    "assignment_status",
    "match_status",
    "unmatched_reason",
    "assignment_objective",
    "assignment_tie_break",
    "pseudo_anchor_reuse_policy",
    "pseudo_segment_overlap_across_pairs_policy",
    "deeptmhmm_source_freeze_version",
    "raw_annotation_freeze_version",
    "tmd_label_policy_freeze_version",
    "tmd_anchor_policy_freeze_version",
    "anchor_table_freeze_version",
    "metagene_window_policy_freeze_version",
    "pseudo_anchor_policy_freeze_version",
    "pseudo_anchor_table_freeze_version",
    "frozen_reference_fasta_sha256",
    "human_features_sha256",
    "normalized_raw_annotation_sha256",
    "primary_anchor_table_sha256",
    "source_job_id",
    "source_batch_id",
    "deeptmhmm_app_version",
    "construction_script_version",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--representativeness-tsv", type=Path, required=True)
    parser.add_argument("--validation-json", type=Path, required=True)
    return parser.parse_args()


def mask_string(mask: Iterable[bool]) -> str:
    return "".join("1" if value else "0" for value in mask)


def linear_quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("Cannot summarize an empty vector.")
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    left = math.floor(position)
    right = math.ceil(position)
    if left == right:
        return ordered[left]
    fraction = position - left
    return ordered[left] * (1.0 - fraction) + ordered[right] * fraction


def numeric_summary(values: list[float]) -> dict[str, float | int]:
    if not values:
        raise ValueError("Numeric summary received no values.")
    return {
        "n": len(values),
        "mean": statistics.fmean(values),
        "population_sd": statistics.pstdev(values),
        "min": min(values),
        "p25": linear_quantile(values, 0.25),
        "median": linear_quantile(values, 0.50),
        "p75": linear_quantile(values, 0.75),
        "p90": linear_quantile(values, 0.90),
        "p95": linear_quantile(values, 0.95),
        "max": max(values),
    }


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def representativeness_rows(pair_rows: list[dict[str, Any]], cds_length_by_transcript: dict[str, int]) -> list[dict[str, Any]]:
    groups = {
        "ALL_REAL": pair_rows,
        "PRIMARY_WINDOW_ELIGIBLE": [row for row in pair_rows if row["real_primary_window_valid"] == "YES"],
        "MATCHED": [row for row in pair_rows if row["match_status"] == "MATCHED"],
        "UNMATCHED_ELIGIBLE": [
            row
            for row in pair_rows
            if row["real_primary_window_valid"] == "YES" and row["match_status"] != "MATCHED"
        ],
    }
    output: list[dict[str, Any]] = []
    categorical = ("TMD_position_class", "pass_class")
    numeric = {
        "TMD_index": lambda row: float(row["TMD_index"]),
        "real_TMD_length_nt": lambda row: float(row["real_TMD_length_nt"]),
        "real_anchor_cds0": lambda row: float(row["real_anchor_cds0"]),
        "normalized_CDS_position": lambda row: float(row["real_anchor_cds0"])
        / cds_length_by_transcript[str(row["transcript"])],
        "CDS_length_nt": lambda row: float(cds_length_by_transcript[str(row["transcript"])]),
        "n_primary_TMDs": lambda row: float(row["n_primary_TMDs"]),
    }
    for group_name, rows in groups.items():
        denominator = len(rows)
        output.append(
            {
                "record_type": "group_size",
                "group": group_name,
                "variable": "rows",
                "level_or_statistic": "n",
                "value": denominator,
                "denominator": denominator,
                "fraction": 1.0,
                "definition": "annotation-only real-TMD rows",
            }
        )
        for variable in categorical:
            counts = Counter(str(row[variable]) for row in rows)
            for level in sorted(counts):
                output.append(
                    {
                        "record_type": "categorical_distribution",
                        "group": group_name,
                        "variable": variable,
                        "level_or_statistic": level,
                        "value": counts[level],
                        "denominator": denominator,
                        "fraction": counts[level] / denominator,
                        "definition": "count and within-group fraction",
                    }
                )
        for variable, getter in numeric.items():
            summary = numeric_summary([getter(row) for row in rows])
            for statistic_name, value in summary.items():
                output.append(
                    {
                        "record_type": "numeric_distribution",
                        "group": group_name,
                        "variable": variable,
                        "level_or_statistic": statistic_name,
                        "value": value,
                        "denominator": denominator,
                        "fraction": NA,
                        "definition": "population SD; linear-interpolated quantiles",
                    }
                )

    eligible = groups["PRIMARY_WINDOW_ELIGIBLE"]
    for variable in ("TMD_position_class", "pass_class"):
        strata: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in eligible:
            strata[str(row[variable])].append(row)
        for level, rows in sorted(strata.items()):
            matched_count = sum(row["match_status"] == "MATCHED" for row in rows)
            output.append(
                {
                    "record_type": "matchability_by_stratum",
                    "group": "PRIMARY_WINDOW_ELIGIBLE",
                    "variable": variable,
                    "level_or_statistic": level,
                    "value": matched_count,
                    "denominator": len(rows),
                    "fraction": matched_count / len(rows),
                    "definition": "matched count and fraction among primary-window-eligible anchors",
                }
            )
    for variable in ("TMD_index", "n_primary_TMDs", "real_TMD_length_nt"):
        strata: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in eligible:
            strata[str(row[variable])].append(row)
        for level, rows in sorted(strata.items(), key=lambda item: int(item[0])):
            matched_count = sum(row["match_status"] == "MATCHED" for row in rows)
            output.append(
                {
                    "record_type": "matchability_by_stratum",
                    "group": "PRIMARY_WINDOW_ELIGIBLE",
                    "variable": variable,
                    "level_or_statistic": level,
                    "value": matched_count,
                    "denominator": len(rows),
                    "fraction": matched_count / len(rows),
                    "definition": "matched count and fraction among primary-window-eligible anchors",
                }
            )
    return output


def main() -> None:
    args = parse_args()
    for path in (args.anchors, args.topology, args.features):
        if not path.is_file():
            raise FileNotFoundError(path)

    anchor_rows = read_tsv(args.anchors)
    topology_rows = read_tsv(args.topology)
    feature_rows = read_tsv(args.features)
    if len(anchor_rows) != EXPECTED_ANCHORS:
        raise RuntimeError(f"Expected {EXPECTED_ANCHORS} anchors, observed {len(anchor_rows)}.")

    anchor_table_sha256 = sha256(args.anchors)
    topology_sha256 = sha256(args.topology)
    features_sha256 = sha256(args.features)
    features: dict[str, tuple[int, int]] = {}
    for row in feature_rows:
        transcript = row["Transcript"]
        if transcript in features:
            raise RuntimeError(f"Duplicate feature transcript: {transcript}")
        features[transcript] = (int(row["UTR5_Length"]), int(row["CDS_Length"]))

    membrane_by_transcript: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for row in topology_rows:
        label = row["raw_deeptmhmm_label"]
        if label in MEMBRANE_FEATURE_LABELS:
            membrane_by_transcript[row["transcript"]].append(
                (int(row["normalized_aa_start0"]), int(row["normalized_aa_end0_exclusive"]), label)
            )

    anchor_by_transcript: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in anchor_rows:
        if row["anchor_table_freeze_version"] != EXPECTED_ANCHOR_FREEZE:
            raise RuntimeError("Unexpected anchor-table freeze version.")
        if row["tmd_label_policy_freeze_version"] != EXPECTED_LABEL_FREEZE:
            raise RuntimeError("Unexpected TMD-label-policy freeze version.")
        anchor_by_transcript[row["transcript"]].append(row)

    primary_index = CENTERS.index(PRIMARY_CENTER)
    output_rows: list[dict[str, Any]] = []
    violations = Counter()
    assignment_total_cost = 0
    assignment_transcripts = 0

    for transcript_number, transcript in enumerate(sorted(anchor_by_transcript), start=1):
        if transcript_number % 500 == 0 or transcript_number == len(anchor_by_transcript):
            print(f"[pairing table] transcript {transcript_number}/{len(anchor_by_transcript)}", file=sys.stderr)
        rows = sorted(anchor_by_transcript[transcript], key=lambda row: int(row["TMD_index"]))
        if [int(row["TMD_index"]) for row in rows] != list(range(1, len(rows) + 1)):
            raise RuntimeError(f"Noncontiguous TMD indices for {transcript}.")
        cds_start_tx0, cds_length = features[transcript]
        protein_lengths = {int(row["protein_length"]) for row in rows}
        if len(protein_lengths) != 1:
            raise RuntimeError(f"Inconsistent protein lengths for {transcript}.")
        protein_length = next(iter(protein_lengths))
        accepted_intervals = [
            (int(row["cds_nt_start0"]), int(row["cds_nt_end0_exclusive"])) for row in rows
        ]
        accepted_starts = {start for start, _end in accepted_intervals}

        blocked = [0] * protein_length
        for start, end, _label in membrane_by_transcript.get(transcript, []):
            if start < 0 or end < start or end > protein_length:
                raise RuntimeError(f"Invalid topology interval for {transcript}.")
            for residue in range(start, end):
                blocked[residue] = 1
        blocked_prefix = [0]
        for value in blocked:
            blocked_prefix.append(blocked_prefix[-1] + value)

        reals: list[RealAnchor] = []
        candidates_by_real: list[list[int]] = []
        candidate_segments: list[dict[int, tuple[int, int]]] = []
        pseudo_mask_cache: dict[int, tuple[bool, ...]] = {}

        for index, row in enumerate(rows):
            aa_start = int(row["aa_start0"])
            aa_end = int(row["aa_end0_exclusive"])
            real_cds_start = int(row["cds_nt_start0"])
            real_cds_end = int(row["cds_nt_end0_exclusive"])
            if real_cds_start != 3 * aa_start or real_cds_end != 3 * aa_end:
                raise RuntimeError(f"Real codon invariant failed for {transcript} TMD {index + 1}.")
            real_mask = build_mask(real_cds_start, cds_length, accepted_intervals, index)
            real = RealAnchor(
                transcript=transcript,
                tmd_index=index + 1,
                aa_start=aa_start,
                aa_end=aa_end,
                cds_start=real_cds_start,
                cds_end=real_cds_end,
                protein_length=protein_length,
                real_mask=real_mask,
            )
            reals.append(real)
            if not real_mask[primary_index]:
                candidates_by_real.append([])
                candidate_segments.append({})
                continue
            length_aa = aa_end - aa_start
            candidates: list[int] = []
            segments: dict[int, tuple[int, int]] = {}
            for pseudo_aa_start in range(protein_length - length_aa + 1):
                pseudo_aa_end = pseudo_aa_start + length_aa
                if blocked_prefix[pseudo_aa_end] - blocked_prefix[pseudo_aa_start] != 0:
                    continue
                pseudo_cds_start = 3 * pseudo_aa_start
                if pseudo_cds_start in accepted_starts:
                    continue
                pseudo_mask = pseudo_mask_cache.get(pseudo_cds_start)
                if pseudo_mask is None:
                    pseudo_mask = build_mask(pseudo_cds_start, cds_length, accepted_intervals, None)
                    pseudo_mask_cache[pseudo_cds_start] = pseudo_mask
                if pseudo_mask != real_mask:
                    continue
                candidates.append(pseudo_cds_start)
                segments[pseudo_cds_start] = (pseudo_aa_start, pseudo_aa_end)
            candidates.sort()
            candidates_by_real.append(candidates)
            candidate_segments.append(segments)

        assignable_indices = [index for index, values in enumerate(candidates_by_real) if values]
        if assignable_indices:
            assignment_transcripts += 1
            partial = deterministic_min_cost_assignment(
                [reals[index] for index in assignable_indices],
                [candidates_by_real[index] for index in assignable_indices],
            )
            assignments = {assignable_indices[index]: coordinate for index, coordinate in partial.items()}
        else:
            assignments = {}

        selected_transcript_coordinates: set[int] = set()
        for index, (row, real) in enumerate(zip(rows, reals)):
            candidate_count = len(candidates_by_real[index])
            if not real.real_mask[primary_index]:
                match_status = "PRIMARY_WINDOW_INVALID"
                assignment_status = "NOT_ELIGIBLE_PRIMARY_WINDOW"
            elif candidate_count == 0:
                match_status = "NO_PSEUDO_CANDIDATE"
                assignment_status = "NO_ELIGIBLE_CANDIDATE"
            elif index not in assignments:
                match_status = "ASSIGNMENT_CONFLICT"
                assignment_status = "UNASSIGNED_ONE_TO_ONE_CONFLICT"
            else:
                match_status = "MATCHED"
                assignment_status = "ASSIGNED_MAX_CARDINALITY_MIN_TOTAL_MISMATCH"

            base: dict[str, Any] = {
                "transcript": transcript,
                "protein_length": row["protein_length"],
                "protein_sha256": row["protein_sha256"],
                "TMD_index": row["TMD_index"],
                "n_primary_TMDs": row["n_primary_TMDs"],
                "pass_class": row["pass_class"],
                "TMD_position_class": row["TMD_position_class"],
                "real_TMD_length_nt": real.cds_end - real.cds_start,
                "real_TMD_aa_start0": real.aa_start,
                "real_TMD_aa_end0_exclusive": real.aa_end,
                "real_TMD_cds_start0": real.cds_start,
                "real_TMD_cds_end0_exclusive": real.cds_end,
                "real_anchor_coordinate": row["primary_anchor_coordinate"],
                "real_anchor_cds0": real.cds_start,
                "CDS_length_nt": cds_length,
                "real_normalized_CDS_position": real.cds_start / cds_length,
                "real_primary_window_valid": "YES" if real.real_mask[primary_index] else "NO",
                "real_scan_validity_mask": mask_string(real.real_mask),
                "pseudo_candidate_count": candidate_count,
                "assignment_status": assignment_status,
                "match_status": match_status,
                "unmatched_reason": NA if match_status == "MATCHED" else match_status,
                "assignment_objective": "MAX_CARDINALITY_THEN_MIN_TOTAL_ABS_CDS_NT_DISPLACEMENT",
                "assignment_tie_break": "STABLE_SSP_SORTED_TMD_INDEX_AND_PSEUDO_COORDINATE_FIRST_DISCOVERED_EQUAL_COST_PATH",
                "pseudo_anchor_reuse_policy": "PROHIBITED",
                "pseudo_segment_overlap_across_pairs_policy": "ALLOWED_AND_AUDITED",
                "deeptmhmm_source_freeze_version": row["deeptmhmm_source_freeze_version"],
                "raw_annotation_freeze_version": row["raw_annotation_freeze_version"],
                "tmd_label_policy_freeze_version": row["tmd_label_policy_freeze_version"],
                "tmd_anchor_policy_freeze_version": row["tmd_anchor_policy_freeze_version"],
                "anchor_table_freeze_version": row["anchor_table_freeze_version"],
                "metagene_window_policy_freeze_version": WINDOW_POLICY_FREEZE,
                "pseudo_anchor_policy_freeze_version": PSEUDO_POLICY_FREEZE,
                "pseudo_anchor_table_freeze_version": TABLE_FREEZE_VERSION,
                "frozen_reference_fasta_sha256": row["frozen_reference_fasta_sha256"],
                "human_features_sha256": row["human_features_sha256"],
                "normalized_raw_annotation_sha256": row["normalized_raw_annotation_sha256"],
                "primary_anchor_table_sha256": anchor_table_sha256,
                "source_job_id": row["source_job_id"],
                "source_batch_id": row["source_batch_id"],
                "deeptmhmm_app_version": row["deeptmhmm_app_version"],
                "construction_script_version": SCRIPT_VERSION,
            }
            pseudo_fields = {
                "pseudo_transcript": NA,
                "pseudo_anchor_coordinate": NA,
                "pseudo_anchor_cds0": NA,
                "pseudo_segment_start0": NA,
                "pseudo_segment_end0_exclusive": NA,
                "pseudo_segment_cds_start0": NA,
                "pseudo_segment_cds_end0_exclusive": NA,
                "pseudo_segment_aa_start0": NA,
                "pseudo_segment_aa_end0_exclusive": NA,
                "pseudo_segment_length_nt": NA,
                "pseudo_primary_window_valid": NA,
                "pseudo_scan_validity_mask": NA,
                "positional_mismatch_nt": NA,
            }
            if match_status == "MATCHED":
                pseudo_cds_start = assignments[index]
                pseudo_aa_start, pseudo_aa_end = candidate_segments[index][pseudo_cds_start]
                pseudo_cds_end = 3 * pseudo_aa_end
                pseudo_tx_start = cds_start_tx0 + pseudo_cds_start
                pseudo_tx_end = cds_start_tx0 + pseudo_cds_end
                pseudo_mask = build_mask(pseudo_cds_start, cds_length, accepted_intervals, None)
                mismatch = abs(pseudo_cds_start - real.cds_start)
                assignment_total_cost += mismatch
                pseudo_fields = {
                    "pseudo_transcript": transcript,
                    "pseudo_anchor_coordinate": pseudo_tx_start,
                    "pseudo_anchor_cds0": pseudo_cds_start,
                    "pseudo_segment_start0": pseudo_tx_start,
                    "pseudo_segment_end0_exclusive": pseudo_tx_end,
                    "pseudo_segment_cds_start0": pseudo_cds_start,
                    "pseudo_segment_cds_end0_exclusive": pseudo_cds_end,
                    "pseudo_segment_aa_start0": pseudo_aa_start,
                    "pseudo_segment_aa_end0_exclusive": pseudo_aa_end,
                    "pseudo_segment_length_nt": pseudo_cds_end - pseudo_cds_start,
                    "pseudo_primary_window_valid": "YES" if pseudo_mask[primary_index] else "NO",
                    "pseudo_scan_validity_mask": mask_string(pseudo_mask),
                    "positional_mismatch_nt": mismatch,
                }
                if pseudo_tx_start in selected_transcript_coordinates:
                    violations["pseudo_anchor_reuse"] += 1
                selected_transcript_coordinates.add(pseudo_tx_start)
                if pseudo_cds_end - pseudo_cds_start != real.cds_end - real.cds_start:
                    violations["tmd_length_match"] += 1
                if pseudo_cds_start % 3 or pseudo_cds_end % 3:
                    violations["codon_alignment"] += 1
                if blocked_prefix[pseudo_aa_end] - blocked_prefix[pseudo_aa_start] != 0:
                    violations["topology_overlap"] += 1
                if pseudo_cds_start in accepted_starts:
                    violations["real_tmd_start_collision"] += 1
                if pseudo_mask != real.real_mask:
                    violations["validity_mask"] += 1
                if not real.real_mask[primary_index] or not pseudo_mask[primary_index]:
                    violations["primary_validity"] += 1
                if mismatch != abs(pseudo_cds_start - real.cds_start):
                    violations["positional_metric"] += 1
            base.update(pseudo_fields)
            output_rows.append(base)

    output_rows.sort(key=lambda row: (str(row["transcript"]), int(row["TMD_index"])))
    status_counts = Counter(str(row["match_status"]) for row in output_rows)
    if len(output_rows) != EXPECTED_ANCHORS:
        raise RuntimeError("Production row count failed.")
    if dict(status_counts) != EXPECTED_COUNTS:
        raise RuntimeError(f"Assignment reproduction failed: {dict(status_counts)}")
    if sum(violations.values()):
        raise RuntimeError(f"Pairing invariants failed: {dict(violations)}")
    if len({(row["transcript"], row["TMD_index"]) for row in output_rows}) != EXPECTED_ANCHORS:
        raise RuntimeError("Duplicate or missing real transcript/TMD_index keys.")

    write_tsv(args.output_table, OUTPUT_FIELDS, output_rows)
    cds_lengths = {transcript: value[1] for transcript, value in features.items()}
    rep_rows = representativeness_rows(output_rows, cds_lengths)
    rep_fields = [
        "record_type", "group", "variable", "level_or_statistic",
        "value", "denominator", "fraction", "definition",
    ]
    write_tsv(args.representativeness_tsv, rep_fields, rep_rows)

    validation = {
        "script_version": SCRIPT_VERSION,
        "table_freeze_version": TABLE_FREEZE_VERSION,
        "scope": "annotation_only_no_riboseq",
        "inputs": {
            "anchor_table": {"path": str(args.anchors), "sha256": anchor_table_sha256},
            "topology": {"path": str(args.topology), "sha256": topology_sha256},
            "features": {"path": str(args.features), "sha256": features_sha256},
        },
        "outputs": {
            "pairing_table": {"path": str(args.output_table), "sha256": sha256(args.output_table)},
            "representativeness": {
                "path": str(args.representativeness_tsv),
                "sha256": sha256(args.representativeness_tsv),
            },
        },
        "counts": {
            "total_real_tmd_rows": len(output_rows),
            "primary_window_eligible": sum(row["real_primary_window_valid"] == "YES" for row in output_rows),
            **dict(status_counts),
            "with_at_least_one_pseudo_candidate": sum(int(row["pseudo_candidate_count"]) > 0 for row in output_rows),
            "assignment_transcripts": assignment_transcripts,
            "total_positional_mismatch_nt": assignment_total_cost,
        },
        "violations": {
            "same_transcript": 0,
            "tmd_length_match": violations["tmd_length_match"],
            "codon_alignment": violations["codon_alignment"],
            "validity_mask": violations["validity_mask"],
            "topology_overlap": violations["topology_overlap"],
            "real_tmd_start_collision": violations["real_tmd_start_collision"],
            "primary_validity": violations["primary_validity"],
            "pseudo_anchor_reuse": violations["pseudo_anchor_reuse"],
            "positional_metric": violations["positional_metric"],
        },
        "invariants": {
            "count_conservation": sum(status_counts.values()) == EXPECTED_ANCHORS,
            "assignment_reproduction": dict(status_counts) == EXPECTED_COUNTS,
            "all_pairing_invariants_zero": sum(violations.values()) == 0,
            "riboseq_used": False,
        },
    }
    args.validation_json.parent.mkdir(parents=True, exist_ok=True)
    with args.validation_json.open("w", encoding="utf-8") as handle:
        json.dump(validation, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
