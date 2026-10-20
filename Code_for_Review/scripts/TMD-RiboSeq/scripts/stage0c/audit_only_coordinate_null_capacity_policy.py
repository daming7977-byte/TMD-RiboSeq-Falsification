#!/usr/bin/env python3
"""Independently audit and classify GSE133111 coordinate-null capacity.

Annotation-only: no Ribo-seq, condition data, null draws, or Module 05 output.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable


SCRIPT_VERSION = "audit_only_coordinate_null_capacity_policy_v1_2026-08-25"
EXPECTED_PRIMARY_VALID_ANCHORS = 15_597
EXPECTED_PRIMARY_VALID_TRANSCRIPTS = 4_465
EXPECTED_ELIGIBLE_ANCHORS = 13_480
EXPECTED_ELIGIBLE_TRANSCRIPTS = 4_103
EXPECTED_INELIGIBLE_ANCHORS = 2_117
EXPECTED_INELIGIBLE_TRANSCRIPTS = 362
TOPOLOGY_LABELS = frozenset({"TMhelix", "signal", "Beta sheet"})
PRIMARY_START = 60
PRIMARY_END_EXCLUSIVE = 91
FREEZE_VERSION = "stage0c_coordinate_null_policy_v1_2026-08-25"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--prior-feasibility", type=Path, required=True)
    parser.add_argument("--universe-tsv", type=Path, required=True)
    parser.add_argument("--representativeness-tsv", type=Path, required=True)
    parser.add_argument("--validation-json", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def overlaps(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return start_a < end_b and start_b < end_a


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return total + current_end - current_start


def linear_quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    left = math.floor(position)
    right = math.ceil(position)
    if left == right:
        return ordered[left]
    fraction = position - left
    return ordered[left] * (1.0 - fraction) + ordered[right] * fraction


def numeric_summary(values: list[float]) -> dict[str, float | int]:
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


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summary_rows(
    anchor_groups: dict[str, list[dict[str, Any]]],
    transcript_groups: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []

    def add_groups(
        analysis_level: str,
        groups: dict[str, list[dict[str, Any]]],
        categorical: list[str],
        numeric: dict[str, Callable[[dict[str, Any]], float]],
    ) -> None:
        for group, rows in groups.items():
            denominator = len(rows)
            output.append({
                "analysis_level": analysis_level,
                "record_type": "group_size",
                "group": group,
                "variable": "rows",
                "level_or_statistic": "n",
                "value": denominator,
                "denominator": denominator,
                "fraction": 1.0,
                "definition": "annotation-only group size",
            })
            for variable in categorical:
                counts = Counter(str(row[variable]) for row in rows)
                for level in sorted(counts):
                    output.append({
                        "analysis_level": analysis_level,
                        "record_type": "categorical_distribution",
                        "group": group,
                        "variable": variable,
                        "level_or_statistic": level,
                        "value": counts[level],
                        "denominator": denominator,
                        "fraction": counts[level] / denominator,
                        "definition": "count and within-group fraction",
                    })
            for variable, getter in numeric.items():
                summary = numeric_summary([getter(row) for row in rows])
                for statistic_name, value in summary.items():
                    output.append({
                        "analysis_level": analysis_level,
                        "record_type": "numeric_distribution",
                        "group": group,
                        "variable": variable,
                        "level_or_statistic": statistic_name,
                        "value": value,
                        "denominator": denominator,
                        "fraction": "NA",
                        "definition": "population SD; linear-interpolated quantiles",
                    })

    add_groups(
        "ANCHOR",
        anchor_groups,
        ["TMD_position_class", "pass_class"],
        {
            "TMD_index": lambda row: float(row["TMD_index"]),
            "n_primary_TMDs": lambda row: float(row["n_primary_TMDs"]),
            "real_TMD_length_nt": lambda row: float(row["real_TMD_length_nt"]),
            "CDS_length_nt": lambda row: float(row["CDS_length_nt"]),
            "real_anchor_cds0": lambda row: float(row["real_anchor_cds0"]),
            "normalized_CDS_position": lambda row: float(row["normalized_CDS_position"]),
            "distance_to_CDS_start_nt": lambda row: float(row["distance_to_CDS_start_nt"]),
            "distance_to_CDS_end_nt": lambda row: float(row["distance_to_CDS_end_nt"]),
            "candidate_count": lambda row: float(row["candidate_count"]),
            "k_per_transcript": lambda row: float(row["k_per_transcript"]),
        },
    )
    add_groups(
        "TRANSCRIPT",
        transcript_groups,
        [],
        {
            "k_analyzed_anchors": lambda row: float(row["k_analyzed_anchors"]),
            "CDS_length_nt": lambda row: float(row["CDS_length_nt"]),
            "candidate_count": lambda row: float(row["candidate_count"]),
            "n_TMhelix_segments": lambda row: float(row["n_TMhelix_segments"]),
            "n_signal_segments": lambda row: float(row["n_signal_segments"]),
            "n_Beta_sheet_segments": lambda row: float(row["n_Beta_sheet_segments"]),
            "n_disallowed_topology_segments": lambda row: float(row["n_disallowed_topology_segments"]),
            "disallowed_topology_union_nt": lambda row: float(row["disallowed_topology_union_nt"]),
            "disallowed_topology_CDS_fraction": lambda row: float(row["disallowed_topology_CDS_fraction"]),
        },
    )
    return output


def main() -> None:
    args = parse_args()
    for path in (args.pairs, args.topology, args.features, args.prior_feasibility):
        if not path.is_file():
            raise FileNotFoundError(path)

    pair_rows = read_tsv(args.pairs)
    topology_rows = read_tsv(args.topology)
    feature_rows = read_tsv(args.features)
    prior_rows = read_tsv(args.prior_feasibility)
    features = {row["Transcript"]: int(row["CDS_Length"]) for row in feature_rows}
    if len(features) != len(feature_rows):
        raise RuntimeError("Duplicate transcript in human_features.tsv")

    primary_valid = [row for row in pair_rows if row["real_primary_window_valid"] == "YES"]
    k_by_transcript = Counter(row["transcript"] for row in primary_valid)
    if len(primary_valid) != EXPECTED_PRIMARY_VALID_ANCHORS or len(k_by_transcript) != EXPECTED_PRIMARY_VALID_TRANSCRIPTS:
        raise RuntimeError("Primary-window-valid source universe does not match expected counts.")

    topology_by_transcript: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    topology_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in topology_rows:
        label = row["raw_deeptmhmm_label"]
        if label in TOPOLOGY_LABELS:
            transcript = row["transcript"]
            start = int(row["cds_nt_start0"])
            end = int(row["cds_nt_end0_exclusive"])
            topology_by_transcript[transcript].append((start, end, label))
            topology_counts[transcript][label] += 1

    transcript_rows: list[dict[str, Any]] = []
    capacity_by_transcript: dict[str, dict[str, Any]] = {}
    for transcript in sorted(k_by_transcript):
        cds_length = features[transcript]
        intervals = sorted(topology_by_transcript.get(transcript, []))
        candidates: list[int] = []
        for coordinate in range(0, cds_length, 3):
            primary_start = coordinate + PRIMARY_START
            primary_end = coordinate + PRIMARY_END_EXCLUSIVE
            if primary_end > cds_length:
                continue
            if any(start <= coordinate < end for start, end, _label in intervals):
                continue
            if any(overlaps(primary_start, primary_end, start, end) for start, end, _label in intervals):
                continue
            candidates.append(coordinate)
        k = k_by_transcript[transcript]
        eligible = len(candidates) >= k
        counts = topology_counts[transcript]
        covered_nt = union_length([(start, end) for start, end, _label in intervals])
        record = {
            "transcript": transcript,
            "k_analyzed_anchors": k,
            "candidate_count": len(candidates),
            "capacity_margin": len(candidates) - k,
            "capacity_status": "COORDINATE_NULL_ELIGIBLE" if eligible else "NULL_CAPACITY_INELIGIBLE",
            "CDS_length_nt": cds_length,
            "n_TMhelix_segments": counts["TMhelix"],
            "n_signal_segments": counts["signal"],
            "n_Beta_sheet_segments": counts["Beta sheet"],
            "n_disallowed_topology_segments": sum(counts.values()),
            "disallowed_topology_union_nt": covered_nt,
            "disallowed_topology_CDS_fraction": covered_nt / cds_length,
        }
        transcript_rows.append(record)
        capacity_by_transcript[transcript] = record

    prior = {row["transcript"]: row for row in prior_rows}
    prior_mismatches = 0
    for row in transcript_rows:
        transcript = row["transcript"]
        if transcript not in prior:
            prior_mismatches += 1
            continue
        if int(prior[transcript]["candidate_count"]) != row["candidate_count"]:
            prior_mismatches += 1
        if int(prior[transcript]["k_primary_window_valid_real_anchors"]) != row["k_analyzed_anchors"]:
            prior_mismatches += 1
    prior_mismatches += len(set(prior) ^ set(capacity_by_transcript))

    universe_rows: list[dict[str, Any]] = []
    for row in sorted(primary_valid, key=lambda value: (value["transcript"], int(value["TMD_index"]))):
        capacity = capacity_by_transcript[row["transcript"]]
        cds_length = capacity["CDS_length_nt"]
        anchor_cds0 = int(row["real_anchor_cds0"])
        universe_rows.append({
            "transcript": row["transcript"],
            "TMD_index": row["TMD_index"],
            "n_primary_TMDs": row["n_primary_TMDs"],
            "pass_class": row["pass_class"],
            "TMD_position_class": row["TMD_position_class"],
            "real_TMD_length_nt": row["real_TMD_length_nt"],
            "real_anchor_coordinate": row["real_anchor_coordinate"],
            "real_anchor_cds0": anchor_cds0,
            "CDS_length_nt": cds_length,
            "normalized_CDS_position": anchor_cds0 / cds_length,
            "distance_to_CDS_start_nt": anchor_cds0,
            "distance_to_CDS_end_nt": cds_length - anchor_cds0,
            "k_per_transcript": capacity["k_analyzed_anchors"],
            "candidate_count": capacity["candidate_count"],
            "capacity_margin": capacity["capacity_margin"],
            "coordinate_null_capacity_status": capacity["capacity_status"],
            "status_scope": "COORDINATE_NULL_ONLY",
            "primary_window_valid": "YES",
            "primary_anchor_table_freeze_version": row["anchor_table_freeze_version"],
            "pseudo_anchor_table_freeze_version": row["pseudo_anchor_table_freeze_version"],
            "coordinate_null_policy_freeze_version": FREEZE_VERSION,
            "candidate_rule": "SAME_CDS_CODON_START;ANCHOR_AND_[+60,+91)_OUTSIDE_TMhelix|signal|Beta_sheet",
            "production_null_generated": "NO",
            "scope": "ANNOTATION_ONLY_NO_RIBOSEQ",
            "script_version": SCRIPT_VERSION,
        })

    eligible_anchors = [row for row in universe_rows if row["coordinate_null_capacity_status"] == "COORDINATE_NULL_ELIGIBLE"]
    ineligible_anchors = [row for row in universe_rows if row["coordinate_null_capacity_status"] == "NULL_CAPACITY_INELIGIBLE"]
    eligible_transcripts = [row for row in transcript_rows if row["capacity_status"] == "COORDINATE_NULL_ELIGIBLE"]
    ineligible_transcripts = [row for row in transcript_rows if row["capacity_status"] == "NULL_CAPACITY_INELIGIBLE"]

    observed = (
        len(universe_rows), len(transcript_rows), len(eligible_anchors), len(eligible_transcripts),
        len(ineligible_anchors), len(ineligible_transcripts),
    )
    expected = (
        EXPECTED_PRIMARY_VALID_ANCHORS, EXPECTED_PRIMARY_VALID_TRANSCRIPTS,
        EXPECTED_ELIGIBLE_ANCHORS, EXPECTED_ELIGIBLE_TRANSCRIPTS,
        EXPECTED_INELIGIBLE_ANCHORS, EXPECTED_INELIGIBLE_TRANSCRIPTS,
    )
    if observed != expected:
        raise RuntimeError(f"Capacity universe mismatch: observed={observed} expected={expected}")
    if prior_mismatches:
        raise RuntimeError(f"Independent candidate recomputation differs from prior feasibility: {prior_mismatches}")
    if any(
        len({row["coordinate_null_capacity_status"] for row in universe_rows if row["transcript"] == transcript}) != 1
        for transcript in k_by_transcript
    ):
        raise RuntimeError("Partial transcript status detected.")

    write_tsv(args.universe_tsv, universe_rows)
    anchor_groups = {
        "PRIMARY_WINDOW_VALID": universe_rows,
        "COORDINATE_NULL_ELIGIBLE": eligible_anchors,
        "NULL_CAPACITY_INELIGIBLE": ineligible_anchors,
    }
    transcript_groups = {
        "PRIMARY_WINDOW_VALID": transcript_rows,
        "COORDINATE_NULL_ELIGIBLE": eligible_transcripts,
        "NULL_CAPACITY_INELIGIBLE": ineligible_transcripts,
    }
    representative_rows = summary_rows(anchor_groups, transcript_groups)
    write_tsv(args.representativeness_tsv, representative_rows)

    validation = {
        "script_version": SCRIPT_VERSION,
        "freeze_version": FREEZE_VERSION,
        "scope": "annotation_only_no_riboseq_no_null_draws",
        "inputs": {
            "pairing_table": {"path": str(args.pairs), "sha256": sha256(args.pairs)},
            "topology": {"path": str(args.topology), "sha256": sha256(args.topology)},
            "features": {"path": str(args.features), "sha256": sha256(args.features)},
            "prior_feasibility": {"path": str(args.prior_feasibility), "sha256": sha256(args.prior_feasibility)},
        },
        "outputs": {
            "universe": {"path": str(args.universe_tsv), "sha256": sha256(args.universe_tsv)},
            "representativeness": {"path": str(args.representativeness_tsv), "sha256": sha256(args.representativeness_tsv)},
        },
        "counts": {
            "primary_window_valid_anchors": len(universe_rows),
            "primary_window_valid_transcripts": len(transcript_rows),
            "coordinate_null_eligible_anchors": len(eligible_anchors),
            "coordinate_null_eligible_transcripts": len(eligible_transcripts),
            "null_capacity_ineligible_anchors": len(ineligible_anchors),
            "null_capacity_ineligible_transcripts": len(ineligible_transcripts),
        },
        "invariants": {
            "anchor_count_conservation": len(eligible_anchors) + len(ineligible_anchors) == len(universe_rows),
            "transcript_count_conservation": len(eligible_transcripts) + len(ineligible_transcripts) == len(transcript_rows),
            "whole_transcript_classification": True,
            "independent_candidate_count_mismatches": prior_mismatches,
            "sampling_with_replacement": False,
            "partial_transcript_dropping": False,
            "riboseq_used": False,
            "production_null_generated": False,
        },
    }
    args.validation_json.parent.mkdir(parents=True, exist_ok=True)
    with args.validation_json.open("w", encoding="utf-8") as handle:
        json.dump(validation, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
