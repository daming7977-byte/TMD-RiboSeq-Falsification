#!/usr/bin/env python3
"""Audit GSE133111 coordinate-null capacity using frozen annotations only.

This script does not read Ribo-seq, draw random coordinates, generate null
realizations, change a frozen policy, or create Module 05.
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
from typing import Any


SCRIPT_VERSION = "audit_only_coordinate_null_feasibility_v1_2026-08-25"
EXPECTED_NULL_ANCHORS = 15_597
EXPECTED_PRIMARY_TMD_ROWS = 19_784
EXPECTED_ANCHOR_FREEZE = "stage0c_primary_tmd_anchors_v1_2026-08-24"
EXPECTED_PSEUDO_TABLE_FREEZE = "stage0c_pseudo_anchor_table_v1_2026-08-25"
TOPOLOGY_LABELS = frozenset({"TMhelix", "signal", "Beta sheet"})
PRIMARY_START = 60
PRIMARY_END_EXCLUSIVE = 91


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
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


def sha256_coordinates(coordinates: list[int]) -> str:
    payload = (",".join(str(value) for value in coordinates) + "\n").encode()
    return hashlib.sha256(payload).hexdigest()


def overlaps(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return start_a < end_b and start_b < end_a


def linear_quantile(values: list[int], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    left = math.floor(position)
    right = math.ceil(position)
    if left == right:
        return float(ordered[left])
    fraction = position - left
    return ordered[left] * (1.0 - fraction) + ordered[right] * fraction


def main() -> None:
    args = parse_args()
    for path in (args.pairs, args.topology, args.features):
        if not path.is_file():
            raise FileNotFoundError(path)

    pair_rows = read_tsv(args.pairs)
    topology_rows = read_tsv(args.topology)
    feature_rows = read_tsv(args.features)
    if len(pair_rows) != EXPECTED_PRIMARY_TMD_ROWS:
        raise RuntimeError(f"Expected {EXPECTED_PRIMARY_TMD_ROWS} pairing rows, found {len(pair_rows)}.")
    if {row["pseudo_anchor_table_freeze_version"] for row in pair_rows} != {EXPECTED_PSEUDO_TABLE_FREEZE}:
        raise RuntimeError("Unexpected pseudo-anchor-table freeze version.")
    if {row["anchor_table_freeze_version"] for row in pair_rows} != {EXPECTED_ANCHOR_FREEZE}:
        raise RuntimeError("Unexpected primary-anchor-table freeze version.")

    features: dict[str, int] = {}
    for row in feature_rows:
        transcript = row["Transcript"]
        if transcript in features:
            raise RuntimeError(f"Duplicate feature transcript: {transcript}")
        features[transcript] = int(row["CDS_Length"])

    k_by_transcript: Counter[str] = Counter()
    source_indices: dict[str, list[int]] = defaultdict(list)
    for row in pair_rows:
        if row["real_primary_window_valid"] == "YES":
            transcript = row["transcript"]
            k_by_transcript[transcript] += 1
            source_indices[transcript].append(int(row["TMD_index"]))
    if sum(k_by_transcript.values()) != EXPECTED_NULL_ANCHORS:
        raise RuntimeError(
            f"Expected {EXPECTED_NULL_ANCHORS} primary-window-valid anchors, "
            f"found {sum(k_by_transcript.values())}."
        )

    topology_by_transcript: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for row in topology_rows:
        label = row["raw_deeptmhmm_label"]
        if label not in TOPOLOGY_LABELS:
            continue
        start = int(row["cds_nt_start0"])
        end = int(row["cds_nt_end0_exclusive"])
        if start < 0 or end <= start:
            raise RuntimeError(f"Invalid topology interval: {row}")
        topology_by_transcript[row["transcript"]].append((start, end, label))
    for transcript in topology_by_transcript:
        topology_by_transcript[transcript].sort(key=lambda value: (value[0], value[1], value[2]))

    output_rows: list[dict[str, Any]] = []
    for transcript in sorted(k_by_transcript):
        cds_length = features[transcript]
        if cds_length <= 0 or cds_length % 3:
            raise RuntimeError(f"Invalid CDS length for null grid: {transcript} {cds_length}")
        intervals = topology_by_transcript.get(transcript, [])
        candidates: list[int] = []
        boundary_excluded = 0
        anchor_topology_excluded = 0
        primary_window_topology_excluded = 0
        for coordinate in range(0, cds_length, 3):
            primary_start = coordinate + PRIMARY_START
            primary_end = coordinate + PRIMARY_END_EXCLUSIVE
            if primary_start < 0 or primary_end > cds_length:
                boundary_excluded += 1
                continue
            if any(start <= coordinate < end for start, end, _label in intervals):
                anchor_topology_excluded += 1
                continue
            if any(overlaps(primary_start, primary_end, start, end) for start, end, _label in intervals):
                primary_window_topology_excluded += 1
                continue
            candidates.append(coordinate)

        k_required = k_by_transcript[transcript]
        feasible = len(candidates) >= k_required
        output_rows.append(
            {
                "transcript": transcript,
                "CDS_length_nt": cds_length,
                "k_primary_window_valid_real_anchors": k_required,
                "source_TMD_indices_in_order": ",".join(str(value) for value in sorted(source_indices[transcript])),
                "total_codon_start_coordinates": cds_length // 3,
                "excluded_primary_window_CDS_boundary": boundary_excluded,
                "excluded_null_anchor_inside_topology": anchor_topology_excluded,
                "excluded_primary_window_topology_overlap": primary_window_topology_excluded,
                "candidate_count": len(candidates),
                "capacity_margin": len(candidates) - k_required,
                "joint_without_replacement_feasible": "YES" if feasible else "NO",
                "real_anchors_affected_if_infeasible": 0 if feasible else k_required,
                "candidate_first_cds0": candidates[0] if candidates else "NA",
                "candidate_last_cds0": candidates[-1] if candidates else "NA",
                "candidate_coordinates_sha256": sha256_coordinates(candidates),
                "candidate_order": "ASCENDING_CDS0",
                "candidate_coordinate_system": "CDS_NT_0BASED_CODON_START",
                "primary_window": "[+60,+91)",
                "topology_exclusions": "TMhelix|signal|Beta sheet",
                "scope": "ANNOTATION_ONLY_NO_RIBOSEQ_NO_NULL_DRAWS",
                "script_version": SCRIPT_VERSION,
            }
        )

    fields = list(output_rows[0])
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    candidate_counts = [int(row["candidate_count"]) for row in output_rows]
    insufficient = [row for row in output_rows if row["joint_without_replacement_feasible"] == "NO"]
    zero_candidate = [row for row in insufficient if int(row["candidate_count"]) == 0]
    positive_insufficient = [row for row in insufficient if int(row["candidate_count"]) > 0]
    exclusion_totals = {
        field: sum(int(row[field]) for row in output_rows)
        for field in (
            "total_codon_start_coordinates",
            "excluded_primary_window_CDS_boundary",
            "excluded_null_anchor_inside_topology",
            "excluded_primary_window_topology_overlap",
            "candidate_count",
        )
    }
    result = {
        "script_version": SCRIPT_VERSION,
        "scope": "annotation_only_no_riboseq_no_null_draws",
        "inputs": {
            "pseudo_anchor_pairing_table": {"path": str(args.pairs), "sha256": sha256(args.pairs)},
            "normalized_topology": {"path": str(args.topology), "sha256": sha256(args.topology)},
            "human_features": {"path": str(args.features), "sha256": sha256(args.features)},
        },
        "output": {"transcript_feasibility_tsv": {"path": str(args.output_tsv), "sha256": sha256(args.output_tsv)}},
        "contract": {
            "null_unit": "joint_within_transcript_codon_start_sampling",
            "candidate_order": "ascending_CDS_0_based_coordinate",
            "sampling": "without_replacement_within_transcript",
            "assignment": "sort_sampled_coordinates_then_assign_to_source_analyzed_anchors_in_TMD_index_order",
            "primary_window": "[+60,+91)",
            "topology_exclusions": sorted(TOPOLOGY_LABELS),
            "exact_11_window_mask_preserved": False,
            "tmd_length_preserved": False,
        },
        "counts": {
            "intended_real_anchor_universe": sum(k_by_transcript.values()),
            "transcripts_in_null_universe": len(output_rows),
            "transcripts_candidate_count_gte_k": len(output_rows) - len(insufficient),
            "transcripts_candidate_count_lt_k": len(insufficient),
            "real_anchors_affected_by_insufficient_capacity": sum(
                int(row["real_anchors_affected_if_infeasible"]) for row in insufficient
            ),
            "zero_candidate_transcripts": len(zero_candidate),
            "real_anchors_in_zero_candidate_transcripts": sum(
                int(row["real_anchors_affected_if_infeasible"]) for row in zero_candidate
            ),
            "positive_but_insufficient_transcripts": len(positive_insufficient),
            "real_anchors_in_positive_but_insufficient_transcripts": sum(
                int(row["real_anchors_affected_if_infeasible"]) for row in positive_insufficient
            ),
            "total_capacity_shortfall_coordinates": sum(
                -int(row["capacity_margin"]) for row in insufficient
            ),
            **exclusion_totals,
        },
        "candidate_count_distribution": {
            "n": len(candidate_counts),
            "mean": statistics.fmean(candidate_counts),
            "min": min(candidate_counts),
            "p25": linear_quantile(candidate_counts, 0.25),
            "median": linear_quantile(candidate_counts, 0.50),
            "p75": linear_quantile(candidate_counts, 0.75),
            "p90": linear_quantile(candidate_counts, 0.90),
            "p95": linear_quantile(candidate_counts, 0.95),
            "max": max(candidate_counts),
        },
        "k_distribution": dict(sorted(Counter(k_by_transcript.values()).items())),
        "invariants": {
            "anchor_count_is_15597": sum(k_by_transcript.values()) == EXPECTED_NULL_ANCHORS,
            "candidate_accounting_per_transcript": all(
                int(row["total_codon_start_coordinates"])
                == int(row["excluded_primary_window_CDS_boundary"])
                + int(row["excluded_null_anchor_inside_topology"])
                + int(row["excluded_primary_window_topology_overlap"])
                + int(row["candidate_count"])
                for row in output_rows
            ),
            "all_candidates_codon_aligned_by_construction": True,
            "candidate_coordinates_unique_by_construction": True,
            "full_universe_joint_sampling_feasible": len(insufficient) == 0,
            "riboseq_used": False,
            "null_realizations_generated": False,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
