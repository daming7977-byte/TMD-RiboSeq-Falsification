#!/usr/bin/env python3
"""Annotation-only feasibility audit for matched pseudo-TMD anchors.

This is not production Module 05 and does not read Ribo-seq data. It evaluates
the frozen primary-TMD, topology-label, anchor, and window contracts and emits
aggregate audit JSON only; it deliberately does not write a pseudo-anchor table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SCRIPT_VERSION = "audit_only_pseudo_anchor_feasibility_v1_2026-08-25"
EXPECTED_ANCHOR_FREEZE = "stage0c_primary_tmd_anchors_v1_2026-08-24"
EXPECTED_LABEL_FREEZE = "stage0c_tmd_label_policy_v1_2026-08-24"
EXPECTED_ANCHORS = 19_784
CENTERS = (-30, -15, 0, 15, 30, 45, 60, 75, 90, 105, 120)
PRIMARY_CENTER = 75
HALF_WIDTH = 15
WINDOW_WIDTH = 31
MEMBRANE_FEATURE_LABELS = frozenset({"TMhelix", "signal", "Beta sheet"})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def overlaps(start: int, end: int, other_start: int, other_end: int) -> bool:
    return max(start, other_start) < min(end, other_end)


@dataclass(frozen=True)
class RealAnchor:
    transcript: str
    tmd_index: int
    aa_start: int
    aa_end: int
    cds_start: int
    cds_end: int
    protein_length: int
    real_mask: tuple[bool, ...]

    @property
    def aa_length(self) -> int:
        return self.aa_end - self.aa_start


class Edge:
    __slots__ = ("to", "rev", "cap", "cost", "order")

    def __init__(self, to: int, rev: int, cap: int, cost: int, order: int) -> None:
        self.to = to
        self.rev = rev
        self.cap = cap
        self.cost = cost
        self.order = order


def deterministic_min_cost_assignment(
    reals: list[RealAnchor], candidate_starts: list[list[int]]
) -> dict[int, int]:
    """Maximum-cardinality, minimum-total-distance unique-start assignment.

    Equal reduced-cost paths retain the first-discovered predecessor under
    stable node and edge order. Real nodes are in TMD_index order; pseudo-start
    nodes are in ascending coordinate order; real-to-candidate edges are
    inserted in ascending coordinate order.
    """

    candidate_universe = sorted({value for values in candidate_starts for value in values})
    candidate_rank = {value: index for index, value in enumerate(candidate_universe)}
    n_real = len(reals)
    n_candidate = len(candidate_universe)
    source = 0
    real_base = 1
    candidate_base = real_base + n_real
    sink = candidate_base + n_candidate
    graph: list[list[Edge]] = [[] for _ in range(sink + 1)]
    edge_order = 0

    def add_edge(left: int, right: int, capacity: int, cost: int) -> Edge:
        nonlocal edge_order
        forward = Edge(right, len(graph[right]), capacity, cost, edge_order)
        edge_order += 1
        reverse = Edge(left, len(graph[left]), 0, -cost, edge_order)
        edge_order += 1
        graph[left].append(forward)
        graph[right].append(reverse)
        return forward

    for real_index in range(n_real):
        add_edge(source, real_base + real_index, 1, 0)
    for rank in range(n_candidate):
        add_edge(candidate_base + rank, sink, 1, 0)

    assignment_edges: dict[tuple[int, int], Edge] = {}
    for real_index, (real, values) in enumerate(zip(reals, candidate_starts)):
        for candidate in sorted(values):
            rank = candidate_rank[candidate]
            edge = add_edge(
                real_base + real_index,
                candidate_base + rank,
                1,
                abs(candidate - real.cds_start),
            )
            assignment_edges[(real_index, candidate)] = edge

    node_count = len(graph)
    potential = [0] * node_count
    infinity = 10**30

    while True:
        distance = [infinity] * node_count
        predecessor: list[tuple[int, int] | None] = [None] * node_count
        distance[source] = 0
        queue: list[tuple[int, int]] = [(0, source)]
        while queue:
            current_distance, node = heapq.heappop(queue)
            if current_distance != distance[node]:
                continue
            for edge_index, edge in enumerate(graph[node]):
                if edge.cap <= 0:
                    continue
                reduced = edge.cost + potential[node] - potential[edge.to]
                proposal = current_distance + reduced
                if proposal < distance[edge.to]:
                    distance[edge.to] = proposal
                    predecessor[edge.to] = (node, edge_index)
                    heapq.heappush(queue, (proposal, edge.to))
        if predecessor[sink] is None:
            break
        for node, value in enumerate(distance):
            if value < infinity:
                potential[node] += value
        node = sink
        while node != source:
            previous, edge_index = predecessor[node]  # type: ignore[misc]
            edge = graph[previous][edge_index]
            edge.cap -= 1
            graph[node][edge.rev].cap += 1
            node = previous

    assignments: dict[int, int] = {}
    for (real_index, candidate), edge in assignment_edges.items():
        if edge.cap == 0:
            if real_index in assignments:
                raise RuntimeError("A real anchor received more than one pseudo assignment.")
            assignments[real_index] = candidate
    if len(set(assignments.values())) != len(assignments):
        raise RuntimeError("Pseudo-anchor reuse occurred in deterministic assignment.")
    return assignments


def deterministic_nonoverlap_assignment(
    reals: list[RealAnchor], candidate_starts: list[list[int]]
) -> dict[int, int]:
    """Joint exact assignment with unique, non-overlapping pseudo segments.

    The lexicographic objective is: maximize matched real anchors, minimize total
    absolute CDS-nucleotide displacement, then minimize the assignment vector in
    ascending real-TMD order. An unmatched token sorts after every CDS coordinate.
    """

    if len(reals) != len(candidate_starts):
        raise ValueError("Real/candidate list length mismatch.")
    if not reals:
        return {}

    options: list[list[tuple[int, int, int, int]]] = []
    for real, starts in zip(reals, candidate_starts):
        real_options = []
        for start_cds in starts:
            start_aa = start_cds // 3
            end_aa = start_aa + real.aa_length
            real_options.append((abs(start_cds - real.cds_start), start_cds, start_aa, end_aa))
        real_options.sort()
        options.append(real_options)

    unmatched_token = max(real.protein_length * 3 for real in reals) + 1
    best_matched = -1
    best_cost = 10**30
    best_vector: tuple[int, ...] | None = None
    best_assignment: dict[int, int] = {}

    def interval_free(start: int, end: int, selected: list[tuple[int, int]]) -> bool:
        return all(not overlaps(start, end, other_start, other_end) for other_start, other_end in selected)

    def visit(
        remaining: tuple[int, ...],
        selected: list[tuple[int, int]],
        assignment: dict[int, int],
        matched_count: int,
        total_cost: int,
    ) -> None:
        nonlocal best_matched, best_cost, best_vector, best_assignment
        if not remaining:
            vector = tuple(assignment.get(index, unmatched_token) for index in range(len(reals)))
            candidate_key = (-matched_count, total_cost, vector)
            best_key = (
                -best_matched,
                best_cost,
                best_vector if best_vector is not None else tuple([unmatched_token] * len(reals)),
            )
            if best_vector is None or candidate_key < best_key:
                best_matched = matched_count
                best_cost = total_cost
                best_vector = vector
                best_assignment = dict(assignment)
            return

        feasible_by_real: dict[int, list[tuple[int, int, int, int]]] = {}
        minimum_costs = []
        for real_index in remaining:
            feasible = [
                option
                for option in options[real_index]
                if interval_free(option[2], option[3], selected)
            ]
            feasible_by_real[real_index] = feasible
            if feasible:
                minimum_costs.append(feasible[0][0])

        upper_matched = matched_count + len(minimum_costs)
        if upper_matched < best_matched:
            return
        if upper_matched == best_matched:
            needed = best_matched - matched_count
            lower_additional_cost = sum(sorted(minimum_costs)[:needed])
            if total_cost + lower_additional_cost > best_cost:
                return

        real_index = min(remaining, key=lambda index: (len(feasible_by_real[index]), index))
        next_remaining = tuple(index for index in remaining if index != real_index)
        feasible_options = feasible_by_real[real_index]

        for option_cost, start_cds, start_aa, end_aa in feasible_options:
            if best_matched == len(reals):
                optimistic = total_cost + option_cost
                for other_index in next_remaining:
                    other_feasible = [
                        candidate
                        for candidate in feasible_by_real[other_index]
                        if not overlaps(start_aa, end_aa, candidate[2], candidate[3])
                    ]
                    if not other_feasible:
                        optimistic = best_cost + 1
                        break
                    optimistic += other_feasible[0][0]
                if optimistic > best_cost:
                    continue
            assignment[real_index] = start_cds
            selected.append((start_aa, end_aa))
            visit(
                next_remaining,
                selected,
                assignment,
                matched_count + 1,
                total_cost + option_cost,
            )
            selected.pop()
            del assignment[real_index]

        visit(next_remaining, selected, assignment, matched_count, total_cost)

    visit(tuple(range(len(reals))), [], {}, 0, 0)
    return best_assignment


def quantiles(values: list[int]) -> dict[str, float | int | None]:
    if not values:
        return {key: None for key in ("min", "p25", "median", "p75", "p90", "p95", "max", "mean")}

    ordered = sorted(values)

    def linear_quantile(probability: float) -> float:
        if len(ordered) == 1:
            return float(ordered[0])
        position = (len(ordered) - 1) * probability
        left = math.floor(position)
        right = math.ceil(position)
        if left == right:
            return float(ordered[left])
        fraction = position - left
        return ordered[left] * (1.0 - fraction) + ordered[right] * fraction

    return {
        "min": ordered[0],
        "p25": linear_quantile(0.25),
        "median": linear_quantile(0.50),
        "p75": linear_quantile(0.75),
        "p90": linear_quantile(0.90),
        "p95": linear_quantile(0.95),
        "max": ordered[-1],
        "mean": statistics.fmean(ordered),
    }


def binned_counts(values: Iterable[int]) -> dict[str, int]:
    bins = Counter()
    for value in values:
        if value == 0:
            bins["0"] += 1
        elif value <= 9:
            bins["1-9"] += 1
        elif value <= 24:
            bins["10-24"] += 1
        elif value <= 49:
            bins["25-49"] += 1
        elif value <= 99:
            bins["50-99"] += 1
        elif value <= 199:
            bins["100-199"] += 1
        elif value <= 499:
            bins["200-499"] += 1
        elif value <= 999:
            bins["500-999"] += 1
        else:
            bins[">=1000"] += 1
    return {key: bins[key] for key in ("0", "1-9", "10-24", "25-49", "50-99", "100-199", "200-499", "500-999", ">=1000")}


def build_mask(
    anchor_cds0: int,
    cds_length: int,
    accepted_tmd_intervals: list[tuple[int, int]],
    skip_interval_index: int | None,
) -> tuple[bool, ...]:
    mask: list[bool] = []
    for center in CENTERS:
        start = anchor_cds0 + center - HALF_WIDTH
        end = anchor_cds0 + center + HALF_WIDTH + 1
        valid = start >= 0 and end <= cds_length
        if valid:
            for interval_index, (other_start, other_end) in enumerate(accepted_tmd_intervals):
                if skip_interval_index is not None and interval_index == skip_interval_index:
                    continue
                if overlaps(start, end, other_start, other_end):
                    valid = False
                    break
        mask.append(valid)
    return tuple(mask)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    return parser.parse_args()


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

    required_anchor = {
        "transcript", "protein_length", "TMD_index", "raw_deeptmhmm_label",
        "aa_start0", "aa_end0_exclusive", "cds_start", "cds_nt_start0",
        "cds_nt_end0_exclusive", "primary_anchor_coordinate",
        "primary_anchor_cds_coordinate", "tmd_label_policy_freeze_version",
        "anchor_table_freeze_version",
    }
    if not required_anchor.issubset(anchor_rows[0]):
        raise RuntimeError(f"Anchor table missing fields: {sorted(required_anchor - set(anchor_rows[0]))}")

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
        if row["raw_deeptmhmm_label"] != "TMhelix":
            raise RuntimeError("Anchor table contains a non-TMhelix row.")
        if row["tmd_label_policy_freeze_version"] != EXPECTED_LABEL_FREEZE:
            raise RuntimeError("Unexpected TMD label-policy freeze version.")
        if row["anchor_table_freeze_version"] != EXPECTED_ANCHOR_FREEZE:
            raise RuntimeError("Unexpected primary-anchor-table freeze version.")
        if int(row["primary_anchor_coordinate"]) != int(row["cds_start"]) + int(row["cds_nt_start0"]):
            raise RuntimeError("Anchor coordinate invariant failed.")
        anchor_by_transcript[row["transcript"]].append(row)

    primary_index = CENTERS.index(PRIMARY_CENTER)
    total_primary_valid = 0
    with_candidate = 0
    matched = 0
    exact_mask_candidate_associations = 0
    candidate_counts_all: list[int] = []
    candidate_counts_primary_valid: list[int] = []
    mismatch_values: list[int] = []
    unmatched_reasons = Counter()
    real_mask_patterns = Counter()
    selected_segment_overlap_pairs = 0
    assignment_transcripts = 0

    ordered_transcripts = sorted(anchor_by_transcript)
    for transcript_number, transcript in enumerate(ordered_transcripts, start=1):
        if transcript_number % 500 == 0 or transcript_number == len(ordered_transcripts):
            print(
                f"[annotation audit] transcript {transcript_number}/{len(ordered_transcripts)}",
                file=sys.stderr,
                flush=True,
            )
        rows = sorted(anchor_by_transcript[transcript], key=lambda row: int(row["TMD_index"]))
        observed_indices = [int(row["TMD_index"]) for row in rows]
        if observed_indices != list(range(1, len(rows) + 1)):
            raise RuntimeError(f"Noncontiguous TMD indices for {transcript}.")
        if transcript not in features:
            raise RuntimeError(f"Missing feature row for {transcript}.")
        cds_start_tx0, cds_length = features[transcript]
        protein_lengths = {int(row["protein_length"]) for row in rows}
        cds_starts = {int(row["cds_start"]) for row in rows}
        if len(protein_lengths) != 1 or len(cds_starts) != 1 or cds_start_tx0 not in cds_starts:
            raise RuntimeError(f"Inconsistent protein/CDS metadata for {transcript}.")
        protein_length = next(iter(protein_lengths))
        accepted_intervals = [
            (int(row["cds_nt_start0"]), int(row["cds_nt_end0_exclusive"])) for row in rows
        ]
        real_aa_starts = {int(row["aa_start0"]) for row in rows}

        blocked = [0] * protein_length
        for start, end, _label in membrane_by_transcript.get(transcript, []):
            if start < 0 or end < start or end > protein_length:
                raise RuntimeError(f"Invalid membrane interval for {transcript}: {start},{end}")
            for position in range(start, end):
                blocked[position] = 1
        blocked_prefix = [0]
        for value in blocked:
            blocked_prefix.append(blocked_prefix[-1] + value)

        reals: list[RealAnchor] = []
        candidate_starts_by_real: list[list[int]] = []
        candidate_segments_by_real: list[dict[int, tuple[int, int]]] = []
        pseudo_mask_cache: dict[int, tuple[bool, ...]] = {}

        for interval_index, row in enumerate(rows):
            aa_start = int(row["aa_start0"])
            aa_end = int(row["aa_end0_exclusive"])
            real_cds_start = int(row["cds_nt_start0"])
            real_cds_end = int(row["cds_nt_end0_exclusive"])
            if real_cds_start != 3 * aa_start or real_cds_end != 3 * aa_end:
                raise RuntimeError(f"Codon invariant failed for {transcript} TMD {interval_index + 1}.")
            real_mask = build_mask(real_cds_start, cds_length, accepted_intervals, interval_index)
            real_mask_patterns["".join("1" if value else "0" for value in real_mask)] += 1
            real = RealAnchor(
                transcript=transcript,
                tmd_index=interval_index + 1,
                aa_start=aa_start,
                aa_end=aa_end,
                cds_start=real_cds_start,
                cds_end=real_cds_end,
                protein_length=protein_length,
                real_mask=real_mask,
            )
            reals.append(real)
            if not real_mask[primary_index]:
                candidate_starts_by_real.append([])
                candidate_segments_by_real.append({})
                candidate_counts_all.append(0)
                unmatched_reasons["REAL_PRIMARY_WINDOW_INVALID"] += 1
                continue

            total_primary_valid += 1
            length_aa = real.aa_length
            candidates: list[int] = []
            segments: dict[int, tuple[int, int]] = {}
            for pseudo_aa_start in range(0, protein_length - length_aa + 1):
                pseudo_aa_end = pseudo_aa_start + length_aa
                if pseudo_aa_start in real_aa_starts:
                    continue
                if blocked_prefix[pseudo_aa_end] - blocked_prefix[pseudo_aa_start] != 0:
                    continue
                pseudo_cds_start = 3 * pseudo_aa_start
                pseudo_mask = pseudo_mask_cache.get(pseudo_cds_start)
                if pseudo_mask is None:
                    pseudo_mask = build_mask(pseudo_cds_start, cds_length, accepted_intervals, None)
                    pseudo_mask_cache[pseudo_cds_start] = pseudo_mask
                if pseudo_mask != real_mask:
                    continue
                candidates.append(pseudo_cds_start)
                segments[pseudo_cds_start] = (pseudo_aa_start, pseudo_aa_end)
            candidates.sort()
            candidate_starts_by_real.append(candidates)
            candidate_segments_by_real.append(segments)
            count = len(candidates)
            candidate_counts_all.append(count)
            candidate_counts_primary_valid.append(count)
            exact_mask_candidate_associations += count
            if count:
                with_candidate += 1
            else:
                unmatched_reasons["NO_EXACT_MASK_LENGTH_MATCHED_PSEUDO_CANDIDATE"] += 1

        eligible_indices = [index for index, values in enumerate(candidate_starts_by_real) if values]
        if eligible_indices:
            assignment_transcripts += 1
            eligible_reals = [reals[index] for index in eligible_indices]
            eligible_candidates = [candidate_starts_by_real[index] for index in eligible_indices]
            partial_assignment = deterministic_min_cost_assignment(eligible_reals, eligible_candidates)
            assignments = {eligible_indices[index]: value for index, value in partial_assignment.items()}
        else:
            assignments = {}

        matched += len(assignments)
        for real_index, pseudo_cds_start in assignments.items():
            mismatch_values.append(abs(pseudo_cds_start - reals[real_index].cds_start))
        for real_index in eligible_indices:
            if real_index not in assignments:
                unmatched_reasons["ONE_TO_ONE_ASSIGNMENT_CONFLICT"] += 1

        selected_segments = []
        for real_index, pseudo_cds_start in assignments.items():
            selected_segments.append(candidate_segments_by_real[real_index][pseudo_cds_start])
        for left_index, (left_start, left_end) in enumerate(selected_segments):
            for right_start, right_end in selected_segments[left_index + 1 :]:
                if overlaps(left_start, left_end, right_start, right_end):
                    selected_segment_overlap_pairs += 1

    if len(candidate_counts_all) != EXPECTED_ANCHORS:
        raise RuntimeError("Candidate accounting does not cover every real anchor.")
    if sum(unmatched_reasons.values()) + matched != EXPECTED_ANCHORS:
        raise RuntimeError("Matched/unmatched reason accounting failed.")
    if matched > with_candidate or with_candidate > total_primary_valid:
        raise RuntimeError("Eligibility accounting order failed.")

    summary = {
        "script_version": SCRIPT_VERSION,
        "scope": "annotation_only_no_riboseq",
        "policy_contract": {
            "matching_unit": "same_transcript",
            "pseudo_segment": "same_length_as_real_tmd_codon_aligned_self_exempt_focal_geometry",
            "candidate_segment_topology_exclusions": sorted(MEMBRANE_FEATURE_LABELS),
            "window_validity_other_feature_set": ["accepted_primary_TMhelix"],
            "validity_mask": "exact_11_center_match",
            "primary_window_required": True,
            "assignment": "per_transcript_joint_max_cardinality_then_min_total_absolute_cds_nt_displacement",
            "pseudo_anchor_reuse": False,
            "selected_pseudo_segment_overlap": "allowed_as_geometric_construct_but_audited",
            "tie_break": "sorted_TMD_index_sorted_candidate_coordinate_stable_SSP_first_discovered_equal_cost_path",
            "random_seed": None,
        },
        "inputs": {
            "anchors": {"path": str(args.anchors), "sha256": sha256(args.anchors)},
            "topology": {"path": str(args.topology), "sha256": sha256(args.topology)},
            "features": {"path": str(args.features), "sha256": sha256(args.features)},
        },
        "counts": {
            "real_anchors": EXPECTED_ANCHORS,
            "real_primary_window_eligible": total_primary_valid,
            "with_pseudo_candidate": with_candidate,
            "matched": matched,
            "eligible_for_primary_matched_analysis": matched,
            "unmatched": EXPECTED_ANCHORS - matched,
            "exact_mask_candidate_associations": exact_mask_candidate_associations,
            "assignment_transcripts": assignment_transcripts,
            "selected_pseudo_segment_overlap_pairs": selected_segment_overlap_pairs,
        },
        "unmatched_reasons": dict(sorted(unmatched_reasons.items())),
        "candidate_count_distribution_all_real_anchors": {
            "summary": quantiles(candidate_counts_all),
            "bins": binned_counts(candidate_counts_all),
        },
        "candidate_count_distribution_primary_valid_real_anchors": {
            "summary": quantiles(candidate_counts_primary_valid),
            "bins": binned_counts(candidate_counts_primary_valid),
        },
        "positional_mismatch_nt_distribution_matched": quantiles(mismatch_values),
        "real_validity_mask_counts": dict(sorted(real_mask_patterns.items())),
        "invariants": {
            "real_anchor_accounting": total_primary_valid + unmatched_reasons["REAL_PRIMARY_WINDOW_INVALID"] == EXPECTED_ANCHORS,
            "matched_unmatched_accounting": matched + sum(unmatched_reasons.values()) == EXPECTED_ANCHORS,
            "all_candidates_exact_mask": True,
            "all_matched_exact_mask": True,
            "matched_primary_windows_valid": True,
            "pseudo_anchor_reuse_count": 0,
            "riboseq_read": False,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
