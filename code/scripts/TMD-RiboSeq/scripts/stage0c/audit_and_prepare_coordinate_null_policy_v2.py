#!/usr/bin/env python3
"""Audit and prepare the annotation-only GSE133111 coordinate-null v2 policy.

This script reads only frozen annotation geometry and v1 coordinate-null
manifests.  It does not read Ribo-seq counts, profiles, sample outcomes,
conditions, benchmark results, or the v1 production realization matrix.

The sole new candidate predicate is that the symmetric point-null background
remaining after the frozen topology mask and [q-30,q+136) exclusion has
positive length.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


SCRIPT_VERSION = "audit_and_prepare_coordinate_null_policy_v2_2026-08-25"
POLICY_VERSION = "stage0c_coordinate_null_policy_v2_2026-08-25"
SUPERSEDED_POLICY_VERSION = "stage0c_coordinate_null_policy_v1_2026-08-25"
SUPERSEDED_REALIZATIONS_VERSION = "stage0c_coordinate_null_realizations_v1_2026-08-25"
TOPOLOGY_LABELS = frozenset({"TMhelix", "signal", "Beta sheet"})
PRIMARY_START = 60
PRIMARY_END_EXCLUSIVE = 91
NULL_EXCLUSION_START = -30
NULL_EXCLUSION_END_EXCLUSIVE = 136
EXPECTED_V1_ANCHORS = 13_480
EXPECTED_V1_TRANSCRIPTS = 4_103
MATRIX_DTYPE = np.dtype("<i4")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v1-anchor-columns", type=Path, required=True)
    parser.add_argument("--v1-candidate-pools", type=Path, required=True)
    parser.add_argument("--v1-candidate-coordinates", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--docs-output-dir", type=Path, required=True)
    parser.add_argument("--reference-output-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_int32(values: np.ndarray) -> str:
    normalized = np.ascontiguousarray(values, dtype=MATRIX_DTYPE)
    return hashlib.sha256(normalized.tobytes(order="C")).hexdigest()


def load_features(path: Path) -> dict[str, int]:
    output: dict[str, int] = {}
    for row in read_tsv(path):
        transcript = row["Transcript"]
        if transcript in output:
            raise RuntimeError(f"Duplicate feature transcript: {transcript}")
        output[transcript] = int(row["CDS_Length"])
    return output


def load_topology(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    output: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for row in read_tsv(path):
        label = row["raw_deeptmhmm_label"]
        if label in TOPOLOGY_LABELS:
            output[row["transcript"]].append(
                (int(row["cds_nt_start0"]), int(row["cds_nt_end0_exclusive"]), label)
            )
    for transcript in output:
        output[transcript].sort(key=lambda value: (value[0], value[1], value[2]))
    return output


def topology_background_mask(cds_length: int, intervals: list[tuple[int, int, str]]) -> np.ndarray:
    """Return the CDS positions retained after inherited Paper 1 topology masking."""
    mask = np.ones(cds_length, dtype=np.bool_)
    for anchor_start, segment_end, _label in intervals:
        start = max(0, anchor_start - 30)
        end = min(cds_length, max(segment_end, anchor_start + 136))
        mask[start:end] = False
    return mask


def real_background_length(base_mask: np.ndarray, real_anchor_cds0: int) -> int:
    mask = base_mask.copy()
    start = max(0, real_anchor_cds0 + PRIMARY_START)
    end = min(len(mask), real_anchor_cds0 + PRIMARY_END_EXCLUSIVE)
    mask[start:end] = False
    return int(np.count_nonzero(mask))


def null_background_length(base_mask: np.ndarray, q: int) -> int:
    mask = base_mask.copy()
    exclusion_start = max(0, q + NULL_EXCLUSION_START)
    exclusion_end = min(len(mask), q + NULL_EXCLUSION_END_EXCLUSIVE)
    mask[exclusion_start:exclusion_end] = False
    target_start = max(0, q + PRIMARY_START)
    target_end = min(len(mask), q + PRIMARY_END_EXCLUSIVE)
    mask[target_start:target_end] = False
    return int(np.count_nonzero(mask))


def main() -> None:
    args = parse_args()
    required = (
        args.v1_anchor_columns,
        args.v1_candidate_pools,
        args.v1_candidate_coordinates,
        args.features,
        args.topology,
    )
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    v1_anchors = read_tsv(args.v1_anchor_columns)
    v1_pools = read_tsv(args.v1_candidate_pools)
    v1_candidates = np.load(args.v1_candidate_coordinates, mmap_mode="r", allow_pickle=False)
    if len(v1_anchors) != EXPECTED_V1_ANCHORS:
        raise RuntimeError(f"v1 anchor count mismatch: {len(v1_anchors)}")
    if len(v1_pools) != EXPECTED_V1_TRANSCRIPTS:
        raise RuntimeError(f"v1 transcript count mismatch: {len(v1_pools)}")
    if v1_candidates.ndim != 1 or v1_candidates.dtype != MATRIX_DTYPE:
        raise RuntimeError(f"v1 candidate-array contract mismatch: {v1_candidates.shape} {v1_candidates.dtype}")

    anchor_indices = [int(row["anchor_column_index0"]) for row in v1_anchors]
    if anchor_indices != list(range(EXPECTED_V1_ANCHORS)):
        raise RuntimeError("v1 source-anchor order is not contiguous and canonical")
    transcripts = [row["transcript"] for row in v1_pools]
    if transcripts != sorted(transcripts) or len(set(transcripts)) != len(transcripts):
        raise RuntimeError("v1 transcript ordering is not unique lexicographic order")

    features = load_features(args.features)
    topology = load_topology(args.topology)
    anchors_by_transcript: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in v1_anchors:
        anchors_by_transcript[row["transcript"]].append(row)

    transcript_rows: list[dict[str, Any]] = []
    anchor_rows: list[dict[str, Any]] = []
    eligible_candidate_chunks: list[np.ndarray] = []
    eligible_pool_rows: list[dict[str, Any]] = []
    eligible_anchor_rows: list[dict[str, Any]] = []
    reason_counts: Counter[str] = Counter()
    candidate_cursor = 0
    anchor_cursor = 0

    for pool in v1_pools:
        transcript = pool["transcript"]
        cds_length = features[transcript]
        base_mask = topology_background_mask(cds_length, topology.get(transcript, []))
        source_rows = sorted(anchors_by_transcript[transcript], key=lambda row: int(row["anchor_column_index0"]))
        k = int(pool["k"])
        if len(source_rows) != k:
            raise RuntimeError(f"k mismatch for {transcript}: {len(source_rows)} != {k}")

        old_start = int(pool["candidate_offset_start0"])
        old_end = int(pool["candidate_offset_end0_exclusive"])
        old_values = np.asarray(v1_candidates[old_start:old_end], dtype=MATRIX_DTYPE)
        old_count = int(pool["candidate_count"])
        if len(old_values) != old_count:
            raise RuntimeError(f"v1 candidate slice mismatch for {transcript}")

        bg_lengths = np.fromiter(
            (null_background_length(base_mask, int(q)) for q in old_values),
            dtype=np.int64,
            count=old_count,
        )
        keep = bg_lengths > 0
        new_values = old_values[keep]
        new_count = len(new_values)
        if new_count > old_count or np.any(
            np.fromiter(
                (null_background_length(base_mask, int(q)) for q in new_values),
                dtype=np.int64,
                count=new_count,
            ) <= 0
        ):
            raise RuntimeError(f"positive-background candidate filter failed for {transcript}")

        real_lengths = [real_background_length(base_mask, int(row["real_anchor_cds0"])) for row in source_rows]
        real_positive = all(value > 0 for value in real_lengths)
        capacity_before = old_count >= k
        capacity_after_candidates = new_count >= k
        capacity_after = real_positive and capacity_after_candidates
        reasons: list[str] = []
        if not real_positive:
            reasons.append("REAL_BACKGROUND_ZERO")
        if not capacity_after_candidates:
            reasons.append("NEW_CANDIDATE_COUNT_LT_K")
        reason = "ELIGIBLE" if not reasons else ";".join(reasons)
        reason_counts[reason] += 1

        transcript_rows.append({
            "transcript_order_index0_v1": int(pool["transcript_order_index0"]),
            "transcript": transcript,
            "old_candidate_count": old_count,
            "new_candidate_count": new_count,
            "removed_zero_background_candidates": old_count - new_count,
            "k": k,
            "capacity_before": "YES" if capacity_before else "NO",
            "all_real_background_positive": "YES" if real_positive else "NO",
            "real_background_zero_anchors": sum(value == 0 for value in real_lengths),
            "capacity_after": "YES" if capacity_after else "NO",
            "v2_status": "COORDINATE_NULL_V2_ELIGIBLE" if capacity_after else "COORDINATE_NULL_V2_INELIGIBLE",
            "exclusion_reason": reason,
            "old_candidate_coordinates_sha256": sha256_int32(old_values),
            "new_candidate_coordinates_sha256": sha256_int32(new_values),
            "candidate_order": "ASCENDING_CDS0",
            "candidate_rule_added_v2": "BACKGROUND_LENGTH_GT_0_UNDER_[q-30,q+136)_SYMMETRIC_BACKGROUND",
            "scope": "ANNOTATION_ONLY_NO_RIBOSEQ_NO_NULL_DRAWS",
            "script_version": SCRIPT_VERSION,
        })

        for row, real_length in zip(source_rows, real_lengths, strict=True):
            out = dict(row)
            out.update({
                "v1_anchor_column_index0": int(row["anchor_column_index0"]),
                "real_background_length_nt": real_length,
                "real_background_positive": "YES" if real_length > 0 else "NO",
                "v2_coordinate_null_status": (
                    "COORDINATE_NULL_V2_ELIGIBLE" if capacity_after else "COORDINATE_NULL_V2_INELIGIBLE"
                ),
                "v2_exclusion_reason": reason,
                "coordinate_null_policy_v2_freeze_version": POLICY_VERSION,
            })
            anchor_rows.append(out)

        if capacity_after:
            expected_tmd_indices = [int(row["TMD_index"]) for row in source_rows]
            eligible_pool_rows.append({
                "transcript_order_index0": len(eligible_pool_rows),
                "transcript_order_index1": len(eligible_pool_rows) + 1,
                "transcript": transcript,
                "k": k,
                "source_TMD_indices_in_order": ",".join(map(str, expected_tmd_indices)),
                "anchor_column_start0": anchor_cursor,
                "anchor_column_end0_exclusive": anchor_cursor + k,
                "candidate_offset_start0": candidate_cursor,
                "candidate_offset_end0_exclusive": candidate_cursor + new_count,
                "old_candidate_count": old_count,
                "candidate_count": new_count,
                "candidate_first_cds0": int(new_values[0]),
                "candidate_last_cds0": int(new_values[-1]),
                "candidate_coordinates_sha256": sha256_int32(new_values),
                "candidate_order": "ASCENDING_CDS0",
                "CDS_length_nt": cds_length,
                "sampling": "WITHOUT_REPLACEMENT",
                "coordinate_null_policy_freeze_version": POLICY_VERSION,
            })
            for source_row in source_rows:
                eligible_anchor_rows.append({
                    "anchor_column_index0": len(eligible_anchor_rows),
                    "anchor_column_index1": len(eligible_anchor_rows) + 1,
                    "v1_anchor_column_index0": int(source_row["anchor_column_index0"]),
                    "transcript": transcript,
                    "TMD_index": source_row["TMD_index"],
                    "real_anchor_coordinate": source_row["real_anchor_coordinate"],
                    "real_anchor_cds0": source_row["real_anchor_cds0"],
                    "TMD_position_class": source_row["TMD_position_class"],
                    "pass_class": source_row["pass_class"],
                    "k_per_transcript": k,
                    "candidate_count": new_count,
                    "coordinate_null_policy_freeze_version": POLICY_VERSION,
                    "matrix_coordinate_system": "CDS_NT_0BASED_CODON_START",
                    "column_order": "V1_RELATIVE_SOURCE_ORDER_AFTER_WHOLE_TRANSCRIPT_EXCLUSIONS",
                })
            eligible_candidate_chunks.append(new_values)
            candidate_cursor += new_count
            anchor_cursor += k

    if len(transcript_rows) != EXPECTED_V1_TRANSCRIPTS or len(anchor_rows) != EXPECTED_V1_ANCHORS:
        raise RuntimeError("v1 count conservation failed")
    if sum(int(row["k"]) for row in transcript_rows) != EXPECTED_V1_ANCHORS:
        raise RuntimeError("transcript k conservation failed")

    eligible_candidates = np.concatenate(eligible_candidate_chunks).astype(MATRIX_DTYPE, copy=False)
    eligible_transcripts = len(eligible_pool_rows)
    eligible_anchors = len(eligible_anchor_rows)
    ineligible_transcripts = EXPECTED_V1_TRANSCRIPTS - eligible_transcripts
    ineligible_anchors = EXPECTED_V1_ANCHORS - eligible_anchors
    real_zero = sum(row["real_background_positive"] == "NO" for row in anchor_rows)
    real_positive = EXPECTED_V1_ANCHORS - real_zero
    affected_transcripts = len({row["transcript"] for row in anchor_rows if row["real_background_positive"] == "NO"})
    removed_candidates = sum(int(row["removed_zero_background_candidates"]) for row in transcript_rows)

    if real_zero != 99:
        raise RuntimeError(f"Mechanical confirmation failed: expected 99 v1 real zero backgrounds, found {real_zero}")
    if any(int(row["candidate_count"]) < int(row["k"]) for row in eligible_pool_rows):
        raise RuntimeError("v2 eligible pool has insufficient capacity")
    if any(row["real_background_positive"] != "YES" for row in anchor_rows if row["v2_coordinate_null_status"] == "COORDINATE_NULL_V2_ELIGIBLE"):
        raise RuntimeError("v2 eligible real anchor has non-positive background")

    docs = args.docs_output_dir
    refs = args.reference_output_dir
    docs.mkdir(parents=True, exist_ok=True)
    refs.mkdir(parents=True, exist_ok=True)
    transcript_path = docs / "STAGE0C_COORDINATE_NULL_V2_TRANSCRIPT_CAPACITY_AUDIT_2026-08-25.tsv"
    anchor_path = docs / "STAGE0C_COORDINATE_NULL_V2_ANCHOR_UNIVERSE_2026-08-25.tsv"
    write_tsv(transcript_path, transcript_rows)
    write_tsv(anchor_path, anchor_rows)

    ref_anchor_path = refs / "coordinate_null_anchor_columns.tsv"
    ref_pool_path = refs / "coordinate_null_candidate_pools.tsv"
    ref_candidate_path = refs / "coordinate_null_candidate_coordinates.npy"
    write_tsv(ref_anchor_path, eligible_anchor_rows)
    write_tsv(ref_pool_path, eligible_pool_rows)
    np.save(ref_candidate_path, eligible_candidates, allow_pickle=False)

    validation = {
        "script_version": SCRIPT_VERSION,
        "policy_freeze_version": POLICY_VERSION,
        "scope": "annotation_only_no_riboseq_no_null_draws",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_change": {
            "candidate_rule_added": "background_length(q) > 0 under symmetric [q-30,q+136) null background",
            "other_v1_candidate_rules_changed": False,
            "null_hypothesis_changed": False,
        },
        "counts": {
            "v1_source_anchors": EXPECTED_V1_ANCHORS,
            "v1_source_transcripts": EXPECTED_V1_TRANSCRIPTS,
            "v1_real_background_positive": real_positive,
            "v1_real_background_zero": real_zero,
            "v1_real_background_zero_affected_transcripts": affected_transcripts,
            "v2_capacity_eligible_transcripts": eligible_transcripts,
            "v2_capacity_eligible_anchors": eligible_anchors,
            "newly_ineligible_transcripts": ineligible_transcripts,
            "newly_ineligible_anchors": ineligible_anchors,
            "zero_background_candidates_removed": removed_candidates,
            "v2_candidate_coordinates": len(eligible_candidates),
        },
        "exclusion_reason_transcript_counts": dict(sorted(reason_counts.items())),
        "invariants": {
            "v1_anchor_count_conservation": eligible_anchors + ineligible_anchors == EXPECTED_V1_ANCHORS,
            "v1_transcript_count_conservation": eligible_transcripts + ineligible_transcripts == EXPECTED_V1_TRANSCRIPTS,
            "whole_transcript_exclusion": True,
            "source_anchor_relative_order_preserved": True,
            "all_v2_real_background_lengths_positive": True,
            "all_v2_candidate_background_lengths_positive": True,
            "fixed_v2_anchor_cardinality_ready": True,
            "realization_specific_filtering_required": False,
            "production_realizations_generated": False,
            "v1_realization_matrix_read": False,
            "riboseq_used": False,
        },
        "randomization_contract_preserved_for_future_v2_production": {
            "realizations": 10_000,
            "seed": 106_330,
            "numpy_version": "2.3.5",
            "prng": "numpy.random.default_rng / PCG64",
            "transcript_order": "lexicographic ascending transcript ID",
            "sampling": "without replacement within transcript",
            "rng_initialization": "once",
            "v1_matrix_reuse_prohibited": True,
        },
        "supersession": {
            "policy": SUPERSEDED_POLICY_VERSION,
            "realizations": SUPERSEDED_REALIZATIONS_VERSION,
            "status": "SUPERSEDED BEFORE BIOLOGICAL UNBLINDING",
            "reason": (
                "The approved symmetric-background log2 density-ratio statistic exposed "
                "structurally undefined zero-length backgrounds in the v1 source/candidate universe."
            ),
        },
        "inputs": {
            "v1_anchor_columns": {"path": str(args.v1_anchor_columns), "sha256": sha256(args.v1_anchor_columns)},
            "v1_candidate_pools": {"path": str(args.v1_candidate_pools), "sha256": sha256(args.v1_candidate_pools)},
            "v1_candidate_coordinates": {"path": str(args.v1_candidate_coordinates), "sha256": sha256(args.v1_candidate_coordinates)},
            "features": {"path": str(args.features), "sha256": sha256(args.features)},
            "topology": {"path": str(args.topology), "sha256": sha256(args.topology)},
        },
        "outputs": {
            "transcript_capacity_audit": {"path": str(transcript_path), "sha256": sha256(transcript_path)},
            "anchor_universe": {"path": str(anchor_path), "sha256": sha256(anchor_path)},
            "v2_anchor_columns": {"path": str(ref_anchor_path), "sha256": sha256(ref_anchor_path)},
            "v2_candidate_pools": {"path": str(ref_pool_path), "sha256": sha256(ref_pool_path)},
            "v2_candidate_coordinates": {
                "path": str(ref_candidate_path),
                "file_sha256": sha256(ref_candidate_path),
                "canonical_int32_content_sha256": sha256_int32(eligible_candidates),
            },
        },
    }
    validation_path = docs / "STAGE0C_COORDINATE_NULL_V2_POLICY_VALIDATION_2026-08-25.json"
    with validation_path.open("w", encoding="utf-8") as handle:
        json.dump(validation, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps({
        "status": "PASS",
        "v1_real_background_positive": real_positive,
        "v1_real_background_zero": real_zero,
        "affected_transcripts": affected_transcripts,
        "v2_capacity_eligible_transcripts": eligible_transcripts,
        "v2_capacity_eligible_anchors": eligible_anchors,
        "newly_ineligible_transcripts": ineligible_transcripts,
        "newly_ineligible_anchors": ineligible_anchors,
        "zero_background_candidates_removed": removed_candidates,
        "v2_candidate_coordinates": len(eligible_candidates),
        "validation_path": str(validation_path),
        "validation_sha256": sha256(validation_path),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
