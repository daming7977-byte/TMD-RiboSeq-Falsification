#!/usr/bin/env python3
"""Generate and exhaustively validate GSE133111 coordinate-null v2 realizations.

This production script consumes only the frozen coordinate-null v2 manifests
and their frozen annotation sources.  It never reads Ribo-seq counts, profile
matrices, sample conditions, benchmark results, or biological statistics.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


SCRIPT_VERSION = "generate_coordinate_null_realizations_v2_2026-08-25"
POLICY_FREEZE_VERSION = "stage0c_coordinate_null_policy_v2_2026-08-25"
PRODUCTION_FREEZE_VERSION = "stage0c_coordinate_null_realizations_v2_2026-08-25"
SUPERSEDED_V1_POLICY = "stage0c_coordinate_null_policy_v1_2026-08-25"
SUPERSEDED_V1_REALIZATIONS = "stage0c_coordinate_null_realizations_v1_2026-08-25"
EXPECTED_NUMPY_VERSION = "2.3.5"
EXPECTED_BIT_GENERATOR = "PCG64"
SEED = 106_330
N_REALIZATIONS = 10_000
BATCH_SIZE = 250
EXPECTED_ANCHORS = 13_060
EXPECTED_TRANSCRIPTS = 3_997
EXPECTED_CANDIDATES = 1_599_733
PRIMARY_START = 60
PRIMARY_END_EXCLUSIVE = 91
NULL_EXCLUSION_START = -30
NULL_EXCLUSION_END_EXCLUSIVE = 136
TOPOLOGY_LABELS = frozenset({"TMhelix", "signal", "Beta sheet"})
MATRIX_DTYPE = np.dtype("<i4")
SENTINEL = np.iinfo(MATRIX_DTYPE).min


class AuditLog:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = path.open("w", encoding="utf-8")

    def write(self, message: str) -> None:
        timestamp = datetime.now(timezone.utc).isoformat()
        line = f"{timestamp}\t{message}"
        print(line, flush=True)
        self.handle.write(line + "\n")
        self.handle.flush()

    def close(self) -> None:
        self.handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-dir", type=Path, required=True)
    parser.add_argument("--policy-validation", type=Path, required=True)
    parser.add_argument("--anchor-universe", type=Path, required=True)
    parser.add_argument("--policy-freeze", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
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


def overlaps(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return start_a < end_b and start_b < end_a


def ensure_numpy_contract() -> dict[str, str]:
    if np.__version__ != EXPECTED_NUMPY_VERSION:
        raise RuntimeError(f"NumPy version mismatch: {np.__version__} != {EXPECTED_NUMPY_VERSION}")
    # Do not initialize the frozen seed here.  Seed 106330 is initialized
    # exactly once for production and exactly once for the independent replay.
    probe = np.random.default_rng()
    bit_generator = type(probe.bit_generator).__name__
    if bit_generator != EXPECTED_BIT_GENERATOR:
        raise RuntimeError(f"Bit generator mismatch: {bit_generator} != {EXPECTED_BIT_GENERATOR}")
    del probe
    return {
        "numpy_version": np.__version__,
        "generator_api": "numpy.random.default_rng",
        "bit_generator": bit_generator,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "byteorder": sys.byteorder,
    }


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
    start = max(0, q + NULL_EXCLUSION_START)
    end = min(len(mask), q + NULL_EXCLUSION_END_EXCLUSIVE)
    mask[start:end] = False
    target_start = max(0, q + PRIMARY_START)
    target_end = min(len(mask), q + PRIMARY_END_EXCLUSIVE)
    mask[target_start:target_end] = False
    return int(np.count_nonzero(mask))


def candidate_geometry(
    candidates: np.ndarray,
    cds_length: int,
    intervals: list[tuple[int, int, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    primary_valid = np.ones(len(candidates), dtype=np.bool_)
    topology_valid = np.ones(len(candidates), dtype=np.bool_)
    codon_valid = candidates % 3 == 0
    base_mask = topology_background_mask(cds_length, intervals)
    backgrounds = np.empty(len(candidates), dtype=np.int64)
    for index, value in enumerate(candidates.tolist()):
        q = int(value)
        primary_start = q + PRIMARY_START
        primary_end = q + PRIMARY_END_EXCLUSIVE
        primary_valid[index] = q >= 0 and primary_end <= cds_length
        topology_valid[index] = not any(start <= q < end for start, end, _label in intervals)
        if topology_valid[index]:
            topology_valid[index] = not any(
                overlaps(primary_start, primary_end, start, end)
                for start, end, _label in intervals
            )
        backgrounds[index] = null_background_length(base_mask, q)
    return primary_valid, topology_valid, codon_valid, backgrounds


def draw_transcript_production(
    rng: np.random.Generator,
    candidates: np.ndarray,
    k: int,
) -> Iterable[tuple[int, int, np.ndarray, np.ndarray]]:
    """Exact frozen v1/predecessor RNG-call semantics for production."""
    candidate_count = len(candidates)
    for start in range(0, N_REALIZATIONS, BATCH_SIZE):
        stop = min(N_REALIZATIONS, start + BATCH_SIZE)
        size = stop - start
        if k == 1:
            chosen = rng.integers(0, candidate_count, size=(size, 1))
        else:
            chosen = np.empty((size, k), dtype=np.int64)
            for row in range(size):
                chosen[row] = rng.choice(candidate_count, size=k, replace=False)
        selected = np.asarray(candidates[chosen], dtype=MATRIX_DTYPE)
        selected.sort(axis=1)
        yield start, stop, chosen, selected


def draw_transcript_replay(
    rng: np.random.Generator,
    candidates: np.ndarray,
    k: int,
) -> Iterable[tuple[int, int, np.ndarray]]:
    """Independent replay implementation of the frozen RNG-call contract."""
    n_candidates = int(candidates.shape[0])
    batch_start = 0
    while batch_start < N_REALIZATIONS:
        batch_stop = min(N_REALIZATIONS, batch_start + BATCH_SIZE)
        batch_size = batch_stop - batch_start
        if k == 1:
            sampled_indices = rng.integers(0, n_candidates, size=(batch_size, 1))
        else:
            sampled_indices = np.empty((batch_size, k), dtype=np.int64)
            for realization_offset in range(batch_size):
                sampled_indices[realization_offset, :] = rng.choice(
                    n_candidates, size=k, replace=False
                )
        replayed = np.asarray(candidates[sampled_indices], dtype=MATRIX_DTYPE)
        replayed.sort(axis=1)
        yield batch_start, batch_stop, replayed
        batch_start = batch_stop


def canonical_matrix_content_sha256(matrix: np.ndarray) -> str:
    digest = hashlib.sha256()
    for start in range(0, N_REALIZATIONS, 50):
        stop = min(N_REALIZATIONS, start + 50)
        block = np.ascontiguousarray(matrix[start:stop, :], dtype=MATRIX_DTYPE)
        digest.update(block.tobytes(order="C"))
    return digest.hexdigest()


def validate_frozen_inputs(
    args: argparse.Namespace,
    policy_validation: dict[str, Any],
    pool_rows: list[dict[str, str]],
    candidates: np.ndarray,
    anchors: list[dict[str, str]],
    anchor_universe: list[dict[str, str]],
    features: dict[str, int],
    topology: dict[str, list[tuple[int, int, str]]],
) -> dict[str, Any]:
    violations = Counter()
    if policy_validation["policy_freeze_version"] != POLICY_FREEZE_VERSION:
        raise RuntimeError("Policy validation version mismatch")
    frozen_counts = policy_validation["counts"]
    if int(frozen_counts["v2_capacity_eligible_transcripts"]) != EXPECTED_TRANSCRIPTS:
        raise RuntimeError("Frozen transcript count mismatch")
    if int(frozen_counts["v2_capacity_eligible_anchors"]) != EXPECTED_ANCHORS:
        raise RuntimeError("Frozen anchor count mismatch")
    if int(frozen_counts["v2_candidate_coordinates"]) != EXPECTED_CANDIDATES:
        raise RuntimeError("Frozen candidate count mismatch")

    expected_outputs = policy_validation["outputs"]
    file_checks = {
        "v2_anchor_columns": args.manifest_dir / "coordinate_null_anchor_columns.tsv",
        "v2_candidate_pools": args.manifest_dir / "coordinate_null_candidate_pools.tsv",
    }
    for key, path in file_checks.items():
        if sha256(path) != expected_outputs[key]["sha256"]:
            raise RuntimeError(f"Frozen v2 manifest hash mismatch: {path}")
    candidate_path = args.manifest_dir / "coordinate_null_candidate_coordinates.npy"
    if sha256(candidate_path) != expected_outputs["v2_candidate_coordinates"]["file_sha256"]:
        raise RuntimeError("Frozen v2 candidate-coordinate file hash mismatch")
    if sha256_int32(candidates) != expected_outputs["v2_candidate_coordinates"]["canonical_int32_content_sha256"]:
        raise RuntimeError("Frozen v2 candidate-coordinate content hash mismatch")
    if sha256(args.anchor_universe) != expected_outputs["anchor_universe"]["sha256"]:
        raise RuntimeError("Frozen v2 anchor-universe hash mismatch")
    for key, path in (("features", args.features), ("topology", args.topology)):
        if sha256(path) != policy_validation["inputs"][key]["sha256"]:
            raise RuntimeError(f"Frozen annotation hash mismatch: {key}")

    if len(pool_rows) != EXPECTED_TRANSCRIPTS or len(anchors) != EXPECTED_ANCHORS:
        raise RuntimeError("Frozen v2 manifest count mismatch")
    if candidates.ndim != 1 or candidates.dtype != MATRIX_DTYPE or len(candidates) != EXPECTED_CANDIDATES:
        raise RuntimeError(f"Candidate array contract mismatch: {candidates.shape} {candidates.dtype}")

    transcripts = [row["transcript"] for row in pool_rows]
    if transcripts != sorted(transcripts) or len(set(transcripts)) != EXPECTED_TRANSCRIPTS:
        violations["transcript_order_or_uniqueness"] += 1
    anchor_indices = [int(row["anchor_column_index0"]) for row in anchors]
    if anchor_indices != list(range(EXPECTED_ANCHORS)):
        violations["anchor_column_order"] += 1
    anchor_keys = [(row["transcript"], row["TMD_index"]) for row in anchors]
    if len(set(anchor_keys)) != EXPECTED_ANCHORS:
        violations["source_anchor_uniqueness"] += 1
    v1_indices = [int(row["v1_anchor_column_index0"]) for row in anchors]
    if any(left >= right for left, right in zip(v1_indices, v1_indices[1:])):
        violations["v1_relative_source_order"] += 1

    eligible_universe_rows = [
        row for row in anchor_universe
        if row["v2_coordinate_null_status"] == "COORDINATE_NULL_V2_ELIGIBLE"
    ]
    eligible_universe = {
        (row["transcript"], row["TMD_index"]): row
        for row in eligible_universe_rows
    }
    if len(eligible_universe_rows) != EXPECTED_ANCHORS:
        violations["eligible_anchor_universe_row_count"] += abs(
            len(eligible_universe_rows) - EXPECTED_ANCHORS
        ) or 1
    if len(eligible_universe) != EXPECTED_ANCHORS:
        violations["eligible_anchor_universe_unique_count"] += abs(
            len(eligible_universe) - EXPECTED_ANCHORS
        ) or 1
    if set(anchor_keys) != set(eligible_universe):
        violations["v1_only_or_missing_anchor"] += len(set(anchor_keys) ^ set(eligible_universe))
    if any(row["real_background_positive"] != "YES" for row in eligible_universe.values()):
        violations["recorded_real_background_nonpositive"] += 1

    expected_anchor_start = 0
    expected_candidate_start = 0
    all_candidate_background_positive = True
    for row in pool_rows:
        transcript = row["transcript"]
        k = int(row["k"])
        a0 = int(row["anchor_column_start0"])
        a1 = int(row["anchor_column_end0_exclusive"])
        c0 = int(row["candidate_offset_start0"])
        c1 = int(row["candidate_offset_end0_exclusive"])
        if a0 != expected_anchor_start or c0 != expected_candidate_start:
            violations["manifest_offset_contiguity"] += 1
        expected_anchor_start = a1
        expected_candidate_start = c1
        if a1 - a0 != k:
            violations["anchor_block_k"] += 1
        values = np.asarray(candidates[c0:c1], dtype=MATRIX_DTYPE)
        if len(values) != int(row["candidate_count"]):
            violations["candidate_count_manifest"] += 1
        if len(values) < k:
            violations["candidate_capacity"] += 1
        if len(values) == 0 or np.any(np.diff(values) <= 0):
            violations["candidate_order_or_uniqueness"] += 1
        if sha256_int32(values) != row["candidate_coordinates_sha256"]:
            violations["candidate_vector_hash"] += 1
        block_anchors = anchors[a0:a1]
        if any(anchor["transcript"] != transcript for anchor in block_anchors):
            violations["anchor_transcript_block"] += 1
        observed_tmds = [int(anchor["TMD_index"]) for anchor in block_anchors]
        expected_tmds = [int(value) for value in row["source_TMD_indices_in_order"].split(",")]
        if observed_tmds != expected_tmds:
            violations["source_anchor_assignment_order"] += 1

        cds_length = features[transcript]
        intervals = topology.get(transcript, [])
        primary_valid, topology_valid, codon_valid, backgrounds = candidate_geometry(
            values, cds_length, intervals
        )
        violations["candidate_primary_window_validity"] += int(np.count_nonzero(~primary_valid))
        violations["candidate_topology_exclusion"] += int(np.count_nonzero(~topology_valid))
        violations["candidate_codon_alignment"] += int(np.count_nonzero(~codon_valid))
        zero_backgrounds = int(np.count_nonzero(backgrounds <= 0))
        violations["candidate_zero_background"] += zero_backgrounds
        all_candidate_background_positive &= zero_backgrounds == 0

        base_mask = topology_background_mask(cds_length, intervals)
        for anchor in block_anchors:
            key = (transcript, anchor["TMD_index"])
            recomputed = real_background_length(base_mask, int(anchor["real_anchor_cds0"]))
            if recomputed <= 0:
                violations["real_zero_background"] += 1
            recorded = int(eligible_universe[key]["real_background_length_nt"])
            if recomputed != recorded:
                violations["real_background_length_mismatch"] += 1

    if expected_anchor_start != EXPECTED_ANCHORS:
        violations["final_anchor_offset"] += 1
    if expected_candidate_start != EXPECTED_CANDIDATES:
        violations["final_candidate_offset"] += 1
    violations_to_report = dict(sorted(violations.items()))
    if sum(violations_to_report.values()):
        raise RuntimeError(f"Frozen v2 input validation failed: {violations_to_report}")
    return {
        "unique_transcripts": len(set(transcripts)),
        "unique_source_anchors": len(set(anchor_keys)),
        "candidate_coordinates": len(candidates),
        "all_real_background_positive": True,
        "all_candidate_background_positive": all_candidate_background_positive,
        "candidate_count_ge_k_all_transcripts": True,
        "v1_only_ineligible_entries": 0,
        "violations": violations_to_report,
    }


def main() -> None:
    args = parse_args()
    manifest_dir = args.manifest_dir
    required = (
        manifest_dir / "coordinate_null_anchor_columns.tsv",
        manifest_dir / "coordinate_null_candidate_pools.tsv",
        manifest_dir / "coordinate_null_candidate_coordinates.npy",
        args.policy_validation,
        args.anchor_universe,
        args.policy_freeze,
        args.topology,
        args.features,
    )
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    final_matrix_path = manifest_dir / "coordinate_null_coordinates.npy"
    inprogress_matrix_path = manifest_dir / "coordinate_null_coordinates.inprogress.npy"
    validation_path = manifest_dir / "coordinate_null_generation_validation.json"
    for path in (final_matrix_path, inprogress_matrix_path, validation_path):
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite existing v2 production artifact: {path}")

    log = AuditLog(manifest_dir / "coordinate_null_generation.log")
    try:
        started = datetime.now(timezone.utc).isoformat()
        runtime = ensure_numpy_contract()
        log.write("V2_PRODUCTION_START")
        policy_validation = json.loads(args.policy_validation.read_text(encoding="utf-8"))
        pool_rows = read_tsv(manifest_dir / "coordinate_null_candidate_pools.tsv")
        anchors = read_tsv(manifest_dir / "coordinate_null_anchor_columns.tsv")
        candidates = np.load(
            manifest_dir / "coordinate_null_candidate_coordinates.npy",
            mmap_mode="r",
            allow_pickle=False,
        )
        anchor_universe = read_tsv(args.anchor_universe)
        features = load_features(args.features)
        topology = load_topology(args.topology)
        input_validation = validate_frozen_inputs(
            args,
            policy_validation,
            pool_rows,
            candidates,
            anchors,
            anchor_universe,
            features,
            topology,
        )
        log.write("V2_INPUT_MANIFEST_VALIDATION_PASS")

        realization_rows = [
            {"realization_row_index0": index, "realization_id": index + 1}
            for index in range(N_REALIZATIONS)
        ]
        realization_rows_path = manifest_dir / "coordinate_null_realization_rows.tsv"
        write_tsv(realization_rows_path, realization_rows)

        randomization_contract = {
            "numpy_version": EXPECTED_NUMPY_VERSION,
            "generator_api": "numpy.random.default_rng",
            "bit_generator": EXPECTED_BIT_GENERATOR,
            "seed": SEED,
            "production_seed_initializations": 1,
            "independent_replay_seed_initializations": 1,
            "realizations": N_REALIZATIONS,
            "realization_order": "1..10000 ascending",
            "batch_size": BATCH_SIZE,
            "transcript_order": "lexicographic ascending transcript ID",
            "candidate_order": "ascending 0-based CDS coordinate",
            "sampling": "without replacement within transcript",
            "source_anchor_order": "frozen source TMD order within transcript",
            "assignment": "sort sampled coordinates then assign to frozen source-anchor order",
            "k_eq_1": "rng.integers(0,candidate_count,size=(batch_size,1))",
            "k_gt_1": "one rng.choice(candidate_count,size=k,replace=False) per realization",
            "candidate_count_eq_k": "still call rng.choice; no shortcut",
        }
        preflight = {
            "script_version": SCRIPT_VERSION,
            "policy_freeze_version": POLICY_FREEZE_VERSION,
            "production_freeze_version": PRODUCTION_FREEZE_VERSION,
            "scope": "annotation_only_no_riboseq_no_biological_statistics",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "runtime": runtime,
            "randomization_contract": randomization_contract,
            "input_validation": input_validation,
            "inputs": {
                "anchor_columns": {
                    "path": str(manifest_dir / "coordinate_null_anchor_columns.tsv"),
                    "sha256": sha256(manifest_dir / "coordinate_null_anchor_columns.tsv"),
                },
                "candidate_pools": {
                    "path": str(manifest_dir / "coordinate_null_candidate_pools.tsv"),
                    "sha256": sha256(manifest_dir / "coordinate_null_candidate_pools.tsv"),
                },
                "candidate_coordinates": {
                    "path": str(manifest_dir / "coordinate_null_candidate_coordinates.npy"),
                    "file_sha256": sha256(manifest_dir / "coordinate_null_candidate_coordinates.npy"),
                    "content_sha256": sha256_int32(candidates),
                },
                "anchor_universe": {"path": str(args.anchor_universe), "sha256": sha256(args.anchor_universe)},
                "policy_validation": {"path": str(args.policy_validation), "sha256": sha256(args.policy_validation)},
                "policy_freeze": {"path": str(args.policy_freeze), "sha256": sha256(args.policy_freeze)},
                "topology": {"path": str(args.topology), "sha256": sha256(args.topology)},
                "features": {"path": str(args.features), "sha256": sha256(args.features)},
            },
            "v1_status": "SUPERSEDED BEFORE BIOLOGICAL UNBLINDING",
            "v1_matrix_read_or_reused": False,
            "production_matrix_generated": False,
            "riboseq_used": False,
        }
        preflight_path = manifest_dir / "coordinate_null_preflight_v2.json"
        with preflight_path.open("w", encoding="utf-8") as handle:
            json.dump(preflight, handle, indent=2, sort_keys=True)
            handle.write("\n")
        log.write(f"V2_PREFLIGHT_PASS sha256={sha256(preflight_path)}")

        matrix = np.lib.format.open_memmap(
            inprogress_matrix_path,
            mode="w+",
            dtype=MATRIX_DTYPE,
            shape=(N_REALIZATIONS, EXPECTED_ANCHORS),
            fortran_order=True,
        )
        matrix[:] = SENTINEL
        matrix.flush()
        rng = np.random.default_rng(SEED)
        production_transcript_hash = hashlib.sha256()
        frequency_rows: list[dict[str, Any]] = []
        runtime_sampling_violations = Counter()
        last_report = time.monotonic()

        for transcript_index, row in enumerate(pool_rows):
            c0 = int(row["candidate_offset_start0"])
            c1 = int(row["candidate_offset_end0_exclusive"])
            a0 = int(row["anchor_column_start0"])
            a1 = int(row["anchor_column_end0_exclusive"])
            k = int(row["k"])
            transcript_candidates = np.asarray(candidates[c0:c1], dtype=MATRIX_DTYPE)
            frequency = np.zeros(len(transcript_candidates), dtype=np.int64)
            for start, stop, chosen, selected in draw_transcript_production(
                rng, transcript_candidates, k
            ):
                matrix[start:stop, a0:a1] = selected
                production_transcript_hash.update(
                    np.ascontiguousarray(selected).tobytes(order="C")
                )
                frequency += np.bincount(chosen.ravel(), minlength=len(transcript_candidates))
                if k > 1:
                    runtime_sampling_violations["within_transcript_reuse"] += int(
                        np.count_nonzero(np.diff(selected, axis=1) == 0)
                    )
                    runtime_sampling_violations["order_assignment"] += int(
                        np.count_nonzero(np.diff(selected, axis=1) <= 0)
                    )
                runtime_sampling_violations["codon_alignment"] += int(
                    np.count_nonzero(selected % 3)
                )
            expected_draws = N_REALIZATIONS * k
            if int(frequency.sum()) != expected_draws:
                runtime_sampling_violations["frequency_total"] += 1
            frequency_rows.append({
                "transcript": row["transcript"],
                "k": k,
                "candidate_count": len(transcript_candidates),
                "total_draws": int(frequency.sum()),
                "expected_total_draws": expected_draws,
                "candidate_frequency_min": int(frequency.min()),
                "candidate_frequency_max": int(frequency.max()),
                "candidate_frequency_mean": float(frequency.mean()),
                "candidate_coordinates_never_sampled": int(np.count_nonzero(frequency == 0)),
                "frequency_vector_sha256": sha256_int32(frequency.astype(MATRIX_DTYPE)),
                "scope": "ANNOTATION_ONLY_RANDOMIZATION_QC",
            })
            now = time.monotonic()
            if now - last_report >= 25 or transcript_index + 1 == len(pool_rows):
                log.write(f"GENERATED_TRANSCRIPTS {transcript_index + 1}/{len(pool_rows)}")
                last_report = now

        matrix.flush()
        del matrix
        matrix = np.load(inprogress_matrix_path, mmap_mode="r", allow_pickle=False)
        if matrix.shape != (N_REALIZATIONS, EXPECTED_ANCHORS) or matrix.dtype != MATRIX_DTYPE:
            raise RuntimeError(f"Matrix contract failed: shape={matrix.shape} dtype={matrix.dtype}")
        if not bool(matrix.flags.f_contiguous):
            raise RuntimeError("Matrix is not stored in frozen Fortran-order NPY layout")

        per_assigned = np.zeros(N_REALIZATIONS, dtype=np.int64)
        per_transcripts = np.zeros(N_REALIZATIONS, dtype=np.int64)
        per_membership = np.zeros(N_REALIZATIONS, dtype=np.int64)
        per_reuse = np.zeros(N_REALIZATIONS, dtype=np.int64)
        per_primary = np.zeros(N_REALIZATIONS, dtype=np.int64)
        per_topology = np.zeros(N_REALIZATIONS, dtype=np.int64)
        per_zero_background = np.zeros(N_REALIZATIONS, dtype=np.int64)
        per_order = np.zeros(N_REALIZATIONS, dtype=np.int64)
        per_codon = np.zeros(N_REALIZATIONS, dtype=np.int64)

        for row in pool_rows:
            c0 = int(row["candidate_offset_start0"])
            c1 = int(row["candidate_offset_end0_exclusive"])
            a0 = int(row["anchor_column_start0"])
            a1 = int(row["anchor_column_end0_exclusive"])
            transcript_candidates = np.asarray(candidates[c0:c1], dtype=MATRIX_DTYPE)
            block = np.asarray(matrix[:, a0:a1])
            present = block != SENTINEL
            per_assigned += np.count_nonzero(present, axis=1)
            per_transcripts += np.all(present, axis=1).astype(np.int64)
            found = np.searchsorted(transcript_candidates, block)
            clipped = np.clip(found, 0, len(transcript_candidates) - 1)
            membership = (found < len(transcript_candidates)) & (
                transcript_candidates[clipped] == block
            )
            per_membership += np.count_nonzero(~membership, axis=1)
            if block.shape[1] > 1:
                differences = np.diff(block, axis=1)
                per_reuse += np.count_nonzero(differences == 0, axis=1)
                per_order += np.count_nonzero(differences <= 0, axis=1)
            per_codon += np.count_nonzero(block % 3 != 0, axis=1)

            cds_length = features[row["transcript"]]
            intervals = topology.get(row["transcript"], [])
            primary_valid, topology_valid, codon_valid, backgrounds = candidate_geometry(
                transcript_candidates, cds_length, intervals
            )
            valid_lookup = membership
            per_primary += np.count_nonzero(valid_lookup & ~primary_valid[clipped], axis=1)
            per_topology += np.count_nonzero(valid_lookup & ~topology_valid[clipped], axis=1)
            per_codon += np.count_nonzero(valid_lookup & ~codon_valid[clipped], axis=1)
            per_zero_background += np.count_nonzero(
                valid_lookup & (backgrounds[clipped] <= 0), axis=1
            )

        per_realization_rows: list[dict[str, Any]] = []
        for index in range(N_REALIZATIONS):
            per_realization_rows.append({
                "realization_row_index0": index,
                "realization_id": index + 1,
                "assigned_coordinates": int(per_assigned[index]),
                "represented_transcripts": int(per_transcripts[index]),
                "candidate_membership_violations": int(per_membership[index]),
                "within_transcript_reuse_violations": int(per_reuse[index]),
                "primary_window_validity_violations": int(per_primary[index]),
                "topology_exclusion_violations": int(per_topology[index]),
                "zero_background_violations": int(per_zero_background[index]),
                "order_assignment_violations": int(per_order[index]),
                "codon_alignment_violations": int(per_codon[index]),
                "status": "PASS" if (
                    per_assigned[index] == EXPECTED_ANCHORS
                    and per_transcripts[index] == EXPECTED_TRANSCRIPTS
                    and per_membership[index] == 0
                    and per_reuse[index] == 0
                    and per_primary[index] == 0
                    and per_topology[index] == 0
                    and per_zero_background[index] == 0
                    and per_order[index] == 0
                    and per_codon[index] == 0
                ) else "FAIL",
            })

        violations = {
            "candidate_membership": int(per_membership.sum()),
            "within_transcript_reuse": int(per_reuse.sum()),
            "primary_window_validity": int(per_primary.sum()),
            "topology_exclusion": int(per_topology.sum()),
            "zero_background": int(per_zero_background.sum()),
            "order_assignment": int(per_order.sum()),
            "codon_alignment": int(per_codon.sum()),
            "assigned_coordinate_count": int(np.count_nonzero(per_assigned != EXPECTED_ANCHORS)),
            "represented_transcript_count": int(np.count_nonzero(per_transcripts != EXPECTED_TRANSCRIPTS)),
            "runtime_within_transcript_reuse": int(runtime_sampling_violations["within_transcript_reuse"]),
            "runtime_order_assignment": int(runtime_sampling_violations["order_assignment"]),
            "runtime_codon_alignment": int(runtime_sampling_violations["codon_alignment"]),
            "runtime_frequency_total": int(runtime_sampling_violations["frequency_total"]),
            "failed_realizations": sum(row["status"] != "PASS" for row in per_realization_rows),
        }
        log.write("EXHAUSTIVE_PER_REALIZATION_VALIDATION_COMPLETE")

        reproduction_mismatches = 0
        replay_transcript_hash = hashlib.sha256()
        replay_rng = np.random.default_rng(SEED)
        last_report = time.monotonic()
        for transcript_index, row in enumerate(pool_rows):
            c0 = int(row["candidate_offset_start0"])
            c1 = int(row["candidate_offset_end0_exclusive"])
            a0 = int(row["anchor_column_start0"])
            a1 = int(row["anchor_column_end0_exclusive"])
            k = int(row["k"])
            transcript_candidates = np.asarray(candidates[c0:c1], dtype=MATRIX_DTYPE)
            for start, stop, replayed in draw_transcript_replay(
                replay_rng, transcript_candidates, k
            ):
                replay_transcript_hash.update(
                    np.ascontiguousarray(replayed).tobytes(order="C")
                )
                reproduction_mismatches += int(
                    np.count_nonzero(matrix[start:stop, a0:a1] != replayed)
                )
            now = time.monotonic()
            if now - last_report >= 25 or transcript_index + 1 == len(pool_rows):
                log.write(f"REPLAYED_TRANSCRIPTS {transcript_index + 1}/{len(pool_rows)}")
                last_report = now

        canonical_content_sha = canonical_matrix_content_sha256(matrix)
        matrix_file_sha = sha256(inprogress_matrix_path)
        production_transcript_sha = production_transcript_hash.hexdigest()
        replay_transcript_sha = replay_transcript_hash.hexdigest()
        if production_transcript_sha != replay_transcript_sha:
            reproduction_mismatches += 1

        frequency_path = manifest_dir / "coordinate_null_randomization_qc.tsv"
        per_realization_path = manifest_dir / "coordinate_null_per_realization_validation.tsv"
        write_tsv(frequency_path, frequency_rows)
        write_tsv(per_realization_path, per_realization_rows)

        if sum(violations.values()) or reproduction_mismatches:
            raise RuntimeError(
                f"V2 production validation failed: violations={violations} "
                f"reproduction_mismatches={reproduction_mismatches}"
            )

        del matrix
        inprogress_matrix_path.rename(final_matrix_path)
        completed = datetime.now(timezone.utc).isoformat()
        validation = {
            "script_version": SCRIPT_VERSION,
            "production_freeze_version": PRODUCTION_FREEZE_VERSION,
            "policy_freeze_version": POLICY_FREEZE_VERSION,
            "scope": "annotation_only_no_riboseq_no_biological_statistics",
            "started_utc": started,
            "completed_utc": completed,
            "runtime": runtime,
            "randomization_contract": randomization_contract,
            "input_validation": input_validation,
            "counts": {
                "transcripts": EXPECTED_TRANSCRIPTS,
                "anchors_per_realization": EXPECTED_ANCHORS,
                "realizations": N_REALIZATIONS,
                "total_null_assignments": N_REALIZATIONS * EXPECTED_ANCHORS,
                "candidate_coordinates": EXPECTED_CANDIDATES,
                "libraries_for_future_statistic": 6,
            },
            "storage": {
                "path": str(final_matrix_path),
                "shape": [N_REALIZATIONS, EXPECTED_ANCHORS],
                "dtype": MATRIX_DTYPE.str,
                "npy_fortran_order": True,
                "matrix_file_sha256": matrix_file_sha,
                "canonical_logical_C_order_int32_content_sha256": canonical_content_sha,
                "production_transcript_major_sha256": production_transcript_sha,
                "independent_replay_transcript_major_sha256": replay_transcript_sha,
            },
            "violations": violations,
            "independent_replay": {
                "status": "PASS",
                "coordinate_assignment_mismatches": reproduction_mismatches,
                "all_assignments_compared": True,
            },
            "fixed_statistic_readiness": {
                "status": "PASS",
                "all_real_source_backgrounds_positive": True,
                "all_null_assignment_backgrounds_positive": True,
                "fixed_anchor_cardinality": EXPECTED_ANCHORS,
                "libraries": 6,
                "realization_specific_filtering_required": False,
                "S_real_calculated": False,
                "S_null_calculated": False,
            },
            "files": {},
            "v1": {
                "policy_version": SUPERSEDED_V1_POLICY,
                "realizations_version": SUPERSEDED_V1_REALIZATIONS,
                "status": "SUPERSEDED BEFORE BIOLOGICAL UNBLINDING",
                "artifacts_modified_or_deleted": False,
                "matrix_read_or_reused": False,
            },
            "invariants": {
                "all_validation_violations_zero": True,
                "all_10000_realizations_pass": True,
                "independent_replay_identical": True,
                "fixed_observation_universe_ready": True,
                "riboseq_used": False,
                "biological_profile_results_inspected": False,
                "biological_statistics_calculated": False,
                "stage1_unblinding_allowed": False,
            },
        }
        log.write(
            "V2_PRODUCTION_PASS "
            f"canonical_content_sha256={canonical_content_sha} "
            f"matrix_file_sha256={matrix_file_sha} reproduction_mismatches=0"
        )
        for name in (
            "coordinate_null_anchor_columns.tsv",
            "coordinate_null_candidate_pools.tsv",
            "coordinate_null_candidate_coordinates.npy",
            "coordinate_null_realization_rows.tsv",
            "coordinate_null_preflight_v2.json",
            "coordinate_null_randomization_qc.tsv",
            "coordinate_null_per_realization_validation.tsv",
            "coordinate_null_generation.log",
            "coordinate_null_coordinates.npy",
        ):
            path = manifest_dir / name
            validation["files"][name] = {"path": str(path), "sha256": sha256(path)}
        with validation_path.open("w", encoding="utf-8") as handle:
            json.dump(validation, handle, indent=2, sort_keys=True)
            handle.write("\n")
    finally:
        log.close()


if __name__ == "__main__":
    main()
