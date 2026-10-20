#!/usr/bin/env python3
"""Frozen Paper 3 Stage 1 counting, inference, and robustness backend.

This module contains no command-line entry point.  It is called only by the
Stage 1 adapter after the canonical provenance preflight succeeds.
"""

from __future__ import annotations

import bisect
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


BACKEND_VERSION = "paper3_stage1_backend_v1.0.0_2026-08-26"
CONTRACT_VERSION = "stage0c_paper3_analysis_contract_v1_2026-08-25"
SAMPLE_ORDER = (
    "GSM3900208",
    "GSM3900209",
    "GSM3900210",
    "GSM3900211",
    "GSM3900212",
    "GSM3900213",
)
SAMPLE_LABELS = {
    "GSM3900208": "NSC",
    "GSM3900209": "NSC",
    "GSM3900210": "NSC",
    "GSM3900211": "siP1_1",
    "GSM3900212": "siP2_1",
    "GSM3900213": "siP2_4",
}
CONDITIONS = {
    gsm: ("Control" if index < 3 else "Knockdown")
    for index, gsm in enumerate(SAMPLE_ORDER)
}
FROZEN_READ_LENGTHS = (27, 28, 29, 30, 32, 33, 34, 35)
MANDATORY_SENSITIVITY_EXCLUDED_LENGTHS = (28, 29)
CENTERS = (-30, -15, 0, 15, 30, 45, 60, 75, 90, 105, 120)
PRIMARY_CENTER = 75
TARGET_LEFT = -15
TARGET_RIGHT = 16
TARGET_LENGTH = 31
MINIMUM_BACKGROUND_COUNT = 20
PRIMARY_STARTING_ANCHORS = 15_597
PSEUDO_STARTING_PAIRS = 9_829
COORDINATE_NULL_VERSION = "V2"
COORDINATE_NULL_ANCHORS = 13_060
COORDINATE_NULL_TRANSCRIPTS = 3_997
COORDINATE_NULL_REALIZATIONS = 10_000
COORDINATE_NULL_SHA256 = (
    "7983d6dc1c382e57521e676669891a03271a4f1ffda0dda8d9ecf101e01071b3"
)
LDB_BAND = (60, 75, 90)
NGS_EPSILON = math.log2(1.20)
TOPOLOGY_MASK_CLASSES = frozenset(("TMhelix", "signal", "Beta sheet"))


@dataclass(frozen=True)
class Feature:
    transcript: str
    utr5_length: int
    cds_length: int
    utr3_length: int


@dataclass(frozen=True)
class Anchor:
    anchor_id: str
    transcript: str
    tmd_index: str
    anchor_cds0: int
    segment_start0: int
    segment_end0: int


@dataclass(frozen=True)
class PseudoPair:
    pair_id: str
    transcript: str
    tmd_index: str
    real_anchor_cds0: int
    real_segment_start0: int
    real_segment_end0: int
    pseudo_anchor_cds0: int
    pseudo_segment_start0: int
    pseudo_segment_end0: int


@dataclass(frozen=True)
class SparseTrack:
    positions: tuple[int, ...]
    prefix_counts: tuple[int, ...]

    @classmethod
    def from_counter(cls, counts: Mapping[int, int]) -> "SparseTrack":
        positions = tuple(sorted(counts))
        running = [0]
        for position in positions:
            running.append(running[-1] + int(counts[position]))
        return cls(positions, tuple(running))

    @property
    def total(self) -> int:
        return self.prefix_counts[-1]

    def count(self, start: int, end: int) -> int:
        if end <= start:
            return 0
        left = bisect.bisect_left(self.positions, start)
        right = bisect.bisect_left(self.positions, end)
        return self.prefix_counts[right] - self.prefix_counts[left]


EMPTY_TRACK = SparseTrack((), (0,))


class IneligibleFeatureGeometry(ValueError):
    """A structurally invalid feature that must be excluded from inference."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_equal(actual: object, expected: object, label: str) -> None:
    if actual != expected:
        raise ValueError(f"{label}: expected {expected!r}, observed {actual!r}")


def merge_intervals(
    intervals: Iterable[tuple[int, int]],
    length: int,
) -> tuple[tuple[int, int], ...]:
    clipped = sorted(
        (max(0, start), min(length, end))
        for start, end in intervals
        if max(0, start) < min(length, end)
    )
    merged: list[tuple[int, int]] = []
    for start, end in clipped:
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return tuple(merged)


def target_interval(anchor_cds0: int, center: int = PRIMARY_CENTER) -> tuple[int, int]:
    return (
        anchor_cds0 + center + TARGET_LEFT,
        anchor_cds0 + center + TARGET_RIGHT,
    )


def focal_interval(start0: int, end0: int, anchor_cds0: int) -> tuple[int, int]:
    return start0 - 30, max(end0, anchor_cds0 + 136)


def condition_design(
    sample_order: Sequence[str] = SAMPLE_ORDER,
) -> tuple[np.ndarray, np.ndarray]:
    if tuple(sample_order) != SAMPLE_ORDER:
        raise ValueError("Paper 3 condition design requires the frozen six-library order")
    design = np.array(
        [[1.0, 0.0] if CONDITIONS[gsm] == "Control" else [1.0, 1.0]
         for gsm in sample_order],
        dtype=float,
    )
    contrast = np.array([0.0, 1.0], dtype=float)
    if np.linalg.matrix_rank(design) != 2:
        raise ValueError("Frozen ~ condition design is rank deficient")
    return design, contrast


def exposure(
    background_count: int,
    background_length: int,
    target_length: int = TARGET_LENGTH,
) -> float:
    if background_count < 0 or background_length <= 0:
        raise ValueError("Invalid background count/length for exposure")
    if target_length != TARGET_LENGTH:
        raise ValueError("Primary target length must remain 31 nt")
    return (background_count + 0.5) * target_length / background_length


def profile_counts(
    track: SparseTrack,
    cds_length: int,
    ordinary_masks: Sequence[tuple[int, int]],
    anchor_cds0: int,
    center: int = PRIMARY_CENTER,
    additional_focal: tuple[int, int] | None = None,
) -> tuple[int, int, int]:
    target_start, target_end = target_interval(anchor_cds0, center)
    if target_start < 0 or target_end > cds_length:
        raise IneligibleFeatureGeometry("Target window is outside the CDS")
    masks = list(ordinary_masks)
    if additional_focal is not None:
        masks.append(additional_focal)
    masks.append((target_start, target_end))
    merged = merge_intervals(masks, cds_length)
    background_length = cds_length - sum(end - start for start, end in merged)
    background_count = track.total - sum(track.count(start, end) for start, end in merged)
    target_count = track.count(target_start, target_end)
    if background_length < 0:
        raise ValueError("Background geometry produced a negative length")
    if background_count < 0:
        raise ValueError("Background subtraction produced a negative count")
    return target_count, background_count, background_length


def strict_six_library_eligible(rows: Sequence[Mapping[str, object]]) -> bool:
    if tuple(str(row["GSM"]) for row in rows) != SAMPLE_ORDER:
        return False
    if any(int(row["background_count"]) < MINIMUM_BACKGROUND_COUNT for row in rows):
        return False
    if any(int(row["background_length"]) <= 0 for row in rows):
        return False
    return sum(int(row["target_count"]) for row in rows) > 0


def eligible_for_samples(
    rows: Sequence[Mapping[str, object]],
    sample_order: Sequence[str],
) -> bool:
    if tuple(str(row["GSM"]) for row in rows) != tuple(sample_order):
        return False
    if any(int(row["background_count"]) < MINIMUM_BACKGROUND_COUNT for row in rows):
        return False
    if any(int(row["background_length"]) <= 0 for row in rows):
        return False
    return sum(int(row["target_count"]) for row in rows) > 0


def paired_eligibility(
    real_rows: Sequence[Mapping[str, object]],
    pseudo_rows: Sequence[Mapping[str, object]],
) -> bool:
    return strict_six_library_eligible(real_rows) and strict_six_library_eligible(
        pseudo_rows
    )


def validate_sensitivity_exclusions(
    exclusions: Mapping[str, Sequence[int]],
) -> None:
    assert_equal(tuple(exclusions), SAMPLE_ORDER, "sensitivity library order")
    for gsm in SAMPLE_ORDER:
        assert_equal(
            tuple(exclusions[gsm]),
            MANDATORY_SENSITIVITY_EXCLUDED_LENGTHS,
            f"mandatory sensitivity exclusions for {gsm}",
        )


def coordinate_log_ratio(
    target_count: np.ndarray | float,
    background_count: np.ndarray | float,
    background_length: np.ndarray | float,
) -> np.ndarray:
    target = np.asarray(target_count, dtype=float)
    background = np.asarray(background_count, dtype=float)
    length = np.asarray(background_length, dtype=float)
    if np.any(target < 0) or np.any(background < 0) or np.any(length <= 0):
        raise ValueError("Invalid coordinate-null count or background length")
    return np.log2(((target + 0.5) / TARGET_LENGTH) / ((background + 0.5) / length))


def coordinate_null_statistic(
    real_log_ratios: np.ndarray,
    null_log_ratios: np.ndarray,
) -> dict[str, object]:
    real = np.asarray(real_log_ratios, dtype=float)
    null = np.asarray(null_log_ratios, dtype=float)
    if real.ndim != 2 or real.shape[0] != 6:
        raise ValueError("Real coordinate-null ratios must be library x anchor")
    if null.ndim != 3 or null.shape[1:] != real.shape:
        raise ValueError("Null ratios must be realization x library x anchor")
    s_real = float(np.mean(np.mean(real, axis=1)))
    s_null = np.mean(np.mean(null, axis=2), axis=1)
    p = float((np.count_nonzero(s_null >= s_real) + 1) / (len(s_null) + 1))
    null_median = float(np.median(s_null))
    return {
        "S_real": s_real,
        "S_null": s_null,
        "null_median": null_median,
        "empirical_p": p,
        "core_pass": bool(s_real > null_median and p < 0.05),
    }


def _average_ranks(values: Sequence[float]) -> tuple[list[float], Counter[float]]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    ties: Counter[float] = Counter(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        rank = ((cursor + 1) + end) / 2.0
        for position in range(cursor, end):
            ranks[order[position]] = rank
        cursor = end
    return ranks, ties


def signed_rank_greater(values: Sequence[float]) -> dict[str, float]:
    nonzero = [float(value) for value in values if float(value) != 0.0]
    if not nonzero:
        raise ValueError("Signed-rank test has no nonzero observations")
    absolute = [abs(value) for value in nonzero]
    ranks, ties = _average_ranks(absolute)
    w_plus = sum(rank for rank, value in zip(ranks, nonzero) if value > 0)
    n = len(nonzero)
    mean_w = n * (n + 1) / 4.0
    variance = n * (n + 1) * (2 * n + 1) / 24.0
    variance -= sum(count**3 - count for count in ties.values()) / 48.0
    if variance <= 0:
        raise ValueError("Signed-rank variance is non-positive")
    z = (w_plus - mean_w - 0.5) / math.sqrt(variance)
    p = 0.5 * math.erfc(z / math.sqrt(2.0))
    return {"W_plus": w_plus, "z": z, "p_greater": p}


def hodges_lehmann_interval(
    values: Sequence[float],
    confidence: float = 0.90,
) -> tuple[float, float, float]:
    x = np.asarray([float(value) for value in values], dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.all(np.isfinite(x)):
        raise ValueError("Hodges-Lehmann interval needs finite one-dimensional data")
    walsh = np.concatenate(
        [(x[index] + x[index:]) / 2.0 for index in range(len(x))]
    )
    walsh.sort()
    n = len(x)
    mean_w = n * (n + 1) / 4.0
    sd_w = math.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    z = 1.6448536269514722 if confidence == 0.90 else None
    if z is None:
        raise ValueError("Only the frozen 90% interval is supported")
    critical = max(0, int(math.floor(mean_w - z * sd_w)))
    lower_index = min(critical, len(walsh) - 1)
    upper_index = max(0, len(walsh) - critical - 1)
    return (
        float(np.median(walsh)),
        float(walsh[lower_index]),
        float(walsh[upper_index]),
    )


def real_pseudo_statistic(
    pair_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    by_transcript: dict[str, list[float]] = defaultdict(list)
    for row in pair_rows:
        by_transcript[str(row["transcript"])].append(float(row["delta"]))
    d_by_transcript = {
        transcript: float(np.mean(values))
        for transcript, values in sorted(by_transcript.items())
    }
    test = signed_rank_greater(list(d_by_transcript.values()))
    location, lower, upper = hodges_lehmann_interval(list(d_by_transcript.values()))
    med = float(np.median(list(d_by_transcript.values())))
    return {
        "D_by_transcript": d_by_transcript,
        "n_pairs": len(pair_rows),
        "n_transcripts": len(d_by_transcript),
        "median_D": med,
        "hodges_lehmann": location,
        "hl_90_lower": lower,
        "hl_90_upper": upper,
        "wilcoxon_p_greater": test["p_greater"],
        "core_pass": bool(med > 0 and test["p_greater"] < 0.05),
    }


def localization_statistics(
    beta_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    by_transcript_center: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in beta_rows:
        by_transcript_center[
            (str(row["transcript"]), int(row["center"]))
        ].append(float(row["beta"]))
    transcripts = sorted({key[0] for key in by_transcript_center})
    r: dict[tuple[str, int], float] = {}
    for transcript in transcripts:
        for center in CENTERS:
            values = by_transcript_center.get((transcript, center))
            if not values:
                raise ValueError("Localization input is not a common all-center universe")
            r[(transcript, center)] = float(np.mean(values))
    off_band = tuple(center for center in CENTERS if center not in LDB_BAND)
    p_values = [float(np.mean([r[(t, c)] for c in LDB_BAND])) for t in transcripts]
    q_values = [float(np.mean([r[(t, c)] for c in off_band])) for t in transcripts]
    differences = [p - q for p, q in zip(p_values, q_values)]
    ldb_test = signed_rank_greater(differences)
    median_p = float(np.median(p_values))
    median_difference = float(np.median(differences))
    _, ngs_lower, ngs_upper = hodges_lehmann_interval(q_values)
    return {
        "n_transcripts": len(transcripts),
        "median_P": median_p,
        "median_P_minus_Q": median_difference,
        "ldb_wilcoxon_p_greater": ldb_test["p_greater"],
        "LDB_pass": bool(
            median_p > 0
            and median_difference > 0
            and ldb_test["p_greater"] < 0.05
        ),
        "NGS_epsilon": NGS_EPSILON,
        "NGS_90_lower": ngs_lower,
        "NGS_90_upper": ngs_upper,
        "NGS_pass": bool(-NGS_EPSILON < ngs_lower and ngs_upper < NGS_EPSILON),
        "P_by_transcript": dict(zip(transcripts, p_values)),
        "Q_by_transcript": dict(zip(transcripts, q_values)),
    }


def load_features(path: Path) -> dict[str, Feature]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"Transcript", "UTR5_Length", "UTR3_Length", "CDS_Length"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Invalid human feature table: {path}")
    features: dict[str, Feature] = {}
    for row in rows:
        transcript = row["Transcript"]
        if transcript in features:
            raise ValueError(f"Duplicate human feature transcript: {transcript}")
        feature = Feature(
            transcript=transcript,
            utr5_length=int(row["UTR5_Length"]),
            cds_length=int(row["CDS_Length"]),
            utr3_length=int(row["UTR3_Length"]),
        )
        if feature.cds_length <= 0:
            raise ValueError(f"Non-positive CDS length: {transcript}")
        features[transcript] = feature
    assert_equal(len(features), 19_192, "frozen human feature transcript count")
    return features


def load_gene_map(path: Path, features: Mapping[str, Feature]) -> dict[str, str]:
    candidates: dict[str, set[str]] = defaultdict(set)
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Invalid APPRIS mapping line {line_number}")
            candidates[fields[2]].add(fields[1])
    mapping: dict[str, str] = {}
    for transcript in features:
        genes = candidates.get(transcript, set())
        if len(genes) != 1:
            raise ValueError(
                f"Transcript-to-gene mapping is not one-to-one for {transcript}: "
                f"{sorted(genes)}"
            )
        mapping[transcript] = next(iter(genes))
    assert_equal(len(mapping), 19_192, "transcript-to-gene mapping coverage")
    return mapping


def load_topology_masks(
    path: Path,
    features: Mapping[str, Feature],
) -> dict[str, tuple[tuple[int, int], ...]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "transcript",
            "raw_deeptmhmm_label",
            "cds_nt_start0",
            "cds_nt_end0_exclusive",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid normalized topology table: {path}")
        for row in reader:
            if row["raw_deeptmhmm_label"] not in TOPOLOGY_MASK_CLASSES:
                continue
            transcript = row["transcript"]
            if transcript not in features:
                raise ValueError(f"Topology transcript lacks frozen features: {transcript}")
            start = int(row["cds_nt_start0"])
            end = int(row["cds_nt_end0_exclusive"])
            intervals[transcript].append(focal_interval(start, end, start))
    return {
        transcript: merge_intervals(
            intervals.get(transcript, ()), feature.cds_length
        )
        for transcript, feature in features.items()
    }


def load_primary_anchors(
    path: Path,
    features: Mapping[str, Feature],
    pseudo_path: Path,
) -> list[Anchor]:
    primary_valid: set[str] = set()
    with pseudo_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"transcript", "TMD_index", "real_primary_window_valid"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Pseudo table lacks primary validity fields: {pseudo_path}")
        for row in reader:
            if row["real_primary_window_valid"] == "YES":
                primary_valid.add(f"{row['transcript']}|{row['TMD_index']}")
    anchors: list[Anchor] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "transcript",
            "TMD_index",
            "cds_nt_start0",
            "cds_nt_end0_exclusive",
            "primary_anchor_cds_coordinate",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid primary anchor table: {path}")
        for row in reader:
            transcript = row["transcript"]
            tmd_index = row["TMD_index"]
            anchor = Anchor(
                anchor_id=f"{transcript}|{tmd_index}",
                transcript=transcript,
                tmd_index=tmd_index,
                anchor_cds0=int(row["primary_anchor_cds_coordinate"]),
                segment_start0=int(row["cds_nt_start0"]),
                segment_end0=int(row["cds_nt_end0_exclusive"]),
            )
            if transcript not in features:
                raise ValueError(f"Primary anchor lacks frozen features: {transcript}")
            start, end = target_interval(anchor.anchor_cds0, PRIMARY_CENTER)
            if (
                anchor.anchor_id in primary_valid
                and 0 <= start < end <= features[transcript].cds_length
            ):
                anchors.append(anchor)
    assert_equal(len(anchors), PRIMARY_STARTING_ANCHORS, "primary starting anchors")
    if len({anchor.anchor_id for anchor in anchors}) != len(anchors):
        raise ValueError("Primary anchor identifiers are not unique")
    return anchors


def load_pseudo_pairs(path: Path) -> list[PseudoPair]:
    pairs: list[PseudoPair] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "transcript",
            "TMD_index",
            "real_anchor_cds0",
            "real_TMD_cds_start0",
            "real_TMD_cds_end0_exclusive",
            "pseudo_anchor_cds0",
            "pseudo_segment_cds_start0",
            "pseudo_segment_cds_end0_exclusive",
            "match_status",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid pseudo-pair table: {path}")
        for row in reader:
            if row["match_status"] != "MATCHED":
                continue
            transcript = row["transcript"]
            tmd_index = row["TMD_index"]
            pairs.append(PseudoPair(
                pair_id=f"{transcript}|{tmd_index}",
                transcript=transcript,
                tmd_index=tmd_index,
                real_anchor_cds0=int(row["real_anchor_cds0"]),
                real_segment_start0=int(row["real_TMD_cds_start0"]),
                real_segment_end0=int(row["real_TMD_cds_end0_exclusive"]),
                pseudo_anchor_cds0=int(row["pseudo_anchor_cds0"]),
                pseudo_segment_start0=int(row["pseudo_segment_cds_start0"]),
                pseudo_segment_end0=int(row["pseudo_segment_cds_end0_exclusive"]),
            ))
    assert_equal(len(pairs), PSEUDO_STARTING_PAIRS, "matched pseudo-pair count")
    if len({pair.pair_id for pair in pairs}) != len(pairs):
        raise ValueError("Pseudo pair identifiers are not unique")
    return pairs


def load_all_tmd_intervals(
    path: Path,
) -> dict[str, tuple[tuple[str, int, int], ...]]:
    intervals: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "transcript",
            "TMD_index",
            "cds_nt_start0",
            "cds_nt_end0_exclusive",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid primary anchor intervals: {path}")
        for row in reader:
            intervals[row["transcript"]].append((
                f"{row['transcript']}|{row['TMD_index']}",
                int(row["cds_nt_start0"]),
                int(row["cds_nt_end0_exclusive"]),
            ))
    return {transcript: tuple(rows) for transcript, rows in intervals.items()}


def load_psite_tracks(
    manifest_rows: Sequence[Mapping[str, str]],
    features: Mapping[str, Feature],
) -> dict[str, dict[str, dict[str, SparseTrack]]]:
    """Read biological P-site rows only from the execution backend."""

    tracks: dict[str, dict[str, dict[str, SparseTrack]]] = {}
    for manifest_row in manifest_rows:
        gsm = str(manifest_row["GSM"])
        sample = manifest_row.get("sample", gsm)
        path = Path(manifest_row["psite_path"]).expanduser().resolve()
        primary_counts: dict[str, Counter[int]] = defaultdict(Counter)
        sensitivity_counts: dict[str, Counter[int]] = defaultdict(Counter)
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"sample", "transcript", "read_length", "relative_psite"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(f"Invalid calibrated P-site table: {path}")
            for row in reader:
                if row["sample"] != sample:
                    raise ValueError(f"P-site sample mismatch in {path}")
                transcript = row["transcript"]
                if transcript not in features:
                    raise ValueError(f"P-site transcript lacks frozen features: {transcript}")
                read_length = int(row["read_length"])
                if read_length not in FROZEN_READ_LENGTHS:
                    raise ValueError(f"Non-frozen read length entered Stage 1: {read_length}")
                position = int(row["relative_psite"])
                if 0 <= position < features[transcript].cds_length:
                    primary_counts[transcript][position] += 1
                    if read_length not in MANDATORY_SENSITIVITY_EXCLUDED_LENGTHS:
                        sensitivity_counts[transcript][position] += 1
        tracks[gsm] = {
            "primary": {
                transcript: SparseTrack.from_counter(counts)
                for transcript, counts in primary_counts.items()
            },
            "exclude_28_29": {
                transcript: SparseTrack.from_counter(counts)
                for transcript, counts in sensitivity_counts.items()
            },
        }
    assert_equal(tuple(tracks), SAMPLE_ORDER, "loaded biological library order")
    return tracks


def build_feature_rows(
    feature_id: str,
    transcript: str,
    anchor_cds0: int,
    features: Mapping[str, Feature],
    masks: Mapping[str, Sequence[tuple[int, int]]],
    tracks: Mapping[str, Mapping[str, Mapping[str, SparseTrack]]],
    variant: str,
    gene_id: str,
    center: int = PRIMARY_CENTER,
    additional_focal: tuple[int, int] | None = None,
) -> list[dict[str, object]]:
    feature = features[transcript]
    rows: list[dict[str, object]] = []
    for gsm in SAMPLE_ORDER:
        track = tracks[gsm][variant].get(transcript, EMPTY_TRACK)
        try:
            target_count, background_count, background_length = profile_counts(
                track,
                feature.cds_length,
                masks[transcript],
                anchor_cds0,
                center,
                additional_focal,
            )
        except IneligibleFeatureGeometry:
            return []
        if background_length == 0:
            return []
        rows.append({
            "GSM": gsm,
            "sample": gsm,
            "condition": CONDITIONS[gsm],
            "feature_id": feature_id,
            "transcript": transcript,
            "gene_id": gene_id,
            "center": center,
            "target_count": target_count,
            "background_count": background_count,
            "target_length": TARGET_LENGTH,
            "background_length": background_length,
            "exposure": exposure(background_count, background_length),
        })
    return rows


def build_primary_and_positional_rows(
    anchors: Sequence[Anchor],
    features: Mapping[str, Feature],
    genes: Mapping[str, str],
    masks: Mapping[str, Sequence[tuple[int, int]]],
    tracks: Mapping[str, Mapping[str, Mapping[str, SparseTrack]]],
    variant: str,
    all_tmd_intervals: Mapping[str, Sequence[tuple[str, int, int]]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    set[str],
]:
    by_anchor_center: dict[tuple[str, int], list[dict[str, object]]] = {}
    common: set[str] = set()
    for anchor in anchors:
        if anchor.transcript not in features:
            raise ValueError(f"Anchor transcript lacks frozen features: {anchor.transcript}")
        all_centers_valid = True
        for center in CENTERS:
            target_start, target_end = target_interval(anchor.anchor_cds0, center)
            if any(
                other_id != anchor.anchor_id
                and max(target_start, other_start)
                < min(target_end, other_end)
                for other_id, other_start, other_end in all_tmd_intervals[
                    anchor.transcript
                ]
            ):
                all_centers_valid = False
                continue
            rows = build_feature_rows(
                anchor.anchor_id,
                anchor.transcript,
                anchor.anchor_cds0,
                features,
                masks,
                tracks,
                variant,
                genes[anchor.transcript],
                center,
            )
            by_anchor_center[(anchor.anchor_id, center)] = rows
            if not strict_six_library_eligible(rows):
                all_centers_valid = False
        if all_centers_valid:
            common.add(anchor.anchor_id)
    primary_rows = [
        row
        for anchor in anchors
        for row in by_anchor_center.get((anchor.anchor_id, PRIMARY_CENTER), ())
        if strict_six_library_eligible(
            by_anchor_center.get((anchor.anchor_id, PRIMARY_CENTER), ())
        )
    ]
    centerwise_rows = [
        row
        for anchor in anchors
        for center in CENTERS
        for row in by_anchor_center.get((anchor.anchor_id, center), ())
        if strict_six_library_eligible(
            by_anchor_center.get((anchor.anchor_id, center), ())
        )
    ]
    common_rows = [
        row
        for anchor in anchors
        if anchor.anchor_id in common
        for center in CENTERS
        for row in by_anchor_center[(anchor.anchor_id, center)]
    ]
    return primary_rows, centerwise_rows, common_rows, common


def build_pseudo_rows(
    pairs: Sequence[PseudoPair],
    features: Mapping[str, Feature],
    genes: Mapping[str, str],
    masks: Mapping[str, Sequence[tuple[int, int]]],
    tracks: Mapping[str, Mapping[str, Mapping[str, SparseTrack]]],
    variant: str,
) -> tuple[list[dict[str, object]], list[str]]:
    model_rows: list[dict[str, object]] = []
    eligible_pairs: list[str] = []
    for pair in pairs:
        real_rows = build_feature_rows(
            f"{pair.pair_id}|real",
            pair.transcript,
            pair.real_anchor_cds0,
            features,
            masks,
            tracks,
            variant,
            genes[pair.transcript],
        )
        pseudo_focal = focal_interval(
            pair.pseudo_segment_start0,
            pair.pseudo_segment_end0,
            pair.pseudo_anchor_cds0,
        )
        pseudo_rows = build_feature_rows(
            f"{pair.pair_id}|pseudo",
            pair.transcript,
            pair.pseudo_anchor_cds0,
            features,
            masks,
            tracks,
            variant,
            genes[pair.transcript],
            additional_focal=pseudo_focal,
        )
        if paired_eligibility(real_rows, pseudo_rows):
            eligible_pairs.append(pair.pair_id)
            model_rows.extend(real_rows)
            model_rows.extend(pseudo_rows)
    return model_rows, eligible_pairs


def build_reagent_omission_rows(
    anchors: Sequence[Anchor],
    features: Mapping[str, Feature],
    genes: Mapping[str, str],
    masks: Mapping[str, Sequence[tuple[int, int]]],
    tracks: Mapping[str, Mapping[str, Mapping[str, SparseTrack]]],
    variant: str,
) -> list[dict[str, object]]:
    plans = {
        "omit_siP1_1": tuple(gsm for gsm in SAMPLE_ORDER if gsm != "GSM3900211"),
        "omit_siP2_1_siP2_4": tuple(
            gsm for gsm in SAMPLE_ORDER if gsm not in {"GSM3900212", "GSM3900213"}
        ),
    }
    output: list[dict[str, object]] = []
    for anchor in anchors:
        full_rows = build_feature_rows(
            anchor.anchor_id,
            anchor.transcript,
            anchor.anchor_cds0,
            features,
            masks,
            tracks,
            variant,
            genes[anchor.transcript],
        )
        if not full_rows:
            continue
        by_gsm = {str(row["GSM"]): row for row in full_rows}
        for family, sample_order in plans.items():
            rows = [by_gsm[gsm] for gsm in sample_order]
            if eligible_for_samples(rows, sample_order):
                for row in rows:
                    item = dict(row)
                    item["analysis_family"] = family
                    item["report_inference"] = "NO"
                    output.append(item)
    return output


def validate_coordinate_null_identity(
    version: str,
    matrix_path: Path,
    expected_sha256: str = COORDINATE_NULL_SHA256,
) -> None:
    assert_equal(version, COORDINATE_NULL_VERSION, "coordinate-null version")
    if "coordinate_null_v2" not in matrix_path.parts:
        raise ValueError("Coordinate-null matrix must come from the authoritative V2 directory")
    assert_equal(sha256(matrix_path), expected_sha256, "coordinate-null V2 matrix SHA256")


def load_coordinate_columns(path: Path) -> list[dict[str, object]]:
    columns: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "anchor_column_index0",
            "transcript",
            "TMD_index",
            "real_anchor_cds0",
            "coordinate_null_policy_freeze_version",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid coordinate-null V2 column table: {path}")
        for row in reader:
            expected_index = len(columns)
            assert_equal(
                int(row["anchor_column_index0"]),
                expected_index,
                "coordinate-null column order",
            )
            if "_v2_" not in row["coordinate_null_policy_freeze_version"]:
                raise ValueError("Coordinate-null V1/superseded policy is prohibited")
            columns.append({
                "anchor_id": f"{row['transcript']}|{row['TMD_index']}",
                "transcript": row["transcript"],
                "real_anchor_cds0": int(row["real_anchor_cds0"]),
            })
    assert_equal(len(columns), COORDINATE_NULL_ANCHORS, "coordinate-null anchor count")
    assert_equal(
        len({str(row["transcript"]) for row in columns}),
        COORDINATE_NULL_TRANSCRIPTS,
        "coordinate-null transcript count",
    )
    return columns


def _coordinate_prefixes(
    track: SparseTrack,
    cds_length: int,
    ordinary_masks: Sequence[tuple[int, int]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    all_counts = np.zeros(cds_length, dtype=np.int64)
    if track.positions:
        all_counts[np.asarray(track.positions, dtype=np.int64)] = np.diff(
            np.asarray(track.prefix_counts, dtype=np.int64)
        )
    ordinary_eligible = np.ones(cds_length, dtype=np.int64)
    for start, end in ordinary_masks:
        ordinary_eligible[start:end] = 0
    eligible_counts = all_counts * ordinary_eligible
    all_prefix = np.concatenate(([0], np.cumsum(all_counts, dtype=np.int64)))
    eligible_count_prefix = np.concatenate(
        ([0], np.cumsum(eligible_counts, dtype=np.int64))
    )
    eligible_length_prefix = np.concatenate(
        ([0], np.cumsum(ordinary_eligible, dtype=np.int64))
    )
    return (
        all_prefix,
        eligible_count_prefix,
        eligible_length_prefix,
        int(eligible_count_prefix[-1]),
        int(eligible_length_prefix[-1]),
    )


def _coordinate_vectors(
    coordinates: np.ndarray,
    cds_length: int,
    prefixes: tuple[np.ndarray, np.ndarray, np.ndarray, int, int],
) -> np.ndarray:
    q = np.asarray(coordinates, dtype=np.int64)
    target_start = q + 60
    target_end = q + 91
    if np.any(target_start < 0) or np.any(target_end > cds_length):
        raise ValueError("Coordinate-null target is outside the CDS")
    focal_start = np.maximum(0, q - 30)
    focal_end = np.minimum(cds_length, q + 136)
    all_prefix, eligible_count_prefix, eligible_length_prefix, base_count, base_length = prefixes
    target = all_prefix[target_end] - all_prefix[target_start]
    removed_count = eligible_count_prefix[focal_end] - eligible_count_prefix[focal_start]
    removed_length = eligible_length_prefix[focal_end] - eligible_length_prefix[focal_start]
    background = base_count - removed_count
    background_length = base_length - removed_length
    return coordinate_log_ratio(target, background, background_length)


def compute_coordinate_null_v2(
    repo: Path,
    features: Mapping[str, Feature],
    masks: Mapping[str, Sequence[tuple[int, int]]],
    tracks: Mapping[str, Mapping[str, Mapping[str, SparseTrack]]],
    variant: str,
    work_dir: Path,
    chunk_size: int = 32,
) -> dict[str, object]:
    base = repo / "docs/stage0c/coordinate_null_v2"
    matrix_path = base / "coordinate_null_coordinates.npy"
    validate_coordinate_null_identity(COORDINATE_NULL_VERSION, matrix_path)
    columns = load_coordinate_columns(base / "coordinate_null_anchor_columns.tsv")
    matrix = np.load(matrix_path, mmap_mode="r")
    assert_equal(
        matrix.shape,
        (COORDINATE_NULL_REALIZATIONS, COORDINATE_NULL_ANCHORS),
        "coordinate-null V2 matrix shape",
    )
    groups: dict[str, np.ndarray] = {}
    for transcript in sorted({str(row["transcript"]) for row in columns}):
        groups[transcript] = np.asarray(
            [index for index, row in enumerate(columns) if row["transcript"] == transcript],
            dtype=np.int64,
        )
    transcripts = tuple(groups)
    transcript_index = {transcript: index for index, transcript in enumerate(transcripts)}
    real_coordinates = np.asarray(
        [int(row["real_anchor_cds0"]) for row in columns], dtype=np.int64
    )
    real_by_sample = np.zeros((len(SAMPLE_ORDER), len(columns)), dtype=float)
    prefixes_by_sample: dict[str, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, int, int]]] = {}
    for sample_index, gsm in enumerate(SAMPLE_ORDER):
        prefixes_by_sample[gsm] = {}
        for transcript, indices in groups.items():
            track = tracks[gsm][variant].get(transcript, EMPTY_TRACK)
            prefixes = _coordinate_prefixes(
                track,
                features[transcript].cds_length,
                masks[transcript],
            )
            prefixes_by_sample[gsm][transcript] = prefixes
            real_by_sample[sample_index, indices] = _coordinate_vectors(
                real_coordinates[indices],
                features[transcript].cds_length,
                prefixes,
            )
    real_anchor_contribution = np.mean(real_by_sample, axis=0)
    global_real_sum = float(np.sum(real_anchor_contribution))
    real_transcript_contribution = np.zeros(len(transcripts), dtype=float)
    transcript_anchor_counts = np.zeros(len(transcripts), dtype=np.int64)
    for transcript, indices in groups.items():
        index = transcript_index[transcript]
        real_transcript_contribution[index] = float(
            np.sum(real_anchor_contribution[indices])
        )
        transcript_anchor_counts[index] = len(indices)

    null_anchor_path = work_dir / f".{variant}.coordinate_null_anchor_contrib.float32"
    null_transcript_path = work_dir / f".{variant}.coordinate_null_transcript_contrib.float32"
    null_anchor = np.memmap(
        null_anchor_path,
        dtype="float32",
        mode="w+",
        shape=(COORDINATE_NULL_REALIZATIONS, COORDINATE_NULL_ANCHORS),
    )
    null_transcript = np.memmap(
        null_transcript_path,
        dtype="float32",
        mode="w+",
        shape=(len(transcripts), COORDINATE_NULL_REALIZATIONS),
    )
    global_null_sum = np.zeros(COORDINATE_NULL_REALIZATIONS, dtype=float)
    for row_start in range(0, COORDINATE_NULL_REALIZATIONS, chunk_size):
        row_end = min(COORDINATE_NULL_REALIZATIONS, row_start + chunk_size)
        chunk_values = np.zeros((row_end - row_start, len(columns)), dtype=float)
        for gsm in SAMPLE_ORDER:
            for transcript, indices in groups.items():
                q = np.asarray(matrix[row_start:row_end, indices], dtype=np.int64)
                chunk_values[:, indices] += _coordinate_vectors(
                    q,
                    features[transcript].cds_length,
                    prefixes_by_sample[gsm][transcript],
                ) / len(SAMPLE_ORDER)
        null_anchor[row_start:row_end, :] = chunk_values.astype("float32")
        global_null_sum[row_start:row_end] = np.sum(chunk_values, axis=1)
        for transcript, indices in groups.items():
            null_transcript[
                transcript_index[transcript], row_start:row_end
            ] = np.sum(chunk_values[:, indices], axis=1).astype("float32")
    null_anchor.flush()
    null_transcript.flush()
    s_real = global_real_sum / len(columns)
    s_null = global_null_sum / len(columns)
    empirical_p = float(
        (np.count_nonzero(s_null >= s_real) + 1)
        / (COORDINATE_NULL_REALIZATIONS + 1)
    )
    summary = {
        "coordinate_null_version": COORDINATE_NULL_VERSION,
        "matrix_sha256": COORDINATE_NULL_SHA256,
        "anchors": len(columns),
        "transcripts": len(transcripts),
        "realizations": len(s_null),
        "S_real": s_real,
        "null_median": float(np.median(s_null)),
        "null_q025": float(np.quantile(s_null, 0.025)),
        "null_q975": float(np.quantile(s_null, 0.975)),
        "empirical_p": empirical_p,
        "core_pass": bool(s_real > np.median(s_null) and empirical_p < 0.05),
    }
    return {
        "summary": summary,
        "columns": columns,
        "transcripts": transcripts,
        "transcript_anchor_counts": transcript_anchor_counts,
        "real_anchor_contribution": real_anchor_contribution,
        "real_transcript_contribution": real_transcript_contribution,
        "global_real_sum": global_real_sum,
        "global_null_sum": global_null_sum,
        "null_anchor": null_anchor,
        "null_transcript": null_transcript,
        "temporary_paths": (null_anchor_path, null_transcript_path),
    }


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty Stage 1 table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_edger(
    repo: Path,
    rows: Sequence[Mapping[str, object]],
    prefix: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    model_input = prefix.with_suffix(".model_input.tsv")
    results_path = prefix.with_suffix(".edger_results.tsv")
    summary_path = prefix.with_suffix(".edger_summary.tsv")
    model_rows = []
    for row in rows:
        model_rows.append({
            "analysis_family": row["analysis_family"],
            "feature_id": row["feature_id"],
            "transcript": row["transcript"],
            "gene_id": row["gene_id"],
            "sample": row["sample"],
            "condition": row["condition"],
            "target_count": row["target_count"],
            "background_count": row["background_count"],
            "target_length": row["target_length"],
            "background_length": row["background_length"],
            "report_inference": row.get("report_inference", "YES"),
        })
    write_tsv(model_input, model_rows)
    command = [
        "/usr/local/bin/Rscript",
        str(repo / "R/fit_paper3_stage1_onefactor.R"),
        str(model_input),
        str(results_path),
        str(summary_path),
    ]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(
            "edgeR Stage 1 fit failed:\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    return read_tsv(results_path), read_tsv(summary_path)


def bh_adjust(p_values: Sequence[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("BH adjustment requires finite probabilities")
    order = np.argsort(p)
    adjusted = np.empty(len(p), dtype=float)
    running = 1.0
    for reverse_rank, index in enumerate(order[::-1], start=1):
        rank = len(p) - reverse_rank + 1
        running = min(running, p[index] * len(p) / rank)
        adjusted[index] = min(1.0, running)
    return adjusted.tolist()


def _rank(values: Sequence[float]) -> np.ndarray:
    ranks, _ = _average_ranks([float(value) for value in values])
    return np.asarray(ranks, dtype=float)


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("Spearman correlation requires paired values")
    rx, ry = _rank(x), _rank(y)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def _coordinate_loo_for_indices(
    coordinate: Mapping[str, object],
    indices: np.ndarray,
) -> dict[str, object]:
    n_total = int(coordinate["summary"]["anchors"])
    if len(indices) == 0:
        return dict(coordinate["summary"])
    remaining = n_total - len(indices)
    if remaining <= 0:
        raise ValueError("Coordinate-null LOO removed the complete universe")
    real_contribution = np.asarray(coordinate["real_anchor_contribution"])[indices]
    null_contribution = np.asarray(coordinate["null_anchor"][:, indices], dtype=float)
    s_real = (
        float(coordinate["global_real_sum"]) - float(np.sum(real_contribution))
    ) / remaining
    s_null = (
        np.asarray(coordinate["global_null_sum"], dtype=float)
        - np.sum(null_contribution, axis=1)
    ) / remaining
    p = float((np.count_nonzero(s_null >= s_real) + 1) / (len(s_null) + 1))
    null_median = float(np.median(s_null))
    return {
        "S_real": s_real,
        "null_median": null_median,
        "empirical_p": p,
        "core_pass": bool(s_real > null_median and p < 0.05),
    }


def cluster_loo(
    real_pseudo: Mapping[str, object],
    coordinate: Mapping[str, object],
    gene_map: Mapping[str, str],
    unit: str,
) -> tuple[list[dict[str, object]], bool]:
    if unit not in {"transcript", "gene"}:
        raise ValueError("LOO unit must be transcript or gene")
    d_by_transcript = {
        str(key): float(value)
        for key, value in real_pseudo["D_by_transcript"].items()
    }
    coordinate_columns = coordinate["columns"]
    all_transcripts = sorted(
        set(d_by_transcript)
        | {str(row["transcript"]) for row in coordinate_columns}
    )
    clusters: dict[str, set[str]] = defaultdict(set)
    for transcript in all_transcripts:
        cluster = transcript if unit == "transcript" else gene_map[transcript]
        clusters[cluster].add(transcript)
    rows: list[dict[str, object]] = []
    for cluster, omitted_transcripts in sorted(clusters.items()):
        remaining_d = [
            value
            for transcript, value in d_by_transcript.items()
            if transcript not in omitted_transcripts
        ]
        rp_test = signed_rank_greater(remaining_d)
        rp_median = float(np.median(remaining_d))
        rp_pass = bool(rp_median > 0 and rp_test["p_greater"] < 0.05)
        indices = np.asarray(
            [
                index
                for index, row in enumerate(coordinate_columns)
                if row["transcript"] in omitted_transcripts
            ],
            dtype=np.int64,
        )
        cn = _coordinate_loo_for_indices(coordinate, indices)
        rows.append({
            "unit": unit,
            "omitted_cluster": cluster,
            "omitted_transcripts": len(omitted_transcripts),
            "real_pseudo_median_D": rp_median,
            "real_pseudo_p": rp_test["p_greater"],
            "real_pseudo_pass": rp_pass,
            "coordinate_S_real": cn["S_real"],
            "coordinate_null_median": cn["null_median"],
            "coordinate_empirical_p": cn["empirical_p"],
            "coordinate_pass": cn["core_pass"],
            "combined_pass": bool(rp_pass and cn["core_pass"]),
        })
    return rows, all(bool(row["combined_pass"]) for row in rows)


def anchor_loo_descriptive(
    real_pseudo_pairs: Sequence[Mapping[str, object]],
    real_pseudo: Mapping[str, object],
    coordinate: Mapping[str, object],
) -> dict[str, object]:
    d_lists: dict[str, list[float]] = defaultdict(list)
    for row in real_pseudo_pairs:
        d_lists[str(row["transcript"])].append(float(row["delta"]))
    full_d = {
        transcript: float(np.mean(values))
        for transcript, values in d_lists.items()
    }
    full_rp_p = float(real_pseudo["wilcoxon_p_greater"])
    max_rp_p_change = 0.0
    max_rp_stat_change = 0.0
    rp_direction_changes = 0
    rp_decision_changes = 0
    for row in real_pseudo_pairs:
        transcript = str(row["transcript"])
        values = d_lists[transcript]
        removed = float(row["delta"])
        remaining_by_transcript = dict(full_d)
        if len(values) == 1:
            del remaining_by_transcript[transcript]
        else:
            remaining_by_transcript[transcript] = (
                sum(values) - removed
            ) / (len(values) - 1)
        remaining = list(remaining_by_transcript.values())
        test = signed_rank_greater(remaining)
        med = float(np.median(remaining))
        max_rp_p_change = max(max_rp_p_change, abs(test["p_greater"] - full_rp_p))
        max_rp_stat_change = max(
            max_rp_stat_change,
            abs(med - float(real_pseudo["median_D"])),
        )
        rp_direction_changes += int((med > 0) != (float(real_pseudo["median_D"]) > 0))
        rp_decision_changes += int(
            (med > 0 and test["p_greater"] < 0.05)
            != bool(real_pseudo["core_pass"])
        )

    full_cn = coordinate["summary"]
    max_cn_p_change = 0.0
    max_cn_stat_change = 0.0
    cn_direction_changes = 0
    cn_decision_changes = 0
    for index in range(int(full_cn["anchors"])):
        result = _coordinate_loo_for_indices(
            coordinate, np.asarray([index], dtype=np.int64)
        )
        max_cn_p_change = max(
            max_cn_p_change,
            abs(float(result["empirical_p"]) - float(full_cn["empirical_p"])),
        )
        max_cn_stat_change = max(
            max_cn_stat_change,
            abs(float(result["S_real"]) - float(full_cn["S_real"])),
        )
        cn_direction_changes += int(
            (float(result["S_real"]) > float(result["null_median"]))
            != (float(full_cn["S_real"]) > float(full_cn["null_median"]))
        )
        cn_decision_changes += int(
            bool(result["core_pass"]) != bool(full_cn["core_pass"])
        )
    return {
        "role": "DESCRIPTIVE_ONLY",
        "real_pseudo_anchor_iterations": len(real_pseudo_pairs),
        "real_pseudo_max_abs_median_change": max_rp_stat_change,
        "real_pseudo_max_abs_p_change": max_rp_p_change,
        "real_pseudo_direction_changes": rp_direction_changes,
        "real_pseudo_decision_changes": rp_decision_changes,
        "coordinate_anchor_iterations": int(full_cn["anchors"]),
        "coordinate_max_abs_S_real_change": max_cn_stat_change,
        "coordinate_max_abs_p_change": max_cn_p_change,
        "coordinate_direction_changes": cn_direction_changes,
        "coordinate_decision_changes": cn_decision_changes,
    }


def _model_rows_with_family(
    rows: Sequence[Mapping[str, object]],
    family: str | None = None,
    report_inference: str = "YES",
) -> list[dict[str, object]]:
    output = []
    for row in rows:
        item = dict(row)
        item["analysis_family"] = family or f"center_{int(row['center']):+d}"
        item["report_inference"] = report_inference
        output.append(item)
    return output


def _pair_effect_rows(
    results: Sequence[Mapping[str, str]],
) -> list[dict[str, object]]:
    effects = {str(row["feature_id"]): float(row["logFC"]) for row in results}
    annotations = {
        str(row["feature_id"]): str(row["transcript"]) for row in results
    }
    pair_ids = sorted(
        feature_id[:-5]
        for feature_id in effects
        if feature_id.endswith("|real")
    )
    rows = []
    for pair_id in pair_ids:
        real_id = f"{pair_id}|real"
        pseudo_id = f"{pair_id}|pseudo"
        if pseudo_id not in effects:
            raise ValueError(f"Eligible real feature lacks matched pseudo: {pair_id}")
        rows.append({
            "pair_id": pair_id,
            "transcript": annotations[real_id],
            "beta_real": effects[real_id],
            "beta_pseudo": effects[pseudo_id],
            "delta": effects[real_id] - effects[pseudo_id],
            "inferential_unit": "TRANSCRIPT_MEAN",
            "pair_level_role": "DESCRIPTIVE_ONLY",
        })
    return rows


def _localization_beta_rows(
    results: Sequence[Mapping[str, str]],
) -> list[dict[str, object]]:
    rows = []
    for row in results:
        family = str(row["analysis_family"])
        if not family.startswith("center_"):
            continue
        rows.append({
            "feature_id": row["feature_id"],
            "transcript": row["transcript"],
            "center": int(family.removeprefix("center_")),
            "beta": float(row["logFC"]),
        })
    return rows


def summarize_reagent_omissions(
    primary_results: Sequence[Mapping[str, str]],
    reagent_results: Sequence[Mapping[str, str]],
) -> dict[str, object]:
    full = {str(row["feature_id"]): float(row["logFC"]) for row in primary_results}
    full_median = float(np.median(list(full.values())))
    summaries = []
    for family in ("omit_siP1_1", "omit_siP2_1_siP2_4"):
        subset = {
            str(row["feature_id"]): float(row["logFC"])
            for row in reagent_results
            if row["analysis_family"] == family
        }
        common = sorted(set(full) & set(subset))
        full_values = [full[key] for key in common]
        subset_values = [subset[key] for key in common]
        subset_median = float(np.median(subset_values))
        sign_concordance = float(np.mean(
            np.sign(full_values) == np.sign(subset_values)
        ))
        summaries.append({
            "omission": family,
            "role": "DESCRIPTIVE_ONLY",
            "tested_anchors": len(subset),
            "common_with_full": len(common),
            "median_logFC": subset_median,
            "spearman_with_full": spearman(full_values, subset_values),
            "sign_concordance": sign_concordance,
            "sign_reversal": bool(
                np.sign(subset_median) != np.sign(full_median)
                and subset_median != 0
                and full_median != 0
            ),
            "reagent_specific_p_fdr_reported": "NO",
        })
    return {
        "full_model_median_logFC": full_median,
        "omissions": summaries,
        "any_sign_reversal": any(row["sign_reversal"] for row in summaries),
    }


def benchmark_decision(
    real_pseudo: Mapping[str, object],
    coordinate: Mapping[str, object],
    localization: Mapping[str, object],
    transcript_loo_pass: bool,
    gene_loo_pass: bool,
    reagent: Mapping[str, object],
) -> dict[str, object]:
    required = {
        "A_real_pseudo": bool(real_pseudo["core_pass"]),
        "B_coordinate_null_v2": bool(coordinate["core_pass"]),
        "C_LDB": bool(localization["LDB_pass"]),
        "D_NGS": bool(localization["NGS_pass"]),
        "E_transcript_LOO": transcript_loo_pass,
        "E_gene_LOO": gene_loo_pass,
    }
    if not all(required.values()):
        decision = "FAIL"
    elif bool(reagent["any_sign_reversal"]):
        decision = "PASS WITH CAUTION"
    else:
        decision = "PASS"
    return {
        "required_intersection": required,
        "reagent_sign_reversal": bool(reagent["any_sign_reversal"]),
        "minimum_primary_fdr_discoveries_required": "NONE",
        "decision": decision,
    }


def _jsonable(value: object) -> object:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.write_text(
        json.dumps(_jsonable(dict(payload)), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def run_analysis_variant(
    repo: Path,
    variant_dir: Path,
    variant: str,
    anchors: Sequence[Anchor],
    pairs: Sequence[PseudoPair],
    features: Mapping[str, Feature],
    genes: Mapping[str, str],
    masks: Mapping[str, Sequence[tuple[int, int]]],
    tracks: Mapping[str, Mapping[str, Mapping[str, SparseTrack]]],
    all_tmd_intervals: Mapping[str, Sequence[tuple[str, int, int]]],
) -> dict[str, object]:
    variant_dir.mkdir()
    primary_rows, positional_rows, common_rows, common = (
        build_primary_and_positional_rows(
        anchors, features, genes, masks, tracks, variant, all_tmd_intervals
        )
    )
    if not primary_rows or not positional_rows or not common_rows:
        raise ValueError(f"No eligible primary/positional rows for {variant}")
    primary_model_rows = _model_rows_with_family(primary_rows, "primary_+75")
    primary_results, primary_fit_summary = run_edger(
        repo, primary_model_rows, variant_dir / "primary"
    )
    positional_model_rows = _model_rows_with_family(positional_rows)
    positional_results, positional_fit_summary = run_edger(
        repo, positional_model_rows, variant_dir / "positional"
    )
    global_fdr = bh_adjust([float(row["PValue"]) for row in positional_results])
    positional_with_global = []
    for row, fdr in zip(positional_results, global_fdr):
        item = dict(row)
        item["complete_anchor_center_BH_FDR"] = fdr
        positional_with_global.append(item)
    write_tsv(
        variant_dir / "positional.edger_results_with_global_fdr.tsv",
        positional_with_global,
    )
    common_results, common_fit_summary = run_edger(
        repo,
        _model_rows_with_family(common_rows),
        variant_dir / "localization_common_all_center",
    )

    pseudo_model_rows, eligible_pairs = build_pseudo_rows(
        pairs, features, genes, masks, tracks, variant
    )
    pseudo_results, pseudo_fit_summary = run_edger(
        repo,
        _model_rows_with_family(
            pseudo_model_rows,
            "real_pseudo",
            report_inference="NO",
        ),
        variant_dir / "real_pseudo",
    )
    pair_effects = _pair_effect_rows(pseudo_results)
    assert_equal(len(pair_effects), len(eligible_pairs), "eligible pair fit count")
    write_tsv(variant_dir / "real_pseudo_pair_deltas.tsv", pair_effects)
    real_pseudo = real_pseudo_statistic(pair_effects)
    write_json(
        variant_dir / "real_pseudo_summary.json",
        {key: value for key, value in real_pseudo.items() if key != "D_by_transcript"},
    )
    write_tsv(
        variant_dir / "real_pseudo_transcript_D.tsv",
        [
            {"transcript": transcript, "D_t": value}
            for transcript, value in real_pseudo["D_by_transcript"].items()
        ],
    )

    localization = localization_statistics(
        _localization_beta_rows(common_results)
    )
    write_json(
        variant_dir / "localization_summary.json",
        {
            key: value
            for key, value in localization.items()
            if key not in {"P_by_transcript", "Q_by_transcript"}
        },
    )
    write_tsv(
        variant_dir / "localization_transcript_statistics.tsv",
        [
            {
                "transcript": transcript,
                "P_t": localization["P_by_transcript"][transcript],
                "Q_t": localization["Q_by_transcript"][transcript],
                "P_minus_Q": (
                    localization["P_by_transcript"][transcript]
                    - localization["Q_by_transcript"][transcript]
                ),
            }
            for transcript in sorted(localization["P_by_transcript"])
        ],
    )

    coordinate = compute_coordinate_null_v2(
        repo, features, masks, tracks, variant, variant_dir
    )
    write_json(variant_dir / "coordinate_null_v2_summary.json", coordinate["summary"])
    write_tsv(
        variant_dir / "coordinate_null_v2_statistics.tsv",
        [
            {"realization_id": index + 1, "S_null": value}
            for index, value in enumerate(
                np.asarray(coordinate["global_null_sum"], dtype=float)
                / COORDINATE_NULL_ANCHORS
            )
        ],
    )

    transcript_loo, transcript_loo_pass = cluster_loo(
        real_pseudo, coordinate, genes, "transcript"
    )
    gene_loo, gene_loo_pass = cluster_loo(
        real_pseudo, coordinate, genes, "gene"
    )
    write_tsv(variant_dir / "transcript_loo.tsv", transcript_loo)
    write_tsv(variant_dir / "gene_loo.tsv", gene_loo)
    anchor_loo = anchor_loo_descriptive(pair_effects, real_pseudo, coordinate)
    write_json(variant_dir / "anchor_loo_descriptive.json", anchor_loo)

    reagent_rows = build_reagent_omission_rows(
        anchors, features, genes, masks, tracks, variant
    )
    reagent_results, reagent_fit_summary = run_edger(
        repo, reagent_rows, variant_dir / "reagent_omission"
    )
    reagent = summarize_reagent_omissions(primary_results, reagent_results)
    write_json(variant_dir / "reagent_omission_descriptive.json", reagent)

    benchmark = benchmark_decision(
        real_pseudo,
        coordinate["summary"],
        localization,
        transcript_loo_pass,
        gene_loo_pass,
        reagent,
    )
    write_json(variant_dir / "benchmark_gate.json", benchmark)
    primary_summary = {
        "starting_anchors": PRIMARY_STARTING_ANCHORS,
        "eligible_anchors": len(primary_results),
        "common_all_center_anchors": len(common),
        "target_center": PRIMARY_CENTER,
        "target_window": "[+60,+91)",
        "target_length": TARGET_LENGTH,
        "background_threshold_all_six": MINIMUM_BACKGROUND_COUNT,
        "design": "~ condition",
        "contrast": "Knockdown-Control",
        "BH_FDR_threshold": 0.05,
        "FDR_discoveries_required_for_gate": "NO",
        "edger_summary": primary_fit_summary,
    }
    write_json(variant_dir / "primary_summary.json", primary_summary)
    del coordinate["null_anchor"]
    del coordinate["null_transcript"]
    for path in coordinate["temporary_paths"]:
        Path(path).unlink(missing_ok=True)
    return {
        "variant": variant,
        "primary": primary_summary,
        "positional_fit_summary": positional_fit_summary,
        "localization_common_fit_summary": common_fit_summary,
        "pseudo_starting_pairs": PSEUDO_STARTING_PAIRS,
        "pseudo_eligible_pairs": len(pair_effects),
        "pseudo_fit_summary": pseudo_fit_summary,
        "real_pseudo_core_pass": real_pseudo["core_pass"],
        "coordinate_null_core_pass": coordinate["summary"]["core_pass"],
        "LDB_pass": localization["LDB_pass"],
        "NGS_pass": localization["NGS_pass"],
        "transcript_loo_pass": transcript_loo_pass,
        "gene_loo_pass": gene_loo_pass,
        "anchor_loo_role": "DESCRIPTIVE_ONLY",
        "reagent_fit_summary": reagent_fit_summary,
        "benchmark": benchmark,
    }


def plan_only_report(
    output_dir: Path | None = None,
) -> dict[str, object]:
    return {
        "backend_version": BACKEND_VERSION,
        "contract_version": CONTRACT_VERSION,
        "mode": "PLAN_ONLY",
        "biological_psite_rows_read": "NO",
        "biological_outputs_generated": "NO",
        "output_directory": str(output_dir.resolve()) if output_dir else None,
        "primary": {
            "libraries": list(SAMPLE_ORDER),
            "sample_labels": SAMPLE_LABELS,
            "design": "~ condition",
            "contrast": "Knockdown-Control",
            "center": PRIMARY_CENTER,
            "target": "[+60,+91)",
            "starting_anchors": PRIMARY_STARTING_ANCHORS,
            "background_count_threshold": MINIMUM_BACKGROUND_COUNT,
            "strict_six_library_intersection": "YES",
            "all_zero_target_excluded": "YES",
            "model": "edgeR robust quasi-likelihood",
            "BH_FDR": "<0.05",
        },
        "real_pseudo": {
            "starting_pairs": PSEUDO_STARTING_PAIRS,
            "identical_pair_eligibility": "YES",
            "inferential_unit": "TRANSCRIPT",
        },
        "coordinate_null": {
            "version": COORDINATE_NULL_VERSION,
            "anchors": COORDINATE_NULL_ANCHORS,
            "transcripts": COORDINATE_NULL_TRANSCRIPTS,
            "realizations": COORDINATE_NULL_REALIZATIONS,
            "matrix_sha256": COORDINATE_NULL_SHA256,
            "v1_prohibited": "YES",
        },
        "localization": {"LDB": "LDB-1", "NGS": "NGS-1"},
        "robustness": {
            "transcript_LOO": "MANDATORY",
            "gene_LOO": "MANDATORY",
            "anchor_LOO": "DESCRIPTIVE_ONLY",
            "reagent_omission": "DESCRIPTIVE_ONLY",
        },
        "technical_sensitivity": {
            "excluded_read_lengths": [28, 29],
            "scope": "ALL_SIX_LIBRARIES",
            "selective_GSM3900210_filtering": "PROHIBITED",
        },
    }


def _output_checksums(root: Path) -> list[dict[str, str]]:
    return [
        {"path": str(path.relative_to(root)), "sha256": sha256(path)}
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name != "output_checksums.json"
    ]


def run_stage1_backend(
    repo: Path,
    stage0c_root: Path,
    manifest_rows: Sequence[Mapping[str, str]],
    output_dir: Path,
    preflight_report: Mapping[str, object],
) -> dict[str, object]:
    if preflight_report.get("input_preflight") != "PASS":
        raise ValueError("Stage 1 run requires a successful canonical preflight")
    if preflight_report.get("biological_psite_rows_read") != "NO":
        raise ValueError("Preflight crossed the biological-row boundary")
    if output_dir.exists():
        state = "nonempty" if output_dir.is_dir() and any(output_dir.iterdir()) else "existing"
        raise FileExistsError(
            f"Refusing to overwrite {state} Stage 1 output path: {output_dir}"
        )
    validate_sensitivity_exclusions({
        gsm: MANDATORY_SENSITIVITY_EXCLUDED_LENGTHS for gsm in SAMPLE_ORDER
    })
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{output_dir.name}.stage1.", dir=output_dir.parent
    ))
    try:
        features_path = stage0c_root / "reference/author_reference/human_features.tsv"
        gene_map_path = (
            stage0c_root
            / "reference/author_reference/human_appris_data.principal.txt"
        )
        topology_path = (
            stage0c_root
            / "reference/deeptmhmm_raw/"
            "deeptmhmm_full_normalized_topology_segments.tsv"
        )
        anchor_path = stage0c_root / "reference/tmd/primary_tmd_anchors.tsv"
        pseudo_path = (
            stage0c_root / "reference/tmd/primary_tmd_pseudo_anchor_pairs.tsv"
        )
        features = load_features(features_path)
        genes = load_gene_map(gene_map_path, features)
        masks = load_topology_masks(topology_path, features)
        anchors = load_primary_anchors(anchor_path, features, pseudo_path)
        all_tmd_intervals = load_all_tmd_intervals(anchor_path)
        pairs = load_pseudo_pairs(pseudo_path)
        tracks = load_psite_tracks(manifest_rows, features)

        variants = []
        variants.append(run_analysis_variant(
            repo,
            temporary / "primary_all_frozen_lengths",
            "primary",
            anchors,
            pairs,
            features,
            genes,
            masks,
            tracks,
            all_tmd_intervals,
        ))
        variants.append(run_analysis_variant(
            repo,
            temporary / "mandatory_sensitivity_exclude_28_29_all_six",
            "exclude_28_29",
            anchors,
            pairs,
            features,
            genes,
            masks,
            tracks,
            all_tmd_intervals,
        ))
        report = {
            "backend_version": BACKEND_VERSION,
            "contract_version": CONTRACT_VERSION,
            "preflight": preflight_report,
            "biological_execution": "COMPLETED",
            "primary_read_lengths": list(FROZEN_READ_LENGTHS),
            "mandatory_sensitivity_excluded_read_lengths": [28, 29],
            "mandatory_sensitivity_scope": "ALL_SIX_LIBRARIES",
            "selective_GSM3900210_filtering": "PROHIBITED",
            "variants": variants,
        }
        write_json(temporary / "run_report.json", report)
        write_json(
            temporary / "output_checksums.json",
            {"files": _output_checksums(temporary)},
        )
        write_json(
            temporary / "_SUCCESS.json",
            {
                "backend_version": BACKEND_VERSION,
                "contract_version": CONTRACT_VERSION,
                "atomic_completion": "PASS",
            },
        )
        os.replace(temporary, output_dir)
        return report
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
