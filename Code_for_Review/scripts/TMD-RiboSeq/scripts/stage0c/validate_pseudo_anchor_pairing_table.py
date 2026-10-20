#!/usr/bin/env python3
"""Independently validate a frozen annotation-only pseudo-anchor pairing table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_only_pseudo_anchor_feasibility import (  # noqa: E402
    CENTERS,
    MEMBRANE_FEATURE_LABELS,
    PRIMARY_CENTER,
    build_mask,
    overlaps,
)


EXPECTED_ROWS = 19_784
EXPECTED_COUNTS = {
    "MATCHED": 9_829,
    "PRIMARY_WINDOW_INVALID": 4_187,
    "NO_PSEUDO_CANDIDATE": 5_736,
    "ASSIGNMENT_CONFLICT": 32,
}
NA = "NA"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--repeat-table", type=Path)
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


def as_int(row: dict[str, str], field: str) -> int:
    return int(row[field])


def mask_string(mask: tuple[bool, ...]) -> str:
    return "".join("1" if value else "0" for value in mask)


def main() -> None:
    args = parse_args()
    table = read_tsv(args.table)
    anchors = read_tsv(args.anchors)
    topology = read_tsv(args.topology)
    features_rows = read_tsv(args.features)

    features = {
        row["Transcript"]: (int(row["UTR5_Length"]), int(row["CDS_Length"]))
        for row in features_rows
    }
    anchor_by_key = {(row["transcript"], row["TMD_index"]): row for row in anchors}
    accepted_by_transcript: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in anchors:
        accepted_by_transcript[row["transcript"]].append(
            (int(row["cds_nt_start0"]), int(row["cds_nt_end0_exclusive"]))
        )
    for transcript in accepted_by_transcript:
        accepted_by_transcript[transcript].sort()

    topology_by_transcript: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for row in topology:
        if row["raw_deeptmhmm_label"] in MEMBRANE_FEATURE_LABELS:
            topology_by_transcript[row["transcript"]].append(
                (
                    int(row["normalized_aa_start0"]),
                    int(row["normalized_aa_end0_exclusive"]),
                    row["raw_deeptmhmm_label"],
                )
            )

    violations = Counter()
    statuses = Counter(row["match_status"] for row in table)
    keys = [(row["transcript"], row["TMD_index"]) for row in table]
    if len(table) != EXPECTED_ROWS:
        violations["row_count"] += abs(len(table) - EXPECTED_ROWS) or 1
    violations["duplicate_or_missing_real_keys"] = EXPECTED_ROWS - len(set(keys))
    if set(keys) != set(anchor_by_key):
        violations["real_key_set"] += len(set(keys) ^ set(anchor_by_key))
    if dict(statuses) != EXPECTED_COUNTS:
        violations["frozen_status_counts"] += sum(
            abs(statuses.get(status, 0) - count) for status, count in EXPECTED_COUNTS.items()
        )

    primary_index = CENTERS.index(PRIMARY_CENTER)
    selected_pseudo_anchors: set[tuple[str, int]] = set()
    pseudo_fields = (
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
    )

    for row in table:
        key = (row["transcript"], row["TMD_index"])
        anchor = anchor_by_key.get(key)
        if anchor is None:
            continue
        transcript = row["transcript"]
        cds_start, cds_length = features[transcript]
        real_index0 = int(row["TMD_index"]) - 1
        accepted = accepted_by_transcript[transcript]
        real_mask = build_mask(as_int(row, "real_anchor_cds0"), cds_length, accepted, real_index0)

        real_equal_fields = {
            "protein_length": "protein_length",
            "protein_sha256": "protein_sha256",
            "n_primary_TMDs": "n_primary_TMDs",
            "pass_class": "pass_class",
            "TMD_position_class": "TMD_position_class",
            "real_TMD_aa_start0": "aa_start0",
            "real_TMD_aa_end0_exclusive": "aa_end0_exclusive",
            "real_TMD_cds_start0": "cds_nt_start0",
            "real_TMD_cds_end0_exclusive": "cds_nt_end0_exclusive",
            "real_anchor_coordinate": "primary_anchor_coordinate",
        }
        for table_field, anchor_field in real_equal_fields.items():
            if row[table_field] != anchor[anchor_field]:
                violations["real_anchor_field_mismatch"] += 1
        if as_int(row, "real_anchor_cds0") != int(anchor["cds_nt_start0"]):
            violations["real_anchor_cds_mismatch"] += 1
        if as_int(row, "real_TMD_length_nt") != int(anchor["cds_nt_end0_exclusive"]) - int(anchor["cds_nt_start0"]):
            violations["real_tmd_length_field"] += 1
        if as_int(row, "CDS_length_nt") != cds_length:
            violations["cds_length_field"] += 1
        if row["real_scan_validity_mask"] != mask_string(real_mask):
            violations["real_mask"] += 1
        expected_real_primary = "YES" if real_mask[primary_index] else "NO"
        if row["real_primary_window_valid"] != expected_real_primary:
            violations["real_primary_validity"] += 1

        if row["match_status"] != "MATCHED":
            if any(row[field] != NA for field in pseudo_fields):
                violations["unmatched_pseudo_not_na"] += 1
            if row["unmatched_reason"] != row["match_status"]:
                violations["unmatched_reason"] += 1
            continue

        if row["unmatched_reason"] != NA:
            violations["matched_unmatched_reason"] += 1
        if row["pseudo_transcript"] != transcript:
            violations["same_transcript"] += 1

        pseudo_tx_start = as_int(row, "pseudo_segment_start0")
        pseudo_tx_end = as_int(row, "pseudo_segment_end0_exclusive")
        pseudo_cds_start = as_int(row, "pseudo_segment_cds_start0")
        pseudo_cds_end = as_int(row, "pseudo_segment_cds_end0_exclusive")
        pseudo_aa_start = as_int(row, "pseudo_segment_aa_start0")
        pseudo_aa_end = as_int(row, "pseudo_segment_aa_end0_exclusive")
        if as_int(row, "pseudo_anchor_coordinate") != pseudo_tx_start:
            violations["pseudo_anchor_start"] += 1
        if as_int(row, "pseudo_anchor_cds0") != pseudo_cds_start:
            violations["pseudo_anchor_cds_start"] += 1
        if pseudo_tx_start != cds_start + pseudo_cds_start or pseudo_tx_end != cds_start + pseudo_cds_end:
            violations["pseudo_transcript_coordinate"] += 1
        if pseudo_cds_start != 3 * pseudo_aa_start or pseudo_cds_end != 3 * pseudo_aa_end:
            violations["codon_alignment"] += 1
        pseudo_length = pseudo_cds_end - pseudo_cds_start
        if pseudo_length != as_int(row, "pseudo_segment_length_nt"):
            violations["pseudo_length_field"] += 1
        if pseudo_length != as_int(row, "real_TMD_length_nt"):
            violations["tmd_length_match"] += 1
        if pseudo_cds_start < 0 or pseudo_cds_end > cds_length or pseudo_cds_start >= pseudo_cds_end:
            violations["pseudo_bounds"] += 1
        if any(overlaps(pseudo_aa_start, pseudo_aa_end, start, end) for start, end, _label in topology_by_transcript[transcript]):
            violations["topology_overlap"] += 1
        if pseudo_cds_start in {start for start, _end in accepted}:
            violations["real_tmd_start_collision"] += 1

        pseudo_mask = build_mask(pseudo_cds_start, cds_length, accepted, None)
        if row["pseudo_scan_validity_mask"] != mask_string(pseudo_mask):
            violations["pseudo_mask"] += 1
        if pseudo_mask != real_mask:
            violations["validity_mask"] += 1
        if not real_mask[primary_index] or not pseudo_mask[primary_index]:
            violations["primary_validity"] += 1
        if row["pseudo_primary_window_valid"] != ("YES" if pseudo_mask[primary_index] else "NO"):
            violations["pseudo_primary_validity_field"] += 1

        expected_metric = abs(pseudo_cds_start - as_int(row, "real_anchor_cds0"))
        if as_int(row, "positional_mismatch_nt") != expected_metric:
            violations["positional_metric"] += 1
        pseudo_key = (transcript, as_int(row, "pseudo_anchor_coordinate"))
        if pseudo_key in selected_pseudo_anchors:
            violations["pseudo_anchor_reuse"] += 1
        selected_pseudo_anchors.add(pseudo_key)

    for name in (
        "same_transcript",
        "tmd_length_match",
        "codon_alignment",
        "validity_mask",
        "topology_overlap",
        "real_tmd_start_collision",
        "primary_validity",
        "pseudo_anchor_reuse",
        "positional_metric",
        "unmatched_pseudo_not_na",
    ):
        violations[name] += 0

    repeat_sha = None
    byte_identical_repeat = None
    if args.repeat_table is not None:
        repeat_sha = sha256(args.repeat_table)
        byte_identical_repeat = args.table.read_bytes() == args.repeat_table.read_bytes()
        if not byte_identical_repeat:
            violations["deterministic_repeat"] += 1

    result = {
        "scope": "annotation_only_no_riboseq",
        "table": {"path": str(args.table), "sha256": sha256(args.table)},
        "repeat_table_sha256": repeat_sha,
        "counts": {
            "rows": len(table),
            "primary_window_eligible": sum(row["real_primary_window_valid"] == "YES" for row in table),
            "with_pseudo_candidate": sum(int(row["pseudo_candidate_count"]) > 0 for row in table),
            **dict(statuses),
        },
        "violations": dict(sorted(violations.items())),
        "invariants": {
            "count_conservation": sum(statuses.values()) == EXPECTED_ROWS,
            "frozen_counts_reproduced": dict(statuses) == EXPECTED_COUNTS,
            "deterministic_byte_identical_repeat": byte_identical_repeat,
            "all_violations_zero": sum(violations.values()) == 0,
            "riboseq_used": False,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if sum(violations.values()) != 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
