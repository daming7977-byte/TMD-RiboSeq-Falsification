#!/usr/bin/env python3
"""Paper 3 Stage 1 adapter frozen against the 2026-08-25 analysis contract.

This program separates code/contract validation from biological execution.  The
``validate`` and ``synthetic-test`` commands never read P-site tables.  The
``preflight`` command verifies a six-file biological input manifest using only
paths, P-site headers/hashes, and sibling technical-provenance tables.  Only
``run`` may parse P-site rows.

Adapter version: stage1_paper3_adapter_v2.0.0_2026-08-26
Contract: stage0c_paper3_analysis_contract_v1_2026-08-25
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from paper3_stage1_backend import plan_only_report, run_stage1_backend

ADAPTER_VERSION = "stage1_paper3_adapter_v2.0.0_2026-08-26"
CONTRACT_VERSION = "stage0c_paper3_analysis_contract_v1_2026-08-25"
EXPECTED_STAGE0C_MODULE_VERSION = "0.2.0"
EXPECTED_PSITE_FREEZE_VERSION = "stage0c_psite_v1_2026-08-24"
EXPECTED_FROZEN_READ_LENGTHS = "27,28,29,30,32,33,34,35"
EXPECTED_UNRESOLVED_READ_LENGTHS = "25,26,31"
EXPECTED_PHASE_WARNING_POLICY = (
    "RECORDED_TECHNICAL_WARNING_VALID_OUTPUTS_RETAINED"
)
BLINDED_AMENDMENT_RELATIVE_PATH = (
    "docs/STAGE0C_MODULE04_BLINDED_TECHNICAL_AMENDMENT_2026-08-26.md"
)
BLINDED_AMENDMENT_SHA256 = (
    "80e9a3a027e17a6278132ca61644d9ceba914df699727d881741eb1c9406f6f2"
)
CENTERS = (-30, -15, 0, 15, 30, 45, 60, 75, 90, 105, 120)
PRIMARY_CENTER = 75
TARGET_LEFT = -15
TARGET_RIGHT = 16
TARGET_LENGTH = 31
BACKGROUND_MIN = 20
TOPOLOGY_CLASSES = frozenset(("TMhelix", "signal", "Beta sheet"))
EXPECTED_GSMS = (
    "GSM3900208", "GSM3900209", "GSM3900210",
    "GSM3900211", "GSM3900212", "GSM3900213",
)
EXPECTED_CONDITIONS = ("Control", "Control", "Control", "Knockdown", "Knockdown", "Knockdown")
EXPECTED_PSITE_SHA256 = {
    "GSM3900208": "ffb3ba4242f9d84837ec0f4bc0e12588c90c77d128e763a173d9938b819cb1e8",
    "GSM3900209": "412d707010605fd3607cbf0a3d48016ec3740d1741ca0cb32a8d21d008f453a5",
    "GSM3900210": "cc8bd02f4e93c88dd90adacd4858ee612cccde43b0161c51f32eccbfd0d20d17",
    "GSM3900211": "cfd418d3c3fcd2d693679c1b0d094d4198eb317a843dca877e03ec47be142340",
    "GSM3900212": "696b449a1e0d4bcc6db75f42241bd3a05b05868118eb7934312e914fcd805e77",
    "GSM3900213": "16bedf4d98fec82c5fec0ef8493cb16321e7ddd77ac11cf9d29abdd38b39cd64",
}
MANDATORY_SENSITIVITY_EXCLUSIONS = {
    gsm: (28, 29) for gsm in EXPECTED_GSMS
}
CONTRACT_SHA = ""  # Contract TSV is validated through its frozen constituent fields below.
DESIGN_SHA = "86679bd767900a7d705867dcbfedb293f4520676c85ea854ddfc8d3b6102821d"
DECISION_SHA = "5009f1e2eccdbefaf7f86f3f909154014d590fb07970afa4682f52a35673c99a"
SESSION_SHA = "fbbdbfc873eba33f89a2f176f55927cc5f7f61bb3486f7401e7574bb81bd8ab4"
NULL_MATRIX_SHA = "7983d6dc1c382e57521e676669891a03271a4f1ffda0dda8d9ecf101e01071b3"
NULL_CONTENT_SHA = "c98f4b899a3f01398000b246b04c22002560b4b6f787f87801673b61e287aaf4"
FROZEN_ARTIFACTS = {
    "reference/author_reference/human_features.tsv": "c54965e9bca3cca852f422402180b8f0811a9ed58e8c9fd7a18b92b0e3099f36",
    "reference/deeptmhmm_raw/deeptmhmm_full_normalized_topology_segments.tsv": "ca60babf576a8547077d68be1f52d29835ce8829a5d2b671c4ef3fd175c94040",
    "reference/tmd/primary_tmd_anchors.tsv": "7257a07c51158ea4235c1c45a01ada899fcd78bcaddba68518bbb9a2146a92df",
    "reference/tmd/primary_tmd_pseudo_anchor_pairs.tsv": "1f6150d9985c82460d9f12ccc565edb416a1539ef70c761149f35b8fefba60a5",
}

REQUIRED_PSITE_COLUMNS = frozenset({
    "sample", "transcript", "psite_coordinate", "cds_start", "relative_psite",
    "freeze_version", "frozen_offset_source",
})
REQUIRED_RUN_STATS = {
    "module_version": EXPECTED_STAGE0C_MODULE_VERSION,
    "accounting_difference": "0",
    "frozen_read_lengths": EXPECTED_FROZEN_READ_LENGTHS,
    "unresolved_read_lengths": EXPECTED_UNRESOLVED_READ_LENGTHS,
    "phase0_not_dominant_policy": EXPECTED_PHASE_WARNING_POLICY,
    "TMD_blind": "YES",
    "benchmark_blind": "YES",
    "biological_result_blind": "YES",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_kv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"field", "value"}:
        raise ValueError(f"Invalid two-column contract table: {path}")
    if len({row['field'] for row in rows}) != len(rows):
        raise ValueError("Contract contains duplicate fields")
    return {row["field"]: row["value"] for row in rows}


def assert_equal(actual: object, expected: object, label: str) -> None:
    if actual != expected:
        raise ValueError(f"{label}: expected {expected!r}, observed {actual!r}")


def validate_contract(repo: Path, stage0c_root: Path | None = None) -> dict[str, object]:
    contract_path = repo / "docs/STAGE0C_PAPER3_ANALYSIS_CONTRACT_FREEZE_2026-08-25.tsv"
    design_path = repo / "docs/STAGE0C_GSE133111_DESIGN_MATRIX_2026-08-25.tsv"
    decision_path = repo / "docs/STAGE0C_PAPER3_STATISTICAL_DECISION_TABLE_2026-08-25.tsv"
    session_path = repo / "docs/STAGE0C_PAPER3_R_SESSIONINFO_2026-08-25.txt"
    values = read_kv(contract_path)
    required = {
        "freeze_version": CONTRACT_VERSION,
        "status": "PASS",
        "decision_table_total": "25",
        "decision_table_resolved": "25",
        "decision_table_unresolved": "0",
        "primary_model_libraries": "6",
        "model_formula": "~ condition",
        "reference_level": "Control",
        "coefficient": "conditionKnockdown",
        "primary_starting_anchors": "15597",
        "primary_center": "+75",
        "primary_target": "[+60,+91)",
        "target_length_nt": "31",
        "ordinary_background_topology_classes": "TMhelix|signal|Beta sheet",
        "ordinary_background_exclusion": "[A-30,max(E,A+136))",
        "background_count_threshold": "20",
        "coverage_intersection": "STRICT_ALL_SIX_RIBO_LIBRARIES",
        "all_zero_target": "EXCLUDE_FROM_INFERENCE_RETAIN_IN_ACCOUNTING",
        "model_family": "edgeR_negative_binomial_robust_QL",
        "offset": "log((background_count+0.5)*31/background_length)",
        "pseudo_annotation_pairs": "9829",
        "pseudo_focal_exclusion": "[A-30,max(E,A+136))",
        "coordinate_null_authoritative_version": "V2",
        "coordinate_null_anchors": "13060",
        "coordinate_null_transcripts": "3997",
        "coordinate_null_realizations": "10000",
        "coordinate_null_reproduction_mismatches": "0",
        "coordinate_null_canonical_content_sha256": NULL_CONTENT_SHA,
        "coordinate_null_matrix_file_sha256": NULL_MATRIX_SHA,
        "coordinate_null_v1": "SUPERSEDED_BEFORE_UNBLINDING_PROHIBITED_FOR_INFERENCE",
        "ldb": "LDB-1",
        "ldb_band": "+60|+75|+90",
        "ngs": "NGS-1",
        "ngs_epsilon": "log2(1.20)",
        "reagent_specific_p_fdr": "PROHIBITED",
        "transcript_loo": "MANDATORY",
        "gene_loo": "MANDATORY_19192_OF_19192_MAPPING",
        "anchor_loo": "DESCRIPTIVE_ONLY",
        "r_version": "4.6.0",
        "edgeR_version": "4.10.3",
        "limma_version": "3.68.5",
        "biological_results_inspected": "NO",
        "stage1_analysis_run": "NO",
    }
    for key, expected in required.items():
        assert_equal(values.get(key), expected, f"contract.{key}")
    assert_equal(sha256(design_path), DESIGN_SHA, "design SHA256")
    assert_equal(sha256(decision_path), DECISION_SHA, "decision-table SHA256")
    assert_equal(sha256(session_path), SESSION_SHA, "sessionInfo SHA256")
    null_matrix = repo / "docs/stage0c/coordinate_null_v2/coordinate_null_coordinates.npy"
    assert_equal(sha256(null_matrix), NULL_MATRIX_SHA, "authoritative v2 matrix SHA256")

    with design_path.open(newline="", encoding="utf-8") as handle:
        design = list(csv.DictReader(handle, delimiter="\t"))
    assert_equal(len(design), 6, "design row count")
    assert_equal(tuple(row["GSM"] for row in design), EXPECTED_GSMS, "design GSM order")
    assert_equal(tuple(row["condition"] for row in design), EXPECTED_CONDITIONS, "design conditions")
    if any(row["assay"] != "Ribo-seq" or row["primary_model_included"] != "YES" for row in design):
        raise ValueError("Design includes a non-Ribo-seq or excluded library")
    artifact_count = 0
    if stage0c_root is not None:
        for relative, expected_sha in FROZEN_ARTIFACTS.items():
            artifact = stage0c_root / relative
            assert_equal(sha256(artifact), expected_sha, f"frozen artifact SHA256 ({artifact})")
            artifact_count += 1
    return {
        "adapter_version": ADAPTER_VERSION,
        "contract_version": CONTRACT_VERSION,
        "contract_validation": "PASS",
        "decision_rows_resolved": "25/25",
        "design_libraries": 6,
        "coordinate_null_authoritative_version": "V2",
        "coordinate_null_v2_matrix_sha256": NULL_MATRIX_SHA,
        "frozen_external_artifacts_validated": artifact_count,
    }


@dataclass(frozen=True)
class Interval:
    start: int
    end: int


def clip_interval(start: int, end: int, length: int) -> Interval | None:
    start, end = max(0, start), min(length, end)
    return Interval(start, end) if start < end else None


def union_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    ordered = sorted(intervals, key=lambda x: (x.start, x.end))
    merged: list[Interval] = []
    for item in ordered:
        if not merged or item.start > merged[-1].end:
            merged.append(item)
        else:
            merged[-1] = Interval(merged[-1].start, max(merged[-1].end, item.end))
    return merged


def target_interval(anchor: int, center: int, cds_length: int) -> Interval | None:
    return clip_interval(anchor + center + TARGET_LEFT, anchor + center + TARGET_RIGHT, cds_length)


def focal_interval(anchor: int, end: int, cds_length: int) -> Interval | None:
    return clip_interval(anchor - 30, max(end, anchor + 136), cds_length)


def background_intervals(cds_length: int, exclusions: Sequence[Interval], target: Interval) -> list[Interval]:
    masked = union_intervals([*exclusions, target])
    result: list[Interval] = []
    cursor = 0
    for item in masked:
        if cursor < item.start:
            result.append(Interval(cursor, item.start))
        cursor = max(cursor, item.end)
    if cursor < cds_length:
        result.append(Interval(cursor, cds_length))
    return result


def count_sorted(positions: Sequence[int], intervals: Sequence[Interval]) -> int:
    import bisect
    return sum(bisect.bisect_left(positions, x.end) - bisect.bisect_left(positions, x.start) for x in intervals)


def synthetic_test(repo: Path, stage0c_root: Path | None = None) -> dict[str, object]:
    report = validate_contract(repo, stage0c_root)
    # Exact half-open primary geometry.
    assert_equal(target_interval(100, 75, 500), Interval(160, 191), "primary target geometry")
    assert_equal(focal_interval(100, 130, 500), Interval(70, 236), "focal geometry")
    # Adjacent/overlapping masks merge, target is removed, and counts respect half-open ends.
    masks = union_intervals([Interval(10, 20), Interval(20, 30), Interval(50, 60)])
    assert_equal(masks, [Interval(10, 30), Interval(50, 60)], "mask union")
    bg = background_intervals(100, masks, Interval(35, 40))
    assert_equal(sum(x.end - x.start for x in bg), 65, "background length")
    positions = [0, 9, 10, 19, 20, 29, 30, 34, 35, 39, 40, 59, 60, 99]
    assert_equal(count_sorted(positions, [Interval(10, 20)]), 2, "half-open count")
    assert_equal(count_sorted(positions, bg), 7, "background count")
    # Contract statistic algebra and empirical +1 correction.
    local = math.log2((((31 + 0.5) / 31) / ((69 + 0.5) / 65)))
    assert math.isfinite(local)
    null = [local - 1, local, local + 1]
    p = (sum(value >= local for value in null) + 1) / (len(null) + 1)
    assert_equal(p, 0.75, "right-tail empirical p")
    report.update({
        "synthetic_geometry": "PASS",
        "synthetic_counting": "PASS",
        "synthetic_statistic_algebra": "PASS",
        "biological_inputs_read": "NO",
    })
    return report


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"GSM", "psite_path", "input_role"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Input manifest requires columns {sorted(required)}")
    return rows


def sibling_provenance_paths(psite_path: Path) -> tuple[Path, Path]:
    suffix = "_calibrated_psites.tsv"
    if not psite_path.name.endswith(suffix):
        raise ValueError(
            f"Canonical P-site filename must end with {suffix}: {psite_path}"
        )
    prefix = psite_path.name[:-len(suffix)]
    return (
        psite_path.with_name(f"{prefix}_calibrated_psite_run_stats.tsv"),
        psite_path.with_name(f"{prefix}_calibrated_psite_phase_check.tsv"),
    )


def read_run_stats(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing sibling Module 04 run stats: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["statistic", "value"]:
            raise ValueError(f"Invalid Module 04 run-stat columns: {path}")
        rows = list(reader)
    if len({row["statistic"] for row in rows}) != len(rows):
        raise ValueError(f"Duplicate Module 04 run statistic: {path}")
    return {row["statistic"]: row["value"] for row in rows}


def validate_run_stats(
    path: Path,
    gsm: str,
    manifest_sample: str,
) -> dict[str, str]:
    values = read_run_stats(path)
    for key, expected in REQUIRED_RUN_STATS.items():
        assert_equal(values.get(key), expected, f"{gsm} run_stats.{key}")
    assert_equal(values.get("sample"), manifest_sample, f"{gsm} run_stats.sample")
    assert_equal(
        values.get("freeze_version"),
        EXPECTED_PSITE_FREEZE_VERSION,
        f"{gsm} run_stats.freeze_version",
    )
    if gsm == "GSM3900210":
        expected_warning_count = "2"
        expected_warning_lengths = "28,29"
    else:
        expected_warning_count = "0"
        expected_warning_lengths = "NONE"
    assert_equal(
        values.get("technical_warning_count"),
        expected_warning_count,
        f"{gsm} run_stats.technical_warning_count",
    )
    assert_equal(
        values.get("technical_warning_read_lengths"),
        expected_warning_lengths,
        f"{gsm} run_stats.technical_warning_read_lengths",
    )
    return values


def validate_gsm3900210_phase_check(path: Path, manifest_sample: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing sibling Module 04 phase check: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample", "read_length", "technical_check_status"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Invalid Module 04 phase-check columns: {path}")
    if any(row["sample"] != manifest_sample for row in rows):
        raise ValueError(f"GSM3900210 phase-check sample identity mismatch: {path}")
    observed_lengths = tuple(int(row["read_length"]) for row in rows)
    assert_equal(
        observed_lengths,
        tuple(int(value) for value in EXPECTED_FROZEN_READ_LENGTHS.split(",")),
        "GSM3900210 phase-check read lengths",
    )
    warning_lengths = tuple(
        int(row["read_length"])
        for row in rows
        if row["technical_check_status"] == "PHASE0_NOT_DOMINANT"
    )
    assert_equal(
        warning_lengths,
        (28, 29),
        "GSM3900210 PHASE0_NOT_DOMINANT read lengths",
    )
    for row in rows:
        expected = (
            "PHASE0_NOT_DOMINANT"
            if int(row["read_length"]) in (28, 29)
            else "PHASE0_DOMINANT"
        )
        assert_equal(
            row["technical_check_status"],
            expected,
            f"GSM3900210 phase-check length {row['read_length']}",
        )


def validate_mandatory_sensitivity_exclusions(
    exclusions: dict[str, Sequence[int]],
) -> None:
    assert_equal(
        tuple(exclusions),
        EXPECTED_GSMS,
        "mandatory sensitivity library order",
    )
    for gsm in EXPECTED_GSMS:
        assert_equal(
            tuple(exclusions[gsm]),
            (28, 29),
            f"mandatory sensitivity exclusions for {gsm}",
        )


def preflight_inputs(repo: Path, manifest: Path, stage0c_root: Path | None = None) -> dict[str, object]:
    report = validate_contract(repo, stage0c_root)
    amendment_path = repo / BLINDED_AMENDMENT_RELATIVE_PATH
    assert_equal(
        sha256(amendment_path),
        BLINDED_AMENDMENT_SHA256,
        "2026-08-26 blinded technical amendment SHA256",
    )
    validate_mandatory_sensitivity_exclusions(MANDATORY_SENSITIVITY_EXCLUSIONS)
    rows = read_manifest(manifest)
    assert_equal(tuple(row["GSM"] for row in rows), EXPECTED_GSMS, "manifest GSM order")
    if any(row["input_role"] != "FULL_FROZEN_CALIBRATED_PSITES" for row in rows):
        raise ValueError("Pilot/partial inputs are prohibited; every input_role must be FULL_FROZEN_CALIBRATED_PSITES")
    provenance: list[dict[str, str]] = []
    for row in rows:
        path = Path(row["psite_path"]).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Missing full calibrated P-site input for {row['GSM']}: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            columns = frozenset(next(csv.reader(handle, delimiter="\t"), ()))
        if not REQUIRED_PSITE_COLUMNS.issubset(columns):
            raise ValueError(
                f"{path} lacks required P-site columns: "
                f"{sorted(REQUIRED_PSITE_COLUMNS - columns)}"
            )
        observed_sha = sha256(path)
        assert_equal(
            observed_sha,
            EXPECTED_PSITE_SHA256[row["GSM"]],
            f"canonical P-site SHA256 ({row['GSM']})",
        )
        run_stats_path, phase_check_path = sibling_provenance_paths(path)
        run_stats = validate_run_stats(
            run_stats_path,
            row["GSM"],
            row.get("sample", row["GSM"]),
        )
        phase_check_provenance = None
        if row["GSM"] == "GSM3900210":
            validate_gsm3900210_phase_check(
                phase_check_path,
                row.get("sample", row["GSM"]),
            )
            phase_check_provenance = {
                "path": str(phase_check_path),
                "sha256": sha256(phase_check_path),
            }
        provenance.append({
            "GSM": row["GSM"],
            "path": str(path),
            "sha256": observed_sha,
            "run_stats_path": str(run_stats_path),
            "run_stats_sha256": sha256(run_stats_path),
            "module_version": run_stats["module_version"],
            "technical_warning_count": run_stats["technical_warning_count"],
            "technical_warning_read_lengths": run_stats[
                "technical_warning_read_lengths"
            ],
            "phase_check": phase_check_provenance,
        })
    report.update({
        "input_preflight": "PASS",
        "inputs": provenance,
        "biological_rows_inspected": "NO",
        "biological_psite_rows_read": "NO",
        "blinded_technical_amendment": {
            "date": "2026-08-26",
            "path": str(amendment_path.resolve()),
            "sha256": BLINDED_AMENDMENT_SHA256,
        },
        "mandatory_sensitivity": {
            "scope": "ALL_SIX_LIBRARIES",
            "excluded_read_lengths": [28, 29],
            "single_library_exclusion_prohibited": "YES",
        },
    })
    return report


def write_json(report: dict[str, object], output: Path | None) -> None:
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(payload, encoding="utf-8")
    sys.stdout.write(payload)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=("validate", "synthetic-test", "preflight", "plan-only", "run"),
    )
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="New, nonexistent atomic production directory required by run",
    )
    parser.add_argument("--stage0c-root", type=Path, default=Path("REPRO_WORKSPACE/GSE133111_stage0C"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.command == "validate":
            report = validate_contract(args.repo, args.stage0c_root)
        elif args.command == "synthetic-test":
            report = synthetic_test(args.repo, args.stage0c_root)
        else:
            if args.manifest is None:
                raise ValueError("--manifest is required for preflight/plan-only/run")
            report = preflight_inputs(args.repo, args.manifest, args.stage0c_root)
            if args.command == "plan-only":
                report = {
                    "preflight": report,
                    "execution_plan": plan_only_report(args.output_dir),
                }
            elif args.command == "run":
                if args.output_dir is None:
                    raise ValueError("--output-dir is required for run")
                report = run_stage1_backend(
                    args.repo,
                    args.stage0c_root,
                    read_manifest(args.manifest),
                    args.output_dir.resolve(),
                    report,
                )
        write_json(report, args.output)
        return 0
    except Exception as exc:
        sys.stderr.write(f"{args.command.upper()}_FAIL: {exc}\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
