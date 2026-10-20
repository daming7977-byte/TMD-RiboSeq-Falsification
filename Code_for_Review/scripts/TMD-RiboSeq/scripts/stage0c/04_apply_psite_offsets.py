#!/usr/bin/env python3
"""
Stage 0C Module 04: apply formally frozen, forward-only P-site offsets.

This module is a downstream integration step. It does not calibrate, select,
or freeze offsets and contains no TMD or benchmark logic.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

from lib.bam import open_bam
from lib.io import read_features


MODULE_VERSION = "0.2.0"
MIN_READ_LENGTH = 25
MAX_READ_LENGTH = 35
VALIDATION_FLANK = 15
EXPECTED_FREEZE_VERSION = "stage0c_psite_v1_2026-08-24"
DEFAULT_FROZEN_OFFSETS = (
    Path(__file__).resolve().parents[2]
    / "docs"
    / "STAGE0C_PSITE_FROZEN_OFFSETS_2026-08-24.tsv"
)

STATUS_CALIBRATED = "CALIBRATED_MATCH"
EXPECTED_UNRESOLVED_STATUSES = {
    25: "UNRESOLVED_AB_DISAGREEMENT",
    26: "UNRESOLVED_AB_DISAGREEMENT",
    31: "VALIDATION_WARNING",
}
EXPECTED_CALIBRATED_LENGTHS = {27, 28, 29, 30, 32, 33, 34, 35}
EXPECTED_FROZEN_OFFSETS = {
    27: 9,
    28: 9,
    29: 9,
    30: 12,
    32: 13,
    33: 13,
    34: 13,
    35: 13,
}

PHASE_DIAGNOSTIC_PROVENANCE = (
    "Full-library phase dominance is diagnostic after independent frozen "
    "calibration; PHASE0_NOT_DOMINANT is recorded as a technical warning and "
    "does not change offsets or exclude samples/read lengths."
)

EXCLUSION_FIELDS = [
    "unmapped_excluded",
    "secondary_excluded",
    "supplementary_excluded",
    "nh_missing_excluded",
    "multimapped_excluded",
    "missing_transcript_features",
    "invalid_read_length",
    "outside_calibration_length_range",
    "reverse_orientation_excluded",
    "unresolved_psite_length_excluded",
    "invalid_raw_five_prime_position",
    "invalid_psite_coordinate",
]

PSITE_FIELDS = [
    "sample",
    "bam_record_index",
    "read_name",
    "transcript",
    "read_length",
    "raw_5prime",
    "frozen_offset_distance",
    "psite_coordinate",
    "cds_start",
    "relative_psite",
    "freeze_version",
    "frozen_offset_source",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Apply the formally frozen forward-only Stage 0C P-site offsets "
            "without recalibration or biological analysis."
        )
    )
    parser.add_argument("--bam", required=True, help="Mapped transcriptome BAM")
    parser.add_argument(
        "--features",
        required=True,
        help="Approved human_features.tsv used in Stage 0C",
    )
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--sample",
        help="Optional sample identifier; otherwise derived from the BAM filename",
    )
    return parser.parse_args()


def derive_sample_name(bam_path):
    sample_name = Path(bam_path).stem
    for suffix in ("_nonrRNA_hisat2", "_hisat2", "_mapped"):
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]
    if not sample_name:
        raise ValueError(f"Could not derive sample name from BAM path: {bam_path}")
    return sample_name


def parse_required_integer(value, field, read_length):
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise RuntimeError(
            f"Length {read_length} has non-integer {field}: {value!r}"
        ) from error
    if parsed < 0:
        raise RuntimeError(
            f"Length {read_length} has negative {field}: {parsed}"
        )
    return parsed


def load_frozen_offsets(path):
    """Load and validate the sole downstream offset source."""

    source_path = Path(path)
    if not source_path.is_file():
        raise RuntimeError(f"Frozen-offset table does not exist: {source_path}")

    required_columns = {
        "read_length",
        "frozen_offset_distance",
        "calibration_status",
        "PILOT_A_offset",
        "PILOT_B_offset",
        "PILOT_A_validation_status",
        "PILOT_B_validation_status",
        "freeze_version",
    }
    rows = {}

    with source_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
            missing = sorted(required_columns - set(reader.fieldnames or []))
            raise RuntimeError(
                f"Frozen-offset table is missing required columns: {missing}"
            )

        for row in reader:
            read_length = parse_required_integer(
                row["read_length"], "read_length", row["read_length"]
            )
            if read_length in rows:
                raise RuntimeError(
                    f"Frozen-offset table contains duplicate length {read_length}"
                )
            rows[read_length] = row

    expected_lengths = set(range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1))
    if len(rows) != 11 or set(rows) != expected_lengths:
        raise RuntimeError(
            "Frozen-offset table must contain lengths 25-35 exactly once; "
            f"observed {sorted(rows)}"
        )

    freeze_versions = {
        (row.get("freeze_version") or "").strip() for row in rows.values()
    }
    if freeze_versions != {EXPECTED_FREEZE_VERSION}:
        raise RuntimeError(
            "Frozen-offset table has a missing or inconsistent freeze_version: "
            f"{sorted(freeze_versions)}"
        )

    calibrated_offsets = {}
    calibrated_lengths = set()

    for read_length in sorted(rows):
        row = rows[read_length]
        status = (row.get("calibration_status") or "").strip()
        frozen_value = (row.get("frozen_offset_distance") or "").strip()

        if status == STATUS_CALIBRATED:
            if frozen_value.upper() == "NA" or not frozen_value:
                raise RuntimeError(
                    f"CALIBRATED_MATCH length {read_length} lacks a numeric offset"
                )
            frozen_offset = parse_required_integer(
                frozen_value, "frozen_offset_distance", read_length
            )
            pilot_a_offset = parse_required_integer(
                row["PILOT_A_offset"], "PILOT_A_offset", read_length
            )
            pilot_b_offset = parse_required_integer(
                row["PILOT_B_offset"], "PILOT_B_offset", read_length
            )
            if frozen_offset != pilot_a_offset or frozen_offset != pilot_b_offset:
                raise RuntimeError(
                    f"Length {read_length} frozen offset does not match both "
                    "formal pilot offsets"
                )
            if (
                row["PILOT_A_validation_status"] != "PHASE0_DOMINANT"
                or row["PILOT_B_validation_status"] != "PHASE0_DOMINANT"
            ):
                raise RuntimeError(
                    f"CALIBRATED_MATCH length {read_length} lacks bilateral "
                    "PHASE0_DOMINANT validation"
                )
            calibrated_offsets[read_length] = frozen_offset
            calibrated_lengths.add(read_length)
        else:
            if frozen_value.upper() != "NA":
                raise RuntimeError(
                    f"Non-CALIBRATED_MATCH length {read_length} has a numeric offset"
                )

    if calibrated_lengths != EXPECTED_CALIBRATED_LENGTHS:
        raise RuntimeError(
            "Calibrated length set differs from the formal freeze: "
            f"{sorted(calibrated_lengths)}"
        )
    if calibrated_offsets != EXPECTED_FROZEN_OFFSETS:
        raise RuntimeError(
            "Frozen offsets differ from the authoritative "
            f"{EXPECTED_FREEZE_VERSION} mapping: {calibrated_offsets}"
        )

    for read_length, expected_status in EXPECTED_UNRESOLVED_STATUSES.items():
        row = rows[read_length]
        if row["calibration_status"] != expected_status:
            raise RuntimeError(
                f"Length {read_length} must remain {expected_status}; observed "
                f"{row['calibration_status']}"
            )
        if (row.get("frozen_offset_distance") or "").strip().upper() != "NA":
            raise RuntimeError(
                f"Unresolved length {read_length} must have frozen offset NA"
            )

    metadata = {
        read_length: {
            "calibration_status": row["calibration_status"],
            "frozen_offset_distance": (
                calibrated_offsets[read_length]
                if read_length in calibrated_offsets
                else "NA"
            ),
        }
        for read_length, row in rows.items()
    }
    return calibrated_offsets, metadata, EXPECTED_FREEZE_VERSION, source_path.resolve()


def get_read_length(read):
    if read.query_sequence is not None:
        return len(read.query_sequence)
    return read.infer_query_length()


def initialize_length_stats(offset_metadata):
    return {
        read_length: {
            "read_length": read_length,
            "calibration_status": offset_metadata[read_length][
                "calibration_status"
            ],
            "frozen_offset_distance": offset_metadata[read_length][
                "frozen_offset_distance"
            ],
            "eligible_before_psite_status": 0,
            "unresolved_psite_length_excluded": 0,
            "reverse_orientation_excluded": 0,
            "invalid_raw_five_prime_position": 0,
            "invalid_psite_coordinate": 0,
            "included_calibrated_psite_reads": 0,
        }
        for read_length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1)
    }


def create_temporary_tsv(output_dir, prefix, fieldnames):
    output_dir.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        newline="",
        encoding="utf-8",
        prefix=prefix,
        suffix=".tmp",
        dir=output_dir,
        delete=False,
    )
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    return handle, writer, Path(handle.name)


def process_bam(
    bam_path,
    features,
    sample_name,
    calibrated_offsets,
    offset_metadata,
    freeze_version,
    frozen_offset_source,
    output_dir,
):
    """Apply frozen offsets with one mutually exclusive outcome per BAM record."""

    exclusion_counts = Counter()
    total_bam_records = 0
    included_count = 0
    length_stats = initialize_length_stats(offset_metadata)
    phase_counts = {length: Counter() for length in calibrated_offsets}
    phase_stats = {
        length: {
            "no_valid_validation_interior": 0,
            "outside_validation_region": 0,
            "n_validation_reads": 0,
        }
        for length in calibrated_offsets
    }
    offset_usage = {length: Counter() for length in calibrated_offsets}

    output_dir = Path(output_dir)
    final_psite_path = output_dir / f"{sample_name}_calibrated_psites.tsv"
    handle, writer, temporary_psite_path = create_temporary_tsv(
        output_dir, f".{sample_name}_calibrated_psites.", PSITE_FIELDS
    )

    bam = None
    try:
        bam = open_bam(bam_path)
        for bam_record_index, read in enumerate(bam.fetch(until_eof=True), start=1):
            total_bam_records += 1

            if read.is_unmapped:
                exclusion_counts["unmapped_excluded"] += 1
                continue
            if read.is_secondary:
                exclusion_counts["secondary_excluded"] += 1
                continue
            if read.is_supplementary:
                exclusion_counts["supplementary_excluded"] += 1
                continue
            try:
                nh_tag = read.get_tag("NH")
            except KeyError:
                exclusion_counts["nh_missing_excluded"] += 1
                continue
            if nh_tag != 1:
                exclusion_counts["multimapped_excluded"] += 1
                continue

            transcript = read.reference_name
            if transcript is None or transcript not in features:
                exclusion_counts["missing_transcript_features"] += 1
                continue

            read_length = get_read_length(read)
            if read_length is None or read_length <= 0:
                exclusion_counts["invalid_read_length"] += 1
                continue
            if read_length < MIN_READ_LENGTH or read_length > MAX_READ_LENGTH:
                exclusion_counts["outside_calibration_length_range"] += 1
                continue

            stats = length_stats[read_length]
            stats["eligible_before_psite_status"] += 1

            if read.is_reverse:
                exclusion_counts["reverse_orientation_excluded"] += 1
                stats["reverse_orientation_excluded"] += 1
                continue
            if read_length not in calibrated_offsets:
                exclusion_counts["unresolved_psite_length_excluded"] += 1
                stats["unresolved_psite_length_excluded"] += 1
                continue

            raw_5prime = read.reference_start
            if raw_5prime is None or raw_5prime < 0:
                exclusion_counts["invalid_raw_five_prime_position"] += 1
                stats["invalid_raw_five_prime_position"] += 1
                continue

            feature = features[transcript]
            cds_start = feature["utr5"]
            cds_end_exclusive = cds_start + feature["cds"]
            transcript_length = cds_end_exclusive + feature["utr3"]
            frozen_offset = calibrated_offsets[read_length]
            psite_coordinate = raw_5prime + frozen_offset

            if psite_coordinate < 0 or psite_coordinate >= transcript_length:
                exclusion_counts["invalid_psite_coordinate"] += 1
                stats["invalid_psite_coordinate"] += 1
                continue

            relative_psite = psite_coordinate - cds_start
            read_name = read.query_name or f"BAM_RECORD_{bam_record_index}"
            writer.writerow({
                "sample": sample_name,
                "bam_record_index": bam_record_index,
                "read_name": read_name,
                "transcript": transcript,
                "read_length": read_length,
                "raw_5prime": raw_5prime,
                "frozen_offset_distance": frozen_offset,
                "psite_coordinate": psite_coordinate,
                "cds_start": cds_start,
                "relative_psite": relative_psite,
                "freeze_version": freeze_version,
                "frozen_offset_source": str(frozen_offset_source),
            })

            included_count += 1
            exclusion_counts["included_calibrated_psite_reads"] += 1
            stats["included_calibrated_psite_reads"] += 1
            offset_usage[read_length][frozen_offset] += 1

            validation_start = cds_start + VALIDATION_FLANK
            validation_end_exclusive = cds_end_exclusive - VALIDATION_FLANK
            if validation_start >= validation_end_exclusive:
                phase_stats[read_length]["no_valid_validation_interior"] += 1
            elif not validation_start <= raw_5prime < validation_end_exclusive:
                phase_stats[read_length]["outside_validation_region"] += 1
            else:
                phase = relative_psite % 3
                phase_counts[read_length][phase] += 1
                phase_stats[read_length]["n_validation_reads"] += 1
    except Exception:
        handle.close()
        temporary_psite_path.unlink(missing_ok=True)
        raise
    finally:
        if bam is not None:
            bam.close()

    handle.close()
    return {
        "total_bam_records": total_bam_records,
        "exclusion_counts": exclusion_counts,
        "included_count": included_count,
        "length_stats": length_stats,
        "phase_counts": phase_counts,
        "phase_stats": phase_stats,
        "offset_usage": offset_usage,
        "temporary_psite_path": temporary_psite_path,
        "final_psite_path": final_psite_path,
    }


def validate_accounting(results, calibrated_offsets):
    exclusion_counts = results["exclusion_counts"]
    total_excluded = sum(exclusion_counts[field] for field in EXCLUSION_FIELDS)
    included = exclusion_counts["included_calibrated_psite_reads"]
    accounted_for = total_excluded + included

    if included != results["included_count"]:
        raise RuntimeError("Included P-site count disagrees with streamed row count")
    if accounted_for != results["total_bam_records"]:
        raise RuntimeError(
            "P-site application count conservation failed: "
            f"{accounted_for:,} accounted != "
            f"{results['total_bam_records']:,} examined"
        )

    length_stats = results["length_stats"]
    for read_length, stats in length_stats.items():
        if read_length not in calibrated_offsets:
            if stats["included_calibrated_psite_reads"] != 0:
                raise RuntimeError(
                    f"Unresolved length {read_length} entered calibrated output"
                )
            if (
                stats["eligible_before_psite_status"]
                != stats["reverse_orientation_excluded"]
                + stats["unresolved_psite_length_excluded"]
            ):
                raise RuntimeError(
                    f"Unresolved length {read_length} accounting is inconsistent"
                )
        else:
            per_length_accounted = (
                stats["reverse_orientation_excluded"]
                + stats["invalid_raw_five_prime_position"]
                + stats["invalid_psite_coordinate"]
                + stats["included_calibrated_psite_reads"]
            )
            if per_length_accounted != stats["eligible_before_psite_status"]:
                raise RuntimeError(
                    f"Frozen length {read_length} accounting is inconsistent"
                )
            expected_usage = Counter()
            if stats["included_calibrated_psite_reads"]:
                expected_usage[calibrated_offsets[read_length]] = stats[
                    "included_calibrated_psite_reads"
                ]
            if results["offset_usage"][read_length] != expected_usage:
                raise RuntimeError(
                    f"Length {read_length} used an offset other than its frozen value"
                )

    if sum(
        stats["included_calibrated_psite_reads"]
        for stats in length_stats.values()
    ) != included:
        raise RuntimeError("Per-length included counts do not match sample total")
    if sum(
        stats["unresolved_psite_length_excluded"]
        for stats in length_stats.values()
    ) != exclusion_counts["unresolved_psite_length_excluded"]:
        raise RuntimeError("Per-length unresolved counts do not match sample total")
    if sum(
        stats["reverse_orientation_excluded"] for stats in length_stats.values()
    ) != exclusion_counts["reverse_orientation_excluded"]:
        raise RuntimeError("Per-length reverse counts do not match sample total")

    return total_excluded, accounted_for


def build_phase_rows(
    sample_name,
    calibrated_offsets,
    offset_metadata,
    phase_counts,
    phase_stats,
    freeze_version,
    frozen_offset_source,
):
    rows = []
    for read_length in sorted(calibrated_offsets):
        counts = phase_counts[read_length]
        n_validation_reads = sum(counts[phase] for phase in (0, 1, 2))
        stats = phase_stats[read_length]
        if n_validation_reads != stats["n_validation_reads"]:
            raise RuntimeError(
                f"Length {read_length} phase counts disagree with validation total"
            )

        if n_validation_reads == 0:
            percentages = ["NA", "NA", "NA"]
            dominant_phase = "NA"
            technical_check_status = "INSUFFICIENT_VALIDATION_EVIDENCE"
        else:
            percentages = [
                round(100.0 * counts[phase] / n_validation_reads, 2)
                for phase in (0, 1, 2)
            ]
            maximum = max(counts[phase] for phase in (0, 1, 2))
            dominant_phases = [
                phase for phase in (0, 1, 2) if counts[phase] == maximum
            ]
            dominant_phase = (
                dominant_phases[0] if len(dominant_phases) == 1 else "TIE"
            )
            technical_check_status = (
                "PHASE0_DOMINANT"
                if counts[0] > counts[1] and counts[0] > counts[2]
                else "PHASE0_NOT_DOMINANT"
            )

        rows.append({
            "sample": sample_name,
            "read_length": read_length,
            "calibration_status": offset_metadata[read_length][
                "calibration_status"
            ],
            "frozen_offset_distance": calibrated_offsets[read_length],
            "no_valid_validation_interior": stats[
                "no_valid_validation_interior"
            ],
            "outside_validation_region": stats["outside_validation_region"],
            "n_validation_reads": n_validation_reads,
            "phase0_count": counts[0],
            "phase1_count": counts[1],
            "phase2_count": counts[2],
            "phase0_pct": percentages[0],
            "phase1_pct": percentages[1],
            "phase2_pct": percentages[2],
            "dominant_phase": dominant_phase,
            "technical_check_status": technical_check_status,
            "freeze_version": freeze_version,
            "frozen_offset_source": str(frozen_offset_source),
        })
    return rows


def validate_phase_rows(phase_rows, calibrated_offsets, length_stats):
    technical_warnings = []
    if {row["read_length"] for row in phase_rows} != set(calibrated_offsets):
        raise RuntimeError("Phase-check rows do not match frozen read lengths")
    for row in phase_rows:
        read_length = row["read_length"]
        if row["frozen_offset_distance"] != calibrated_offsets[read_length]:
            raise RuntimeError(
                f"Length {read_length} phase check used a non-frozen offset"
            )
        phase_sum = row["phase0_count"] + row["phase1_count"] + row["phase2_count"]
        if phase_sum != row["n_validation_reads"]:
            raise RuntimeError(f"Length {read_length} phase accounting failed")
        validation_accounted = (
            row["no_valid_validation_interior"]
            + row["outside_validation_region"]
            + row["n_validation_reads"]
        )
        if (
            validation_accounted
            != length_stats[read_length]["included_calibrated_psite_reads"]
        ):
            raise RuntimeError(
                f"Length {read_length} validation-region accounting failed"
            )
        if row["technical_check_status"] == "PHASE0_NOT_DOMINANT":
            technical_warnings.append({
                "sample": row["sample"],
                "read_length": read_length,
                "technical_check_status": row["technical_check_status"],
            })
        elif row["technical_check_status"] != "PHASE0_DOMINANT":
            raise RuntimeError(
                f"Length {read_length} post-application phase check discrepancy: "
                f"{row['technical_check_status']}; frozen offsets were not changed"
            )
    return technical_warnings


def write_tsv_atomic(path, fieldnames, rows):
    path = Path(path)
    handle, writer, temporary_path = create_temporary_tsv(
        path.parent, f".{path.stem}.", fieldnames
    )
    try:
        writer.writerows(rows)
        handle.close()
        os.replace(temporary_path, path)
    except Exception:
        handle.close()
        temporary_path.unlink(missing_ok=True)
        raise


def finalize_outputs(
    results,
    sample_name,
    bam_path,
    features_path,
    offset_metadata,
    calibrated_offsets,
    freeze_version,
    frozen_offset_source,
    output_dir,
):
    total_excluded, accounted_for = validate_accounting(
        results, calibrated_offsets
    )
    phase_rows = build_phase_rows(
        sample_name,
        calibrated_offsets,
        offset_metadata,
        results["phase_counts"],
        results["phase_stats"],
        freeze_version,
        frozen_offset_source,
    )
    technical_warnings = validate_phase_rows(
        phase_rows,
        calibrated_offsets,
        results["length_stats"],
    )
    results["technical_warnings"] = technical_warnings

    output_dir = Path(output_dir)
    length_rows = []
    for read_length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1):
        row = dict(results["length_stats"][read_length])
        row["freeze_version"] = freeze_version
        row["frozen_offset_source"] = str(frozen_offset_source)
        length_rows.append(row)
    run_rows = [
        {"statistic": "sample", "value": sample_name},
        {"statistic": "module_version", "value": MODULE_VERSION},
        {"statistic": "bam", "value": str(Path(bam_path).resolve())},
        {"statistic": "features", "value": str(Path(features_path).resolve())},
        {"statistic": "freeze_version", "value": freeze_version},
        {
            "statistic": "frozen_offset_source",
            "value": str(frozen_offset_source),
        },
        {
            "statistic": "filter_order",
            "value": ",".join(EXCLUSION_FIELDS)
            + ",included_calibrated_psite_reads",
        },
        {
            "statistic": "total_bam_records_examined",
            "value": results["total_bam_records"],
        },
    ]
    for field in EXCLUSION_FIELDS:
        run_rows.append({"statistic": field, "value": results["exclusion_counts"][field]})
    run_rows.extend([
        {
            "statistic": "included_calibrated_psite_reads",
            "value": results["exclusion_counts"][
                "included_calibrated_psite_reads"
            ],
        },
        {"statistic": "total_excluded", "value": total_excluded},
        {"statistic": "accounted_for", "value": accounted_for},
        {
            "statistic": "accounting_difference",
            "value": results["total_bam_records"] - accounted_for,
        },
        {
            "statistic": "frozen_read_lengths",
            "value": ",".join(str(length) for length in sorted(calibrated_offsets)),
        },
        {
            "statistic": "unresolved_read_lengths",
            "value": ",".join(
                str(length)
                for length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1)
                if length not in calibrated_offsets
            ),
        },
        {
            "statistic": "full_library_phase_dominance_provenance",
            "value": PHASE_DIAGNOSTIC_PROVENANCE,
        },
        {
            "statistic": "phase0_not_dominant_policy",
            "value": "RECORDED_TECHNICAL_WARNING_VALID_OUTPUTS_RETAINED",
        },
        {
            "statistic": "technical_warning_count",
            "value": len(technical_warnings),
        },
        {
            "statistic": "technical_warning_read_lengths",
            "value": ",".join(
                str(warning["read_length"])
                for warning in technical_warnings
            ) or "NONE",
        },
        {"statistic": "TMD_blind", "value": "YES"},
        {"statistic": "benchmark_blind", "value": "YES"},
        {"statistic": "biological_result_blind", "value": "YES"},
    ])

    length_fields = list(length_rows[0].keys())
    phase_fields = list(phase_rows[0].keys())
    write_tsv_atomic(
        output_dir / f"{sample_name}_calibrated_psite_run_stats.tsv",
        ["statistic", "value"],
        run_rows,
    )
    write_tsv_atomic(
        output_dir / f"{sample_name}_calibrated_psite_length_stats.tsv",
        length_fields,
        length_rows,
    )
    write_tsv_atomic(
        output_dir / f"{sample_name}_calibrated_psite_phase_check.tsv",
        phase_fields,
        phase_rows,
    )
    os.replace(results["temporary_psite_path"], results["final_psite_path"])

    return [
        results["final_psite_path"],
        output_dir / f"{sample_name}_calibrated_psite_run_stats.tsv",
        output_dir / f"{sample_name}_calibrated_psite_length_stats.tsv",
        output_dir / f"{sample_name}_calibrated_psite_phase_check.tsv",
    ]


def main():
    args = parse_args()
    sample_name = args.sample or derive_sample_name(args.bam)

    print("=" * 60)
    print("Stage0C Module 04")
    print(f"Frozen P-Site Offset Application v{MODULE_VERSION}")
    print("=" * 60)
    print()
    print("TMD-BLIND: YES")
    print("BENCHMARK-BLIND: YES")
    print("BIOLOGICAL-RESULT-BLIND: YES")
    print("FORWARD-ONLY P-SITE APPLICATION: YES")
    print()

    print("[1/5] Loading and validating the authoritative frozen-offset table...")
    (
        calibrated_offsets,
        offset_metadata,
        freeze_version,
        frozen_offset_source,
    ) = load_frozen_offsets(DEFAULT_FROZEN_OFFSETS)

    print("[2/5] Loading approved human transcript features...")
    features = read_features(args.features)
    if not features:
        raise RuntimeError("Approved human feature table contains no transcripts")

    print("[3/5] Applying frozen offsets to eligible forward reads...")
    results = process_bam(
        args.bam,
        features,
        sample_name,
        calibrated_offsets,
        offset_metadata,
        freeze_version,
        frozen_offset_source,
        args.output,
    )

    print("[4/5] Validating accounting, offsets, and CDS-interior phase...")
    try:
        output_paths = finalize_outputs(
            results,
            sample_name,
            args.bam,
            args.features,
            offset_metadata,
            calibrated_offsets,
            freeze_version,
            frozen_offset_source,
            args.output,
        )
    except Exception:
        results["temporary_psite_path"].unlink(missing_ok=True)
        raise

    print("[5/5] Validated outputs written:")
    for path in output_paths:
        print(path.resolve())
    print()
    print(f"Total BAM records examined: {results['total_bam_records']:,}")
    print(f"Included calibrated P-site reads: {results['included_count']:,}")
    print("COUNT CONSERVATION: PASS")
    technical_warnings = results["technical_warnings"]
    if technical_warnings:
        for warning in technical_warnings:
            print(
                "TECHNICAL WARNING: full-library length "
                f"{warning['read_length']} is PHASE0_NOT_DOMINANT; "
                "valid outputs retained and frozen offsets unchanged"
            )
        print("POST-APPLICATION PHASE CHECK: PASS WITH TECHNICAL WARNING")
    else:
        print("POST-APPLICATION PHASE CHECK: PASS")
    print("FROZEN OFFSETS MODIFIED: NO")
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
