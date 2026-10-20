#!/usr/bin/env python3
"""
Stage 0C

Module 03

Read-Length-Specific P-Site Calibration Candidate Audit

Version:
0.1.0

This module exposes provisional start-codon candidate evidence and independent
CDS-interior validation. It does not freeze final P-site offsets.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

from lib.bam import open_bam
from lib.io import read_features


MIN_READ_LENGTH = 25
MAX_READ_LENGTH = 35
RAW_PROFILE_MIN = -60
RAW_PROFILE_MAX = 120
CANDIDATE_OFFSET_MIN = 6
CANDIDATE_OFFSET_MAX = 18
VALIDATION_FLANK = 15

STATUS_SELECTED = "SELECTED"
STATUS_AMBIGUOUS = "AMBIGUOUS"
STATUS_INSUFFICIENT = "INSUFFICIENT_EVIDENCE"
STATUS_NOT_EVALUATED = "NOT_EVALUATED"


# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Generate a forward-only, read-length-specific P-site candidate "
            "evidence audit without freezing final offsets."
        )
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Mapped BAM file",
    )

    parser.add_argument(
        "--features",
        required=True,
        help="Approved human_features.tsv",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output directory",
    )

    return parser.parse_args()


# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------

def derive_sample_name(bam_path):
    """Derive a stable sample identifier from the BAM filename."""

    sample_name = Path(bam_path).stem

    for suffix in ("_nonrRNA_hisat2", "_hisat2", "_mapped"):
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]

    if not sample_name:
        raise ValueError(f"Could not derive sample name from BAM path: {bam_path}")

    return sample_name


def get_read_length(read):
    """Return query-sequence length, falling back to CIGAR inference."""

    if read.query_sequence is not None:
        return len(read.query_sequence)

    return read.infer_query_length()


def get_raw_five_prime_position(read):
    """Return the frozen raw 5' coordinate used for eligibility auditing."""

    if read.is_reverse:
        if read.reference_end is None:
            return None
        return read.reference_end - 1

    return read.reference_start


def classify_basic_eligible_read(read, features, exclusion_counts=None):
    """
    Apply mutually exclusive basic eligibility filters.

    Returns
    -------
    tuple | None
        (read_length, transcript, cds_start, cds_end_exclusive,
         read_5prime, is_reverse) for a basic-eligible record.
    """

    def exclude(category):
        if exclusion_counts is not None:
            exclusion_counts[category] += 1
        return None

    if read.is_unmapped:
        return exclude("unmapped_excluded")

    if read.is_secondary:
        return exclude("secondary_excluded")

    if read.is_supplementary:
        return exclude("supplementary_excluded")

    try:
        nh_tag = read.get_tag("NH")
    except KeyError:
        return exclude("nh_missing_excluded")

    if nh_tag != 1:
        return exclude("multimapped_excluded")

    transcript = read.reference_name
    if transcript is None or transcript not in features:
        return exclude("missing_transcript_features")

    read_length = get_read_length(read)
    if read_length is None or read_length <= 0:
        return exclude("invalid_read_length")

    if read_length < MIN_READ_LENGTH or read_length > MAX_READ_LENGTH:
        return exclude("outside_read_length_range")

    read_5prime = get_raw_five_prime_position(read)
    if read_5prime is None or read_5prime < 0:
        return exclude("invalid_five_prime_position")

    cds_start = features[transcript]["utr5"]
    cds_end_exclusive = cds_start + features[transcript]["cds"]

    return (
        read_length,
        transcript,
        cds_start,
        cds_end_exclusive,
        read_5prime,
        read.is_reverse,
    )


def select_candidate(candidate_counts):
    """Apply zero-evidence, unique-maximum, and exact-tie rules."""

    offsets = list(range(CANDIDATE_OFFSET_MIN, CANDIDATE_OFFSET_MAX + 1))
    candidate_window_total_count = sum(candidate_counts[offset] for offset in offsets)
    top_score = max(candidate_counts[offset] for offset in offsets)
    top_offsets = [
        offset for offset in offsets if candidate_counts[offset] == top_score
    ]

    if candidate_window_total_count == 0:
        return {
            "candidate_window_total_count": 0,
            "top1_offset": "NA",
            "top1_score": 0,
            "top2_offset": "NA",
            "top2_score": 0,
            "top1_minus_top2": 0,
            "top1_to_top2_ratio": "NA",
            "selected_offset_distance": "NA",
            "calibration_status": STATUS_INSUFFICIENT,
        }

    if len(top_offsets) > 1:
        return {
            "candidate_window_total_count": candidate_window_total_count,
            "top1_offset": "NA",
            "top1_score": top_score,
            "top2_offset": "NA",
            "top2_score": top_score,
            "top1_minus_top2": 0,
            "top1_to_top2_ratio": 1.0,
            "selected_offset_distance": "NA",
            "calibration_status": STATUS_AMBIGUOUS,
        }

    top1_offset = top_offsets[0]
    remaining_offsets = [offset for offset in offsets if offset != top1_offset]
    top2_score = max(candidate_counts[offset] for offset in remaining_offsets)
    top2_offsets = [
        offset
        for offset in remaining_offsets
        if candidate_counts[offset] == top2_score
    ]
    top2_offset = top2_offsets[0] if len(top2_offsets) == 1 else "NA"
    top1_to_top2_ratio = (
        round(top_score / top2_score, 6) if top2_score > 0 else "NA"
    )

    return {
        "candidate_window_total_count": candidate_window_total_count,
        "top1_offset": top1_offset,
        "top1_score": top_score,
        "top2_offset": top2_offset,
        "top2_score": top2_score,
        "top1_minus_top2": top_score - top2_score,
        "top1_to_top2_ratio": top1_to_top2_ratio,
        "selected_offset_distance": top1_offset,
        "calibration_status": STATUS_SELECTED,
    }


def write_tsv(path, fieldnames, rows):
    """Write deterministic TSV output."""

    with path.open("w", newline="", encoding="utf-8") as output_handle:
        writer = csv.DictWriter(
            output_handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------
# Processing
# ---------------------------------------------------------

def build_raw_profile(bam_path, features):
    """Construct the direct-BAM, forward-only raw start profile."""

    exclusion_counts = Counter()
    total_bam_records = 0
    raw_profile_counts = Counter()
    length_stats = {
        length: {
            "read_length": length,
            "basic_eligible_before_orientation": 0,
            "forward_included": 0,
            "reverse_orientation_excluded": 0,
            "outside_raw_profile_range": 0,
            "raw_profile_count": 0,
        }
        for length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1)
    }

    bam = open_bam(bam_path)

    try:
        for read in bam.fetch(until_eof=True):
            total_bam_records += 1
            eligible = classify_basic_eligible_read(read, features, exclusion_counts)

            if eligible is None:
                continue

            read_length, _, cds_start, _, read_5prime, is_reverse = eligible
            length_stats[read_length]["basic_eligible_before_orientation"] += 1

            if is_reverse:
                exclusion_counts["reverse_orientation_excluded"] += 1
                length_stats[read_length]["reverse_orientation_excluded"] += 1
                continue

            length_stats[read_length]["forward_included"] += 1
            relative_5prime = read_5prime - cds_start

            if not RAW_PROFILE_MIN <= relative_5prime <= RAW_PROFILE_MAX:
                exclusion_counts["outside_raw_profile_range"] += 1
                length_stats[read_length]["outside_raw_profile_range"] += 1
                continue

            # Deliberately aggregate raw eligible read counts across approved
            # human transcripts; do not normalize or weight transcripts.
            raw_profile_counts[(read_length, relative_5prime)] += 1
            exclusion_counts["included_in_raw_calibration_profile"] += 1
            length_stats[read_length]["raw_profile_count"] += 1

    finally:
        bam.close()

    return total_bam_records, exclusion_counts, raw_profile_counts, length_stats


def build_candidate_outputs(sample_name, raw_profile_counts):
    """Build candidate and sample/read-length selection evidence rows."""

    candidate_rows = []
    summary_rows = []
    selections = {}

    for read_length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1):
        candidate_counts = {
            offset: raw_profile_counts[(read_length, -offset)]
            for offset in range(CANDIDATE_OFFSET_MIN, CANDIDATE_OFFSET_MAX + 1)
        }
        selection = select_candidate(candidate_counts)
        selections[read_length] = selection

        for offset in range(CANDIDATE_OFFSET_MIN, CANDIDATE_OFFSET_MAX + 1):
            candidate_count = candidate_counts[offset]
            candidate_rank = 1 + sum(
                count > candidate_count for count in candidate_counts.values()
            )
            candidate_rows.append({
                "sample": sample_name,
                "read_length": read_length,
                "candidate_offset_distance": offset,
                "relative_5prime_position": -offset,
                "candidate_count": candidate_count,
                "candidate_rank": candidate_rank,
                "selected_offset_distance": selection["selected_offset_distance"],
                "calibration_status": selection["calibration_status"],
            })

        summary_rows.append({
            "sample": sample_name,
            "read_length": read_length,
            **selection,
        })

    return candidate_rows, summary_rows, selections


def build_validation_rows(bam_path, features, sample_name, selections):
    """Build independent CDS-interior validation evidence for selected lengths."""

    selected_offsets = {
        length: selection["selected_offset_distance"]
        for length, selection in selections.items()
        if selection["calibration_status"] == STATUS_SELECTED
    }
    phase_counts = {length: Counter() for length in selected_offsets}
    validation_stats = {
        length: {
            "validation_basic_eligible": 0,
            "no_valid_validation_interior": 0,
            "outside_validation_region": 0,
            "n_validation_reads": 0,
        }
        for length in selected_offsets
    }

    if selected_offsets:
        bam = open_bam(bam_path)

        try:
            for read in bam.fetch(until_eof=True):
                eligible = classify_basic_eligible_read(read, features)

                if eligible is None:
                    continue

                (
                    read_length,
                    _,
                    cds_start,
                    cds_end_exclusive,
                    read_5prime,
                    is_reverse,
                ) = eligible

                if is_reverse or read_length not in selected_offsets:
                    continue

                stats = validation_stats[read_length]
                stats["validation_basic_eligible"] += 1
                validation_start = cds_start + VALIDATION_FLANK
                validation_end_exclusive = cds_end_exclusive - VALIDATION_FLANK

                if validation_start >= validation_end_exclusive:
                    stats["no_valid_validation_interior"] += 1
                    continue

                if not validation_start <= read_5prime < validation_end_exclusive:
                    stats["outside_validation_region"] += 1
                    continue

                relative_psite = (
                    read_5prime - cds_start + selected_offsets[read_length]
                )
                calibrated_phase = relative_psite % 3
                phase_counts[read_length][calibrated_phase] += 1
                stats["n_validation_reads"] += 1

        finally:
            bam.close()

    validation_rows = []

    for read_length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1):
        selection = selections[read_length]

        if selection["calibration_status"] != STATUS_SELECTED:
            validation_rows.append({
                "sample": sample_name,
                "read_length": read_length,
                "selected_offset_distance": "NA",
                "n_validation_reads": "NA",
                "phase0_count": "NA",
                "phase1_count": "NA",
                "phase2_count": "NA",
                "phase0_pct": "NA",
                "phase1_pct": "NA",
                "phase2_pct": "NA",
                "dominant_calibrated_phase": "NA",
                "dominant_calibrated_phase_pct": "NA",
                "validation_status": STATUS_NOT_EVALUATED,
            })
            continue

        counts = phase_counts[read_length]
        phase0_count = counts[0]
        phase1_count = counts[1]
        phase2_count = counts[2]
        n_validation_reads = phase0_count + phase1_count + phase2_count

        if n_validation_reads == 0:
            phase0_pct = "NA"
            phase1_pct = "NA"
            phase2_pct = "NA"
            dominant_phase = "NA"
            dominant_phase_pct = "NA"
            validation_status = "INSUFFICIENT_VALIDATION_EVIDENCE"
        else:
            phase0_pct = round(100.0 * phase0_count / n_validation_reads, 2)
            phase1_pct = round(100.0 * phase1_count / n_validation_reads, 2)
            phase2_pct = round(100.0 * phase2_count / n_validation_reads, 2)
            counts_by_phase = {
                0: phase0_count,
                1: phase1_count,
                2: phase2_count,
            }
            max_count = max(counts_by_phase.values())
            dominant_phases = [
                phase for phase, count in counts_by_phase.items() if count == max_count
            ]
            dominant_phase = (
                dominant_phases[0] if len(dominant_phases) == 1 else "TIE"
            )
            dominant_phase_pct = round(100.0 * max_count / n_validation_reads, 2)
            validation_status = (
                "PHASE0_DOMINANT"
                if phase0_count > phase1_count and phase0_count > phase2_count
                else "PHASE0_NOT_DOMINANT"
            )

        validation_rows.append({
            "sample": sample_name,
            "read_length": read_length,
            "selected_offset_distance": selection["selected_offset_distance"],
            "n_validation_reads": n_validation_reads,
            "phase0_count": phase0_count,
            "phase1_count": phase1_count,
            "phase2_count": phase2_count,
            "phase0_pct": phase0_pct,
            "phase1_pct": phase1_pct,
            "phase2_pct": phase2_pct,
            "dominant_calibrated_phase": dominant_phase,
            "dominant_calibrated_phase_pct": dominant_phase_pct,
            "validation_status": validation_status,
        })

    return validation_rows, validation_stats


# ---------------------------------------------------------
# Validation
# ---------------------------------------------------------

def validate_raw_accounting(
    total_bam_records,
    exclusion_counts,
    raw_profile_rows,
    length_stats,
):
    """Validate record-level, profile, and per-length conservation."""

    exclusion_fields = [
        "unmapped_excluded",
        "secondary_excluded",
        "supplementary_excluded",
        "nh_missing_excluded",
        "multimapped_excluded",
        "missing_transcript_features",
        "invalid_read_length",
        "outside_read_length_range",
        "invalid_five_prime_position",
        "reverse_orientation_excluded",
        "outside_raw_profile_range",
    ]
    total_excluded = sum(exclusion_counts[field] for field in exclusion_fields)
    included = exclusion_counts["included_in_raw_calibration_profile"]
    accounted_for = total_excluded + included
    accounting_difference = total_bam_records - accounted_for

    if accounting_difference != 0:
        raise RuntimeError(
            "Raw-profile count conservation failed: "
            f"accounted_for ({accounted_for:,}) != "
            f"total_bam_records ({total_bam_records:,}); "
            f"difference = {accounting_difference:,}"
        )

    if len(raw_profile_rows) != 11 * 181:
        raise RuntimeError(
            f"Raw profile must contain 1,991 rows; observed {len(raw_profile_rows):,}"
        )

    raw_profile_sum = sum(row["count"] for row in raw_profile_rows)
    if raw_profile_sum != included:
        raise RuntimeError(
            f"Raw-profile count sum ({raw_profile_sum:,}) != included count "
            f"({included:,})"
        )

    for read_length, stats in length_stats.items():
        if stats["forward_included"] != (
            stats["outside_raw_profile_range"] + stats["raw_profile_count"]
        ):
            raise RuntimeError(
                f"Length {read_length} forward accounting failed: "
                f"forward_included ({stats['forward_included']:,}) != "
                "outside_raw_profile_range + raw_profile_count "
                f"({stats['outside_raw_profile_range']:,} + "
                f"{stats['raw_profile_count']:,})"
            )

    if sum(stats["raw_profile_count"] for stats in length_stats.values()) != included:
        raise RuntimeError(
            "Per-length raw_profile_count sum does not match "
            "included_in_raw_calibration_profile"
        )

    if sum(
        stats["reverse_orientation_excluded"] for stats in length_stats.values()
    ) != exclusion_counts["reverse_orientation_excluded"]:
        raise RuntimeError(
            "Per-length reverse_orientation_excluded sum does not match "
            "the sample-level count"
        )

    if sum(
        stats["outside_raw_profile_range"] for stats in length_stats.values()
    ) != exclusion_counts["outside_raw_profile_range"]:
        raise RuntimeError(
            "Per-length outside_raw_profile_range sum does not match "
            "the sample-level count"
        )

    for read_length, stats in length_stats.items():
        if stats["basic_eligible_before_orientation"] != (
            stats["forward_included"] + stats["reverse_orientation_excluded"]
        ):
            raise RuntimeError(
                f"Length {read_length} orientation accounting failed: "
                "basic_eligible_before_orientation != "
                "forward_included + reverse_orientation_excluded"
            )

        profile_length_sum = sum(
            row["count"]
            for row in raw_profile_rows
            if row["read_length"] == read_length
        )
        if profile_length_sum != stats["raw_profile_count"]:
            raise RuntimeError(
                f"Length {read_length} profile count sum ({profile_length_sum:,}) "
                f"!= raw_profile_count ({stats['raw_profile_count']:,})"
            )

    return total_excluded, accounted_for, accounting_difference


def validate_candidate_outputs(candidate_rows, summary_rows, raw_profile_rows):
    """Validate candidate grid, scores, ranks, and provisional selections."""

    if len(candidate_rows) != 11 * 13:
        raise RuntimeError(
            f"Candidate table must contain 143 rows; observed {len(candidate_rows):,}"
        )

    if len(summary_rows) != 11:
        raise RuntimeError(
            f"Selection summary must contain 11 rows; observed {len(summary_rows):,}"
        )

    summary_by_length = {row["read_length"]: row for row in summary_rows}

    for read_length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1):
        rows = [
            row for row in candidate_rows if row["read_length"] == read_length
        ]
        expected_offsets = list(
            range(CANDIDATE_OFFSET_MIN, CANDIDATE_OFFSET_MAX + 1)
        )
        observed_offsets = [row["candidate_offset_distance"] for row in rows]

        if len(rows) != 13 or sorted(observed_offsets) != expected_offsets:
            raise RuntimeError(
                f"Length {read_length} does not contain candidate offsets 6-18 "
                "exactly once"
            )

        for row in rows:
            if row["relative_5prime_position"] != -row["candidate_offset_distance"]:
                raise RuntimeError(
                    f"Length {read_length}, offset "
                    f"{row['candidate_offset_distance']} has inconsistent raw position"
                )

        summary = summary_by_length[read_length]
        score_sum = sum(row["candidate_count"] for row in rows)
        candidate_window_total_count = summary["candidate_window_total_count"]
        if score_sum != candidate_window_total_count:
            raise RuntimeError(
                f"Length {read_length} candidate score sum ({score_sum:,}) != "
                "candidate_window_total_count "
                f"({candidate_window_total_count:,})"
            )

        score_positions = {
            -offset
            for offset in range(CANDIDATE_OFFSET_MIN, CANDIDATE_OFFSET_MAX + 1)
        }
        length_profile_rows = [
            row
            for row in raw_profile_rows
            if row["read_length"] == read_length
        ]
        raw_score_region_count = sum(
            row["count"]
            for row in length_profile_rows
            if row["relative_position"] in score_positions
        )
        counts_outside_score_region = sum(
            row["count"]
            for row in length_profile_rows
            if row["relative_position"] not in score_positions
        )
        raw_profile_count = sum(row["count"] for row in length_profile_rows)

        if candidate_window_total_count != raw_score_region_count:
            raise RuntimeError(
                f"Length {read_length} candidate_window_total_count "
                f"({candidate_window_total_count:,}) != raw-profile score-region count "
                f"({raw_score_region_count:,})"
            )
        if (
            candidate_window_total_count + counts_outside_score_region
            != raw_profile_count
        ):
            raise RuntimeError(
                f"Length {read_length} score-region complement check failed: "
                f"candidate_window_total_count ({candidate_window_total_count:,}) + "
                f"counts_outside_score_region "
                f"({counts_outside_score_region:,}) != raw_profile_count "
                f"({raw_profile_count:,})"
            )

        candidate_counts = {
            row["candidate_offset_distance"]: row["candidate_count"]
            for row in rows
        }
        expected_selection = select_candidate(candidate_counts)
        for field, expected_value in expected_selection.items():
            if summary[field] != expected_value:
                raise RuntimeError(
                    f"Length {read_length} has inconsistent {field}: "
                    f"{summary[field]} != {expected_value}"
                )

        max_count = max(row["candidate_count"] for row in rows)
        max_rows = [row for row in rows if row["candidate_count"] == max_count]
        status = summary["calibration_status"]
        selected = summary["selected_offset_distance"]

        if score_sum == 0:
            if status != STATUS_INSUFFICIENT or selected != "NA":
                raise RuntimeError(
                    f"Length {read_length} violates the zero-evidence rule"
                )
        elif len(max_rows) > 1:
            if status != STATUS_AMBIGUOUS or selected != "NA":
                raise RuntimeError(
                    f"Length {read_length} violates the exact-tie rule"
                )
        else:
            unique_max_offset = max_rows[0]["candidate_offset_distance"]
            if status != STATUS_SELECTED or selected != unique_max_offset:
                raise RuntimeError(
                    f"Length {read_length} selected offset is not the unique maximum"
                )

        for row in rows:
            expected_rank = 1 + sum(
                other["candidate_count"] > row["candidate_count"] for other in rows
            )
            if row["candidate_rank"] != expected_rank:
                raise RuntimeError(
                    f"Length {read_length}, offset "
                    f"{row['candidate_offset_distance']} has incorrect rank"
                )
            if (
                row["selected_offset_distance"] != selected
                or row["calibration_status"] != status
            ):
                raise RuntimeError(
                    f"Length {read_length} candidate rows disagree with summary"
                )


def validate_validation_outputs(validation_rows, validation_stats, selections):
    """Validate CDS-interior phase accounting and status assignment."""

    if len(validation_rows) != 11:
        raise RuntimeError(
            f"Validation table must contain 11 rows; observed {len(validation_rows):,}"
        )

    for row in validation_rows:
        read_length = row["read_length"]
        selection = selections[read_length]

        if selection["calibration_status"] != STATUS_SELECTED:
            na_fields = [
                "selected_offset_distance",
                "n_validation_reads",
                "phase0_count",
                "phase1_count",
                "phase2_count",
                "phase0_pct",
                "phase1_pct",
                "phase2_pct",
                "dominant_calibrated_phase",
                "dominant_calibrated_phase_pct",
            ]
            if any(row[field] != "NA" for field in na_fields):
                raise RuntimeError(
                    f"Length {read_length} has fabricated validation values"
                )
            if row["validation_status"] != STATUS_NOT_EVALUATED:
                raise RuntimeError(
                    f"Length {read_length} must have NOT_EVALUATED validation status"
                )
            continue

        phase_sum = row["phase0_count"] + row["phase1_count"] + row["phase2_count"]
        if (
            row["selected_offset_distance"]
            != selection["selected_offset_distance"]
        ):
            raise RuntimeError(
                f"Length {read_length} validation offset disagrees with selection"
            )
        if phase_sum != row["n_validation_reads"]:
            raise RuntimeError(
                f"Length {read_length} validation phase counts ({phase_sum:,}) != "
                f"n_validation_reads ({row['n_validation_reads']:,})"
            )

        stats = validation_stats[read_length]
        accounted = (
            stats["no_valid_validation_interior"]
            + stats["outside_validation_region"]
            + stats["n_validation_reads"]
        )
        if accounted != stats["validation_basic_eligible"]:
            raise RuntimeError(
                f"Length {read_length} validation accounting failed: "
                f"{accounted:,} != {stats['validation_basic_eligible']:,}"
            )

        if row["n_validation_reads"] == 0:
            if row["validation_status"] != "INSUFFICIENT_VALIDATION_EVIDENCE":
                raise RuntimeError(
                    f"Length {read_length} has incorrect zero-validation status"
                )
            for field in [
                "phase0_pct",
                "phase1_pct",
                "phase2_pct",
                "dominant_calibrated_phase",
                "dominant_calibrated_phase_pct",
            ]:
                if row[field] != "NA":
                    raise RuntimeError(
                        f"Length {read_length} must report NA for {field}"
                    )
            continue

        pct_sum = row["phase0_pct"] + row["phase1_pct"] + row["phase2_pct"]
        if abs(pct_sum - 100.0) > 0.02 + 1e-9:
            raise RuntimeError(
                f"Length {read_length} rounded phase percentages sum to "
                f"{pct_sum:.2f}, not approximately 100"
            )

        expected_status = (
            "PHASE0_DOMINANT"
            if row["phase0_count"] > row["phase1_count"]
            and row["phase0_count"] > row["phase2_count"]
            else "PHASE0_NOT_DOMINANT"
        )
        if row["validation_status"] != expected_status:
            raise RuntimeError(
                f"Length {read_length} has inconsistent validation status"
            )

        phase_counts = {
            0: row["phase0_count"],
            1: row["phase1_count"],
            2: row["phase2_count"],
        }
        max_count = max(phase_counts.values())
        dominant_phases = [
            phase for phase, count in phase_counts.items() if count == max_count
        ]
        expected_dominant = (
            dominant_phases[0] if len(dominant_phases) == 1 else "TIE"
        )
        if row["dominant_calibrated_phase"] != expected_dominant:
            raise RuntimeError(
                f"Length {read_length} has inconsistent dominant phase"
            )

        expected_dominant_pct = round(
            100.0 * max_count / row["n_validation_reads"],
            2,
        )
        if row["dominant_calibrated_phase_pct"] != expected_dominant_pct:
            raise RuntimeError(
                f"Length {read_length} has inconsistent dominant phase percentage"
            )


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    args = parse_args()
    sample_name = derive_sample_name(args.bam)

    print("=" * 60)
    print("Stage0C Module 03")
    print("Read-Length-Specific P-Site Candidate Audit v0.1.0")
    print("=" * 60)
    print()
    print("TMD-BLIND: YES")
    print("BENCHMARK-BLIND: YES")
    print("BIOLOGICAL-RESULT-BLIND: YES")
    print("FINAL OFFSET FREEZE: NO")
    print()

    print("[1/7] Loading approved human transcript features...")
    features = read_features(args.features)
    if not features:
        raise RuntimeError("Approved human feature table contains no transcripts")
    print(f"Loaded {len(features):,} transcripts")

    print()
    print("[2/7] Constructing direct-BAM forward-only raw profile...")
    (
        total_bam_records,
        exclusion_counts,
        raw_profile_counts,
        length_stats,
    ) = build_raw_profile(args.bam, features)

    raw_profile_rows = [
        {
            "sample": sample_name,
            "read_length": read_length,
            "relative_position": relative_position,
            "count": raw_profile_counts[(read_length, relative_position)],
        }
        for read_length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1)
        for relative_position in range(RAW_PROFILE_MIN, RAW_PROFILE_MAX + 1)
    ]

    print()
    print("[3/7] Validating raw-profile accounting...")
    total_excluded, accounted_for, accounting_difference = validate_raw_accounting(
        total_bam_records,
        exclusion_counts,
        raw_profile_rows,
        length_stats,
    )

    print()
    print("[4/7] Evaluating provisional candidate offsets...")
    candidate_rows, summary_rows, selections = build_candidate_outputs(
        sample_name,
        raw_profile_counts,
    )

    print()
    print("[5/7] Validating candidate outputs and score-region complements...")
    validate_candidate_outputs(candidate_rows, summary_rows, raw_profile_rows)

    print()
    print("[6/7] Computing and validating CDS-interior validation evidence...")
    validation_rows, validation_stats = build_validation_rows(
        args.bam,
        features,
        sample_name,
        selections,
    )
    validate_validation_outputs(validation_rows, validation_stats, selections)

    length_stats_rows = []
    for read_length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1):
        row = dict(length_stats[read_length])
        if read_length in validation_stats:
            row.update(validation_stats[read_length])
        else:
            row.update({
                "validation_basic_eligible": "NA",
                "no_valid_validation_interior": "NA",
                "outside_validation_region": "NA",
                "n_validation_reads": "NA",
            })
        length_stats_rows.append(row)

    run_stats_rows = [
        {"statistic": "sample", "value": sample_name},
        {
            "statistic": "total_bam_records_examined",
            "value": total_bam_records,
        },
        {
            "statistic": "basic_eligible_before_orientation",
            "value": sum(
                stats["basic_eligible_before_orientation"]
                for stats in length_stats.values()
            ),
        },
        {
            "statistic": "forward_included",
            "value": sum(
                stats["forward_included"] for stats in length_stats.values()
            ),
        },
    ]
    for statistic in [
        "unmapped_excluded",
        "secondary_excluded",
        "supplementary_excluded",
        "nh_missing_excluded",
        "multimapped_excluded",
        "missing_transcript_features",
        "invalid_read_length",
        "outside_read_length_range",
        "invalid_five_prime_position",
        "reverse_orientation_excluded",
        "outside_raw_profile_range",
        "included_in_raw_calibration_profile",
    ]:
        run_stats_rows.append({
            "statistic": statistic,
            "value": exclusion_counts[statistic],
        })
    run_stats_rows.extend([
        {"statistic": "total_excluded", "value": total_excluded},
        {"statistic": "accounted_for", "value": accounted_for},
        {"statistic": "accounting_difference", "value": accounting_difference},
        {"statistic": "raw_profile_min", "value": RAW_PROFILE_MIN},
        {"statistic": "raw_profile_max", "value": RAW_PROFILE_MAX},
        {
            "statistic": "candidate_offset_min",
            "value": CANDIDATE_OFFSET_MIN,
        },
        {
            "statistic": "candidate_offset_max",
            "value": CANDIDATE_OFFSET_MAX,
        },
        {"statistic": "raw_profile_rows", "value": len(raw_profile_rows)},
        {"statistic": "candidate_rows", "value": len(candidate_rows)},
        {"statistic": "selection_summary_rows", "value": len(summary_rows)},
        {"statistic": "validation_rows", "value": len(validation_rows)},
        {
            "statistic": "selected_lengths",
            "value": sum(
                selection["calibration_status"] == STATUS_SELECTED
                for selection in selections.values()
            ),
        },
        {
            "statistic": "ambiguous_lengths",
            "value": sum(
                selection["calibration_status"] == STATUS_AMBIGUOUS
                for selection in selections.values()
            ),
        },
        {
            "statistic": "insufficient_evidence_lengths",
            "value": sum(
                selection["calibration_status"] == STATUS_INSUFFICIENT
                for selection in selections.values()
            ),
        },
        {"statistic": "final_offset_freeze", "value": "NO"},
    ])

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_profile_path = output_dir / f"{sample_name}_psite_raw_profile.tsv"
    candidates_path = output_dir / f"{sample_name}_psite_candidates.tsv"
    summary_path = output_dir / f"{sample_name}_psite_selection_summary.tsv"
    validation_path = output_dir / f"{sample_name}_psite_validation.tsv"
    run_stats_path = output_dir / f"{sample_name}_psite_run_stats.tsv"
    length_stats_path = output_dir / f"{sample_name}_psite_length_stats.tsv"

    print()
    print("[7/7] Writing validated audit outputs...")
    write_tsv(
        raw_profile_path,
        ["sample", "read_length", "relative_position", "count"],
        raw_profile_rows,
    )
    write_tsv(
        candidates_path,
        [
            "sample",
            "read_length",
            "candidate_offset_distance",
            "relative_5prime_position",
            "candidate_count",
            "candidate_rank",
            "selected_offset_distance",
            "calibration_status",
        ],
        candidate_rows,
    )
    write_tsv(
        summary_path,
        [
            "sample",
            "read_length",
            "candidate_window_total_count",
            "top1_offset",
            "top1_score",
            "top2_offset",
            "top2_score",
            "top1_minus_top2",
            "top1_to_top2_ratio",
            "selected_offset_distance",
            "calibration_status",
        ],
        summary_rows,
    )
    write_tsv(
        validation_path,
        [
            "sample",
            "read_length",
            "selected_offset_distance",
            "n_validation_reads",
            "phase0_count",
            "phase1_count",
            "phase2_count",
            "phase0_pct",
            "phase1_pct",
            "phase2_pct",
            "dominant_calibrated_phase",
            "dominant_calibrated_phase_pct",
            "validation_status",
        ],
        validation_rows,
    )
    write_tsv(
        run_stats_path,
        ["statistic", "value"],
        run_stats_rows,
    )
    write_tsv(
        length_stats_path,
        [
            "read_length",
            "basic_eligible_before_orientation",
            "forward_included",
            "reverse_orientation_excluded",
            "outside_raw_profile_range",
            "raw_profile_count",
            "validation_basic_eligible",
            "no_valid_validation_interior",
            "outside_validation_region",
            "n_validation_reads",
        ],
        length_stats_rows,
    )

    print()
    print("Run statistics")
    print(f"Sample: {sample_name}")
    print(f"Total BAM records examined: {total_bam_records:,}")
    print(f"Reverse orientation excluded: "
          f"{exclusion_counts['reverse_orientation_excluded']:,}")
    print(f"Outside raw-profile range: "
          f"{exclusion_counts['outside_raw_profile_range']:,}")
    print(f"Included in raw calibration profile: "
          f"{exclusion_counts['included_in_raw_calibration_profile']:,}")
    print(f"Accounting difference: {accounting_difference:,}")

    print()
    print("Provisional candidate evidence by read length")
    for row in summary_rows:
        print(
            f"Length {row['read_length']:2d}: "
            f"window_n={row['candidate_window_total_count']:>7,}, "
            f"top1={row['top1_offset']} ({row['top1_score']}), "
            f"top2={row['top2_offset']} ({row['top2_score']}), "
            f"selected={row['selected_offset_distance']}, "
            f"status={row['calibration_status']}"
        )

    print()
    for label, path in [
        ("Raw profile", raw_profile_path),
        ("Candidates", candidates_path),
        ("Selection summary", summary_path),
        ("Validation", validation_path),
        ("Run statistics", run_stats_path),
        ("Length statistics", length_stats_path),
    ]:
        print(f"{label}: {path.resolve()}")

    print()
    print("CANDIDATE AUDIT MODE: YES")
    print("FINAL OFFSET FREEZE: NO")
    print("PASS")
    print()
    print("P-site calibration candidate audit completed.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
