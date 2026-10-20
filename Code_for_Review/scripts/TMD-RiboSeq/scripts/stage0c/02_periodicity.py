#!/usr/bin/env python3
"""
Stage 0C

Module 02

Read-Length-Specific 3-nt Periodicity QC

Version:
0.1.0
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from lib.io import read_features
from lib.bam import open_bam


MIN_READ_LENGTH = 25
MAX_READ_LENGTH = 35
PERIODICITY_FLANK = 15  # Exclude first and last 5 codons = 15 nt each


# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description="Generate Read-Length-Specific 3-nt Periodicity QC."
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Unique mapped BAM file"
    )

    parser.add_argument(
        "--features",
        required=True,
        help="human_features.tsv"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output directory"
    )

    return parser.parse_args()


# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------

def get_read_length(read):
    """Return query-sequence length, falling back to CIGAR inference."""

    if read.query_sequence is not None:
        return len(read.query_sequence)

    return read.infer_query_length()


def get_five_prime_position(read):
    """Return the zero-based mapped 5' end in transcript coordinates."""

    if read.is_reverse:
        if read.reference_end is None:
            return None
        return read.reference_end - 1

    return read.reference_start


def get_cds_region(features, transcript):
    """
    Compute CDS start and end (exclusive) for a transcript.

    Returns (cds_start, cds_end_exclusive) or (None, None) if transcript not in features.
    """
    if transcript not in features:
        return None, None

    utr5 = features[transcript]["utr5"]
    cds_len = features[transcript]["cds"]
    cds_end_exclusive = utr5 + cds_len
    cds_start = utr5

    return cds_start, cds_end_exclusive


def get_periodicity_region(features, transcript):
    """
    Compute periodicity analysis region for a transcript.

    Returns (periodicity_start, periodicity_end_exclusive) or (None, None) if:
    - transcript not in features
    - transcript has no valid interior (periodicity_start >= periodicity_end_exclusive)

    The periodicity region excludes first and last 5 codons (15 nt each).
    """
    cds_start, cds_end_exclusive = get_cds_region(features, transcript)

    if cds_start is None:
        return None, None

    periodicity_start = cds_start + PERIODICITY_FLANK
    periodicity_end_exclusive = cds_end_exclusive - PERIODICITY_FLANK

    if periodicity_start >= periodicity_end_exclusive:
        return None, None

    return periodicity_start, periodicity_end_exclusive


def compute_phase(read_5prime, cds_start):
    """
    Compute phase relative to true CDS start.

    Phase = (relative_position) % 3
    where relative_position = read_5prime - cds_start

    Uses Python modulo behavior which handles negative numbers correctly.
    """
    relative_position = read_5prime - cds_start
    return relative_position % 3


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    args = parse_args()

    print("=" * 60)
    print("Stage0C Module 02")
    print("Read-Length-Specific 3-nt Periodicity QC v0.1.0")
    print("=" * 60)

    print()

    print("[1/3] Loading transcript features...")

    features = read_features(args.features)

    print(f"Loaded {len(features):,} transcripts")

    print()

    print("[2/3] Opening BAM...")

    bam = open_bam(args.bam)

    print("BAM opened successfully")

    print()

    print("[3/3] Processing reads and computing periodicity...")

    # Initialize counters for run-level accounting
    total_bam_records = 0
    unmapped_excluded = 0
    secondary_excluded = 0
    supplementary_excluded = 0
    nh_missing_excluded = 0
    multimapped_excluded = 0
    outside_read_length_range = 0
    invalid_read_length = 0
    invalid_five_prime_position = 0
    missing_transcript_features = 0
    outside_analysis_region = 0
    no_valid_periodicity_interior = 0
    included_in_periodicity = 0

    # Per-read-length technical accounting
    # Key: read_length -> dict of counts
    length_stats = {}
    for length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1):
        length_stats[length] = {
            "read_length": length,
            "outside_analysis_region": 0,
            "no_valid_periodicity_interior": 0,
            "included_in_periodicity": 0,
            "phase0_count": 0,
            "phase1_count": 0,
            "phase2_count": 0,
        }

    try:
        for read in bam.fetch(until_eof=True):

            total_bam_records += 1

            # Filter 1: Unmapped
            if read.is_unmapped:
                unmapped_excluded += 1
                continue

            # Filter 2: Secondary alignment
            if read.is_secondary:
                secondary_excluded += 1
                continue

            # Filter 3: Supplementary alignment
            if read.is_supplementary:
                supplementary_excluded += 1
                continue

            # Filter 4: Check NH tag for unique mapping
            try:
                nh_tag = read.get_tag("NH")
            except KeyError:
                nh_missing_excluded += 1
                continue

            if nh_tag != 1:
                multimapped_excluded += 1
                continue

            # Get transcript name
            transcript = read.reference_name
            if transcript is None:
                # This should be extremely rare, but handle it
                missing_transcript_features += 1
                continue

            # Filter 5: Check if transcript is in approved human features table
            if transcript not in features:
                missing_transcript_features += 1
                continue

            # Get read length
            length = get_read_length(read)
            if length is None or length <= 0:
                # Should not happen for properly aligned reads, but handle gracefully
                invalid_read_length += 1
                continue

            # Filter 6: Length eligibility (25-35 nt)
            if length < MIN_READ_LENGTH or length > MAX_READ_LENGTH:
                outside_read_length_range += 1
                continue

            # Get 5' position
            read_5prime = get_five_prime_position(read)
            if read_5prime is None or read_5prime < 0:
                invalid_five_prime_position += 1
                continue

            # Get periodicity region for this transcript
            periodicity_start, periodicity_end_exclusive = get_periodicity_region(
                features, transcript
            )

            if periodicity_start is None:
                # Feature-table membership was already checked above, so this
                # category means only that the CDS has no valid interior.
                no_valid_periodicity_interior += 1
                length_stats[length]["no_valid_periodicity_interior"] += 1
                continue

            # Filter 7: Check if read falls within periodicity region
            if not (periodicity_start <= read_5prime < periodicity_end_exclusive):
                outside_analysis_region += 1
                length_stats[length]["outside_analysis_region"] += 1
                continue

            # Compute phase relative to TRUE CDS start (not periodicity_start)
            cds_start, _ = get_cds_region(features, transcript)
            phase = compute_phase(read_5prime, cds_start)

            # Record phase count for this read length
            length_stats[length]["included_in_periodicity"] += 1
            length_stats[length][f"phase{phase}_count"] += 1
            included_in_periodicity += 1

    finally:
        bam.close()

    # Determine sample name from BAM filename
    bam_path = Path(args.bam)
    sample_name = bam_path.stem
    for suffix in ["_nonrRNA_hisat2", "_hisat2", "_mapped"]:
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-(len(suffix))]

    # Build output rows for main TSV
    output_rows = []
    for length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1):
        phase0_count = length_stats[length]["phase0_count"]
        phase1_count = length_stats[length]["phase1_count"]
        phase2_count = length_stats[length]["phase2_count"]

        n_reads = phase0_count + phase1_count + phase2_count

        if n_reads > 0:
            phase0_pct = round(100.0 * phase0_count / n_reads, 2)
            phase1_pct = round(100.0 * phase1_count / n_reads, 2)
            phase2_pct = round(100.0 * phase2_count / n_reads, 2)

            # Determine dominant phase
            phase_counts = [(0, phase0_count), (1, phase1_count), (2, phase2_count)]
            max_count = max(phase_counts, key=lambda x: x[1])[1]
            # Check for ties
            phases_with_max = [p for p, c in phase_counts if c == max_count]
            if len(phases_with_max) > 1:
                dominant_phase = "TIE"
            else:
                dominant_phase = phases_with_max[0]
            dominant_phase_pct = round(100.0 * max_count / n_reads, 2)
        else:
            phase0_pct = 0.0
            phase1_pct = 0.0
            phase2_pct = 0.0
            dominant_phase = "NA"
            dominant_phase_pct = "NA"

        output_rows.append({
            "sample": sample_name,
            "read_length": length,
            "n_reads": n_reads,
            "phase0_count": phase0_count,
            "phase1_count": phase1_count,
            "phase2_count": phase2_count,
            "phase0_pct": phase0_pct,
            "phase1_pct": phase1_pct,
            "phase2_pct": phase2_pct,
            "dominant_phase": dominant_phase,
            "dominant_phase_pct": dominant_phase_pct,
        })

    total_excluded = (
        unmapped_excluded +
        secondary_excluded +
        supplementary_excluded +
        nh_missing_excluded +
        multimapped_excluded +
        outside_read_length_range +
        invalid_read_length +
        invalid_five_prime_position +
        missing_transcript_features +
        outside_analysis_region +
        no_valid_periodicity_interior
    )

    accounted_for = total_excluded + included_in_periodicity
    accounting_difference = total_bam_records - accounted_for

    # Count conservation check 1: all mutually exclusive outcomes must account
    # for every BAM record examined.
    if accounting_difference != 0:
        raise RuntimeError(
            f"Count conservation failed: total_excluded ({total_excluded:,}) + "
            f"included_in_periodicity ({included_in_periodicity:,}) = "
            f"{accounted_for:,} != total_bam_records ({total_bam_records:,}); "
            f"accounting_difference = {accounting_difference:,}"
        )

    # Count conservation check 2: sum of n_reads across lengths must equal included_in_periodicity
    sum_n_reads = sum(row["n_reads"] for row in output_rows)
    if sum_n_reads != included_in_periodicity:
        raise RuntimeError(
            f"Count conservation failed: sum of n_reads ({sum_n_reads:,}) != "
            f"included_in_periodicity ({included_in_periodicity:,})"
        )

    # Count conservation check 3: phase0 + phase1 + phase2 must equal n_reads for each length
    for row in output_rows:
        phase_sum = row["phase0_count"] + row["phase1_count"] + row["phase2_count"]
        if phase_sum != row["n_reads"]:
            raise RuntimeError(
                f"Count conservation failed for length {row['read_length']}: "
                f"phase counts ({row['phase0_count']} + {row['phase1_count']} + "
                f"{row['phase2_count']} = {phase_sum}) != n_reads ({row['n_reads']})"
            )
        length = row["read_length"]
        if length_stats[length]["included_in_periodicity"] != row["n_reads"]:
            raise RuntimeError(
                f"Count conservation failed for length {length}: "
                f"length_stats included_in_periodicity "
                f"({length_stats[length]['included_in_periodicity']}) != "
                f"n_reads ({row['n_reads']})"
            )

    # Per-length accounting must reconcile with the corresponding global counts.
    length_outside_total = sum(
        stats["outside_analysis_region"] for stats in length_stats.values()
    )
    if length_outside_total != outside_analysis_region:
        raise RuntimeError(
            "Per-length outside_analysis_region counts do not match the global "
            f"count: {length_outside_total:,} != {outside_analysis_region:,}"
        )

    length_no_interior_total = sum(
        stats["no_valid_periodicity_interior"] for stats in length_stats.values()
    )
    if length_no_interior_total != no_valid_periodicity_interior:
        raise RuntimeError(
            "Per-length no_valid_periodicity_interior counts do not match the "
            f"global count: {length_no_interior_total:,} != "
            f"{no_valid_periodicity_interior:,}"
        )

    length_included_total = sum(
        stats["included_in_periodicity"] for stats in length_stats.values()
    )
    if length_included_total != included_in_periodicity:
        raise RuntimeError(
            "Per-length included_in_periodicity counts do not match the global "
            f"count: {length_included_total:,} != {included_in_periodicity:,}"
        )

    run_stats_rows = [
        {
            "statistic": "total_bam_records_examined",
            "value": total_bam_records,
        },
        {
            "statistic": "unmapped_excluded",
            "value": unmapped_excluded,
        },
        {
            "statistic": "secondary_excluded",
            "value": secondary_excluded,
        },
        {
            "statistic": "supplementary_excluded",
            "value": supplementary_excluded,
        },
        {
            "statistic": "nh_missing_excluded",
            "value": nh_missing_excluded,
        },
        {
            "statistic": "multimapped_excluded",
            "value": multimapped_excluded,
        },
        {
            "statistic": "outside_read_length_range",
            "value": outside_read_length_range,
        },
        {
            "statistic": "invalid_read_length",
            "value": invalid_read_length,
        },
        {
            "statistic": "invalid_five_prime_position",
            "value": invalid_five_prime_position,
        },
        {
            "statistic": "missing_transcript_features",
            "value": missing_transcript_features,
        },
        {
            "statistic": "outside_analysis_region",
            "value": outside_analysis_region,
        },
        {
            "statistic": "no_valid_periodicity_interior",
            "value": no_valid_periodicity_interior,
        },
        {
            "statistic": "included_in_periodicity",
            "value": included_in_periodicity,
        },
        {
            "statistic": "total_excluded",
            "value": total_excluded,
        },
        {
            "statistic": "accounted_for",
            "value": accounted_for,
        },
        {
            "statistic": "accounting_difference",
            "value": accounting_difference,
        },
    ]

    # Write outputs only after all validations pass.
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    periodicity_path = output_dir / f"{sample_name}_periodicity.tsv"
    run_stats_path = output_dir / f"{sample_name}_periodicity_run_stats.tsv"
    length_stats_path = output_dir / f"{sample_name}_periodicity_length_stats.tsv"

    with periodicity_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "sample",
                "read_length",
                "n_reads",
                "phase0_count",
                "phase1_count",
                "phase2_count",
                "phase0_pct",
                "phase1_pct",
                "phase2_pct",
                "dominant_phase",
                "dominant_phase_pct",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    with run_stats_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["statistic", "value"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(run_stats_rows)

    with length_stats_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "read_length",
                "outside_analysis_region",
                "no_valid_periodicity_interior",
                "included_in_periodicity",
                "phase0_count",
                "phase1_count",
                "phase2_count",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            length_stats[length]
            for length in range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1)
        )

    # Print summary
    print()
    print("Run statistics")
    print(f"Total BAM records examined: {total_bam_records:,}")
    print(f"Unmapped excluded: {unmapped_excluded:,}")
    print(f"Secondary alignment excluded: {secondary_excluded:,}")
    print(f"Supplementary alignment excluded: {supplementary_excluded:,}")
    print(f"NH tag missing excluded: {nh_missing_excluded:,}")
    print(f"Non-unique (NH != 1) excluded: {multimapped_excluded:,}")
    print(f"Outside read length range (25-35) excluded: {outside_read_length_range:,}")
    print(f"Invalid read length excluded: {invalid_read_length:,}")
    print(f"Invalid 5' position excluded: {invalid_five_prime_position:,}")
    print(f"Missing transcript features excluded: {missing_transcript_features:,}")
    print(f"Outside analysis region excluded: {outside_analysis_region:,}")
    print(f"No valid periodicity interior excluded: {no_valid_periodicity_interior:,}")
    print(f"Included in periodicity analysis: {included_in_periodicity:,}")

    print()
    print("Accounting verification")
    print(f"Total excluded: {total_excluded:,}")
    print(f"Accounted for (excluded + included): {accounted_for:,}")
    print(f"Original total: {total_bam_records:,}")
    print(f"Difference: {accounting_difference:,}")

    # Print periodicity summary
    print()
    print("Periodicity summary by read length")
    for row in output_rows:
        dominant_phase = row["dominant_phase"]
        dominant_pct = row["dominant_phase_pct"]
        if dominant_phase == "NA":
            dominant_str = "dominant=NA"
        elif dominant_phase == "TIE":
            dominant_str = f"dominant=TIE ({dominant_pct:5.2f}%)"
        else:
            dominant_str = f"dominant=phase{dominant_phase} ({dominant_pct:5.2f}%)"
        print(
            f"Length {row['read_length']:2d}: "
            f"n={row['n_reads']:>6,}, "
            f"phase0={row['phase0_pct']:5.2f}%, "
            f"phase1={row['phase1_pct']:5.2f}%, "
            f"phase2={row['phase2_pct']:5.2f}%, "
            f"{dominant_str}"
        )

    print()
    print(f"Output path (periodicity TSV): {periodicity_path.resolve()}")
    print(f"Output path (run statistics TSV): {run_stats_path.resolve()}")
    print(f"Output path (length statistics TSV): {length_stats_path.resolve()}")

    print()
    print("PASS")

    print()
    print("Read-length-specific 3-nt periodicity QC completed.")

    return 0


if __name__ == "__main__":

    sys.exit(main())
