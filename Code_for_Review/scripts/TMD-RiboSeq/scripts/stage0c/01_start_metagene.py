#!/usr/bin/env python3
"""
Stage 0C

Module 01

Start Codon Metagene

Version:
0.2.0
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

from lib.io import read_features
from lib.bam import open_bam


START_WINDOW_MIN = -60
START_WINDOW_MAX = 120
MIN_READ_LENGTH = 25
MAX_READ_LENGTH = 35
OUTPUT_FILENAME = "start_metagene.tsv"


# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description="Generate Start Codon Metagene."
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Unique mapped BAM"
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

    parser.add_argument(
        "--window-min",
        type=int,
        default=START_WINDOW_MIN,
        help=f"Minimum relative position (default: {START_WINDOW_MIN})"
    )

    parser.add_argument(
        "--window-max",
        type=int,
        default=START_WINDOW_MAX,
        help=f"Maximum relative position (default: {START_WINDOW_MAX})"
    )

    return parser.parse_args()


# ---------------------------------------------------------
# Metagene helpers
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


def write_metagene(counts, read_lengths, window_min, window_max, output_path):
    """Write a deterministic, complete read-length/position grid."""

    output_rows = 0
    output_count = 0

    with output_path.open("w", newline="", encoding="utf-8") as output_handle:

        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "read_length",
            "relative_position",
            "count",
        ])

        for length in sorted(read_lengths):

            for position in range(window_min, window_max + 1):

                count = counts[(length, position)]
                writer.writerow([length, position, count])
                output_rows += 1
                output_count += count

    return output_rows, output_count


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    args = parse_args()

    if args.window_min > args.window_max:
        print(
            "ERROR: --window-min must be less than or equal to --window-max.",
            file=sys.stderr,
        )
        return 1

    print("=" * 60)
    print("Stage0C Module 01")
    print("Start Codon Metagene v0.2")
    print("=" * 60)

    print()

    print("[1/3] Loading transcript features...")

    features = read_features(args.features)

    print(
        f"Loaded {len(features):,} transcripts"
    )

    print()

    print("[2/3] Opening BAM...")

    bam = open_bam(args.bam)

    print("BAM opened successfully")

    print()

    print("[3/3] Calculating raw 5'-end start-codon metagene...")

    counts = Counter()
    read_lengths = set(range(MIN_READ_LENGTH, MAX_READ_LENGTH + 1))
    skipped = Counter()

    total_records = 0
    primary_mapped = 0
    assigned_reads = 0
    reads_inside_window = 0

    try:

        for read in bam.fetch(until_eof=True):

            total_records += 1

            if read.is_unmapped:
                skipped["unmapped"] += 1
                continue

            if read.is_secondary:
                skipped["secondary_alignment"] += 1
                continue

            if read.is_supplementary:
                skipped["supplementary_alignment"] += 1
                continue

            primary_mapped += 1

            transcript = read.reference_name

            if transcript is None:
                skipped["missing_reference_name"] += 1
                continue

            if transcript not in features:
                skipped["reference_not_in_features"] += 1
                continue

            length = get_read_length(read)

            if length is None or length <= 0:
                skipped["missing_or_invalid_read_length"] += 1
                continue

            if length < MIN_READ_LENGTH or length > MAX_READ_LENGTH:
                skipped["read_length_outside_25_35"] += 1
                continue

            read_5prime = get_five_prime_position(read)

            if read_5prime is None or read_5prime < 0:
                skipped["missing_or_invalid_5prime_position"] += 1
                continue

            cds_start = features[transcript]["utr5"]
            relative_position = read_5prime - cds_start

            assigned_reads += 1

            if args.window_min <= relative_position <= args.window_max:
                counts[(length, relative_position)] += 1
                reads_inside_window += 1

    finally:

        bam.close()

    output_dir = Path(args.output)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path = output_dir / OUTPUT_FILENAME

    output_rows, output_count = write_metagene(
        counts,
        read_lengths,
        args.window_min,
        args.window_max,
        output_path,
    )

    if output_count != reads_inside_window:
        print(
            "ERROR: Internal consistency check failed: "
            f"output count sum ({output_count:,}) != reads inside window "
            f"({reads_inside_window:,}).",
            file=sys.stderr,
        )
        return 1

    print()
    print("Run statistics")
    print(f"Total BAM records examined: {total_records:,}")
    print(f"Primary mapped reads: {primary_mapped:,}")
    print(f"Reads assigned to transcript/CDS coordinates: {assigned_reads:,}")
    print(f"Reads inside metagene window: {reads_inside_window:,}")
    print(f"Read lengths observed: {len(read_lengths):,}")
    print(f"Observed read lengths: {','.join(map(str, sorted(read_lengths)))}")
    print(f"Output rows: {output_rows:,}")
    print(f"Output path: {output_path.resolve()}")

    print()
    print("Skipped records")

    if skipped:
        for reason in sorted(skipped):
            print(f"{reason}: {skipped[reason]:,}")
    else:
        print("None")

    print()

    print("PASS")

    print()

    print("Raw 5'-end start-codon metagene completed.")

    return 0


if __name__ == "__main__":

    sys.exit(main())
