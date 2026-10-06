#!/usr/bin/env python3
"""Exact-family UMI sensitivity with 5-nt occupancy collision diagnostics."""

import argparse
import csv
import heapq
import json
import math
import os
import re
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata" / "gsm_manifest.tsv"
SOURCE = ROOT / "derived" / "alignments"
OUT = ROOT / "derived" / "umi_deduplicated"
QC = ROOT / "qc"
UMI_RE = re.compile(r"_UMI_([ACGTN]{5})$")
UMI_SPACE = 4**5


def family_size_bin(size):
    if size <= 5:
        return str(size)
    if size <= 10:
        return "6-10"
    if size <= 20:
        return "11-20"
    if size <= 50:
        return "21-50"
    return "51+"


def ranking(read):
    nm = read.get_tag("NM") if read.has_tag("NM") else 10**9
    return (-int(read.mapping_quality), int(nm), read.query_name)


def process(sample, threads):
    gsm = sample["gsm"]
    source = SOURCE / f"{gsm}.sorted.bam"
    final = OUT / f"{gsm}.umi_dedup.sorted.bam"
    metrics_path = OUT / f"{gsm}.umi_metrics.json"
    if final.is_file() and Path(str(final) + ".bai").is_file() and metrics_path.is_file():
        metrics = json.loads(metrics_path.read_text())
        if metrics.get("algorithm_schema") == 1:
            return metrics

    OUT.mkdir(parents=True, exist_ok=True)
    unsorted = OUT / f"{gsm}.umi_dedup.unsorted.bam"
    before_lengths = Counter()
    after_lengths = Counter()
    size_distribution = Counter()
    metrics = Counter()
    collision_adjusted = 0.0
    saturated_groups = 0

    with pysam.AlignmentFile(str(source), "rb") as bam, pysam.AlignmentFile(
        str(unsorted), "wb", header=bam.header
    ) as output:
        groups = defaultdict(dict)
        heap = []
        queued = set()
        contig = None

        def flush_group(group_key):
            nonlocal collision_adjusted, saturated_groups
            families = groups.pop(group_key, None)
            queued.discard(group_key)
            if families is None:
                return
            reads_in_group = sum(item[1] for item in families.values())
            observed_umis = len(families)
            if observed_umis >= UMI_SPACE:
                adjusted = float(reads_in_group)
                saturated_groups += 1
            else:
                adjusted = min(
                    float(reads_in_group),
                    -UMI_SPACE * math.log1p(-observed_umis / UMI_SPACE),
                )
            collision_adjusted += adjusted
            metrics["positional_groups"] += 1
            metrics["observed_families"] += observed_umis
            for representative, size in families.values():
                output.write(representative)
                after_lengths[int(representative.query_length or 0)] += 1
                size_distribution[family_size_bin(size)] += 1

        def flush_before(reference_start):
            while heap and heap[0][0] < reference_start:
                _, group_key = heapq.heappop(heap)
                if group_key in queued:
                    flush_group(group_key)

        def flush_all():
            for group_key in list(groups):
                flush_group(group_key)
            heap.clear()

        for read in bam.fetch(until_eof=True):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.mapping_quality < 20
            ):
                continue
            if contig is None or read.reference_name != contig:
                if contig is not None:
                    flush_all()
                contig = read.reference_name
            flush_before(int(read.reference_start))
            match = UMI_RE.search(read.query_name)
            metrics["eligible_before"] += 1
            before_lengths[int(read.query_length or 0)] += 1
            if not match:
                metrics["umi_missing"] += 1
                continue
            metrics["umi_recognized"] += 1
            aligned = read.get_reference_positions(full_length=False)
            if not aligned:
                metrics["no_aligned_query_base"] += 1
                continue
            five_prime = int(aligned[-1] if read.is_reverse else aligned[0])
            strand = "-" if read.is_reverse else "+"
            length = int(read.query_length or 0)
            group_key = (five_prime, strand, length)
            if group_key not in queued:
                heapq.heappush(heap, (five_prime, group_key))
                queued.add(group_key)
            umi = match.group(1)
            previous = groups[group_key].get(umi)
            if previous is None:
                groups[group_key][umi] = [read, 1]
            else:
                previous[1] += 1
                if ranking(read) < ranking(previous[0]):
                    previous[0] = read
        flush_all()

    sorted_tmp = OUT / f"{gsm}.umi_dedup.sorted.bam.tmp"
    subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", str(sorted_tmp), str(unsorted)],
        check=True,
    )
    os.replace(sorted_tmp, final)
    unsorted.unlink()
    subprocess.run(["samtools", "index", "-@", str(threads), str(final)], check=True)

    before = int(metrics["umi_recognized"])
    after = int(metrics["observed_families"])
    result = {
        "algorithm_schema": 1,
        "sample_id": gsm,
        "series_id": sample["series_id"],
        "time_min": sample["auxin_time_min"],
        "fraction": sample["fraction"],
        "reads_before_dedup": before,
        "reads_after_exact_family_dedup": after,
        "apparent_duplicate_reads": before - after,
        "apparent_duplicate_fraction": (before - after) / before if before else 0.0,
        "exact_family_retained_fraction": after / before if before else 0.0,
        "collision_adjusted_molecules": collision_adjusted,
        "collision_adjusted_retained_fraction": collision_adjusted / before if before else 0.0,
        "collision_adjustment_over_observed_families": collision_adjusted / after if after else 0.0,
        "positional_groups": int(metrics["positional_groups"]),
        "saturated_positional_groups": saturated_groups,
        "umi_missing": int(metrics["umi_missing"]),
        "umi_recognition_fraction": before / int(metrics["eligible_before"]),
        "before_length_counts": dict(sorted(before_lengths.items())),
        "after_length_counts": dict(sorted(after_lengths.items())),
        "family_size_distribution": dict(size_distribution),
        "bam": str(final),
    }
    metrics_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--threads-per-worker", type=int, default=2)
    parser.add_argument("--gsm")
    args = parser.parse_args()
    samples = [
        row for row in csv.DictReader(META.open(), delimiter="\t")
        if row["linker_type"] == "barcoded"
    ]
    if args.gsm:
        samples = [row for row in samples if row["gsm"] == args.gsm]
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(process, sample, args.threads_per_worker): sample
            for sample in samples
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(
                result["sample_id"], result["reads_before_dedup"],
                result["reads_after_exact_family_dedup"],
                result["collision_adjusted_molecules"], flush=True,
            )
    if args.gsm:
        return
    fields = [
        "sample_id", "series_id", "time_min", "fraction", "reads_before_dedup",
        "reads_after_exact_family_dedup", "apparent_duplicate_reads",
        "apparent_duplicate_fraction", "exact_family_retained_fraction",
        "collision_adjusted_molecules", "collision_adjusted_retained_fraction",
        "collision_adjustment_over_observed_families", "positional_groups",
        "saturated_positional_groups", "umi_missing", "umi_recognition_fraction", "bam",
    ]
    with (QC / "umi_dedup_metrics_by_gsm.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(sorted(results, key=lambda x: int(x["sample_id"][3:])))
    with (QC / "umi_family_size_distribution.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "family_size_bin", "observed_families"])
        for result in sorted(results, key=lambda x: int(x["sample_id"][3:])):
            for size_bin in ["1", "2", "3", "4", "5", "6-10", "11-20", "21-50", "51+"]:
                writer.writerow([result["sample_id"], size_bin, result["family_size_distribution"].get(size_bin, 0)])
    with (QC / "umi_read_length_before_after.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "read_length", "before", "after"])
        for result in sorted(results, key=lambda x: int(x["sample_id"][3:])):
            lengths = sorted(set(result["before_length_counts"]) | set(result["after_length_counts"]), key=int)
            for length in lengths:
                writer.writerow([
                    result["sample_id"], length,
                    result["before_length_counts"].get(length, 0),
                    result["after_length_counts"].get(length, 0),
                ])


if __name__ == "__main__":
    main()
