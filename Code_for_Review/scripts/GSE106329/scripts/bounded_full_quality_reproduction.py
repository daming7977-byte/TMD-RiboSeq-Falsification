#!/usr/bin/env python3
"""Bounded full-quality versus SRA Lite reproduction on three locked runs."""

import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]
FROZEN_SRC = ROOT / "frozen/v1.0.0_owner_review_candidate/wheel_extracted"
sys.path.insert(0, str(FROZEN_SRC))
from tmdriboseq.coordinates import load_gtf_transcripts  # noqa: E402
from tmdriboseq.psite import psite_reference_position0  # noqa: E402
from tmdriboseq.qc import AMBIGUOUS_FRAME, build_cds_frame_lookup  # noqa: E402


RUNS = {
    "Sec63_noCHX": "SRR6234866",
    "Sec63_CHX": "SRR6234818",
    "Ubc6_CHX": "SRR6234842",
}
N = 1_000_000
META = {r["gsm"]: r for r in csv.DictReader((ROOT / "metadata/gsm_manifest.tsv").open(), delimiter="\t")}
SRR = {r["srr"]: r for r in csv.DictReader((ROOT / "metadata/srr_manifest.tsv").open(), delimiter="\t")}
OUT = ROOT / "bounded_full_quality"
CONTAM = ROOT / "references/genome/indices/contaminants/rrna_trna_adapters"
GENOME = ROOT / "references/genome/indices/genome/R64-1-1"
GTF = ROOT / "references/genome/Saccharomyces_cerevisiae.R64-1-1.115.gtf"


def full_path(srr):
    row = SRR[srr]
    pilot = "pilot_a" if row["pilot"] == "A" else "pilot_b"
    return ROOT / "raw" / pilot / row["gsm"] / f"{srr}.fastq.gz"


def lite_path(srr):
    directory = ROOT / "raw/sra_lite_cache" / srr
    for suffix in ("sralite", "sra"):
        path = directory / f"{srr}.{suffix}"
        if path.is_file():
            return path
    raise FileNotFoundError(srr)


def extract_full(source, destination):
    count = 0
    with gzip.open(source, "rt") as inp, gzip.open(destination, "wt", compresslevel=4) as out:
        while count < N:
            lines = [inp.readline() for _ in range(4)]
            if not all(lines):
                break
            out.writelines(lines)
            count += 1
    if count != N:
        raise RuntimeError(f"{source} supplied only {count} complete reads")
    return count


def extract_lite(source, destination):
    dump = subprocess.Popen(
        ["fastq-dump", "-X", str(N), "-Z", str(source)], stdout=subprocess.PIPE, text=True
    )
    count = 0
    with gzip.open(destination, "wt", compresslevel=4) as out:
        while True:
            block = dump.stdout.read(1024 * 1024)
            if not block:
                break
            count += block.count("\n")
            out.write(block)
    if dump.wait() or count // 4 != N:
        raise RuntimeError(f"SRA Lite export failed: {source}, lines={count}")
    return N


def sequence_digest(path):
    digest = hashlib.sha256()
    count = 0
    quality = Counter()
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            handle.readline()
            qualities = handle.readline().rstrip("\n")
            digest.update(sequence.encode()); digest.update(b"\n")
            quality.update(ord(x) - 33 for x in qualities)
            count += 1
    total_q = sum(quality.values())
    return count, digest.hexdigest(), sum(q * n for q, n in quality.items()) / total_q, sum(n for q, n in quality.items() if q < 20) / total_q


def structural_trim(raw, destination, sample, work):
    adapter_out = work / (destination.stem + ".adapter.fastq.gz")
    report = work / (destination.stem + ".cutadapt.json")
    subprocess.run([
        "cutadapt", "-a", sample["three_prime_adapter"], "-e", "0.10", "-O", "8",
        "-m", "24", "--discard-untrimmed", "--json", str(report),
        "-o", str(adapter_out), str(raw),
    ], check=True, stdout=subprocess.DEVNULL)
    retained = mismatch = 0
    with gzip.open(adapter_out, "rt") as inp, gzip.open(destination, "wt", compresslevel=4) as out:
        while True:
            header = inp.readline()
            if not header:
                break
            seq = inp.readline().rstrip("\n"); plus = inp.readline(); qual = inp.readline().rstrip("\n")
            if sample["linker_type"] == "barcoded":
                barcode = sample["sample_barcode"]
                if len(seq) < 10 or not seq.endswith(barcode):
                    mismatch += 1; continue
                umi = seq[-9:-4]
                seq, qual = seq[1:-9], qual[1:-9]
                header = header.rstrip("\n") + "_UMI_" + umi + "\n"
            else:
                seq, qual = seq[1:], qual[1:]
            if not seq:
                continue
            out.write(header); out.write(seq + "\n"); out.write(plus); out.write(qual + "\n")
            retained += 1
    cut = json.loads(report.read_text())
    return int(cut["read_counts"]["output"]), mismatch, retained


def quality_trim(source, destination):
    report = Path(str(destination) + ".json")
    subprocess.run([
        "cutadapt", "-q", "20", "-m", "20", "--json", str(report),
        "-o", str(destination), str(source),
    ], check=True, stdout=subprocess.DEVNULL)
    return int(json.loads(report.read_text())["read_counts"]["output"])


def bowtie_counts(text):
    patterns = {
        "total": r"^\s*([0-9,]+) reads; of these:",
        "zero": r"^\s*([0-9,]+) \([^\n]+\) aligned 0 times",
        "once": r"^\s*([0-9,]+) \([^\n]+\) aligned exactly 1 time",
        "multi": r"^\s*([0-9,]+) \([^\n]+\) aligned >1 times",
    }
    return {key: int(re.search(pattern, text, re.M).group(1).replace(",", "")) for key, pattern in patterns.items()}


def align(source, prefix):
    filtered = Path(str(prefix) + ".noncontam.fastq.gz")
    contam_log = Path(str(prefix) + ".contam.log")
    with contam_log.open("w") as log:
        subprocess.run([
            "bowtie2", "--end-to-end", "--very-sensitive", "-N", "0", "-L", "15",
            "-x", str(CONTAM), "-U", str(source), "--un-gz", str(filtered),
            "--no-unal", "-S", "/dev/null",
        ], check=True, stdout=subprocess.DEVNULL, stderr=log)
    bam = Path(str(prefix) + ".sorted.bam")
    genome_log = Path(str(prefix) + ".genome.log")
    with genome_log.open("w") as log:
        bowtie = subprocess.Popen([
            "bowtie2", "--end-to-end", "--very-sensitive", "-N", "0", "-L", "20",
            "-x", str(GENOME), "-U", str(filtered), "--no-unal",
        ], stdout=subprocess.PIPE, stderr=log)
        view = subprocess.Popen(["samtools", "view", "-b", "-"], stdin=bowtie.stdout, stdout=subprocess.PIPE)
        bowtie.stdout.close()
        sort = subprocess.run(["samtools", "sort", "-o", str(bam), "-"], stdin=view.stdout)
        view.stdout.close(); view_rc = view.wait(); bowtie_rc = bowtie.wait()
    if sort.returncode or view_rc or bowtie_rc:
        raise RuntimeError(f"alignment failed for {source}")
    subprocess.run(["samtools", "index", str(bam)], check=True)
    contam = bowtie_counts(contam_log.read_text()); genome = bowtie_counts(genome_log.read_text())
    unique = int(subprocess.run(["samtools", "view", "-c", "-q", "20", "-F", "260", str(bam)], capture_output=True, text=True, check=True).stdout)
    return bam, contam, genome, unique


transcripts = load_gtf_transcripts(GTF)
frame_lookup, _ = build_cds_frame_lookup(transcripts)


def calibrate(bam):
    lengths = Counter(); frames = defaultdict(Counter)
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < 20:
                continue
            length = int(read.query_length or 0); lengths[length] += 1
            if length not in range(26, 33):
                continue
            strand = "-" if read.is_reverse else "+"
            for offset in range(10, 18):
                position = psite_reference_position0(read, offset)
                frame = frame_lookup.get((read.reference_name, strand, position)) if position is not None else None
                if frame is not None and frame != AMBIGUOUS_FRAME:
                    frames[(length, offset)][frame] += 1
    offsets = {}; per_length = []
    for length in range(26, 33):
        candidates = []
        for phase, canonical in [(0, 12), (1, 13), (2, 14)]:
            aggregate = Counter()
            for offset in range(10, 18):
                if offset % 3 == phase:
                    aggregate.update(frames[(length, offset)])
            total = sum(aggregate.values()); fraction = aggregate[0] / total if total else 0.0
            candidates.append((fraction, total, canonical))
        fraction, total, offset = max(candidates); offsets[length] = offset
        per_length.append((length, offset, fraction, total))
    aggregate = Counter()
    for length, offset in offsets.items():
        aggregate.update(frames[(length, offset)])
    cds = sum(aggregate.values()); eligible = sum(lengths[x] for x in range(26, 33))
    return {
        "eligible": sum(lengths.values()), "usable": eligible,
        "usable_fraction": eligible / sum(lengths.values()) if lengths else 0.0,
        "cds": cds, "cds_fraction": cds / eligible if eligible else 0.0,
        "frame0_fraction": aggregate[0] / cds if cds else 0.0,
        "per_length": per_length,
    }


def main():
    OUT.mkdir(exist_ok=True)
    metrics_rows = []; offset_rows = []
    for label, srr in RUNS.items():
        work = OUT / srr; work.mkdir(exist_ok=True)
        sample = META[SRR[srr]["gsm"]]
        raw_paths = {"full_quality": work / "full.raw.fastq.gz", "sra_lite": work / "lite.raw.fastq.gz"}
        extract_full(full_path(srr), raw_paths["full_quality"])
        extract_lite(lite_path(srr), raw_paths["sra_lite"])
        raw_stats = {source: sequence_digest(path) for source, path in raw_paths.items()}
        if raw_stats["full_quality"][1] != raw_stats["sra_lite"][1]:
            raise RuntimeError(f"{srr}: first {N} raw sequences differ")
        for source, raw in raw_paths.items():
            structural = work / f"{source}.structural.fastq.gz"
            adapter_out, mismatch, retained = structural_trim(raw, structural, sample, work)
            for mode in ["structural", "q20_after_structural"]:
                reads = structural
                mode_retained = retained
                if mode == "q20_after_structural":
                    reads = work / f"{source}.q20.fastq.gz"
                    mode_retained = quality_trim(structural, reads)
                prefix = work / f"{source}.{mode}"
                bam, contam, genome, unique = align(reads, prefix)
                qc = calibrate(bam)
                metrics_rows.append({
                    "class": label, "srr": srr, "source": source, "mode": mode,
                    "raw_subset_reads": raw_stats[source][0],
                    "sequence_sha256": raw_stats[source][1],
                    "mean_raw_phred": raw_stats[source][2],
                    "raw_bases_below_q20_fraction": raw_stats[source][3],
                    "adapter_trimmed_reads": adapter_out,
                    "barcode_mismatch": mismatch,
                    "retained_reads": mode_retained,
                    "retained_fraction": mode_retained / N,
                    "contaminant_fraction": (contam["total"] - contam["zero"]) / contam["total"],
                    "unique_mapq20": unique,
                    "unique_fraction_of_retained": unique / mode_retained if mode_retained else 0.0,
                    "usable_26_32_fraction": qc["usable_fraction"],
                    "cds_assignment_fraction": qc["cds_fraction"],
                    "frame0_fraction": qc["frame0_fraction"],
                })
                for length, offset, frame0, cds_reads in qc["per_length"]:
                    offset_rows.append({
                        "class": label, "srr": srr, "source": source, "mode": mode,
                        "read_length": length, "offset": offset,
                        "frame0_fraction": frame0, "cds_psites_across_phase": cds_reads,
                    })
                print(label, source, mode, mode_retained, unique, qc["frame0_fraction"], flush=True)
    for path, rows in [
        (OUT / "bounded_reproduction_metrics.tsv", metrics_rows),
        (OUT / "bounded_reproduction_offsets.tsv", offset_rows),
    ]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()
