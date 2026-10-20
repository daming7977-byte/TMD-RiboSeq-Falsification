#!/usr/bin/env python3
"""Filter contaminants and align trimmed Pilot GSM libraries to R64-1-1."""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SAMPLES = ROOT / "metadata" / "gsm_manifest.tsv"
TRIM = ROOT / "derived" / "trimmed" / "gsm"
FILTERED = ROOT / "derived" / "filtered"
ALIGN = ROOT / "derived" / "alignments"
QC = ROOT / "qc"
CONTAM = ROOT / "references" / "genome" / "indices" / "contaminants" / "rrna_trna_adapters"
GENOME = ROOT / "references" / "genome" / "indices" / "genome" / "R64-1-1"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bowtie_counts(text):
    total = re.search(r"^\s*([0-9,]+) reads; of these:", text, re.M)
    zero = re.search(r"^\s*([0-9,]+) \([^\n]+\) aligned 0 times", text, re.M)
    once = re.search(r"^\s*([0-9,]+) \([^\n]+\) aligned exactly 1 time", text, re.M)
    multi = re.search(r"^\s*([0-9,]+) \([^\n]+\) aligned >1 times", text, re.M)
    if not all((total, zero, once, multi)):
        raise ValueError("Could not parse Bowtie2 summary:\n" + text[-4000:])
    return {k: int(v.group(1).replace(",", "")) for k, v in
            (("total", total), ("zero", zero), ("once", once), ("multi", multi))}


def count_bam(path, args):
    completed = subprocess.run(["samtools", "view", "-c"] + args + [str(path)],
                               capture_output=True, text=True, check=True)
    return int(completed.stdout.strip())


def process(sample, threads):
    gsm = sample["gsm"]
    FILTERED.mkdir(parents=True, exist_ok=True)
    ALIGN.mkdir(parents=True, exist_ok=True)
    filtered = FILTERED / f"{gsm}.noncontaminant.fastq.gz"
    bam = ALIGN / f"{gsm}.sorted.bam"
    metrics_path = ALIGN / f"{gsm}.metrics.json"
    contam_log = ALIGN / f"{gsm}.contaminant.log"
    genome_log = ALIGN / f"{gsm}.genome.log"
    if bam.is_file() and Path(str(bam) + ".bai").is_file() and metrics_path.is_file():
        metrics = json.loads(metrics_path.read_text())
        if metrics.get("status") == "PASS":
            return metrics

    filtered_tmp = Path(str(filtered) + ".tmp.gz")
    with contam_log.open("w") as log:
        completed = subprocess.run([
            "bowtie2", "--end-to-end", "--very-sensitive", "-N", "0", "-L", "15",
            "-p", str(threads), "-x", str(CONTAM), "-U", str(TRIM / f"{gsm}.trimmed.fastq.gz"),
            "--un-gz", str(filtered_tmp), "--no-unal", "-S", "/dev/null",
        ], stdout=subprocess.DEVNULL, stderr=log)
    if completed.returncode:
        raise RuntimeError(f"{gsm}: contaminant Bowtie2 failed")
    os.replace(filtered_tmp, filtered)
    contam_text = contam_log.read_text()
    contaminant = bowtie_counts(contam_text)

    bam_tmp = Path(str(bam) + ".tmp")
    with genome_log.open("w") as log:
        bowtie = subprocess.Popen([
            "bowtie2", "--end-to-end", "--very-sensitive", "-N", "0", "-L", "20",
            "-p", str(threads), "--rg-id", gsm, "--rg", f"SM:{gsm}",
            "-x", str(GENOME), "-U", str(filtered), "--no-unal",
        ], stdout=subprocess.PIPE, stderr=log)
        view = subprocess.Popen(["samtools", "view", "-b", "-"], stdin=bowtie.stdout,
                                stdout=subprocess.PIPE, stderr=log)
        bowtie.stdout.close()
        sort = subprocess.run(["samtools", "sort", "-@", str(threads), "-o", str(bam_tmp), "-"],
                              stdin=view.stdout, stderr=log)
        view.stdout.close()
        view_rc = view.wait()
        bowtie_rc = bowtie.wait()
    if any((sort.returncode, view_rc, bowtie_rc)):
        raise RuntimeError(f"{gsm}: genome alignment pipeline failed")
    os.replace(bam_tmp, bam)
    subprocess.run(["samtools", "index", "-@", str(threads), str(bam)], check=True)
    genome = bowtie_counts(genome_log.read_text())
    aligned_primary = count_bam(bam, ["-F", "260"])
    unique = count_bam(bam, ["-F", "260", "-q", "20"])
    unique_nm0 = count_bam(bam, ["-F", "260", "-q", "20", "-e", "[NM] == 0"])
    metrics = {
        "status": "PASS", "pilot": sample["pilot"], "gsm": gsm,
        "series_id": sample["series_id"], "fraction": sample["fraction"],
        "trimmed_reads": contaminant["total"],
        "contaminant_reads": contaminant["total"] - contaminant["zero"],
        "noncontaminant_reads": contaminant["zero"],
        "genome_input_reads": genome["total"], "genome_unmapped_reads": genome["zero"],
        "genome_aligned_primary": aligned_primary, "unique_mapq20": unique,
        "unique_mapq20_nm0": unique_nm0, "bam": str(bam),
        "bam_bytes": bam.stat().st_size, "bam_sha256": sha256(bam),
    }
    if metrics["noncontaminant_reads"] != metrics["genome_input_reads"]:
        raise RuntimeError(f"{gsm}: filter/alignment handoff count mismatch")
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    return metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", choices=["pilot", "dataset"], default="pilot")
    parser.add_argument("--gsm")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--threads-per-worker", type=int, default=2)
    args = parser.parse_args()
    samples = list(csv.DictReader(SAMPLES.open(), delimiter="\t"))
    if args.scope == "pilot":
        samples = [r for r in samples if r["pilot"] in {"A", "B"}]
    if args.gsm:
        samples = [r for r in samples if r["gsm"] == args.gsm]
    metrics = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process, sample, args.threads_per_worker): sample for sample in samples}
        for future in as_completed(futures):
            result = future.result(); metrics.append(result)
            print(f"PASS\t{result['pilot']}\t{result['gsm']}\t{result['unique_mapq20']}", flush=True)
    if not args.gsm:
        fields = ["pilot", "gsm", "series_id", "fraction", "trimmed_reads", "contaminant_reads",
                  "noncontaminant_reads", "genome_input_reads", "genome_unmapped_reads",
                  "genome_aligned_primary", "unique_mapq20", "unique_mapq20_nm0",
                  "bam", "bam_bytes", "bam_sha256"]
        QC.mkdir(exist_ok=True)
        with (QC / "mapping_metrics_by_gsm.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader(); writer.writerows(sorted(metrics, key=lambda x: int(x["gsm"][3:])))
        print(f"PASS: all {len(samples)} selected GSM libraries aligned", flush=True)


if __name__ == "__main__":
    main()
