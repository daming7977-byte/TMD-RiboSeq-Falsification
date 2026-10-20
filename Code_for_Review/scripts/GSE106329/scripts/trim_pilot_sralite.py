#!/usr/bin/env python3
"""Trim all validated Pilot A/B SRA Lite runs and merge technical runs by GSM."""

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "metadata" / "srr_manifest.tsv"
SAMPLES = ROOT / "metadata" / "gsm_manifest.tsv"
CACHE = ROOT / "raw" / "sra_lite_cache"
OUT = ROOT / "derived" / "trimmed"
QC = ROOT / "qc"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sra_path(srr):
    directory = CACHE / srr
    for suffix in ("sralite", "sra"):
        path = directory / f"{srr}.{suffix}"
        if path.is_file():
            return path
    raise FileNotFoundError(srr)


def transform(row, sample, cores):
    run_dir = OUT / "runs" / row["srr"]
    run_dir.mkdir(parents=True, exist_ok=True)
    final = run_dir / f"{row['srr']}.trimmed.fastq.gz"
    metrics_path = run_dir / "metrics.json"
    if final.is_file() and metrics_path.is_file():
        metrics = json.loads(metrics_path.read_text())
        if metrics.get("status") == "PASS" and metrics.get("sha256") == sha256(final):
            return metrics

    temporary = run_dir / f"{row['srr']}.trimmed.fastq.gz.tmp"
    cutadapt_json = run_dir / "cutadapt.json"
    log_path = run_dir / "processing.log"
    lengths = Counter()
    retained = barcode_mismatch = malformed = 0
    with log_path.open("w") as log, gzip.open(temporary, "wt", compresslevel=4) as output:
        dump = subprocess.Popen(
            ["fastq-dump", "-Z", str(sra_path(row["srr"]))],
            stdout=subprocess.PIPE, stderr=log, text=True,
        )
        cut = subprocess.Popen(
            [
                "cutadapt", "-a", sample["three_prime_adapter"],
                "-e", "0.10", "-O", "8", "-m", "24",
                "--discard-untrimmed", "--cores", str(cores),
                "--json", str(cutadapt_json), "-o", "-", "-",
            ],
            stdin=dump.stdout, stdout=subprocess.PIPE, stderr=log, text=True,
        )
        dump.stdout.close()
        handle = cut.stdout
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip("\n")
            plus = handle.readline()
            qual = handle.readline().rstrip("\n")
            if not plus or len(seq) != len(qual):
                malformed += 1
                continue
            if sample["linker_type"] == "barcoded":
                barcode = sample["sample_barcode"]
                if len(seq) < 10 or not seq.endswith(barcode):
                    barcode_mismatch += 1
                    continue
                umi = seq[-9:-4]
                new_seq, new_qual = seq[1:-9], qual[1:-9]
                first, *rest = header.rstrip("\n").split(" ", 1)
                header = first + "_UMI_" + umi + ((" " + rest[0]) if rest else "") + "\n"
            else:
                new_seq, new_qual = seq[1:], qual[1:]
            if not new_seq:
                malformed += 1
                continue
            output.write(header)
            output.write(new_seq + "\n")
            output.write(plus)
            output.write(new_qual + "\n")
            lengths[len(new_seq)] += 1
            retained += 1
        handle.close()
        cut_rc = cut.wait()
        dump_rc = dump.wait()
    if cut_rc or dump_rc or malformed:
        raise RuntimeError(f"{row['srr']}: dump={dump_rc}, cutadapt={cut_rc}, malformed={malformed}")
    os.replace(temporary, final)
    report = json.loads(cutadapt_json.read_text())
    metrics = {
        "status": "PASS", "pilot": row["pilot"], "gsm": row["gsm"],
        "srr": row["srr"], "expected_raw_spots": int(row["read_count"]),
        "cutadapt_input_reads": int(report["read_counts"]["input"]),
        "cutadapt_output_reads": int(report["read_counts"]["output"]),
        "barcode_mismatch_after_clipping": barcode_mismatch,
        "retained_reads": retained,
        "length_counts": dict(sorted(lengths.items())),
        "path": str(final), "bytes": final.stat().st_size, "sha256": sha256(final),
    }
    if metrics["cutadapt_input_reads"] != metrics["expected_raw_spots"]:
        raise RuntimeError(f"{row['srr']}: exported spot count disagrees with manifest")
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    return metrics


def merge_gsm(gsm, rows):
    directory = OUT / "gsm"
    directory.mkdir(parents=True, exist_ok=True)
    final = directory / f"{gsm}.trimmed.fastq.gz"
    temporary = directory / f"{gsm}.trimmed.fastq.gz.tmp"
    with temporary.open("wb") as output:
        for row in sorted(rows, key=lambda x: int(x["technical_run_index"])):
            source = OUT / "runs" / row["srr"] / f"{row['srr']}.trimmed.fastq.gz"
            with source.open("rb") as handle:
                shutil.copyfileobj(handle, output, 8 * 1024 * 1024)
    os.replace(temporary, final)
    return final


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", choices=["pilot", "dataset"], default="pilot")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--cores-per-worker", type=int, default=2)
    args = parser.parse_args()
    samples = {r["gsm"]: r for r in csv.DictReader(SAMPLES.open(), delimiter="\t")}
    rows = list(csv.DictReader(RUNS.open(), delimiter="\t"))
    if args.scope == "pilot":
        rows = [r for r in rows if r["pilot"] in {"A", "B"}]
    OUT.mkdir(parents=True, exist_ok=True)
    metrics = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(transform, row, samples[row["gsm"]], args.cores_per_worker): row for row in rows}
        for future in as_completed(futures):
            result = future.result()
            metrics.append(result)
            print(f"PASS\t{result['pilot']}\t{result['gsm']}\t{result['srr']}\t{result['retained_reads']}", flush=True)

    by_gsm = defaultdict(list)
    for row in rows:
        by_gsm[row["gsm"]].append(row)
    gsm_paths = {gsm: merge_gsm(gsm, members) for gsm, members in sorted(by_gsm.items())}
    fields = ["pilot", "gsm", "srr", "expected_raw_spots", "cutadapt_input_reads",
              "cutadapt_output_reads", "barcode_mismatch_after_clipping", "retained_reads",
              "bytes", "sha256", "path"]
    QC.mkdir(exist_ok=True)
    with (QC / "trimming_metrics_by_srr.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(sorted(metrics, key=lambda x: int(x["srr"][3:])))
    with (QC / "trimmed_length_counts_by_srr.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["pilot", "gsm", "srr", "length", "count"])
        for x in sorted(metrics, key=lambda z: int(z["srr"][3:])):
            for length, count in x["length_counts"].items():
                writer.writerow([x["pilot"], x["gsm"], x["srr"], length, count])
    with (QC / "trimmed_gsm_files.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["gsm", "technical_run_count", "path", "bytes", "sha256"])
        for gsm, path in gsm_paths.items():
            writer.writerow([gsm, len(by_gsm[gsm]), path, path.stat().st_size, sha256(path)])
    print(f"PASS: all runs trimmed and technical runs merged into {len(gsm_paths)} GSM files", flush=True)


if __name__ == "__main__":
    main()
