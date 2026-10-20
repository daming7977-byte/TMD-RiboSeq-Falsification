#!/usr/bin/env python3
"""Create deterministic genome-wide calibration BAMs for all 30 GSM."""

import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAPPING = ROOT / "qc/mapping_metrics_by_gsm.tsv"
META = ROOT / "metadata/gsm_manifest.tsv"
SOURCE = ROOT / "derived/alignments"
OUT = ROOT / "derived/stage1b_calibration_bams"
SAMPLE_OUT = ROOT / "metadata/stage1b_calibration_samples.tsv"
TARGET = 5_000_000
SEED = 731


mapping = {r["gsm"]: r for r in csv.DictReader(MAPPING.open(), delimiter="\t")}
samples = list(csv.DictReader(META.open(), delimiter="\t"))
OUT.mkdir(parents=True, exist_ok=True)
rows = []
for sample in samples:
    gsm = sample["gsm"]
    unique = int(mapping[gsm]["unique_mapq20"])
    fraction = min(1.0, TARGET / max(1, unique))
    source = SOURCE / f"{gsm}.sorted.bam"
    destination = OUT / f"{gsm}.calibration.bam"
    command = ["samtools", "view", "-bh", "-q", "20", "-F", "260"]
    if fraction < 1.0:
        digits = f"{fraction:.12f}".split(".", 1)[1]
        command += ["-s", f"{SEED}.{digits}"]
    command += ["-o", str(destination), str(source)]
    subprocess.run(command, check=True)
    subprocess.run(["samtools", "index", str(destination)], check=True)
    observed = int(subprocess.run(
        ["samtools", "view", "-c", str(destination)], capture_output=True,
        text=True, check=True,
    ).stdout)
    rows.append({
        "sample_id": gsm, "bam": destination, "replicate": gsm, "role": "qc",
        "series": sample["series_id"], "fraction": sample["fraction"],
        "source_unique_reads": unique, "sampling_fraction": fraction,
        "calibration_reads": observed,
    })
    print(gsm, observed, flush=True)

with SAMPLE_OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader(); writer.writerows(rows)

