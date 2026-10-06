#!/usr/bin/env python3
"""Freeze TMD-blind series-consensus and prespecified P-site sensitivities."""

import csv
import hashlib
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CAL = ROOT / "qc/stage1b_calibration/qc"
PERIODICITY = CAL / "periodicity_by_sample_length_offset.tsv"
SELECTED = CAL / "psite_offsets.tsv"
META = ROOT / "metadata/gsm_manifest.tsv"
OUT = ROOT / "frozen/psite"
PRIMARY = OUT / "psite_series_consensus_frozen.tsv"
SENSITIVITY = OUT / "psite_sensitivity_schemes_frozen.tsv"
CHECKSUM = OUT / "SHA256SUMS.txt"


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


meta = {r["gsm"]: r for r in csv.DictReader(META.open(), delimiter="\t")}
periodicity = list(csv.DictReader(PERIODICITY.open(), delimiter="\t"))
by_key = {(r["sample_id"], int(r["read_length"]), int(r["offset"])): r for r in periodicity}
sample_offsets = {
    (r["sample_id"], int(r["read_length"])): int(r["offset"])
    for r in csv.DictReader(SELECTED.open(), delimiter="\t")
}
series_samples = defaultdict(list)
for sample, row in meta.items():
    series_samples[row["series_id"]].append(sample)

primary_rows = []; consensus = {}
for series, samples in sorted(series_samples.items()):
    for length in range(26, 33):
        phase_candidates = []
        for phase, canonical in [(0, 12), (1, 13), (2, 14)]:
            counts = Counter()
            for sample in samples:
                for offset in range(10, 18):
                    if offset % 3 != phase:
                        continue
                    row = by_key[(sample, length, offset)]
                    for frame in range(3):
                        counts[frame] += int(row[f"frame{frame}"])
            total = sum(counts.values())
            phase_candidates.append((counts[0] / total if total else 0.0, total, canonical))
        _, _, offset = max(phase_candidates)
        consensus[(series, length)] = offset
        support = Counter()
        sample_support = []
        disagreement = Counter()
        for sample in samples:
            row = by_key[(sample, length, offset)]
            counts = Counter({frame: int(row[f"frame{frame}"]) for frame in range(3)})
            support.update(counts)
            total = sum(counts.values())
            sample_support.append(counts[0] / total if total else 0.0)
            disagreement[sample_offsets[(sample, length)]] += 1
        total = sum(support.values()); frame0 = support[0] / total if total else 0.0
        status = "PASS" if total >= 10000 and frame0 >= 0.50 else "WARN"
        primary_rows.append({
            "series_id": series, "read_length": length, "selected_offset_5prime_0based": offset,
            "pooled_frame0": support[0], "pooled_frame1": support[1], "pooled_frame2": support[2],
            "pooled_cds_psites": total, "pooled_frame0_fraction": frame0,
            "sample_frame0_min": min(sample_support), "sample_frame0_max": max(sample_support),
            "independent_sample_offset_counts": ";".join(f"{k}:{v}" for k, v in sorted(disagreement.items())),
            "samples_matching_consensus": disagreement[offset], "samples_total": len(samples),
            "agreement_fraction": disagreement[offset] / len(samples), "qc_status": status,
            "selection_evidence": "genome-wide CDS frame/start-stop only; no TMD coordinates or phenotypes",
        })

sensitivity_rows = []
for sample, row in sorted(meta.items()):
    series = row["series_id"]
    for length in range(26, 33):
        for scheme, offset in [
            ("primary_series_consensus", consensus[(series, length)]),
            ("sample_specific", sample_offsets[(sample, length)]),
            ("author_L_minus_15", length - 15),
            ("adjacent_L_minus_16", length - 16),
        ]:
            support = by_key[(sample, length, offset)]
            sensitivity_rows.append({
                "sample_id": sample, "series_id": series, "time_min": row["auxin_time_min"],
                "fraction": row["fraction"], "read_length": length, "scheme": scheme,
                "offset_5prime_0based": offset, "frame0_fraction": support["frame0_fraction"],
                "cds_psites": support["cds_psites"],
            })

OUT.mkdir(parents=True, exist_ok=True)
write(PRIMARY, primary_rows); write(SENSITIVITY, sensitivity_rows)
with CHECKSUM.open("w") as handle:
    for path in [PRIMARY, SENSITIVITY]:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        handle.write(f"{digest}  {path.name}\n")
print(PRIMARY); print(SENSITIVITY); print(CHECKSUM)
