#!/usr/bin/env python3
"""Apply frozen P-site offsets and prespecified held-out Ubc6 eligibility."""

import csv
import sys
from collections import Counter
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "frozen/v1.0.0_owner_review_candidate/wheel_extracted"))
from tmdriboseq.coordinates import load_gtf_transcripts  # noqa: E402
from tmdriboseq.psite import psite_reference_position0  # noqa: E402
from tmdriboseq.qc import AMBIGUOUS_FRAME, build_cds_frame_lookup  # noqa: E402


META = ROOT / "metadata/gsm_manifest.tsv"
MAPPING = ROOT / "qc/mapping_metrics_by_gsm.tsv"
BAMS = ROOT / "derived/stage1b_calibration_bams"
OFFSETS = ROOT / "frozen/psite/psite_series_consensus_frozen.tsv"
GTF = ROOT / "references/genome/Saccharomyces_cerevisiae.R64-1-1.115.gtf"
QC_OUT = ROOT / "qc/stage1b_qc_by_gsm.tsv"
ELIG_OUT = ROOT / "qc/ubc6_validation_eligibility_by_gsm.tsv"
PAIR_OUT = ROOT / "qc/ubc6_validation_eligibility_by_pair.tsv"
THRESHOLDS = {
    "unique": 5_000_000, "usable": 0.75, "cds": 0.70,
    "frame0": 0.45, "assignable": 0.99, "contamination": 0.80,
}


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


meta = {r["gsm"]: r for r in csv.DictReader(META.open(), delimiter="\t")}
mapping = {r["gsm"]: r for r in csv.DictReader(MAPPING.open(), delimiter="\t")}
offsets = {
    (r["series_id"], int(r["read_length"])): int(r["selected_offset_5prime_0based"])
    for r in csv.DictReader(OFFSETS.open(), delimiter="\t")
}
frame_lookup, conflicts = build_cds_frame_lookup(load_gtf_transcripts(GTF))
qc_rows = []
for sample, m in sorted(meta.items()):
    counts = Counter(); frames = Counter()
    with pysam.AlignmentFile(str(BAMS / f"{sample}.calibration.bam"), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            counts["eligible"] += 1
            length = int(read.query_length or 0)
            if length not in range(26, 33):
                continue
            counts["usable"] += 1
            position = psite_reference_position0(read, offsets[(m["series_id"], length)])
            if position is None:
                continue
            counts["assignable"] += 1
            strand = "-" if read.is_reverse else "+"
            frame = frame_lookup.get((read.reference_name, strand, position))
            if frame is not None and frame != AMBIGUOUS_FRAME:
                counts["cds"] += 1; frames[frame] += 1
    mp = mapping[sample]
    trimmed = int(mp["trimmed_reads"]); contaminant = int(mp["contaminant_reads"])
    qc_rows.append({
        "sample_id": sample, "series_id": m["series_id"], "time_min": m["auxin_time_min"],
        "fraction": m["fraction"], "pair_id": m["pair_id"],
        "unique_mapq20_reads": int(mp["unique_mapq20"]),
        "contaminant_fraction": contaminant / trimmed,
        "calibration_alignments": counts["eligible"], "usable_26_32": counts["usable"],
        "usable_26_32_fraction": counts["usable"] / counts["eligible"],
        "psite_assignable": counts["assignable"],
        "psite_assignability_fraction": counts["assignable"] / counts["usable"],
        "cds_assigned": counts["cds"],
        "cds_assignment_fraction": counts["cds"] / counts["usable"],
        "frame0": frames[0], "frame1": frames[1], "frame2": frames[2],
        "frame0_fraction": frames[0] / sum(frames.values()) if frames else 0.0,
        "ambiguous_frame_loci_masked": conflicts,
    })
write(QC_OUT, qc_rows)

eligibility = []
for row in qc_rows:
    if row["series_id"] != "Ubc6_S288C_CHX":
        continue
    tests = {
        "unique_pass": row["unique_mapq20_reads"] >= THRESHOLDS["unique"],
        "usable_pass": row["usable_26_32_fraction"] >= THRESHOLDS["usable"],
        "cds_pass": row["cds_assignment_fraction"] >= THRESHOLDS["cds"],
        "frame0_pass": row["frame0_fraction"] >= THRESHOLDS["frame0"],
        "assignability_pass": row["psite_assignability_fraction"] >= THRESHOLDS["assignable"],
        "contamination_pass": row["contaminant_fraction"] <= THRESHOLDS["contamination"],
    }
    eligibility.append({**row, **tests, "GSM_eligible": all(tests.values())})
write(ELIG_OUT, eligibility)

pair_rows = []
by_pair = {}
for row in eligibility:
    by_pair.setdefault(row["pair_id"], {})[row["fraction"]] = row
for pair_id, members in sorted(by_pair.items()):
    input_row = members.get("input"); pulldown_row = members.get("pulldown")
    complete = input_row is not None and pulldown_row is not None
    pair_pass = complete and input_row["GSM_eligible"] and pulldown_row["GSM_eligible"]
    pair_rows.append({
        "pair_id": pair_id, "time_min": input_row["time_min"] if input_row else pulldown_row["time_min"],
        "input_GSM": input_row["sample_id"] if input_row else "",
        "pulldown_GSM": pulldown_row["sample_id"] if pulldown_row else "",
        "pair_complete": complete,
        "input_eligible": input_row["GSM_eligible"] if input_row else False,
        "pulldown_eligible": pulldown_row["GSM_eligible"] if pulldown_row else False,
        "pair_eligible": pair_pass,
        "formal_0_to_30_validation_eligible": False,
    })
formal = {r["time_min"]: r["pair_eligible"] for r in pair_rows}.get("0", False) and {r["time_min"]: r["pair_eligible"] for r in pair_rows}.get("30", False)
for row in pair_rows:
    row["formal_0_to_30_validation_eligible"] = formal
write(PAIR_OUT, pair_rows)
print(QC_OUT); print(ELIG_OUT); print(PAIR_OUT); print("formal_0_to_30", formal)
