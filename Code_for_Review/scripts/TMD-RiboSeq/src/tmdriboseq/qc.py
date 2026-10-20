"""Replicate-level Ribo-seq QC and independent P-site calibration."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Mapping

import pysam

from .coordinates import TranscriptModel, load_gtf_transcripts
from .errors import ConfigurationError
from .profiles import SampleRecord, load_samples
from .psite import ReadFilter, psite_reference_position0, read_is_eligible


AMBIGUOUS_FRAME = -1


def build_cds_frame_lookup(
    transcripts: Mapping[str, TranscriptModel],
) -> tuple[dict[tuple[str, str, int], int], int]:
    """Build a strand-aware CDS frame lookup, excluding conflicting isoforms."""

    lookup: dict[tuple[str, str, int], int] = {}
    conflicts = 0
    for model in transcripts.values():
        for index, genomic_position in enumerate(model.genomic_positions0):
            key = (model.contig, model.strand, genomic_position)
            frame = index % 3
            previous = lookup.get(key)
            if previous is None:
                lookup[key] = frame
            elif previous != frame and previous != AMBIGUOUS_FRAME:
                lookup[key] = AMBIGUOUS_FRAME
                conflicts += 1
    return lookup, conflicts


def build_cds_position_lookup(
    transcripts: Mapping[str, TranscriptModel],
) -> tuple[dict[tuple[str, str, int], int], int]:
    """Map genomic CDS loci to 0-based CDS positions, masking isoform conflicts."""

    lookup: dict[tuple[str, str, int], int] = {}
    conflicts = 0
    for model in transcripts.values():
        for index, genomic_position in enumerate(model.genomic_positions0):
            key = (model.contig, model.strand, genomic_position)
            previous = lookup.get(key)
            if previous is None:
                lookup[key] = index
            elif previous != index and previous != AMBIGUOUS_FRAME:
                lookup[key] = AMBIGUOUS_FRAME
                conflicts += 1
    return lookup, conflicts


def calibrate_sample(
    sample: SampleRecord,
    frame_lookup: Mapping[tuple[str, str, int], int],
    read_lengths: list[int],
    candidate_offsets: list[int],
    policy: ReadFilter,
    max_alignments: int | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    lengths: Counter[int] = Counter()
    frames: dict[tuple[int, int], Counter[int]] = defaultdict(Counter)
    metrics: Counter[str] = Counter()
    with pysam.AlignmentFile(str(sample.bam), "rb") as bam:
        metrics["bam_mapped_index_count"] = int(bam.mapped)
        metrics["bam_unmapped_index_count"] = int(bam.unmapped)
        for read in bam.fetch(until_eof=True):
            metrics["alignments_seen"] += 1
            if max_alignments is not None and metrics["alignments_seen"] > max_alignments:
                break
            if not read_is_eligible(read, policy):
                metrics["alignments_filtered"] += 1
                continue
            length = read.query_length
            if length is None:
                continue
            length = int(length)
            lengths[length] += 1
            metrics["eligible_alignments"] += 1
            if length not in read_lengths:
                continue
            strand = "-" if read.is_reverse else "+"
            for offset in candidate_offsets:
                position = psite_reference_position0(read, offset)
                if position is None:
                    continue
                frame = frame_lookup.get((read.reference_name, strand, position))
                if frame is None or frame == AMBIGUOUS_FRAME:
                    continue
                frames[(length, offset)][frame] += 1

    length_rows = [
        {
            "sample_id": sample.sample_id,
            "read_length": length,
            "eligible_read_count": count,
            "fraction_of_eligible_alignments": (
                count / metrics["eligible_alignments"]
                if metrics["eligible_alignments"]
                else 0.0
            ),
        }
        for length, count in sorted(lengths.items())
    ]
    periodicity_rows: list[dict[str, Any]] = []
    for length in read_lengths:
        for offset in candidate_offsets:
            counts = frames[(length, offset)]
            total = sum(counts.values())
            periodicity_rows.append(
                {
                    "sample_id": sample.sample_id,
                    "read_length": length,
                    "offset": offset,
                    "frame0": counts[0],
                    "frame1": counts[1],
                    "frame2": counts[2],
                    "cds_psites": total,
                    "frame0_fraction": counts[0] / total if total else None,
                }
            )
    return length_rows, periodicity_rows, {"sample_id": sample.sample_id, **metrics}


def recommend_offsets(
    periodicity_rows: list[dict[str, Any]],
    minimum_reads: int,
    warning_frame_fraction: float,
    canonical_by_phase: Mapping[int, int],
) -> list[dict[str, Any]]:
    by_sample_length_phase: dict[tuple[str, int, int], list[dict[str, Any]]] = defaultdict(list)
    for row in periodicity_rows:
        by_sample_length_phase[
            (str(row["sample_id"]), int(row["read_length"]), int(row["offset"]) % 3)
        ].append(row)

    by_sample_length: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for (sample, length, phase), rows in by_sample_length_phase.items():
        valid = [row for row in rows if row["frame0_fraction"] is not None]
        if not valid:
            continue
        by_sample_length[(sample, length)].append(
            {
                "phase": phase,
                "mean_frame0_fraction": sum(float(row["frame0_fraction"]) for row in valid)
                / len(valid),
                "maximum_cds_psites": max(int(row["cds_psites"]) for row in valid),
            }
        )

    recommendations: list[dict[str, Any]] = []
    for (sample, length), candidates in sorted(by_sample_length.items()):
        best = max(candidates, key=lambda row: row["mean_frame0_fraction"])
        reads = int(best["maximum_cds_psites"])
        fraction = float(best["mean_frame0_fraction"])
        status = "PASS"
        reasons: list[str] = []
        if reads < minimum_reads:
            status = "WARN"
            reasons.append(f"CDS P-sites {reads} < {minimum_reads}")
        if fraction < warning_frame_fraction:
            status = "WARN"
            reasons.append(
                f"frame-0 fraction {fraction:.3f} < {warning_frame_fraction:.3f}"
            )
        phase = int(best["phase"])
        recommendations.append(
            {
                "sample_id": sample,
                "read_length": length,
                "offset_phase": phase,
                "offset": int(canonical_by_phase[phase]),
                "frame0_fraction": fraction,
                "cds_psites": reads,
                "qc_status": status,
                "qc_reason": "; ".join(reasons),
            }
        )
    return recommendations


def cds_start_metagene_for_sample(
    sample: SampleRecord,
    cds_position_lookup: Mapping[tuple[str, str, int], int],
    offsets: Mapping[int, int],
    policy: ReadFilter,
    start_nt: int,
    end_nt: int,
    max_alignments: int | None,
) -> list[dict[str, Any]]:
    """Count calibrated P-sites by CDS-relative position near translation start."""

    if start_nt < 0 or end_nt < start_nt:
        raise ConfigurationError(
            "The validated CDS-start metagene range must satisfy 0 <= start <= end."
        )
    counts: Counter[int] = Counter()
    alignments_seen = 0
    with pysam.AlignmentFile(str(sample.bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            alignments_seen += 1
            if max_alignments is not None and alignments_seen > max_alignments:
                break
            if not read_is_eligible(read, policy):
                continue
            length = read.query_length
            if length is None or int(length) not in offsets:
                continue
            position = psite_reference_position0(read, offsets[int(length)])
            if position is None:
                continue
            strand = "-" if read.is_reverse else "+"
            cds_position = cds_position_lookup.get(
                (read.reference_name, strand, position)
            )
            if (
                cds_position is None
                or cds_position == AMBIGUOUS_FRAME
                or cds_position < start_nt
                or cds_position > end_nt
            ):
                continue
            counts[int(cds_position)] += 1
    total = sum(counts.values())
    return [
        {
            "sample_id": sample.sample_id,
            "cds_position_nt_0based": position,
            "psite_count": counts[position],
            "fraction_within_reported_range": (
                counts[position] / total if total else 0.0
            ),
        }
        for position in range(start_nt, end_nt + 1)
    ]


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Cannot write empty QC table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_calibration_stage(config: dict[str, Any]) -> dict[str, Any]:
    psite = config["psite"]
    read_lengths = [int(value) for value in config["reads"]["allowed_lengths"]]
    candidate_offsets = [int(value) for value in psite.get("candidate_offsets", range(10, 18))]
    phase_config = psite.get("canonical_offset_by_phase", {0: 12, 1: 13, 2: 14})
    canonical_by_phase = {int(key): int(value) for key, value in phase_config.items()}
    if set(canonical_by_phase) != {0, 1, 2}:
        raise ConfigurationError("psite.canonical_offset_by_phase must define phases 0, 1, and 2.")

    policy = ReadFilter(
        minimum_mapping_quality=int(config["reads"].get("minimum_mapping_quality", 0)),
        include_secondary=bool(config["reads"].get("include_secondary", False)),
        include_supplementary=bool(config["reads"].get("include_supplementary", False)),
        include_duplicates=bool(config["reads"].get("include_duplicates", False)),
    )
    transcripts = load_gtf_transcripts(config["reference"]["annotation"])
    frame_lookup, frame_conflicts = build_cds_frame_lookup(transcripts)
    cds_position_lookup, position_conflicts = build_cds_position_lookup(transcripts)
    samples = load_samples(config)
    all_length_rows: list[dict[str, Any]] = []
    all_periodicity_rows: list[dict[str, Any]] = []
    all_mapping_rows: list[dict[str, Any]] = []
    max_alignments_raw = psite.get("maximum_alignments_per_sample")
    max_alignments = int(max_alignments_raw) if max_alignments_raw is not None else None
    for sample in samples:
        length_rows, periodicity_rows, mapping = calibrate_sample(
            sample,
            frame_lookup,
            read_lengths,
            candidate_offsets,
            policy,
            max_alignments,
        )
        all_length_rows.extend(length_rows)
        all_periodicity_rows.extend(periodicity_rows)
        all_mapping_rows.append(mapping)

    recommendations = recommend_offsets(
        all_periodicity_rows,
        minimum_reads=int(psite.get("minimum_reads_per_length", 1000)),
        warning_frame_fraction=float(psite.get("warning_frame0_fraction", 0.5)),
        canonical_by_phase=canonical_by_phase,
    )
    selected_by_sample: dict[str, dict[int, int]] = defaultdict(dict)
    for row in recommendations:
        selected_by_sample[str(row["sample_id"])][int(row["read_length"])] = int(
            row["offset"]
        )
    metagene_start = int(psite.get("start_codon_metagene_start_nt", 0))
    metagene_end = int(psite.get("start_codon_metagene_end_nt", 120))
    metagene_rows: list[dict[str, Any]] = []
    for sample in samples:
        metagene_rows.extend(
            cds_start_metagene_for_sample(
                sample,
                cds_position_lookup,
                selected_by_sample.get(sample.sample_id, {}),
                policy,
                metagene_start,
                metagene_end,
                max_alignments,
            )
        )
    output = Path(config["output"]["directory"]) / "qc"
    _write_tsv(output / "read_length_distribution.tsv", all_length_rows)
    _write_tsv(output / "periodicity_by_sample_length_offset.tsv", all_periodicity_rows)
    _write_tsv(output / "mapping_summary.tsv", all_mapping_rows)
    _write_tsv(output / "psite_offsets.tsv", recommendations)
    _write_tsv(output / "cds_start_metagene.tsv", metagene_rows)
    return {
        "samples": len(samples),
        "transcripts": len(transcripts),
        "ambiguous_frame_loci": frame_conflicts,
        "ambiguous_cds_position_loci": position_conflicts,
        "recommended_sample_length_offsets": len(recommendations),
        "warning_offsets": sum(row["qc_status"] == "WARN" for row in recommendations),
        "offset_table": str(output / "psite_offsets.tsv"),
        "cds_start_metagene": str(output / "cds_start_metagene.tsv"),
    }
