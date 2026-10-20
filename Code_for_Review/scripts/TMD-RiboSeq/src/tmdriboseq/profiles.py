"""TMD-centered target and background counting for canonical inputs."""

from __future__ import annotations

import csv
import gzip
import math
import os
import statistics
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import pysam

from .coordinates import TranscriptModel, load_gtf_transcripts
from .coverage import summarize_coverage, write_tsv
from .errors import ConfigurationError, ValidationFailed
from .psite import ReadFilter, psite_for_read, read_is_eligible
from .tmds import TMDRecord, load_canonical_tmds, load_paper_legacy_tmds


@dataclass(frozen=True)
class SampleRecord:
    sample_id: str
    bam: Path
    replicate: str
    role: str
    factors: dict[str, str]


def _table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def load_samples(config: dict[str, Any]) -> list[SampleRecord]:
    block = config["samples"]
    path = Path(block["metadata"])
    mapping = block.get("columns", {})
    factors = config["design"]["factors"]
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=_table_delimiter(path))
        samples: list[SampleRecord] = []
        for row in reader:
            bam_text = os.path.expandvars(row[mapping.get("bam", "bam")].strip())
            if "${" in bam_text:
                raise ConfigurationError(
                    f"Unresolved environment variable in BAM path for "
                    f"{row[mapping.get('sample_id', 'sample_id')]}"
                )
            bam = Path(bam_text).expanduser()
            if not bam.is_absolute():
                bam = path.parent / bam
            samples.append(
                SampleRecord(
                    sample_id=row[mapping.get("sample_id", "sample_id")].strip(),
                    bam=bam.resolve(),
                    replicate=row[mapping.get("replicate", "replicate")].strip(),
                    role=row[mapping.get("role", "role")].strip(),
                    factors={
                        factor: row[mapping.get(factor, factor)].strip() for factor in factors
                    },
                )
            )
    return samples


def load_offsets(path: str | Path) -> dict[str, dict[int, int]]:
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=_table_delimiter(path))
        fields = set(reader.fieldnames or [])
        sample_column = "sample_id" if "sample_id" in fields else "accession"
        offset_column = (
            "recommended_offset" if "recommended_offset" in fields else "offset"
        )
        required = {sample_column, "read_length", offset_column}
        if not required.issubset(fields):
            raise ValidationFailed(
                f"P-site offset table must contain {', '.join(sorted(required))}."
            )
        offsets: dict[str, dict[int, int]] = defaultdict(dict)
        for row in reader:
            offsets[row[sample_column]][int(row["read_length"])] = int(row[offset_column])
    return dict(offsets)


def count_psite_loci(
    bam_path: str | Path,
    offsets: Mapping[int, int],
    policy: ReadFilter,
) -> tuple[Counter[tuple[str, str, int]], dict[str, int]]:
    """Scan a BAM once and count eligible P sites by contig, strand, and position."""

    counts: Counter[tuple[str, str, int]] = Counter()
    metrics: Counter[str] = Counter()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            metrics["alignments_seen"] += 1
            if not read_is_eligible(read, policy):
                metrics["alignments_filtered"] += 1
                continue
            metrics["eligible_alignments"] += 1
            psite = psite_for_read(read, offsets)
            if psite is None:
                metrics["unsupported_length_or_offset"] += 1
                continue
            strand = "-" if read.is_reverse else "+"
            counts[(read.reference_name, strand, psite)] += 1
            metrics["assigned_psites"] += 1
    return counts, dict(metrics)


def scan_centers(config: dict[str, Any]) -> list[int]:
    profiles = config["profiles"]
    start = int(profiles["scan_start_nt"])
    end = int(profiles["scan_end_nt"])
    step = int(profiles["step_nt"])
    if step <= 0 or end < start:
        raise ConfigurationError("Invalid positional scan range or step.")
    centers = list(range(start, end + 1, step))
    primary = int(profiles["primary_center_nt"])
    if primary not in centers:
        centers.append(primary)
        centers.sort()
    return centers


def _half_width(width: int) -> int:
    if width <= 0 or width % 2 == 0:
        raise ConfigurationError("profiles.window_width_nt must be a positive odd integer.")
    return width // 2


def _overlaps_other_tmd(
    tmd: TMDRecord,
    start_cds: int,
    end_cds: int,
    by_transcript: Mapping[str, list[TMDRecord]],
) -> bool:
    for other in by_transcript[tmd.transcript_id]:
        if other.tmd_id == tmd.tmd_id:
            continue
        if max(start_cds, other.start_cds_nt) <= min(end_cds, other.end_cds_nt):
            return True
    return False


def build_background_masks(
    transcripts: Mapping[str, TranscriptModel],
    tmds: Iterable[TMDRecord],
    config: dict[str, Any],
) -> dict[str, tuple[bool, ...]]:
    """Build transcript-level gene-wide non-TMD background masks."""

    strategy = config["background"]["strategy"]
    if strategy != "gene_wide_non_tmd_cds":
        raise ConfigurationError(
            "The canonical 0.1 implementation currently supports "
            "background.strategy=gene_wide_non_tmd_cds only."
        )
    profiles = config["profiles"]
    half = _half_width(int(profiles["window_width_nt"]))
    exclusion = config["background"].get("exclusion")
    if isinstance(exclusion, dict):
        scan_left = -int(exclusion["left_of_tmd_start_nt"])
        scan_right = int(exclusion["right_of_tmd_start_nt"])
    else:
        scan_left = int(profiles["scan_start_nt"]) - half
        scan_right = int(profiles["scan_end_nt"]) + half
    by_transcript: dict[str, list[TMDRecord]] = defaultdict(list)
    for tmd in tmds:
        by_transcript[tmd.transcript_id].append(tmd)

    masks: dict[str, tuple[bool, ...]] = {}
    for transcript_id, transcript_tmds in by_transcript.items():
        model = transcripts[transcript_id]
        eligible = [True] * model.cds_length_nt
        for tmd in transcript_tmds:
            left = max(1, tmd.start_cds_nt + scan_left)
            right = min(
                model.cds_length_nt,
                max(tmd.end_cds_nt, tmd.start_cds_nt + scan_right),
            )
            for index in range(left - 1, right):
                eligible[index] = False
        masks[transcript_id] = tuple(eligible)
    return masks


def _count_cds_interval(
    model: TranscriptModel,
    counts: Mapping[tuple[str, str, int], int],
    start_cds: int,
    end_cds: int,
) -> int:
    return sum(
        counts.get((model.contig, model.strand, position), 0)
        for position in model.genomic_interval0(start_cds, end_cds)
    )


def _background_excluding_target(
    model: TranscriptModel,
    counts: Mapping[tuple[str, str, int], int],
    mask: tuple[bool, ...],
    background_count: int,
    background_length: int,
    target_start_cds: int,
    target_end_cds: int,
) -> tuple[int, int]:
    """Return a background that is structurally disjoint from one target window."""

    overlap_start = max(1, target_start_cds)
    overlap_end = min(model.cds_length_nt, target_end_cds)
    if overlap_start > overlap_end:
        return background_count, background_length
    included_target_positions = [
        cds_position
        for cds_position in range(overlap_start, overlap_end + 1)
        if mask[cds_position - 1]
    ]
    if not included_target_positions:
        return background_count, background_length
    overlapping_count = sum(
        counts.get((model.contig, model.strand, genomic), 0)
        for cds_position in included_target_positions
        for genomic in model.genomic_interval0(cds_position, cds_position)
    )
    return (
        background_count - overlapping_count,
        background_length - len(included_target_positions),
    )


def count_profiles_for_sample(
    sample: SampleRecord,
    offsets: Mapping[int, int],
    transcripts: Mapping[str, TranscriptModel],
    tmds: list[TMDRecord],
    config: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    reads = config.get("reads", {})
    policy = ReadFilter(
        minimum_mapping_quality=int(reads.get("minimum_mapping_quality", 0)),
        include_secondary=bool(reads.get("include_secondary", False)),
        include_supplementary=bool(reads.get("include_supplementary", False)),
        include_duplicates=bool(reads.get("include_duplicates", False)),
    )
    psite_counts, metrics = count_psite_loci(sample.bam, offsets, policy)
    centers = scan_centers(config)
    width = int(config["profiles"]["window_width_nt"])
    half = _half_width(width)
    minimum_background = int(config["coverage"]["minimum_background_count"])
    masks = build_background_masks(transcripts, tmds, config)
    by_transcript: dict[str, list[TMDRecord]] = defaultdict(list)
    for tmd in tmds:
        by_transcript[tmd.transcript_id].append(tmd)

    background_counts: dict[str, tuple[int, int]] = {}
    for transcript_id, mask in masks.items():
        model = transcripts[transcript_id]
        background_length = sum(mask)
        background_count = sum(
            psite_counts.get((model.contig, model.strand, genomic), 0)
            for genomic, include in zip(model.genomic_positions0, mask)
            if include
        )
        background_counts[transcript_id] = (background_count, background_length)

    rows: list[dict[str, Any]] = []
    for tmd in tmds:
        model = transcripts[tmd.transcript_id]
        mask = masks[tmd.transcript_id]
        transcript_background_count, transcript_background_length = background_counts[
            tmd.transcript_id
        ]
        for center in centers:
            target_start = tmd.start_cds_nt + center - half
            target_end = tmd.start_cds_nt + center + half
            background_count, background_length = _background_excluding_target(
                model,
                psite_counts,
                mask,
                transcript_background_count,
                transcript_background_length,
                target_start,
                target_end,
            )
            target_within_cds = (
                target_start >= 1 and target_end <= model.cds_length_nt
            )
            target_overlaps_other = False
            target_count: int | None = None
            if target_within_cds:
                target_overlaps_other = _overlaps_other_tmd(
                    tmd, target_start, target_end, by_transcript
                )
                target_count = _count_cds_interval(
                    model, psite_counts, target_start, target_end
                )
            structurally_valid = (
                target_within_cds
                and not target_overlaps_other
                and background_length > 0
            )
            coverage_pass = structurally_valid and background_count >= minimum_background
            signal = None
            if structurally_valid and target_count is not None:
                signal = math.log2(
                    ((target_count + 0.5) / width)
                    / ((background_count + 0.5) / background_length)
                )
            row: dict[str, Any] = {
                "sample_id": sample.sample_id,
                "replicate": sample.replicate,
                "role": sample.role,
                **sample.factors,
                "gene_id": tmd.gene_id,
                "transcript_id": tmd.transcript_id,
                "tmd_id": tmd.tmd_id,
                "center_nt": center,
                "target_start_cds_nt": target_start,
                "target_end_cds_nt": target_end,
                "target_length_nt": width,
                "target_count": target_count,
                "background_count": background_count,
                "background_length_nt": background_length,
                "local_log2_occupancy": signal,
                "target_within_cds": target_within_cds,
                "target_overlaps_other_tmd": target_overlaps_other,
                "structurally_valid": structurally_valid,
                "coverage_pass": coverage_pass,
            }
            rows.append(row)
    return rows, metrics


def write_tsv_gz(path: str | Path, rows: list[dict[str, Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("Cannot write an empty result table.")
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sample_signal_correlations(
    rows: Iterable[Mapping[str, Any]], primary_center_nt: int
) -> list[dict[str, Any]]:
    """Return pairwise Pearson correlations over common primary-window TMDs."""

    by_sample: dict[str, dict[str, float]] = defaultdict(dict)
    for row in rows:
        if (
            int(row["center_nt"]) != primary_center_nt
            or not bool(row["coverage_pass"])
            or row["local_log2_occupancy"] is None
        ):
            continue
        by_sample[str(row["sample_id"])][str(row["tmd_id"])] = float(
            row["local_log2_occupancy"]
        )
    output: list[dict[str, Any]] = []
    sample_ids = sorted(by_sample)
    for left in sample_ids:
        for right in sample_ids:
            common = sorted(set(by_sample[left]) & set(by_sample[right]))
            left_values = [by_sample[left][key] for key in common]
            right_values = [by_sample[right][key] for key in common]
            correlation: float | None = None
            if len(common) >= 2:
                try:
                    correlation = float(statistics.correlation(left_values, right_values))
                except statistics.StatisticsError:
                    correlation = None
            output.append(
                {
                    "sample_a": left,
                    "sample_b": right,
                    "n_common_tmds": len(common),
                    "pearson_r": correlation,
                }
            )
    return output


def baseline_equal_weight_metagene(
    rows: Iterable[Mapping[str, Any]], baseline_role: str, baseline_filters: Mapping[str, str]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Summarize per-TMD-normalized signals with equal TMD and sample weights."""

    values: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in rows:
        if (
            str(row["role"]) != baseline_role
            or any(str(row.get(name, "")) != str(value) for name, value in baseline_filters.items())
            or not bool(row["coverage_pass"])
            or row["local_log2_occupancy"] is None
        ):
            continue
        values[(str(row["sample_id"]), int(row["center_nt"]))].append(
            float(row["local_log2_occupancy"])
        )
    sample_rows: list[dict[str, Any]] = []
    for (sample_id, center), scores in sorted(values.items()):
        sample_rows.append(
            {
                "sample_id": sample_id,
                "center_nt": center,
                "n_equal_weight_tmds": len(scores),
                "mean_normalized_signal": statistics.mean(scores),
                "median_normalized_signal": statistics.median(scores),
            }
        )
    by_center: dict[int, list[float]] = defaultdict(list)
    by_center_n: dict[int, list[int]] = defaultdict(list)
    for row in sample_rows:
        center = int(row["center_nt"])
        by_center[center].append(float(row["mean_normalized_signal"]))
        by_center_n[center].append(int(row["n_equal_weight_tmds"]))
    summary_rows: list[dict[str, Any]] = []
    for center, means in sorted(by_center.items()):
        summary_rows.append(
            {
                "center_nt": center,
                "n_baseline_samples": len(means),
                "minimum_tmds_per_sample": min(by_center_n[center]),
                "mean_of_sample_means": statistics.mean(means),
                "sd_across_sample_means": (
                    statistics.stdev(means) if len(means) > 1 else 0.0
                ),
            }
        )
    return sample_rows, summary_rows


def run_profile_stage(config: dict[str, Any]) -> dict[str, Any]:
    """Run canonical profile counting and write machine-readable audit outputs."""

    transcripts = load_gtf_transcripts(config["reference"]["annotation"])
    coordinate_system = config["tmd_annotation"]["coordinate_system"]
    if coordinate_system == "cds_nt_1based_inclusive":
        tmds = load_canonical_tmds(
            config["tmd_annotation"]["path"], config["tmd_annotation"].get("columns")
        )
    elif (
        coordinate_system == "paper_legacy_genomic"
        and config["project"].get("mode") == "paper_reproduction"
    ):
        tmds = load_paper_legacy_tmds(
            config["tmd_annotation"]["path"], config["reference"]["annotation"]
        )
    else:
        raise ConfigurationError(
            f"Unsupported TMD coordinate system for profile stage: {coordinate_system}"
        )
    for tmd in tmds:
        model = transcripts.get(tmd.transcript_id)
        if model is None:
            raise ValidationFailed(f"TMD {tmd.tmd_id} transcript is absent from GTF.")
        if tmd.end_cds_nt > model.cds_length_nt:
            raise ValidationFailed(f"TMD {tmd.tmd_id} extends outside its CDS.")

    psite = config["psite"]
    if psite.get("mode") != "fixed_table" or not psite.get("offset_table"):
        raise ConfigurationError(
            "Profile stage currently requires psite.mode=fixed_table. "
            "Calibration will be migrated as the next scientific stage."
        )
    offsets_by_sample = load_offsets(psite["offset_table"])
    samples = load_samples(config)
    output = Path(config["output"]["directory"])
    profile_dir = output / "profiles"
    qc_rows: list[dict[str, Any]] = []
    all_rows: list[dict[str, Any]] = []
    for sample in samples:
        if sample.sample_id not in offsets_by_sample:
            raise ValidationFailed(
                f"No P-site offsets found for sample {sample.sample_id}."
            )
        rows, metrics = count_profiles_for_sample(
            sample,
            offsets_by_sample[sample.sample_id],
            transcripts,
            tmds,
            config,
        )
        write_tsv_gz(profile_dir / f"{sample.sample_id}.tmd_profiles.tsv.gz", rows)
        all_rows.extend(rows)
        qc_rows.append({"sample_id": sample.sample_id, **metrics})

    write_tsv_gz(output / "tmd_profiles.tsv.gz", all_rows)
    qc_path = output / "qc" / "psite_assignment_summary.tsv"
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    with qc_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qc_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_rows)
    coverage_flow, coverage_membership = summarize_coverage(
        all_rows,
        primary_center_nt=int(config["profiles"]["primary_center_nt"]),
        model_role=str(config["samples"].get("model_role", "model")),
        factor_names=[str(value) for value in config["design"]["factors"]],
    )
    coverage_dir = output / "coverage"
    coverage_flow_path = coverage_dir / "coverage_flow.tsv"
    coverage_membership_path = coverage_dir / "coverage_membership.tsv"
    write_tsv(coverage_flow_path, coverage_flow)
    write_tsv(coverage_membership_path, coverage_membership)
    correlation_path = output / "qc" / "sample_primary_signal_correlations.tsv"
    write_tsv(
        correlation_path,
        sample_signal_correlations(
            all_rows, int(config["profiles"]["primary_center_nt"])
        ),
    )
    metagene_rows, metagene_summary = baseline_equal_weight_metagene(
        all_rows,
        str(config["samples"].get("baseline_role", "baseline")),
        {
            str(key): str(value)
            for key, value in config["samples"].get("baseline_filters", {}).items()
        },
    )
    metagene_path = output / "qc" / "baseline_equal_weight_metagene.tsv"
    metagene_summary_path = (
        output / "qc" / "baseline_equal_weight_metagene_summary.tsv"
    )
    write_tsv(metagene_path, metagene_rows)
    write_tsv(metagene_summary_path, metagene_summary)
    return {
        "samples": len(samples),
        "tmds": len(tmds),
        "centers": scan_centers(config),
        "profile_rows": len(all_rows),
        "profile_table": str(output / "tmd_profiles.tsv.gz"),
        "assignment_summary": str(qc_path),
        "coverage_flow": str(coverage_flow_path),
        "coverage_membership": str(coverage_membership_path),
        "sample_correlations": str(correlation_path),
        "baseline_metagene": str(metagene_path),
        "baseline_metagene_summary": str(metagene_summary_path),
        "strict_primary_eligible_tmds": int(
            coverage_flow[-1]["eligible_tmds"]
        ),
    }
