"""Input validation for supported v1 configurations."""

from __future__ import annotations

import csv
import json
import os
import re
from dataclasses import asdict, dataclass, field
from itertools import product
from pathlib import Path
from typing import Any, Iterable

import pysam

from .config import unresolved_environment_tokens


@dataclass(frozen=True)
class Diagnostic:
    severity: str
    code: str
    message: str
    context: dict[str, Any] = field(default_factory=dict)


@dataclass
class ValidationReport:
    diagnostics: list[Diagnostic] = field(default_factory=list)
    summary: dict[str, Any] = field(default_factory=dict)

    def add(
        self,
        severity: str,
        code: str,
        message: str,
        **context: Any,
    ) -> None:
        self.diagnostics.append(Diagnostic(severity, code, message, context))

    @property
    def errors(self) -> list[Diagnostic]:
        return [item for item in self.diagnostics if item.severity == "ERROR"]

    @property
    def warnings(self) -> list[Diagnostic]:
        return [item for item in self.diagnostics if item.severity == "WARN"]

    @property
    def passed(self) -> bool:
        return not self.errors

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "error_count": len(self.errors),
            "warning_count": len(self.warnings),
            "summary": self.summary,
            "diagnostics": [asdict(item) for item in self.diagnostics],
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2)


def _require_mapping(
    config: dict[str, Any], section: str, report: ValidationReport
) -> dict[str, Any]:
    value = config.get(section)
    if not isinstance(value, dict):
        report.add("ERROR", "CONFIG_SECTION", f"Missing mapping: {section}")
        return {}
    return value


def _check_file(
    value: Any,
    label: str,
    report: ValidationReport,
    required: bool = True,
) -> Path | None:
    if value in (None, ""):
        if required:
            report.add("ERROR", "MISSING_PATH", f"Missing required path: {label}")
        return None
    path = Path(str(value))
    if not path.is_file():
        report.add("ERROR", "FILE_NOT_FOUND", f"File not found for {label}: {path}")
        return None
    return path


def _delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def _read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=_delimiter(path))
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def _parse_attributes(text: str, annotation_format: str) -> dict[str, str]:
    if annotation_format == "gff3":
        attrs: dict[str, str] = {}
        for item in text.strip().split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attrs[key.strip()] = value.strip()
        return attrs
    return {
        key: value
        for key, value in re.findall(r'(\S+)\s+"([^"]*)"', text)
    }


def _annotation_index(
    path: Path, annotation_format: str, report: ValidationReport
) -> tuple[set[str], dict[str, int], dict[str, str]]:
    contigs: set[str] = set()
    cds_lengths: dict[str, int] = {}
    transcript_strands: dict[str, str] = {}
    malformed = 0
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                malformed += 1
                continue
            contig, _, feature, start, end, _, strand, _, attributes = fields
            contigs.add(contig)
            if feature != "CDS":
                continue
            attrs = _parse_attributes(attributes, annotation_format)
            transcript_id = (
                attrs.get("transcript_id")
                or attrs.get("Parent")
                or attrs.get("protein_id")
            )
            if not transcript_id:
                continue
            transcript_id = transcript_id.split(",")[0]
            try:
                length = int(end) - int(start) + 1
            except ValueError:
                malformed += 1
                continue
            cds_lengths[transcript_id] = cds_lengths.get(transcript_id, 0) + length
            previous = transcript_strands.setdefault(transcript_id, strand)
            if previous != strand:
                report.add(
                    "ERROR",
                    "TRANSCRIPT_STRAND_CONFLICT",
                    f"Transcript {transcript_id} has CDS features on both strands.",
                    line=line_number,
                )
    if malformed:
        report.add(
            "WARN",
            "ANNOTATION_MALFORMED_LINES",
            f"Ignored {malformed} malformed annotation line(s).",
        )
    if not cds_lengths:
        report.add(
            "ERROR",
            "NO_TRANSCRIPT_CDS",
            "No transcript-linked CDS features were parsed from the annotation.",
        )
    return contigs, cds_lengths, transcript_strands


def _resolve_metadata_path(metadata_path: Path, value: str) -> Path:
    candidate = Path(os.path.expandvars(value)).expanduser()
    if not candidate.is_absolute():
        candidate = metadata_path.parent / candidate
    return candidate.resolve(strict=False)


def _validate_metadata(
    block: dict[str, Any],
    design: dict[str, Any],
    report: ValidationReport,
) -> tuple[list[dict[str, str]], set[str]]:
    path = _check_file(block.get("metadata"), "samples.metadata", report)
    if path is None:
        return [], set()

    columns, rows = _read_table(path)
    mapping = block.get("columns", {})
    if not isinstance(mapping, dict):
        report.add("ERROR", "SAMPLE_COLUMN_MAPPING", "samples.columns must be a mapping.")
        return [], set()

    logical_required = ["sample_id", "bam", "replicate", "role"]
    factors = design.get("factors", [])
    if not isinstance(factors, list):
        factors = []
    logical_required.extend(str(item) for item in factors)

    physical_required: list[str] = []
    for logical in logical_required:
        physical = mapping.get(logical, logical)
        physical_required.append(physical)
        if physical not in columns:
            report.add(
                "ERROR",
                "SAMPLE_COLUMN_MISSING",
                f"Sample metadata is missing column {physical!r} for {logical!r}.",
            )

    if any(item not in columns for item in physical_required):
        return rows, set()

    sample_col = mapping.get("sample_id", "sample_id")
    bam_col = mapping.get("bam", "bam")
    role_col = mapping.get("role", "role")
    replicate_col = mapping.get("replicate", "replicate")
    model_role = str(block.get("model_role", "model"))

    seen: set[str] = set()
    duplicate_ids: set[str] = set()
    bam_contigs: set[str] = set()
    cell_counts: dict[tuple[str, ...], int] = {}
    cell_replicates: dict[tuple[str, ...], set[str]] = {}
    factor_levels: dict[str, set[str]] = {str(factor): set() for factor in factors}
    bam_count = 0
    for row_number, row in enumerate(rows, start=2):
        sample_id = row.get(sample_col, "").strip()
        if not sample_id:
            report.add(
                "ERROR", "EMPTY_SAMPLE_ID", "Sample ID is empty.", row=row_number
            )
            continue
        if sample_id in seen:
            duplicate_ids.add(sample_id)
        seen.add(sample_id)

        replicate = row.get(replicate_col, "").strip()
        if not replicate:
            report.add(
                "ERROR",
                "EMPTY_REPLICATE",
                f"Sample {sample_id} has no replicate value.",
            )

        role = row.get(role_col, "").strip()
        if role == model_role:
            cell: list[str] = []
            for factor in factors:
                physical = mapping.get(factor, factor)
                value = row.get(physical, "").strip()
                if not value:
                    report.add(
                        "ERROR",
                        "EMPTY_FACTOR",
                        f"Model sample {sample_id} has no value for factor {factor}.",
                    )
                else:
                    factor_levels[str(factor)].add(value)
                cell.append(value)
            cell_key = tuple(cell)
            cell_counts[cell_key] = cell_counts.get(cell_key, 0) + 1
            seen_replicates = cell_replicates.setdefault(cell_key, set())
            if replicate and replicate in seen_replicates:
                report.add(
                    "ERROR",
                    "DUPLICATE_REPLICATE",
                    f"Design cell {cell_key} repeats replicate label {replicate!r}.",
                )
            seen_replicates.add(replicate)

        bam_value = row.get(bam_col, "").strip()
        if not bam_value:
            report.add(
                "ERROR", "EMPTY_BAM", f"Sample {sample_id} has no BAM path."
            )
            continue
        expanded_bam_value = os.path.expandvars(bam_value)
        if "${" in expanded_bam_value:
            report.add(
                "ERROR",
                "UNRESOLVED_ENVIRONMENT",
                f"Sample {sample_id} has an unresolved environment variable in its BAM path.",
                value=bam_value,
            )
            continue

        bam_path = _resolve_metadata_path(path, expanded_bam_value)
        if not bam_path.is_file():
            report.add(
                "ERROR",
                "BAM_NOT_FOUND",
                f"BAM file not found for sample {sample_id}: {bam_path}",
            )
            continue

        bai_candidates = [Path(f"{bam_path}.bai"), bam_path.with_suffix(".bai")]
        if not any(candidate.is_file() for candidate in bai_candidates):
            report.add(
                "ERROR",
                "BAI_NOT_FOUND",
                f"BAM index not found for sample {sample_id}: {bam_path}",
            )
            continue

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                bam_contigs.update(bam.references)
        except (ValueError, OSError) as exc:
            report.add(
                "ERROR",
                "BAM_OPEN_FAILED",
                f"Could not open BAM for sample {sample_id}: {exc}",
            )
            continue
        bam_count += 1

    for sample_id in sorted(duplicate_ids):
        report.add(
            "ERROR", "DUPLICATE_SAMPLE_ID", f"Duplicate sample ID: {sample_id}"
        )

    minimum_replicates = int(design.get("minimum_replicates_per_cell", 2))
    for cell, count in sorted(cell_counts.items()):
        if count < minimum_replicates:
            report.add(
                "ERROR",
                "INSUFFICIENT_REPLICATES",
                f"Design cell {cell} has {count} replicate(s); "
                f"at least {minimum_replicates} are required.",
            )

    if design.get("type") == "factorial_2x2" and len(factors) == 2:
        for factor in factors:
            levels = factor_levels.get(str(factor), set())
            if len(levels) != 2:
                report.add(
                    "ERROR",
                    "FACTOR_LEVELS",
                    f"Factor {factor!r} must have exactly two non-empty levels; "
                    f"found {sorted(levels)}.",
                )
        if all(len(factor_levels.get(str(factor), set())) == 2 for factor in factors):
            expected_cells = set(
                product(*(sorted(factor_levels[str(factor)]) for factor in factors))
            )
            missing_cells = sorted(expected_cells - set(cell_counts))
            if missing_cells:
                report.add(
                    "ERROR",
                    "MISSING_DESIGN_CELLS",
                    f"The 2 × 2 design is missing cell(s): {missing_cells}.",
                )
        references = design.get("reference_levels", {})
        if not isinstance(references, dict):
            report.add(
                "ERROR",
                "REFERENCE_LEVELS",
                "design.reference_levels must be a mapping.",
            )
        else:
            for factor in factors:
                reference = references.get(factor)
                if reference not in factor_levels.get(str(factor), set()):
                    report.add(
                        "ERROR",
                        "REFERENCE_LEVEL_MISSING",
                        f"Reference level {reference!r} is not present for factor {factor!r}.",
                    )

    report.summary.update(
        {
            "sample_rows": len(rows),
            "readable_indexed_bams": bam_count,
            "model_cells": {" | ".join(cell): count for cell, count in cell_counts.items()},
        }
    )
    return rows, bam_contigs


def _integer_parameter(
    block: dict[str, Any],
    key: str,
    report: ValidationReport,
    *,
    minimum: int | None = None,
) -> int | None:
    value = block.get(key)
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        report.add("ERROR", "INTEGER_PARAMETER", f"{key} must be an integer.")
        return None
    if minimum is not None and parsed < minimum:
        report.add(
            "ERROR",
            "INTEGER_PARAMETER_RANGE",
            f"{key} must be at least {minimum}; received {parsed}.",
        )
    return parsed


def _validate_analysis_parameters(
    config: dict[str, Any], report: ValidationReport
) -> None:
    reads = _require_mapping(config, "reads", report)
    psite = _require_mapping(config, "psite", report)
    profiles = _require_mapping(config, "profiles", report)
    background = _require_mapping(config, "background", report)
    coverage = _require_mapping(config, "coverage", report)
    statistics = _require_mapping(config, "statistics", report)
    qc = config.get("qc", {})
    if qc is not None and not isinstance(qc, dict):
        report.add("ERROR", "QC_CONFIG", "qc must be a mapping.")
        qc = {}

    lengths = reads.get("allowed_lengths")
    if not isinstance(lengths, list) or not lengths:
        report.add(
            "ERROR", "READ_LENGTHS", "reads.allowed_lengths must be a non-empty list."
        )
    else:
        try:
            parsed_lengths = [int(value) for value in lengths]
        except (TypeError, ValueError):
            report.add("ERROR", "READ_LENGTHS", "All allowed read lengths must be integers.")
        else:
            if any(value <= 0 for value in parsed_lengths) or len(set(parsed_lengths)) != len(parsed_lengths):
                report.add(
                    "ERROR",
                    "READ_LENGTHS",
                    "Allowed read lengths must be unique positive integers.",
                )
    mapq = _integer_parameter(reads, "minimum_mapping_quality", report, minimum=0)
    if mapq is not None and mapq > 255:
        report.add("ERROR", "MAPQ_RANGE", "reads.minimum_mapping_quality cannot exceed 255.")

    mode = psite.get("mode")
    if mode not in {"calibrate", "fixed_table"}:
        report.add("ERROR", "PSITE_MODE", "psite.mode must be calibrate or fixed_table.")
    if mode == "fixed_table":
        _check_file(psite.get("offset_table"), "psite.offset_table", report)
    candidates = psite.get("candidate_offsets")
    if not isinstance(candidates, list) or not candidates:
        report.add(
            "ERROR", "PSITE_CANDIDATES", "psite.candidate_offsets must be a non-empty list."
        )
    else:
        try:
            parsed_candidates = [int(value) for value in candidates]
        except (TypeError, ValueError):
            report.add("ERROR", "PSITE_CANDIDATES", "P-site offsets must be integers.")
        else:
            if any(value < 0 for value in parsed_candidates):
                report.add("ERROR", "PSITE_CANDIDATES", "P-site offsets cannot be negative.")
    _integer_parameter(psite, "minimum_reads_per_length", report, minimum=1)
    metagene_start = _integer_parameter(
        psite, "start_codon_metagene_start_nt", report, minimum=0
    )
    metagene_end = _integer_parameter(
        psite, "start_codon_metagene_end_nt", report, minimum=0
    )
    if (
        metagene_start is not None
        and metagene_end is not None
        and metagene_end < metagene_start
    ):
        report.add(
            "ERROR",
            "CDS_START_METAGENE_RANGE",
            "psite.start_codon_metagene_end_nt must be >= start_nt.",
        )

    scan_start = _integer_parameter(profiles, "scan_start_nt", report)
    scan_end = _integer_parameter(profiles, "scan_end_nt", report)
    step = _integer_parameter(profiles, "step_nt", report, minimum=1)
    width = _integer_parameter(profiles, "window_width_nt", report, minimum=1)
    primary = _integer_parameter(profiles, "primary_center_nt", report)
    if width is not None and width % 2 == 0:
        report.add(
            "ERROR", "WINDOW_WIDTH", "profiles.window_width_nt must be an odd integer."
        )
    if scan_start is not None and scan_end is not None and scan_end < scan_start:
        report.add("ERROR", "SCAN_RANGE", "profiles.scan_end_nt must be >= scan_start_nt.")
    if all(value is not None for value in (scan_start, scan_end, step, primary)):
        centers = list(range(scan_start, scan_end + 1, step))
        if primary not in centers:
            report.add(
                "ERROR",
                "PRIMARY_CENTER",
                "profiles.primary_center_nt must be one of the configured scan centers.",
            )

    if background.get("strategy") != "gene_wide_non_tmd_cds":
        report.add(
            "ERROR",
            "BACKGROUND_STRATEGY",
            "v0.1 supports background.strategy=gene_wide_non_tmd_cds only.",
        )
    if background.get("exclude_neighboring_tmds") is not True:
        report.add(
            "ERROR",
            "BACKGROUND_NEIGHBORS",
            "v0.1 requires background.exclude_neighboring_tmds=true.",
        )
    if background.get("exclude_target_scan_region") is not True:
        report.add(
            "ERROR",
            "BACKGROUND_SCAN_EXCLUSION",
            "v0.1 requires background.exclude_target_scan_region=true.",
        )

    _integer_parameter(coverage, "minimum_background_count", report, minimum=0)
    if coverage.get("eligibility") != "strict_all_model_samples":
        report.add(
            "ERROR",
            "COVERAGE_ELIGIBILITY",
            "v0.1 supports coverage.eligibility=strict_all_model_samples only.",
        )

    if statistics.get("backend") != "edgeR_ql":
        report.add("ERROR", "STATISTICS_BACKEND", "statistics.backend must be edgeR_ql.")
    if statistics.get("fdr_method") != "BH":
        report.add("ERROR", "FDR_METHOD", "statistics.fdr_method must be BH in v0.1.")
    try:
        fdr = float(statistics.get("fdr_threshold"))
        if not 0 < fdr <= 1:
            raise ValueError
    except (TypeError, ValueError):
        report.add(
            "ERROR", "FDR_THRESHOLD", "statistics.fdr_threshold must be in (0, 1]."
        )

    for key in ("minimum_frame0_fraction_warn", "minimum_replicate_correlation_warn"):
        if key not in qc:
            continue
        try:
            value = float(qc[key])
            if not -1 <= value <= 1:
                raise ValueError
        except (TypeError, ValueError):
            report.add(
                "ERROR", "QC_THRESHOLD", f"qc.{key} must be a number in [-1, 1]."
            )
    if "minimum_common_tmds_warn" in qc:
        _integer_parameter(qc, "minimum_common_tmds_warn", report, minimum=1)


def _validate_tmds(
    block: dict[str, Any],
    project_mode: str,
    cds_lengths: dict[str, int],
    report: ValidationReport,
) -> None:
    path = _check_file(block.get("path"), "tmd_annotation.path", report)
    if path is None:
        return
    columns, rows = _read_table(path)
    coordinate_system = block.get("coordinate_system")
    mapping = block.get("columns", {})
    if not isinstance(mapping, dict):
        report.add("ERROR", "TMD_COLUMN_MAPPING", "tmd_annotation.columns must be a mapping.")
        return

    if coordinate_system == "paper_legacy_genomic":
        if project_mode != "paper_reproduction":
            report.add(
                "ERROR",
                "LEGACY_COORDINATE_MODE",
                "paper_legacy_genomic is allowed only in project.mode=paper_reproduction.",
            )
            return
        required = ["gene", "tmd_index", "tmd_start_nt", "tmd_end_nt", "chrom", "strand"]
        missing = [name for name in required if name not in columns]
        if missing:
            report.add(
                "ERROR",
                "LEGACY_TMD_COLUMNS",
                f"Legacy paper TMD table is missing: {', '.join(missing)}",
            )
        else:
            report.add(
                "WARN",
                "LEGACY_PAPER_ADAPTER",
                "Using the paper-only yeast genomic TMD schema; this is not an external-data input format.",
            )
        report.summary["tmd_rows"] = len(rows)
        return

    if coordinate_system != "cds_nt_1based_inclusive":
        report.add(
            "ERROR",
            "TMD_COORDINATE_SYSTEM",
            "External v1 input requires coordinate_system=cds_nt_1based_inclusive.",
            received=coordinate_system,
        )
        return

    required_logical = [
        "gene_id",
        "transcript_id",
        "tmd_id",
        "start",
        "end",
    ]
    physical = {name: mapping.get(name, name) for name in required_logical}
    for logical, column in physical.items():
        if column not in columns:
            report.add(
                "ERROR",
                "TMD_COLUMN_MISSING",
                f"TMD annotation is missing column {column!r} for {logical!r}.",
            )
    if any(column not in columns for column in physical.values()):
        return

    seen_ids: set[str] = set()
    duplicate_ids: set[str] = set()
    missing_transcripts = 0
    outside_cds = 0
    invalid_coordinates = 0
    for row in rows:
        tmd_id = row[physical["tmd_id"]].strip()
        transcript_id = row[physical["transcript_id"]].strip()
        if tmd_id in seen_ids:
            duplicate_ids.add(tmd_id)
        seen_ids.add(tmd_id)
        try:
            start = int(row[physical["start"]])
            end = int(row[physical["end"]])
        except ValueError:
            invalid_coordinates += 1
            continue
        if start < 1 or end < start:
            invalid_coordinates += 1
            continue
        cds_length = cds_lengths.get(transcript_id)
        if cds_length is None:
            missing_transcripts += 1
        elif end > cds_length:
            outside_cds += 1

    if duplicate_ids:
        report.add(
            "ERROR",
            "DUPLICATE_TMD_ID",
            f"Found {len(duplicate_ids)} duplicate TMD ID(s).",
            examples=sorted(duplicate_ids)[:10],
        )
    if invalid_coordinates:
        report.add(
            "ERROR",
            "INVALID_TMD_COORDINATES",
            f"Found {invalid_coordinates} TMD row(s) with invalid coordinates.",
        )
    if missing_transcripts:
        report.add(
            "ERROR",
            "TMD_TRANSCRIPT_NOT_FOUND",
            f"Found {missing_transcripts} TMD row(s) whose transcript_id is absent from annotated CDS features.",
        )
    if outside_cds:
        report.add(
            "ERROR",
            "TMD_OUTSIDE_CDS",
            f"Found {outside_cds} TMD row(s) extending outside annotated CDS coordinates. "
            "Check whether the TMD file uses amino-acid or nucleotide coordinates.",
        )
    report.summary["tmd_rows"] = len(rows)


def _validate_optional_modules(
    config: dict[str, Any], report: ValidationReport
) -> None:
    optional = config.get("optional_modules", {})
    if optional in (None, {}):
        return
    if not isinstance(optional, dict):
        report.add(
            "ERROR", "OPTIONAL_MODULES", "optional_modules must be a mapping."
        )
        return
    topology = optional.get("topology", {})
    topology_enabled = isinstance(topology, dict) and bool(topology.get("enabled", False))

    sequence = optional.get("sequence_features", {})
    if isinstance(sequence, dict) and bool(sequence.get("enabled", False)):
        if not topology_enabled:
            report.add(
                "ERROR",
                "SEQUENCE_REQUIRES_TOPOLOGY",
                "sequence_features requires optional_modules.topology.enabled=true.",
            )
        _check_file(
            config.get("reference", {}).get("peptide_fasta"),
            "reference.peptide_fasta",
            report,
        )
        mapping = sequence.get("protein_mapping")
        if not isinstance(mapping, dict):
            report.add(
                "ERROR",
                "PROTEIN_MAPPING",
                "sequence_features.protein_mapping must be a mapping.",
            )
        else:
            path = _check_file(
                mapping.get("path"), "sequence_features.protein_mapping.path", report
            )
            if mapping.get("key") not in {"gene_id", "transcript_id"}:
                report.add(
                    "ERROR",
                    "PROTEIN_MAPPING_KEY",
                    "protein_mapping.key must be gene_id or transcript_id.",
                )
            if path is not None:
                columns, _ = _read_table(path)
                required = [
                    str(mapping.get("key_column", mapping.get("key", "transcript_id"))),
                    str(mapping.get("protein_id_column", "protein_id")),
                ]
                missing = [column for column in required if column not in columns]
                if missing:
                    report.add(
                        "ERROR",
                        "PROTEIN_MAPPING_COLUMNS",
                        f"Protein mapping table lacks columns: {missing}.",
                    )

    overlap = optional.get("overlap", {})
    if isinstance(overlap, dict) and bool(overlap.get("enabled", False)):
        if not topology_enabled:
            report.add(
                "ERROR",
                "OVERLAP_REQUIRES_TOPOLOGY",
                "overlap analysis requires optional_modules.topology.enabled=true.",
            )

    enrichment = optional.get("enrichment", {})
    if isinstance(enrichment, dict) and bool(enrichment.get("enabled", False)):
        mapping = enrichment.get("mapping")
        if not isinstance(mapping, dict):
            report.add(
                "ERROR",
                "ENRICHMENT_MAPPING",
                "enrichment.mapping must be a mapping.",
            )
        else:
            path = _check_file(
                mapping.get("path"), "enrichment.mapping.path", report
            )
            if path is not None:
                columns, _ = _read_table(path)
                required = [
                    str(mapping.get("gene_column", "gene_id")),
                    str(mapping.get("term_column", "term_id")),
                ]
                for name in ("ontology_column", "term_name_column"):
                    if mapping.get(name):
                        required.append(str(mapping[name]))
                missing = [column for column in required if column not in columns]
                if missing:
                    report.add(
                        "ERROR",
                        "ENRICHMENT_MAPPING_COLUMNS",
                        f"Enrichment mapping table lacks columns: {missing}.",
                    )


def validate_config(config: dict[str, Any]) -> ValidationReport:
    """Validate config structure and currently supported external inputs."""

    report = ValidationReport()
    for location in unresolved_environment_tokens(config):
        report.add(
            "ERROR",
            "UNRESOLVED_ENVIRONMENT",
            f"Unresolved environment variable at {location}.",
        )

    project = _require_mapping(config, "project", report)
    reference = _require_mapping(config, "reference", report)
    tmd_block = _require_mapping(config, "tmd_annotation", report)
    samples = _require_mapping(config, "samples", report)
    design = _require_mapping(config, "design", report)
    _require_mapping(config, "profiles", report)
    _require_mapping(config, "background", report)
    _require_mapping(config, "coverage", report)
    _require_mapping(config, "statistics", report)
    output = _require_mapping(config, "output", report)

    _validate_analysis_parameters(config, report)

    project_mode = str(project.get("mode", "external"))
    if project_mode not in {"external", "paper_reproduction"}:
        report.add(
            "ERROR",
            "PROJECT_MODE",
            "project.mode must be external or paper_reproduction.",
        )

    design_type = design.get("type")
    if design_type != "factorial_2x2":
        report.add(
            "ERROR",
            "DESIGN_TYPE",
            "The validated v0.1 backend supports design.type=factorial_2x2 only.",
            received=design_type,
        )
    factors = design.get("factors", [])
    expected_factor_count = 2
    if not isinstance(factors, list) or len(factors) != expected_factor_count:
        report.add(
            "ERROR",
            "DESIGN_FACTORS",
            f"{design_type} requires exactly {expected_factor_count} factor name(s).",
        )

    annotation_format = str(reference.get("annotation_format", "gtf")).lower()
    if annotation_format not in {"gtf", "gff3"}:
        report.add(
            "ERROR",
            "ANNOTATION_FORMAT",
            "reference.annotation_format must be gtf or gff3.",
        )
    if annotation_format == "gff3" and project_mode == "external":
        report.add(
            "ERROR",
            "GFF3_NOT_VALIDATED",
            "GFF3 parsing is not yet part of the validated v1 support contract; use GTF.",
        )

    annotation_path = _check_file(
        reference.get("annotation"), "reference.annotation", report
    )
    annotation_contigs: set[str] = set()
    cds_lengths: dict[str, int] = {}
    if annotation_path is not None:
        annotation_contigs, cds_lengths, _ = _annotation_index(
            annotation_path, annotation_format, report
        )

    _, bam_contigs = _validate_metadata(samples, design, report)
    if annotation_contigs and bam_contigs:
        overlap = annotation_contigs & bam_contigs
        if not overlap:
            report.add(
                "ERROR",
                "CONTIG_MISMATCH",
                "BAM and annotation have no chromosome/contig names in common.",
                annotation_examples=sorted(annotation_contigs)[:10],
                bam_examples=sorted(bam_contigs)[:10],
            )
        else:
            report.summary["shared_annotation_bam_contigs"] = len(overlap)

    _validate_tmds(tmd_block, project_mode, cds_lengths, report)
    _validate_optional_modules(config, report)

    output_directory = output.get("directory")
    if output_directory in (None, ""):
        report.add("ERROR", "OUTPUT_DIRECTORY", "output.directory is required.")

    report.add(
        "INFO",
        "COVERAGE_DEFERRED",
        "Coverage sufficiency is evaluated after read counting, not during structural input validation.",
    )
    return report


def format_report(report: ValidationReport) -> str:
    status = "PASS" if report.passed else "FAIL"
    lines = [
        f"Validation: {status} ({len(report.errors)} error(s), {len(report.warnings)} warning(s))"
    ]
    for item in report.diagnostics:
        lines.append(f"{item.severity}: [{item.code}] {item.message}")
    if report.summary:
        lines.append("Summary:")
        for key, value in sorted(report.summary.items()):
            lines.append(f"  {key}: {value}")
    return "\n".join(lines)
