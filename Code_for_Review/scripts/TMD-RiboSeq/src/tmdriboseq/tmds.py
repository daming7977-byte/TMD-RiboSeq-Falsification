"""Canonical TMD annotation records."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .errors import ValidationFailed


@dataclass(frozen=True)
class TMDRecord:
    gene_id: str
    transcript_id: str
    tmd_id: str
    start_cds_nt: int
    end_cds_nt: int

    @property
    def length_nt(self) -> int:
        return self.end_cds_nt - self.start_cds_nt + 1


def load_canonical_tmds(
    path: str | Path,
    column_mapping: dict[str, str] | None = None,
) -> list[TMDRecord]:
    path = Path(path)
    mapping = {
        "gene_id": "gene_id",
        "transcript_id": "transcript_id",
        "tmd_id": "tmd_id",
        "start": "tmd_start_cds_nt",
        "end": "tmd_end_cds_nt",
    }
    if column_mapping:
        mapping.update(column_mapping)
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    records: list[TMDRecord] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        missing = [column for column in mapping.values() if column not in (reader.fieldnames or [])]
        if missing:
            raise ValidationFailed(f"TMD table is missing columns: {', '.join(missing)}")
        for row_number, row in enumerate(reader, start=2):
            try:
                start = int(row[mapping["start"]])
                end = int(row[mapping["end"]])
            except ValueError as exc:
                raise ValidationFailed(
                    f"TMD table row {row_number} has non-integer coordinates."
                ) from exc
            if start < 1 or end < start:
                raise ValidationFailed(
                    f"TMD table row {row_number} has invalid coordinates {start}–{end}."
                )
            records.append(
                TMDRecord(
                    gene_id=row[mapping["gene_id"]].strip(),
                    transcript_id=row[mapping["transcript_id"]].strip(),
                    tmd_id=row[mapping["tmd_id"]].strip(),
                    start_cds_nt=start,
                    end_cds_nt=end,
                )
            )
    return records


def _gtf_gene_name_to_transcript(path: str | Path) -> dict[str, str]:
    associations: dict[str, set[str]] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = {
                key: value
                for key, value in re.findall(r'(\S+)\s+"([^"]*)"', fields[8])
            }
            gene_name = attrs.get("gene_name")
            transcript_id = attrs.get("transcript_id")
            if gene_name and transcript_id:
                associations.setdefault(gene_name, set()).add(transcript_id)
    ambiguous = {gene: values for gene, values in associations.items() if len(values) != 1}
    if ambiguous:
        examples = ", ".join(sorted(ambiguous)[:5])
        raise ValidationFailed(
            "The paper legacy adapter requires exactly one transcript per gene; "
            f"ambiguous examples: {examples}."
        )
    return {gene: next(iter(values)) for gene, values in associations.items()}


def load_paper_legacy_tmds(
    path: str | Path,
    gtf_path: str | Path,
) -> list[TMDRecord]:
    """Convert the frozen yeast paper table into canonical in-memory records.

    This adapter is deliberately unavailable for external configurations.
    """

    gene_to_transcript = _gtf_gene_name_to_transcript(gtf_path)
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "gene",
            "tmd_index",
            "tmd_start_nt",
            "tmd_end_nt",
            "tmd_start_genome",
            "chrom",
            "strand",
        }
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValidationFailed(
                f"Paper legacy TMD table is missing: {', '.join(missing)}"
            )
        records: list[TMDRecord] = []
        for row_number, row in enumerate(reader, start=2):
            gene = row["gene"].strip()
            transcript_id = gene_to_transcript.get(gene)
            if transcript_id is None:
                raise ValidationFailed(
                    f"Paper TMD row {row_number}: gene {gene} has no unique GTF transcript."
                )
            records.append(
                TMDRecord(
                    gene_id=gene,
                    transcript_id=transcript_id,
                    tmd_id=f"{gene}_TMD{int(row['tmd_index'])}",
                    start_cds_nt=int(row["tmd_start_nt"]),
                    end_cds_nt=int(row["tmd_end_nt"]),
                )
            )
    return records
