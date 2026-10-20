"""Transcript/CDS coordinate models used by all spatial analyses."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .errors import ValidationFailed


def parse_gtf_attributes(text: str) -> dict[str, str]:
    """Parse standard GTF ``key "value";`` attributes."""

    return {key: value for key, value in re.findall(r'(\S+)\s+"([^"]*)"', text)}


@dataclass(frozen=True)
class CDSSegment:
    contig: str
    start0: int
    end0: int
    strand: str
    phase: int | None

    @property
    def length(self) -> int:
        return self.end0 - self.start0


@dataclass(frozen=True)
class TranscriptModel:
    gene_id: str
    transcript_id: str
    contig: str
    strand: str
    genomic_positions0: tuple[int, ...]

    @property
    def cds_length_nt(self) -> int:
        return len(self.genomic_positions0)

    def genomic_position0(self, cds_position_1based: int) -> int:
        if cds_position_1based < 1 or cds_position_1based > self.cds_length_nt:
            raise IndexError(
                f"CDS position {cds_position_1based} is outside transcript "
                f"{self.transcript_id} (length {self.cds_length_nt})."
            )
        return self.genomic_positions0[cds_position_1based - 1]

    def genomic_interval0(
        self, start_cds_1based: int, end_cds_1based: int
    ) -> tuple[int, ...]:
        if end_cds_1based < start_cds_1based:
            raise ValueError("CDS interval end precedes start.")
        self.genomic_position0(start_cds_1based)
        self.genomic_position0(end_cds_1based)
        return self.genomic_positions0[start_cds_1based - 1 : end_cds_1based]

    def cds_index_by_genomic_position0(self) -> dict[int, int]:
        """Map genomic positions to 1-based CDS positions for this transcript."""

        return {position: index for index, position in enumerate(self.genomic_positions0, 1)}


def _build_model(
    gene_id: str,
    transcript_id: str,
    segments: Iterable[CDSSegment],
) -> TranscriptModel:
    items = list(segments)
    contigs = {item.contig for item in items}
    strands = {item.strand for item in items}
    if len(contigs) != 1 or len(strands) != 1:
        raise ValidationFailed(
            f"Transcript {transcript_id} has CDS segments on conflicting contigs/strands."
        )
    strand = next(iter(strands))
    if strand not in {"+", "-"}:
        raise ValidationFailed(f"Transcript {transcript_id} has invalid strand {strand!r}.")
    ordered = sorted(items, key=lambda item: item.start0, reverse=(strand == "-"))
    positions: list[int] = []
    for segment in ordered:
        if segment.length <= 0:
            raise ValidationFailed(f"Transcript {transcript_id} has an empty CDS segment.")
        if strand == "+":
            positions.extend(range(segment.start0, segment.end0))
        else:
            positions.extend(range(segment.end0 - 1, segment.start0 - 1, -1))
    if len(set(positions)) != len(positions):
        raise ValidationFailed(f"Transcript {transcript_id} has overlapping CDS segments.")
    return TranscriptModel(
        gene_id=gene_id,
        transcript_id=transcript_id,
        contig=next(iter(contigs)),
        strand=strand,
        genomic_positions0=tuple(positions),
    )


def load_gtf_transcripts(path: str | Path) -> dict[str, TranscriptModel]:
    """Load transcript-oriented CDS coordinate maps from a GTF file.

    GTF coordinates are converted from 1-based closed intervals to 0-based
    genomic positions. The returned CDS coordinates always increase in the
    biological translation direction, including minus-strand transcripts.
    """

    grouped: dict[str, list[CDSSegment]] = {}
    genes: dict[str, str] = {}
    path = Path(path)
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id") or attrs.get("gene_name")
            if not transcript_id or not gene_id:
                raise ValidationFailed(
                    f"CDS line {line_number} lacks transcript_id or gene_id."
                )
            try:
                start0 = int(fields[3]) - 1
                end0 = int(fields[4])
            except ValueError as exc:
                raise ValidationFailed(
                    f"Invalid CDS coordinates at line {line_number}."
                ) from exc
            phase = None if fields[7] == "." else int(fields[7])
            grouped.setdefault(transcript_id, []).append(
                CDSSegment(fields[0], start0, end0, fields[6], phase)
            )
            previous_gene = genes.setdefault(transcript_id, gene_id)
            if previous_gene != gene_id:
                raise ValidationFailed(
                    f"Transcript {transcript_id} is assigned to multiple genes."
                )
    if not grouped:
        raise ValidationFailed(f"No transcript-linked CDS features found in {path}.")
    return {
        transcript_id: _build_model(genes[transcript_id], transcript_id, segments)
        for transcript_id, segments in grouped.items()
    }

