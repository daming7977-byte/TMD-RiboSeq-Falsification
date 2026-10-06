"""CIGAR-aware P-site mapping and read filters."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import pysam


@dataclass(frozen=True)
class ReadFilter:
    minimum_mapping_quality: int = 0
    include_secondary: bool = False
    include_supplementary: bool = False
    include_duplicates: bool = False


def read_is_eligible(read: pysam.AlignedSegment, policy: ReadFilter) -> bool:
    if read.is_unmapped:
        return False
    if read.mapping_quality < policy.minimum_mapping_quality:
        return False
    if read.is_secondary and not policy.include_secondary:
        return False
    if read.is_supplementary and not policy.include_supplementary:
        return False
    if read.is_duplicate and not policy.include_duplicates:
        return False
    return True


def aligned_five_prime_reference_positions(
    read: pysam.AlignedSegment,
) -> list[int]:
    """Return aligned query-base positions from the footprint's 5' end.

    Soft-clipped and inserted query bases do not have reference coordinates and
    are excluded. Deletions and splice gaps are traversed correctly because the
    list is derived from the CIGAR rather than arithmetic on alignment bounds.
    """

    positions = list(read.get_reference_positions(full_length=False))
    if read.is_reverse:
        positions.reverse()
    return positions


def psite_reference_position0(
    read: pysam.AlignedSegment, offset: int
) -> int | None:
    """Map a 0-based P-site offset from the aligned 5' end to the reference."""

    if offset < 0:
        raise ValueError("P-site offset must be non-negative.")
    positions = aligned_five_prime_reference_positions(read)
    if offset >= len(positions):
        return None
    return positions[offset]


def psite_for_read(
    read: pysam.AlignedSegment,
    offsets_by_read_length: Mapping[int, int],
) -> int | None:
    read_length = read.query_length
    if read_length is None:
        return None
    offset = offsets_by_read_length.get(int(read_length))
    if offset is None:
        return None
    return psite_reference_position0(read, int(offset))

