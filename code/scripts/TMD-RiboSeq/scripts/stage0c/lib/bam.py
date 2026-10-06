"""
Stage0C BAM utilities.
"""

from pathlib import Path
import pysam


def open_bam(path):
    """
    Open BAM file.

    Parameters
    ----------
    path : str | Path

    Returns
    -------
    pysam.AlignmentFile
    """

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(path)

    return pysam.AlignmentFile(path, "rb")


def iter_primary_reads(bam):
    """
    Iterate over primary mapped reads only.

    Yields
    ------
    pysam.AlignedSegment
    """

    for read in bam.fetch(until_eof=True):

        if read.is_unmapped:
            continue

        if read.is_secondary:
            continue

        if read.is_supplementary:
            continue

        yield read


def five_prime_position(read):
    """
    Return transcript-coordinate 5' position.

    Assumes transcriptome alignment.
    """

    return read.reference_start


def read_length(read):
    """
    Return query length.
    """

    return read.query_length
