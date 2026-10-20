#!/usr/bin/env python3
"""Build an auditable rRNA/tRNA contaminant FASTA from the frozen genome/GTF."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "references" / "genome"
FASTA = REF / "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
GTF = REF / "Saccharomyces_cerevisiae.R64-1-1.115.gtf"
OUT = REF / "contaminants"


def revcomp(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def wrap(sequence, width=80):
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


genome = {}
name = None
parts = []
with FASTA.open() as handle:
    for line in handle:
        if line.startswith(">"):
            if name is not None:
                genome[name] = "".join(parts).upper()
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if name is not None:
        genome[name] = "".join(parts).upper()

records = []
with GTF.open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields[2] != "gene":
            continue
        attrs = fields[8]
        match = re.search(r'gene_biotype "([^"]+)"', attrs)
        if not match or match.group(1) not in {"rRNA", "tRNA"}:
            continue
        gene_id = re.search(r'gene_id "([^"]+)"', attrs).group(1)
        chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
        sequence = genome[chrom][start - 1:end]
        if strand == "-":
            sequence = revcomp(sequence)
        records.append((f"{match.group(1)}|{gene_id}|{chrom}:{start}-{end}:{strand}", sequence))

# Known adapter sequences are included as synthetic controls. The original
# cloning-oligo sequences are not supplied in the accessible GEO metadata.
adapters = {
    "adapter|barcoded": "AGATCGGAAGAGCACACGTCTGAACTC",
    "adapter|unbarcoded": "CTGTAGGCACCATCAAT",
}
for label, sequence in adapters.items():
    records.append((label, sequence))
    records.append((label + "|revcomp", revcomp(sequence)))

OUT.mkdir(parents=True, exist_ok=True)
with (OUT / "rrna_trna_adapters.fa").open("w") as handle:
    for header, sequence in records:
        handle.write(f">{header}\n{wrap(sequence)}\n")
with (OUT / "README.md").open("w") as handle:
    handle.write(
        "# Contaminant reference\n\n"
        "Generated from gene features annotated as `rRNA` or `tRNA` in "
        "Ensembl R64-1-1 release 115, plus the two GEO-documented adapter "
        "sequences and their reverse complements. Original cloning-oligo "
        "sequences remain unavailable and are not guessed.\n"
    )
print(f"Wrote {len(records)} records to {OUT / 'rrna_trna_adapters.fa'}")
