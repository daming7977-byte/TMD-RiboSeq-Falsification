# P3 eLife Upgrade — U2B1 Exact-Family Duplicate-Collapse Contract v1

Status: FROZEN BEFORE U2B1 DUPLICATE-COLLAPSE OUTPUTS

## Purpose
Construct a conservative duplicate-collapsed version of the 12 frozen GSE300977 formal counting input BAMs using the recoverable 5-nt UMI field retained in BAM query names.

U2B1 is an engineering/technical stage only. It does not calculate TMD-centered effects, topology interactions, or any biological result.

## Frozen input BAMs
Exactly:
`unique_mapping_blinded/bam/B01.NH1.unique.bam` through `B12.NH1.unique.bam`.

These are the same `NH:i:1` unique BAMs used by the formal P3 CenterMap counting workflow.

Frozen expected total input-read count across all 12 BAMs:
`215999661`

## Recoverable UMI field
Read identifier parser:
`UMI5=([ACGTN]{5})`

The sample barcode field (`BC4`) is not part of molecule identity and is not used for filtering.

## Exact-family key
For reads with a valid A/C/G/T-only UMI5, duplicate families are defined within each library by:
- reference contig
- strand
- CIGAR-aware biological 5' genomic coordinate
  - forward: `reference_start`
  - reverse: `reference_end - 1`
- `reference_start`
- `reference_end`
- exact CIGAR string
- full query length
- aligned query length
- exact 5-nt UMI sequence

The explicit reference start/end and CIGAR fields make this an intentionally conservative exact-alignment family definition and prevent collapsing distinct splice/alignment geometries.

## Ambiguous or missing UMI5
If UMI5 contains `N`, or is absent/malformed, retain the read unchanged as a singleton and never collapse it with another read.

## Representative rule
Within an exact family, retain the first occurrence in the frozen coordinate-sorted input BAM.

No MAPQ/NM-based rescue or result-dependent representative selection is allowed.

## Explicitly prohibited
- no Hamming-distance UMI correction
- no adjacency/directional UMI correction
- no barcode filtering
- no fuzzy coordinate matching
- no collapsing across libraries
- no collapsing N-containing UMIs
- no remapping
- no re-alignment
- no pseudo-anchor rematching
- no TMD-centered biological calculation in U2B1
- no change to the v6 baseline

## Output
Create 12 coordinate-preserving BAMs under a new U2B1 directory and index them.

Report per library:
- input reads
- reads with A/C/G/T-only UMI5
- ambiguous/malformed UMI reads
- exact valid-UMI families
- duplicate families
- maximum family size
- collapsed reads
- output reads
- retention fraction
- collapsed fraction
- unexpected unmapped/secondary/supplementary reads

## Gates
U2B1_INPUT_COUNT_GATE = PASS only if total input reads = 215999661.

U2B1_QNAME_UMI_GATE = PASS only if every mapped primary read contains a parseable 5-character `UMI5=` field consisting of A/C/G/T/N.

U2B1_OUTPUT_ACCOUNTING_GATE = PASS only if, for every library:
`output_reads + collapsed_reads = input_reads`
and output BAM indexing succeeds.

## Interpretation boundary
U2B1 output is a technical sensitivity input. It is not a replacement for the frozen primary BAMs and does not alter any existing conclusion.

A downstream U2B2 biological sensitivity may be run only after U2B1 completes and its accounting is reviewed, with U2B2 rules frozen before viewing deduplicated TMD-centered results.
