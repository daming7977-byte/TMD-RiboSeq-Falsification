# B2 matched pseudo-anchor contract v0

Status: annotation-only definitions frozen before B2 positional outcome access. Seed **6101207**. The machine-readable rules in annotation_rules.json were SHA-frozen before eligibility enumeration. The final pseudo table is authoritative for later execution; do not regenerate from a changing annotation service.

For each eligible real anchor a in a protein of L translated amino acids, enumerate integer CDS codon positions q. The start codon has index 0; the terminal stop is excluded from L. Both real and pseudo positions require the complete inclusive visualization window [-60,+100] inside translated CDS. Every nucleotide in that window must belong to exactly one same-strand retained canonical coding model. The same support rule applies on both strands and across spliced coding exons.

For every source TMD [s,e) in the same gene, exclude pseudo q in inclusive [s-60,e-1+60]. Thus the buffer is **60 codons beyond each complete TMD interval**, not merely beyond its start. This also makes the primary pseudo window [-30,+60] free of annotated TMD residues. The broader visualization window can intersect a TMD and is not a formal non-TMD-region estimand. Non-TMD means absent from this frozen annotation, not proof of no membrane interaction.

Relative CDS position is q/(L-1). Require |q-a|/(L-1) <= **0.20**, implemented exactly as 5*|q-a| <= L-1. Do not enlarge this tolerance for short or TMD-dense genes. It bounds residual position mismatch; it does not establish position equivalence or remove all CDS-position confounding.

Choose exactly **five** distinct positions per real anchor. Rank valid candidates by:

1. Absolute codon distance |q-a| ascending.
2. SHA256 of UTF-8 `6101207|anchor_id|q`, interpreted as lowercase hexadecimal, ascending.
3. Numeric q ascending in the event of a hash collision.

Here anchor_id is ASCII `gene_locus|source_start_0based|source_end_exclusive`, and q is an unpadded base-10 integer. Take the first five. This is deterministic nearest-position matching with seeded hash randomization of distance ties, not uniform sampling of the CDS. It uses no read counts. Anchors are ordered by ASCII locus, numeric start, numeric end. Identical source intervals collapse with all source-row provenance retained. No duplicate intervals were found in this source.

Fewer than five valid candidates excludes the real anchor. There is no fallback, replacement within a matched set, tolerance expansion, alternative seed or second matching run. Pseudo positions may be reused across different real anchors, and real/pseudo windows may overlap across sets. Reuse and spatial correlation are preserved by resampling whole genes. Pseudo positions are not independent observations or biological replicates.

The frozen result is **1,068 real anchors, 338 genes and 5,340 pseudo definitions**. Use the identical table for 28-nt primary and 32-nt sensitivity. The signed and relative distances and seeded tie hashes are recorded per match. No later count/coverage filter may change the anchor population.

Sources: RULES, SUPP2, UNIPROT and TASK6_REF in B2_BENCHMARK_DESIGN_SOURCE_LOG_v0.tsv. The numerical buffer, tolerance and matching multiplicity are investigator design choices, not claims that the source study used these parameters.
