# B2 Ubc6 benchmark design freeze v0

Status: **DESIGN FROZEN; BENCHMARK NOT EXECUTED; SCIENTIFIC ADJUDICATION REQUIRED.** Created 2026-09-10 in an isolated post-submission directory. This contract incorporates B2_PSEUDO_ANCHOR_CONTRACT_v0.md and B2_PRIMARY_ESTIMAND_CONTRACT_v0.md; their explicit rules are normative. SHA256SUMS.txt identifies the exact frozen files. Any failure of a required file/hash or coordinate check stops execution, without discretionary repair.

## Technical scope

Primary: **28 nt / A-site +15 nt** from mapped 5-prime end in coding direction. Separate sensitivity: **32 nt / +16 nt**, retaining the Task 6 GSM-pooling qualification. Exclude 26,27,29,30 nt for failed periodicity and 31 nt for ambiguous calibration. No reintroduction or combined optimized population.

Input GSM1495513/SRX691037 combines SRR1562895 and SRR1562896. Pulldown GSM1495514/SRX691038 combines SRR1562897 and SRR1562898. No additional samples, UMI deduction, duplicate collapse or biological replication claim.

Use exactly the Task 6 technical contract, SHA256 `9bee7f2fd74b253e32dc3d5a128f6f96a419107e071d066f3f1f5e9dbc5bb8a2`, and its frozen references and alignments. GENOME_BUILD=sacCer3, SGD April 2011 sequence. ANNOTATION_RELEASE=archived UCSC sgdGene/sgdOther 2011-10-06 snapshots. This is an independent reconstruction reference, not the historical author annotation/exporter.

Key hashes:

- genome.fa: `92634e39d9d7013b52e82de31606fad36da15fed9704cdf1f487e3b29322388b`
- coding_models.json: `1a8539be61e88ed3f94305e5eb31233d8dce11380ff749be748cf73cf9e93a76`
- sgdGene.txt.gz: `52f8b74d483367dd835e3128e64631e8840f3eecfbe07504967953515d5e9f5f`
- sgdOther.txt.gz: `f26b57dc042d742b11b8272d89c355aa7f889533ced1efffe081f36859224cb1`

Full technical references and index hashes remain in the Task 6 manifest, referenced in the source log. Existing raw and submitted P3 files were not modified.

## Exact TMD annotation source

Use **Shurtleff et al. 2018, eLife 7:e37018, Supplementary file 2**, DOI **10.7554/eLife.37018.013**, publisher asset **elife-37018-supp2-v2.xlsx**, retrieved 2026-09-10. SHA256 `3439f6eff1026215b85f0b4a159d521211fa853ff40aaa15114b5fbbffad11d3`. The workbook contains 5,760 interval rows in sheet `Yeast TMDs`. This is the final source-study annotation list; it is not an outcome-filtered client list. The source Methods merge SGD TMHMM annotations dated 2017-08-23 with reviewed UniProt intervals, preferring UniProt for overlaps and excluding mitochondrial genes, dubious ORFs and pseudogenes. The historical UniProt release and per-row originating method are not completely specified. Preserve this limitation rather than inventing a source release.

The deposited start/end fields are used as zero-based half-open amino-acid intervals [s,e), checked directly against the supplied TMD peptide and upstream/downstream sequences at those exact positions. The one-based inclusive amino-acid interval is [s+1,e]. No coordinate shift, local sequence alignment, profile fitting or boundary inference is allowed. A mismatch excludes the whole gene to prevent an incomplete pseudo exclusion mask. These exact sequence checks support the interval semantics for the retained population; the workbook does not explicitly label its index convention in its header.

Independently require the complete translated frozen CDS, excluding terminal stop, to equal the canonical reviewed UniProt S288C sequence from **release 2026_03, released 2026-09-02**. Retain the sequence-service payload, HTTP release headers, entry version, sequence version and entry modification date. The sequence payload SHA256 is `f8e28e2c5ecc9ba9f7c7f20d8971efd825dab0f3e8ec6d6941ec6f3f3433b961`. UniProt is used here for exact protein identity only; no new TMD features are merged into the source-study annotation.

Use the supplement accession when present; if absent, require exactly one reviewed canonical UniProt entry with the identical ordered yeast locus name. Verify that locus and full sequence. The original accession field is retained separately from the resolved accession. This identity bridge does not establish the complete historical source-study protein release, but each retained TMD peptide and flank is also checked directly against the frozen sequence.

## Eligibility and anchors

Anchor is the first amino acid of each annotated TMD: codon index s, CDS nucleotide index 3*s. Genomic coordinates are obtained from the frozen spliced coding-strand exon map, with minus-strand order reversed. All eligible TMDs enter together. FIRST/INTERNAL labels are metadata only.

Require one exact nuclear coding model, exact canonical protein/locus identity, all source TMD intervals and supplied flanks matching the frozen sequence, no overlapping source TMD intervals, and full inclusive [-60,+100] visualization-window support. Every base of this window must have unambiguous same-strand canonical CDS membership. Opposite-strand CDSs do not create ambiguity for a strand-specific assignment. The translated protein excludes stop, so neither start-edge padding nor stop-containing windows are allowed. No read coverage, expression, enrichment, profile quality or published-effect filter is applied.

The pseudo contract additionally requires five valid same-CDS positions. Every source row, source coordinate, exact-identity check and exclusion is recorded in the two annotation TSVs. Failure reasons are staged: annotation failures prevent later window/matching checks; no unperformed downstream check is called passed. Gene-level annotation failures apply to all its anchors.

Frozen annotation-only accounting:

| Category | Unique real anchors |
|---|---:|
| Source intervals | 5,760 |
| Eligible after all annotation and pseudo rules | 1,068 |
| Excluded: insufficient full window | 1,795 |
| Excluded: fewer than five valid pseudo positions | 2,865 |
| Excluded: ambiguous same-strand CDS window | 12 |
| Excluded: gene-level protein/peptide/flank discrepancies | 20 |

The retained anchors span **338 genes** and have **5,340 pseudo definitions**. The 20 gene-annotation-excluded rows can have multiple discrepancy labels; those labels must not be summed as disjoint counts. Excluded genes/anchors are not described as biologically negative. Long CDSs with sufficient non-TMD regions are preferentially retained by design; conclusions are conditional on this subset, not all membrane proteins.

## Locked benchmark choices

- Primary window [-30,+60] codons; baseline [-30,0] and emergence [+30,+60]. Visualization only [-60,+100], unsmoothed.
- Five nearest valid same-CDS pseudo anchors, 60-codon buffer beyond all TMD intervals, relative-position tolerance 0.20, deterministic seeded SHA256 distance-tie ordering, seed 6101207.
- Fixed count-depth normalization, 0.5 raw-read pseudocount, log2 pulldown/input ratio, no smoothing or outcome filtering.
- Primary statistic: equal-gene mean of real minus mean-five-pseudo downstream-minus-baseline log-ratio contrasts.
- Whole-gene paired bootstrap, 10,000 draws, PCG64 seed 6101208; one-sided 95% lower percentile bound and its prespecified >0 directional-support rule. This quantifies conditional gene variation, not biological replication.
- The 32-nt sensitivity repeats the exact anchors, matches, windows, formula and resampling indices. It cannot replace or rescue the 28-nt primary.

Exact equations, zero handling, aggregation order, inference limitations and future technical stop conditions are specified in B2_PRIMARY_ESTIMAND_CONTRACT_v0.md. No biological statistic or profile has been computed. The design is executable without further discretionary analytical choices once separately adjudicated.

## Audit of implementation corrections and boundaries

The initial annotation-only implementation mistakenly queried UniProt alternative ORF names instead of ordered locus names and treated blank supplement accession cells as conflicting identities. The corrected lookup adds the ordered-locus field and resolves blank accessions uniquely by locus. No reference sequence, TMD coordinate, numerical eligibility/window/pseudo parameter, seed or biological formula changed in response. The initial and corrected sequence payloads are both retained. An unavailable optional sequence library was replaced with explicit standard nuclear-code translation; this did not change the scientific rules. These are transparent implementation corrections before design freeze, not biological-outcome tuning.

The original Downloads handoff path was unavailable; the relocated authoritative copy was found and read in full. Submitted P3 work remains unchanged and closed. This new work is development/method evaluation with known literature exposure. Do not execute TMD profiles, biological contrasts, bootstrap tests or FIRST/INTERNAL analyses until scientific adjudication.
