# Route A benchmark design freeze v0

DESIGN FROZEN; BENCHMARK NOT EXECUTED. Scientific adjudication required. Normative components are this file, ROUTE_A_SUPPORT_AND_MISSINGNESS_CONTRACT_v0.md, ROUTE_A_PRIMARY_FALSIFICATION_ESTIMAND_v0.md, ROUTE_A_DOSE_TREND_SECONDARY_v0.md, annotation_rules.json and the final exact real/pseudo tables. SHA256SUMS.txt freezes these files and their source/implementation assets. Any hash/reference mismatch stops execution, without discretionary repair.

## Technical input

Use the unchanged Task12 independent sacCer3/S288C reference and archived UCSC sgdGene/sgdOther2011-10-06 annotation. Task12 contract SHA256 716439411b875322332bd81792cce36bcb7777a58163dea7dec10687ec2e0dc3; common reconstruction contract SHA256 2c4d660ca3b6e64c5e8dc0209ad2848f2d52f7743369232728385880f4e5a7af. Genome SHA256 92634e39d9d7013b52e82de31606fad36da15fed9704cdf1f487e3b29322388b; coding models SHA256 1a8539be61e88ed3f94305e5eb31233d8dce11380ff749be748cf73cf9e93a76. All source references and eight aligned files are identified by the unchanged Task12 manifest. This is not claimed to be the exact historical author reference.

Primary27nt/A-site+14; separate sensitivity28nt/+15. No pooling lengths, excluded-length reintroduction, new trimming, UMI/coordinate deduplication, or condition-specific calibration. Identical Task12 N exclusion, contaminants, exact end-to-end mapping and operational uniqueness. Same-strand unique frame0 CDS A-site assignment, startcodon0, explicit stopslot. This design's gene-wide denominator includes that stopslot; real/pseudo windows exclude stop.

| GSM | Run | Condition |
|---|---|---|
| GSM1439584 | SRR1520311 | noCHX |
| GSM1439585 | SRR1520312 | 1/64x |
| GSM1439586 | SRR1520313 | 1/16x |
| GSM1439587 | SRR1520314 | 1/4x |
| GSM1439588 | SRR1520315 | 1x |
| GSM1439589 | SRR1520316 | 8x |
| GSM1439590 | SRR1520317 | 100x rep1 |
| GSM1439591 | SRR1520318 | 100x rep2 |

## Annotation and coordinate identity

Shurtleff et al.2018, eLife7:e37018, Supplementary file2 DOI10.7554/eLife.37018.013, publisher asset elife-37018-supp2-v2.xlsx, SHA256 3439f6eff1026215b85f0b4a159d521211fa853ff40aaa15114b5fbbffad11d3. Exact copied source payload, not B2 result data. Rechecked independently against Route-A genome/CDS. Its5760 intervals in Yeast TMDs are interpreted as zero-based half-open [s,e) only when supplied peptides/flanks agree at the declared positions. One-based amino-acid coordinates are[s+1,e]; anchor at the first residue of each TMD is codon s/CDS nucleotide3s. No boundary shifts, local fitting or inferred TMDs. Full canonical protein identity must agree with reviewed UniProt/S288C release2026_03 (2026-09-02), exact retained payload SHA256 f8e28e2c5ecc9ba9f7c7f20d8971efd825dab0f3e8ec6d6941ec6f3f3433b961; this is an identity bridge, not a new TMD feature source. Retain per-entry versions and deposited/resolved accessions. Missing supplement accessions require a unique reviewed ordered-locus match. Historical per-row annotation provenance remains incompletely documented; no release is invented.

Every interval/peptide/flank in a gene must pass; any failed interval excludes all anchors in that gene, preserving a complete pseudo exclusion mask. Require unique nuclear canonical CDS, no overlapping source TMD intervals, and full inclusive[-60,+100] translated-CDS support with no ambiguous same-strand coding base. No FIRST/INTERNAL metadata or subgroup analysis is used. All eligible TMDs enter together.

Annotation accounting:5760source intervals;1068eligible real anchors in338genes;1795excluded for incomplete full window,2865for fewer than5pseudo candidates,12for ambiguous sense-CDS windows,20for gene annotation discrepancies (labels can overlap within those20). Exact row reasons are in the provenance/eligibility tables. Long-CDS and available non-TMD-region selection is a design limitation; this population is not all membrane proteins.

## Exact pseudo definitions and geometry

For each real anchor a in a protein of L translated residues excluding stop, enumerate integer q with full[-60,+100] support and unambiguous same-strand windows. Exclude q inside inclusive[s-60,e-1+60] around EVERY annotated TMD[s,e). Thus the nearest allowed residue is beyond the60-codon inclusive buffer. Require5*abs(q-a)<=L-1, exactly relative-position distance<=0.20. Relative coordinate is q/(L-1). Sort by abs(q-a), then lowercase hexadecimal SHA256 of UTF8 '5957313|anchor_id|q', then numeric q. anchor_id is ASCII gene_locus|source_start_0based|source_end_exclusive; q is unpadded decimal. Take5distinct candidates, without replacement within a set; reuse across real anchors allowed. Order anchors by ASCII locus,numeric start,end. Collapse exact source-interval duplicates retaining source-row provenance. No fallback or rematching. Final table contains5340pseudo definitions and is authoritative; do not regenerate against a changing service.

Formal window[-30,+60]; baseline[-30,0]; downstream[+30,+60], both31codons. Visualization only[-60,+100]. These reuse the pre-adjudicated geometry and cotranslational emergence rationale; exact bands are fixed design choices, not outcome-fitted regions. Pseudo formal windows are TMD-free under this annotation; broader visualization windows may intersect TMDs. No scan, smoothing or optimization.

## Estimands, support and seeds

Normalized codon density K*C/N, zero-gene total non-evaluable, no pseudocount/log. Downstream-minus-baseline real contrast minus mean5pseudo contrasts, equal anchors within gene, equal genes. Primary averages the two100x experiments equally per gene then subtracts noCHX, two-sided percentile bootstrap CI. Secondary is mean per-gene OLS slope over ranks0..6 with100x equally collapsed at gene level; never replaces primary. Exact equations and degeneracy rules are in the normative contracts.

Annotation population1068anchors/338genes. Primary common support: 304genes/946anchors; seven-dose support: 280genes/879anchors. Both footprints use the identical support set for each estimand, requiring positive gene-wide counts in both lengths. Membership and per-condition evaluable-anchor totals are frozen in SUPPORT_MEMBERSHIP_FROZEN.json and the real-anchor table. This task inspected only gene-wide technical support totals, never TMD-relative counts, local occupancies, contrasts or bootstrap outcomes.

Matching seed5957313; primary bootstrap5957314; ordinal bootstrap5957315. Each bootstrap10000whole-gene paired draws using NumPy PCG64 and linear2.5/97.5percentiles; same draws reused across footprints for each estimand. No biological replication claim. The fixed first execution is retained, including null/negative or discordant results. No relationship with B2 may motivate changes. No discretionary analytical choices remain for the declared executable support populations. Stop for adjudication before execution.
