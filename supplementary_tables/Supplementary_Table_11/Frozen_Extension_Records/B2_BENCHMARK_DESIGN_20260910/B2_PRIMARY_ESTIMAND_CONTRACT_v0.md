# B2 primary specificity estimand v0

Status: formula and inference frozen; **not executed**. This is one directional local-structure estimand, not an omnibus test for every possible TMD-centered shape.

## Coordinate and read reconstruction

Primary is exactly 28 nt with A-site offset +15 nt from the mapped 5-prime end in coding direction. The separate sensitivity is exactly 32 nt with +16 nt. No other lengths enter either population and they are never combined. Apply Task 6 exact alignment, trimming, contaminant, multimapper and duplicate rules and SHA-frozen references. Use one reported exact alignment, full-length M CIGAR and NM=0; reject multiple reported placements. Determine the A-site genomic base by adding the offset to the plus-strand leftmost base or subtracting it from the minus-strand rightmost base. Require exactly one retained nuclear sense CDS assignment at frame 0. Off-frame, ambiguous and outside-CDS reads are not rounded or allocated fractionally.

CDS bases follow spliced coding-strand order. CDS nt 0 is the first start-codon base, codon 0 the start. For L amino acids, stop is codon L; record it separately and exclude it from normalization and analysis windows. No clipping, padding, flank inference or hidden start/stop deletion. All fixed windows fit actual translated CDS by annotation eligibility.

Sum run counts within GSM before normalization: SRR1562895+SRR1562896 form input GSM1495513, and SRR1562897+SRR1562898 form pulldown GSM1495514. Retain sequence and coordinate duplicates. UMI_STATUS=NOT_DOCUMENTED; UMI_DEDUPLICATION=OFF. There are two fractions and no established biological replication.

## One profile formula

For each length population separately, let C_l(g,j) be the integer GSM count at coding codon j. Let N_l be total assigned in-frame counts over all unambiguous retained canonical nuclear translated CDS codons (including start, excluding stop), including genes without TMDs. Count each read once. N is not a total over overlapping anchor windows. Record N for both libraries and lengths. If either library has N=0, that length is TECHNICALLY_NOT_EVALUABLE; do not substitute another population or alter the anchors.

Define depth-normalized pseudocount density and positional log ratio:

`D_l(g,j) = 1,000,000 * (C_l(g,j) + 0.5) / N_l`

`E(g,j) = log2(D_pulldown(g,j) / D_input(g,j))`.

The fixed pseudocount is **0.5 raw reads per codon**, added before depth normalization. This is a regularized density, not a claim of molecule abundance. No smoothing, winsorization, gene-wise rescaling, expression cutoff or exclusion of zero windows. Zero/zero codons give the library-size constant log2(N_input/N_pulldown), which cancels in the within-anchor band difference. Genes with no counts remain present; all-zero paired windows contribute zero structure. Do not treat unavailable/corrupt data as zeros. Any missing required technical input stops that population.

## Primary window and scalar statistic

Anchor is the first annotated TMD amino acid. Primary support is inclusive codons **[-30,+60]** relative to the anchor. Baseline B=[-30,0] and downstream emergence band D=[+30,+60], each 31 codons. Positions +1..+29 are deliberately outside the contrast; there is no scan or strongest-center selection.

For any real or pseudo position a:

`Z(g,a) = mean_{r=30..60} E(g,a+r) - mean_{r=-30..0} E(g,a+r)`.

For real anchor a and its five frozen pseudo positions p(a,k):

`d(g,a) = Z(g,a) - (1/5) * sum_{k=1..5} Z(g,p(a,k))`.

Let A_g be all eligible real anchors in gene g. The sole primary statistic is:

`T = (1/G) * sum_g [ (1/|A_g|) * sum_{a in A_g} d(g,a) ]`.

Direction is **T>0**: downstream-minus-baseline positional log-ratio structure is greater at real TMD coordinates than at their matched non-TMD positions. Negative values must be retained and reported. This does not test for an arbitrary peak, total ER recruitment, causal TMD action or FIRST-versus-INTERNAL differences.

Aggregation order is codon log ratios -> band means -> paired real-minus-five-pseudo contrast -> within-gene anchor mean -> equal-gene mean. Equal weighting of genes prevents genes with many TMDs dominating by anchor count. Neither anchors nor genes are biological replicates.

## Conditional gene resampling

After separately authorized execution, calculate 10,000 whole-gene paired bootstrap draws. Sort the 338 contributing genes by ASCII locus. Use NumPy **2.3.5**, `Generator(PCG64(6101208))`, drawing `integers(0,G,size=(10000,G),dtype=int64)` in row-major order. Each draw keeps every selected gene's real anchors, five pseudo matches, input/pulldown pair and all codons together. Densities and N remain fixed at the observed library values; only gene contributions are resampled. Do not independently resample anchors, reads or technical SRRs.

The primary inferential summary is the one-sided 95% lower percentile bound `quantile(T_boot,0.05,method='linear')`. The frozen directional-support rule is lower bound >0. This is conditional gene-resampling support in these libraries, not population-level biological replication, causal randomization inference, or a GO/HOLD/NO-GO classification. No permutation p-value is claimed because real and matched pseudo labels were not experimentally randomized. Fewer than 20 eligible genes or any nonfinite bootstrap statistic triggers technical non-evaluability with no redesign. No bootstrap or statistic has been calculated in Task 7.

## Visualization and sensitivity

The only broader visualization window is inclusive [-60,+100] codons, unsmoothed. At each r, average E across real anchors within gene, then equally across genes. For pseudo, first average the five matches per real anchor, then anchors within gene, then genes. Display real and pseudo curves together. No peak finding, center scan, pointwise significance bands or FIRST/INTERNAL outcome panels.

Repeat the exact population, pseudo table, formula and bootstrap gene-index draws for the 32-nt sensitivity. Do not drop anchors because 32-nt counts are low or zero. If technically unavailable, report NOT_EXECUTED/TECHNICALLY_NOT_EVALUABLE and why. Agreement or disagreement is descriptive technical robustness; 32 nt cannot rescue the primary or supply an additional biological replicate.

## Rationale and limitations

Codon resolution is supported by Task 6 frame/site calibration. Exit-tunnel geometry commonly occludes roughly 30 residues and can vary with conformation, motivating a broad fixed +30..+60 interval after the first TMD residue rather than a fitted peak. Cotranslational NMR work discusses this geometry (GEOMETRY). The original proximity-Ribo study aligns to the first codon of a hydrophobic feature (JAN2014, Figure 4 caption), providing an anchor precedent, not the numerical window or expected B2 effect. The exact bands are independent design choices and do not assume every TMD has completed emergence at a particular codon.

Nearest same-CDS matching controls gene identity and bounds position differences. It does not match sequence, mappability, local elongation, TMD length or every nuisance factor. Annotation incompleteness, fixed pseudocount sensitivity and window-overlap correlation remain limitations. Any future redesign requires a separately labeled post-outcome contract; it cannot replace this primary.
