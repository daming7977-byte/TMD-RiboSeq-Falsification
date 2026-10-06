# GSE300977 post-reveal region-level topology interaction analysis v1

Date frozen: 2026-10-02

## Status and evidential role

This is an outcome-informed, post-reveal region-level robustness analysis.
The +90/+105/+120 pattern was known before this contract was written. The
analysis is not confirmatory, is not an independent validation, and does not
replace the prespecified +75-nt primary analysis or the prespecified
FIRST-versus-INTERNAL secondary analysis.

## Purpose

1. Summarize the correlated +90, +105, and +120 pointwise topology-by-
   compartment estimates as one observed-plateau statistic.
2. Test sensitivity to a broader +60, +75, +90, +105, and +120 summary.
3. Address dependence among multiple TMD units assigned to the same gene by
   resampling whole genes.
4. Quantify single-gene influence without treating TMD units as biological
   replicates.

## Frozen regions

- `PLATEAU_90_120`: centers +90, +105, +120 nt. This region was selected after
  the pointwise results were known.
- `BROAD_60_120`: centers +60, +75, +90, +105, +120 nt. This is a broader
  sensitivity region and does not replace the plateau summary.

The centers refer to the existing 31-nt center-specific windows. A region
statistic is the unweighted mean of the already defined center-specific
estimands. It is not a newly counted continuous nucleotide window.

## Analytic population

- Reuse the frozen structural, Crossmap, real-anchor, pseudo-anchor, masking,
  read-length, library, and center definitions.
- Restrict to units for which every contributing target transcript has frozen
  `pass_class = MULTI_PASS`.
- Retain only consensus `FIRST` and `INTERNAL` units.
- Require a unique gene assignment for gene-clustered resampling.
- Within each region, require the same canonical geometry-resolved TMD unit to
  be eligible at every included center. No read-count threshold is added.
- A canonical unit is identified by chromosome, strand, real TMD anchor,
  frozen target-transcript set, unique gene, position class, and TMD-index set.

## Estimand preservation

The frozen pointwise estimand first aggregates real and matched-pseudo
densities within a TMD class and then computes:

`D = log2(real_PD / real_input) - log2(pseudo_PD / pseudo_input)`.

The region analysis preserves that order of operations. It does not compute a
per-TMD log ratio, because doing so would require a new pseudocount or
library-dependent exclusion for units with zero density and would therefore
change the estimand.

For center `c`, location `l`, and biological replicate `r`:

`FI(c,l,r) = D_FIRST(c,l,r) - D_INTERNAL(c,l,r)`.

For a frozen region `R`:

`FI(R,l,r) = mean_c_in_R FI(c,l,r)`.

The region-level interaction is:

`I(R,r) = FI(R,mitochondria,r) - FI(R,ER,r)`.

The displayed overall point estimate is the arithmetic mean of the two
replicate-specific interactions. Cytosol is not part of the frozen primary
region interaction.

## Whole-gene bootstrap

- Sampling unit: unique gene.
- Number of resamples: 10,000.
- Seeds: 30097890 for `PLATEAU_90_120`; 30097860 for `BROAD_60_120`.
- Sample the observed gene universe with replacement.
- A sampled gene multiplicity is applied jointly to all of its FIRST and
  INTERNAL TMD units, all region centers, ER and mitochondria, and both
  biological replicates.
- Recompute the complete density-aggregation estimand in every resample.
- Report percentile 95% intervals. No confirmatory P value is assigned.

Gene-clustered resampling addresses within-gene dependence but does not create
additional biological replication.

## Influence analysis

For each region, delete one gene at a time and recompute the region-level
interaction. Report:

- absolute and relative change from the full-data point estimate;
- sign reversal status;
- genes with the largest absolute deletion changes;
- each gene's retained FIRST, INTERNAL, and total TMD-unit counts;
- Spearman correlation between total TMD-unit count and absolute deletion
  change;
- distribution quantiles of signed and absolute deletion changes;
- fraction of total absolute deletion change contributed by the top decile of
  genes ranked by TMD-unit count.

Deletion-change summaries are influence diagnostics, not independent effects
or biological-replicate estimates.

## Interpretation constraints

Allowed: `post-reveal`, `region-level robustness`, `replicate-concordant` (only
if both replicate-specific interactions have the same sign), `internally
stable` (only if supported by both regions and influence diagnostics), and
`gene-clustered interval`.

Not allowed: `validated`, `independently replicated`, `confirmatory`,
`mechanistic`, `causal`, `unbiased confirmation`, or `reproducible topology
program`.

