# Region-composite analysis amendment v1.1

Date recorded: 2026-10-02, after execution of the v1 complete-case analysis.

The v1 contract required a TMD unit to be eligible at every center in a frozen
region. After reveal, the complete-case restriction was found to alter the
point estimates at +60 and +75 materially relative to the previously reported
center-specific populations. The complete-case results remain preserved and
will be reported as a population sensitivity; they are not discarded.

To distinguish a region-average effect from a complete-case population effect,
an additional `AVAILABLE_CASE` analysis will be run. It will:

- retain the original frozen eligible population separately at each center;
- use the union of uniquely assigned genes across those center populations as
  the gene-resampling universe;
- apply each bootstrap gene multiplicity jointly to every available TMD record
  for that gene across classes, centers, compartments, and replicates;
- calculate the same center-specific FIRST-minus-INTERNAL and
  MITO-minus-ER estimands before averaging over the frozen region;
- use the same 10,000 resamples and seeds as v1;
- repeat the same leave-one-gene-out influence diagnostics.

This amendment was motivated by the revealed v1 population shift and is itself
post-reveal. It does not rescue, validate, or confirm the topology association.
For manuscript interpretation, the available-case analysis is the more direct
region summary of the already reported pointwise estimands; the complete-case
analysis tests sensitivity to a fixed TMD population across centers.

