# Pre-result implementation conventions

Recorded before the first new sensitivity calculation; supplements the unchanged hashed protocol without changing its design.

Variant indices (seed = 202610160000 + 100*region_index + variant_index, except original frozen seeds): 0 FROZEN_UNIT_ALL; 1 EQUAL_GENE_ALL; 2 UNIT_SHARED; 3 EQUAL_GENE_SHARED; 4 UNIT_POSITIVE_SHARED; 5 EQUAL_GENE_POSITIVE_SHARED; 6 WITHIN_GENE_MEAN; 7 WITHIN_GENE_MEDIAN; 8 WITHIN_GENE_COMMON_TRANSCRIPT; 9 BALANCE_CORE_RESTRICTED_UNIT; 10 BALANCE_CORE_RESTRICTED_EQUAL_GENE; 11 BALANCE_CORE_WEIGHTED; 12 BALANCE_EXTENDED_RESTRICTED_UNIT; 13 BALANCE_EXTENDED_RESTRICTED_EQUAL_GENE; 14 BALANCE_EXTENDED_WEIGHTED. Region index 1 plateau and 2 broad. FROZEN_UNIT_ALL uses 30097890/30097860.

If pooled weighted variance is zero, SMD=0 for equal means; SMD is infinite/unsupported for unequal means. Variances use population weighted second moments. Positive-support and common-transcript n<30 remain sparse descriptive evidence. An empty branch is not estimable. Core count bin [2,3) includes median 2.5; do not label it exact count two. Baseline input-density balance uses the predeclared mean across eight real/pseudo input entries, without library-depth adjustment; the limitations of such coverage and its outcome-denominator coupling are reported.

Within-gene median column `interaction_mean` denotes median of per-gene replicate-averaged interactions. Its side columns contain marginal medians for compartments/replicates and need not obey mean-estimand identities. The mean and all aggregate-density variants do obey those identities. No estimated biological-replicate interval is added.
