# Supplementary Table S14

Post-reveal annotation-overlap exclusion and operational UMI-like collapse sensitivities. Scientific version 1.0.1-rc.2; delivery 1.0.1-rc.2+delivery.20261018.5. Twelve rows preserve full-precision saved results; no estimates, intervals or tests were recalculated during documentation integration.

`SUPPLEMENTARY_TABLE_S14.tsv` contains five frozen-baseline reconstruction rows (two region interactions and three primary compartments), five overlap-exclusion rows, and two operational-collapse rows. Frozen baseline reconstructions were validated against the original frozen values within absolute 1e-10; insignificant trailing digits can differ from historical summaries. All effects are log2 units. Regional contrasts are mitochondrial-minus-ER differences in FIRST-minus-INTERNAL real-minus-pseudo enrichment. Primary rows are compartment-specific +75-nt enrichment. Conditional 95% intervals resample genes for regions and paired TMD units for primary; all use 10,000 percentile draws. The source tables and S17 document seeds and populations.

`source_output` resolves relative to THIS table directory as `../../Code_for_Review/nargab_targeted_robustness/expected/<source_output>`.

`verification_source_output` retains the original verification A/B path for traceability; `source_row_identifier` identifies the saved source row. The code directory's PACKAGED_SOURCE_IDENTITY.tsv records original paths and source SHA256s. The inventory occurs once at `../../Source_Data/Targeted_Robustness_S17/P3_OVERLAP_INVENTORY.tsv`; supporting masks, population flags, compatibility records and executable reconstruction are in Code for Review.

Primary owner-gene totals are unions of frozen owner IDs (including multi-gene units), not resampling counts. The primary bootstrap uses 9,098 frozen or 8,750 retained paired units. Regional center-specific pair counts differ from union TMD-unit counts because availability varies by center. Exclusion removes whole matched pairs, including any pseudo-overlap case, without partial remasking or redefining the estimand. Different populations answer a diagnostic sensitivity question and do not replace frozen inference. The molecular identity of the UMI-like field and a comparable read-length sensitivity remain unverified.

Primary PD/input contrasts are ER B01/B02 and B11/B12; mitochondria B09/B10 and B07/B08; cytosol B05/B06 and B03/B04. Cross-compartment replicate labels do not establish shared cultures or batches.
