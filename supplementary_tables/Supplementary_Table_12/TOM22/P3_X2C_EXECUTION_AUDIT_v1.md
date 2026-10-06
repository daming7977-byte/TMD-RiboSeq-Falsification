# X2C final frozen biological execution audit

X2C_STATUS = COMPLETE
FROZEN_VERDICT = X2_NO_DIRECTIONAL_SUPPORT
V4_SHA256_VERIFIED = YES
FINAL_PREBIOLOGY_HASHES_VERIFIED = YES
RESULT_DEPENDENT_MODIFICATION = NO
ALTERNATIVE_NORMALIZATION_RUN = NO
CENTER_SELECTION = NO
POPULATION_CHANGED_AFTER_REVEAL = NO
CURRENT_MANUSCRIPT_MODIFIED = NO

## Chronology and authorization
1. v1 failed paired transcript-accession equality before biology.
2. v2 amended the technical gate to six-library exact-protein coordinate compatibility before biology.
3. X2B stopped before count inspection because sparse-position semantics were ambiguous.
4. v3 clarified SPARSE_ZERO_IMPLICIT from author reader code; writer/Script2 recovery limitations remain documented.
5. v4 fixed the common BOTH-class gene population across all 11 centers before reveal. This contract is final; no v5 was created.
6. The present user request first authorized count inspection and biological execution. Pre-execution checks completed at 2026-09-08T04:15:42.964321+00:00.
7. First count-inspection timestamp (UTC, recorded immediately before the first H5 count read): **2026-09-08T04:17:29.420751+00:00**.
8. Biological computation completed at 2026-09-08T04:17:34.213571+00:00.

Two dependency-import attempts stopped before any H5 access or timestamp creation: bundled Python lacked h5py, then SciPy. The existing `/private/tmp/p3_x2_deps` h5py installation was used. Before count access, the Student-t df=2 critical value implementation was changed from unavailable SciPy to the exact analytic inverse CDF sqrt(2)*0.95/sqrt(1-0.95^2); the statistical rule is identical. No biological execution or count access occurred during those import failures. The successful runner was hashed before count access and was unchanged afterward.

## Frozen inputs
The v4 contract was read in full. SHA256 = `59af0746e41920ce72e888598e6bbed2af01dfcecc5ec9cbc48a44c4e920a981`.
Final freeze manifest SHA256 = `b438cea59ee6fbbb5b3c4b085f6bc766f5bd5f86441f0cce15ce61c0b6238155`.
The following hashes were verified before count access and again after execution:
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/01_BASELINE_H5/GSM8076696_Translatome_rep1.h5`: `3dbb124a4a251fe30db67d7c10f2cf758d269e2091568bfe59dd4bc4d5057293` (PASS)
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/01_BASELINE_H5/GSM8076697_TOM_interactome_rep1.h5`: `aaa3e4c36e14a31ab2ea139290105203ef0ceaebd2cc4d9cb61b442299aa8645` (PASS)
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/01_BASELINE_H5/GSM8076698_Translatome_rep2.h5`: `c9d91634f228a8d015ee71a1fec796fb36cd7ee4377cd028c989c6860691b1ab` (PASS)
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/01_BASELINE_H5/GSM8076699_TOM_interactome_rep2.h5`: `46a835f6ba671afd318354f7c7ccb3d2e4a2575bbf49cf54aa724ee4fc0dd3d8` (PASS)
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/01_BASELINE_H5/GSM8076700_Translatome_rep3.h5`: `443223767e3d94aff06f5df9b6f526cc83c7143ad0918db1fb53d83930e468ff` (PASS)
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/01_BASELINE_H5/GSM8076701_TOM_interactome_rep3.h5`: `10b488035a02c3ac8206dc288feba7bc7f9c82119417507179ea75bf67ab9d4f` (PASS)
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/03_CONTRACT/P3_X2_GSE255657_DIRECTIONAL_VALIDATION_CONTRACT_v4.md`: `59af0746e41920ce72e888598e6bbed2af01dfcecc5ec9cbc48a44c4e920a981` (PASS)
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/02_TECHNICAL_AUDIT/P3_X2B1_FIXED_SUPPORT_GENE_MAP_v1.tsv`: `a9349ecacaffb32426cdd89cf6026c0e223f84eba4432dc5ed4606e3fed61c32` (PASS)
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/02_TECHNICAL_AUDIT/P3_X2B1_FIXED_SUPPORT_SUMMARY_v1.tsv`: `b0ddf089088305c46feb0d8e552ed483b5813c61e86c6934e7eec4b73cf9185d` (PASS)
- `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/02_TECHNICAL_AUDIT/P3_X2B1_FIXED_SUPPORT_CENTER_ANCHOR_KEYS_v1.tsv`: `75dbbb95e67d06d73c68dca3469be2d5257c779dda4e2ad55c738bd487b8b41c` (PASS)

## Technical count audit
All six required libraries passed. Counts were real numeric, finite and nonnegative; positions were finite integers, unique and within deposited CDS bounds. Shapes were two-row datasets; position/count lengths matched. No malformed gene or missing frozen gene was found. A full CDS-length zero vector was initialized for each valid gene; stored counts were assigned at position minus one. No clipping, interpolation, smoothing, codon aggregation or pseudocount was used. The entire-library total included all deposited gene datasets, not only analysis genes.

| Library | Genes | Stored entries | Total assigned P-sites | Audit |
|---|---:|---:|---:|---|
| GSM8076696_Translatome_rep1.h5 | 13737 | 1718563 | 3181050 | PASS |
| GSM8076697_TOM_interactome_rep1.h5 | 14062 | 730106 | 1442231 | PASS |
| GSM8076698_Translatome_rep2.h5 | 14105 | 2701197 | 6070916 | PASS |
| GSM8076699_TOM_interactome_rep2.h5 | 13560 | 822455 | 1760884 | PASS |
| GSM8076700_Translatome_rep3.h5 | 14222 | 3089999 | 7431435 | PASS |
| GSM8076701_TOM_interactome_rep3.h5 | 13753 | 1052753 | 2392455 | PASS |

Detailed technical errors table is empty (header only). The audit preserved all 131 genes irrespective of local zero counts. Required observed aggregate means were positive.

## Population, aggregation and normalization
Exactly 131 common-support genes, 131 inherited FIRST anchors and 797 inherited INTERNAL anchors. All 11 centers contain both classes for exactly the same 131 genes. Every original valid center/anchor row was taken from the v4-hashed key table. The other 634 structural candidates were not used. No remapping, annotation change, center-specific gene addition or count-based exclusion occurred.

Each valid 31-nt window was summed over its frozen zero-based half-open boundaries after sparse-zero expansion. Within-gene class anchor sums were averaged, divided by the fixed whole-library total and 31, then averaged equally across genes. Applying the shared linear normalization after within-gene averaging is algebraically identical to normalization per anchor. Library totals were computed once and held fixed during resampling. No alternative normalization was used.

The bootstrap sampled 131 genes with replacement, with all libraries, both classes and all center masks together. Deterministic within-gene class averages were cached once and carried with each sampled gene; recomputing these same averages for each occurrence would give the same values. Across-gene means, D, FI, replicate broad means and the three-replicate mean were recalculated for every resample. NumPy 2.3.5 default_rng PCG64, seed 255657300; 10,000 iterations. NumPy linear percentile endpoints. These implementation choices were recorded in the runner before count inspection. The frozen v4 provides no numerical invalid-fraction tolerance; none was invented. All 10,000 iterations were valid, so the unresolved need for a nonzero tolerance never arose.

## Frozen results
- Replicate 1 broad FI: -0.219734761212424
- Replicate 2 broad FI: 0.0977471452095004
- Replicate 3 broad FI: 0.0455137969965111
- Mean broad FI: -0.0254912730021376
- Gene-bootstrap percentile 95% CI: [-0.300147204376851, 0.155451839850487]
- Valid iterations: 10000; invalid: 0; invalid fraction: 0.
- Descriptive replicate-only Student-t 95% reference (n=3, df=2): [-0.448378638587901, 0.397396092583626]
- Verdict: **X2_NO_DIRECTIONAL_SUPPORT**. Strong support fails; partial support also fails because the mean is negative, despite two positive replicate estimates.

The full 11-center profile and all 10,000 bootstrap draws are retained. No downstream-only statistic, new CI class, P value or maxT inference was calculated. The t reference did not alter the verdict.

## Execution and output verification
Runner SHA256: `0f90c9c74bb70005b93d376f195cd6ef50c62998571aafd51cf8084757664e4a`.
Figure build script SHA256: `ab00a90e2604202a7602979ce2d8c4986b6184acdf3ba7e3ebbeca3e7171c3ac`.
Read-only checks confirmed 11 ordered centers, consistency of saved profile means with saved broad statistics, and 10,000 explicitly valid bootstrap rows. All 69 previously snapshotted non-H5 project files plus manuscript (count includes manuscript) matched their pre-execution hashes; all six input H5 hashes also matched. Manuscript SHA256 remained `cf1eed6a74db873506e7339066275ebd4740c9d6615bf265a650ca4123399522`.

## Figure and visual QA
A single-page five-panel PDF and matching SVG/PNG were produced from saved results only. The rendered PDF/PNG was visually inspected: A-E labels, three replicate markers plus mean, all 11 centers, complete profiles, distinct bootstrap/t uncertainty panels and common D/E scale are readable with no clipping or overlap. No center is privileged and no mechanism cartoon or significance star is shown. SVG and PDF share the exact same ReportLab drawing. PDF text extraction confirmed all panel titles and the no-support verdict. Initial Poppler rendering had a local fontconfig cache problem; rerendering with a temporary font configuration completed successfully and the final image was inspected. No result changed during figure work.

All new artifacts are installed as new files under `HISTORICAL_PROJECT/P3_POSTREVEAL_STRENGTHENING_v1/03_EXTERNAL_VALIDATION/GSE255657/04_RESULTS`. Existing contracts, audits, manuscript, figures and frozen inputs are not overwritten. The output SHA256 manifest lists every new deliverable except itself. This execution is explicitly post-reveal independent TOM22 directional evaluation, not replication of the full MITO_MINUS_ER contrast and not a causal finding.
