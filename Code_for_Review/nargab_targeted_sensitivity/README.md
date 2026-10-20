> Current release package: **v1.0.1**, preserving scientific version **1.0.1-rc.2**. Read the package-root README first. Historical delivery/version, original-byte, publication-state and path statements below describe the preserved reviewer history. The release identity manifest separates historical candidate hashes from current release hashes. Preparation did not publish files or rerun scientific analyses.

# Targeted within-gene and composition sensitivity: reviewer rerun

Scientific candidate 1.0.1-rc.2; unpublished delivery 1.0.1-rc.2+sensitivity.20261016.1. This extension corresponds to Supplementary Table S13. It preserves the frozen equal-TMD-unit region estimates and all original figures. The alternative estimands are post-reveal exploratory sensitivity analyses, not independent biological replication, confirmatory tests, or causal mechanism evidence.

## Inputs, software and execution

All seventeen required processed inputs are included in this Code_for_Review archive. Sixteen are reused at their existing paths; one annotation-derived genomic anchor table is added. No other archive, raw sequencing download, alignment, network connection, or original P3 statistical function is required. INPUT_MANIFEST.tsv records exact paths, sizes and checksums; INPUT_PROVENANCE_AND_RIGHTS.md states the annotation and licensing boundaries.

Reference environment: Python 3.12.14, NumPy 2.3.5, pandas 2.2.3, macOS arm64 Accelerate. No SciPy is needed. requirements.txt pins the two Python packages. Creating this environment on another platform does not establish byte identity or biological inferential validity. The runner sets single-thread limits for the calculations. From the extracted Code_for_Review directory, with that environment active:

```bash
python nargab_targeted_sensitivity/run_all.py --output-dir ../targeted_sensitivity_reproduction
```

Choose a new, nonexistent output directory. The command refuses overwrites, checks input/extension hashes, reconstructs 29,927 eligible center-units and their annotations, and executes the 15 prespecified variants in each of two frozen regions (10,000 whole-gene draws per variant). It writes only to the explicit output directory. It never invokes the original statistical functions or historical figure runner. It preserves all seeds, mathematical functions, windows, exclusions, balance strata, positivity rules and influence tie rules of the independent audit. The only adaptations are portable file locations, output-directory controls and packaged checksum/result verification.

The expected results, including all compartment/replicate estimates and percentile/basic intervals, are in expected/calculation_results.tsv. expected/reported_results.tsv and expected/reported_eligibility_and_balance.tsv additionally retain the completed audit's descriptive reporting flags and balance-membership columns. Raw eligibility identities are compared to expected/eligibility_and_balance.tsv; the full weight table is compared separately. Every calculation result field is compared with absolute tolerance 1e-10 and all text/identity fields exactly. Gene eligibility and eleven diagnostic/weight/influence tables are also compared. Metadata-output hashes must match; all 30 saved bootstrap arrays have reference hashes, whose byte identity is reported separately. This tolerance comparison is not a cross-platform test. Success produces REVIEWER_REPRODUCTION_CHECK.json, detailed logs and OUTPUT_SHA256.tsv. A mismatch fails without adjusting thresholds or recalculating with a selected specification. To repeat verification of a completed output directory without new calculation, add --verify-existing; that mode is explicitly recorded in its verification report.

## Interpretation and historical protocol

The original protocol and freeze metadata are retained as explicitly path-normalized documentary copies, with original freeze hashes retained for provenance; the executable calculation path has no dependency on those locations. Protocol SHA-256: e573192d8ea67a56da558420efa6086dbabf13f8a20212e2bd505730eedf6961. It was recorded before the new sensitivity calculations, but after the original manuscript and frozen region outcomes were known. The pre-result implementation conventions preserve the exact seed and boundary definitions.

The frozen estimates are 0.9325467977 and 0.6518202468 log2 units. Equal-gene, shared-population, within-gene and composition-standardized quantities answer distinct questions. Extended standardization gives approximately 0.590465 and 0.408269, reductions of 36.7% and 37.4%. Direct within-gene log contrasts require positive densities in every required library/class/center; only the same 17 genes qualify. Do not generalize that sparse population to all 1,750/1,848 original genes. The within-gene median's interaction_mean field is the median of gene-specific replicate-averaged interactions, not an algebraic combination of its marginal compartment medians.

Conditional whole-gene intervals condition on the observed libraries, eligible population and fixed balancing design. They do not create biological replicates or establish biological confidence coverage. Residual continuous composition imbalance remains, baseline input coverage is coupled to outcome denominators, and matching replicate labels does not establish culture/batch pairing. No new search, significance-driven choice or raw-to-processed validation is performed.

## Licensing

Original new scripts fall within the recorded author-direction MIT code scope and are listed individually in the package FILE_LEVEL_CODE_LICENSE_MANIFEST.tsv. Annotation/data/outputs, original manuscript materials, third-party code and models are not blanket-relicensed. Preserve LICENSE_SCOPE.md and existing third-party notices. Release preparation did not assign a DOI; the v1.0.0 DOI does not identify this v1.0.1 package.
