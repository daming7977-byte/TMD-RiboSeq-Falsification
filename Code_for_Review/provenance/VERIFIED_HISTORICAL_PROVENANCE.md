# Verified historical runner and freeze provenance

This file documents a bounded comparison for local package 1.0.1-rc.2+public-candidate.20261019.2. It does not certify the entire current candidate as identical to public v1.0.0. The current candidate includes later S16/S13 and S17/S14 materials and other recorded revision changes.

## Five Route A/B2 runner and helper files

The five files in `../historical_extension_runners/` were compared directly with surviving historical source files, their original SHA256SUMS records, and the v1.0.0 Git blobs at commit `b0fbe8bf091ab45635096e63c34319ef0d2fb1a3`. All five current copies are byte-identical to their v1.0.0 counterparts. Four are also byte-identical to their original sources. The sole difference in `ROUTE_A_EXECUTION_20260910/run_benchmark.py` replaces its machine-specific absolute Node executable locator with `NODE='node'`. All other source bytes are unchanged; parsed Python structure is identical after normalizing only that string constant. No computational function, scientific parameter, threshold, eligibility rule, estimand, numerical operation, seed or output definition differs in this five-file comparison.

`HISTORICAL_RUNNER_THREE_WAY_SHA256.tsv` records original, v1.0.0 and current hashes. Original source locators retain the historical directory/filename, without exposing personal machine paths. The Route A execution input manifest records the original run_benchmark.py digest; both B2 numerical freeze JSONs record the execute.py digest. Original execution manifests also bind the helper files. These are retained project records, not externally timestamped pre-execution Git attestations.

## Table 11 contracts, source records and numerical freezes

All 30 rows of `../../Supplementary_Tables/Supplementary_Table_11/SOURCE_SHA_MAP.tsv` were independently reconciled to surviving originals, original-directory SHA256SUMS entries, v1.0.0 Git blobs, and current packaged files. Every source/recorded/public/current digest agrees. No mapped original is missing. The source map itself is byte-identical to v1.0.0. See `HISTORICAL_TABLE11_RECONCILIATION.tsv` for each binding, including the GSE106329 reconstruction-freeze completion record.

This finding applies to those 30 records. It must not be generalized to all historical documentary copies in the release: other S16/S17 provenance records have separately documented privacy normalization, and their original freeze digests identify originals rather than the sanitized documentary bytes.

## Distinction from the five S17 privacy-normalized scripts

The five Route A/B2 files above are not the five S17 files named in the root README. The S17 files are `primary_exclusion.py`, `region_exclusion.py`, `overlap_density_contributions.py`, `build_overlap_inventory.py` and `umi_bridge.py`, under `../nargab_targeted_robustness/original_scripts/`. Their existing source identities and deployment adapters have a separate provenance record in that directory. No Route A/B2 equivalence finding is transferred to S17, and no S17 source or adapter was changed in this documentation update.

## Limits of the evidence

Historical pre-execution Git commits were not independently established. The v1.0.0 public repository starts with a clean release commit, so it binds publication bytes but does not supply the earlier development or execution history. Historical attribution is supported by surviving source files, saved manifests, embedded code digests and execution records. Hash equality alone does not establish the time of outcome exposure.

The Node locator substitution is VERIFIED PATH-ONLY DIFFERENCE. It changes executable resolution to PATH; runtime equivalence under arbitrary Node installations is UNVERIFIED. This review did not run historical scripts or scientific calculations, validate external runtime installations, or establish complete upstream/raw-read reproducibility. Table 11 identity and the five-script comparison do not expand the supported reproduction boundary.

## Manifest scope

Current file bytes are covered by the root `RELEASE_FILE_MANIFEST.tsv`, which excludes only itself. `RELEASE_SOURCE_IDENTITY.tsv` preserves reviewer-source digests and adds predecessor-candidate digests; it excludes itself and the root manifest to avoid a checksum cycle, and retains the prior explicitly omitted stale deposit-metadata row. Active nested manifests are refreshed for documentation changes. Historical source hashes and historical-scope manifests are preserved and are not redefined as checksums of this successor candidate.
