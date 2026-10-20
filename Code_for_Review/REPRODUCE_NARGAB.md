> Current release package: **v1.0.1**, preserving scientific version **1.0.1-rc.2**. Read the package-root README first. Historical delivery/version, original-byte, publication-state and path statements below describe the preserved reviewer history. The release identity manifest separates historical candidate hashes from current release hashes. Preparation did not publish files or rerun scientific analyses.

# Reproduce the controlled NARGAB revision

## Scope

This controlled candidate retains the post-reveal region-composite calculation, with the two ancillary influence-summary fields revised under the Phase 2 policy; it verifies the primary +75-nt statistics, and rebuilds main Figures 1, 2 and 4. This candidate contains the historical Figure 3 SVG input and a separate terminology-only presentation helper; this command neither regenerates nor copies Figure 3, and no standalone final Figure 3 export is included here. The primary and benchmark statistics are read/checked as saved outputs, not recomputed by this command. The package does not perform discovery, choose windows or reinterpret P values.

The included 90-MB `frozen_intermediates/region_composite_project/` directory contains the exact processed inputs required by the region-composite script. Raw FASTQ-to-alignment workflows remain documented historical workflows and are not silently represented as a one-command raw-data reproduction.

Known inherited-code limitation: `scripts/TMD-RiboSeq/R/fit_target_background_sensitivity.R` and `scripts/TMD-RiboSeq/src/tmdriboseq/r_scripts/fit_target_background_sensitivity.R` are retained unchanged for provenance. Their `*_target_minus_background_log2_effect` columns contain natural-log-scale effects, despite the log2 labels. This routine did not generate the reported P3 results; GSE133111 used `fit_paper3_stage1_onefactor.R`. The legacy `tmdriboseq sensitivity` command and general `tmdriboseq run` workflow can invoke it, with the latter enabling sensitivity by default. Do not interpret these legacy effect columns as log2 values. The provenance manifests record inclusion, not execution.

## Environment

- Python 3.12.14
- NumPy 2.3.5
- pandas 2.2.3
- matplotlib 3.11.2 for figure rendering
- macOS Accelerate BLAS/LAPACK (reference environment for byte-identical regenerated TSV files)

Create the pinned environment with either:

```bash
conda env create -f environment-nargab.yml
conda activate p3-nargab-2026
```

or install `requirements-nargab-lock.txt` into a Python 3.12 virtual environment. The conda specification is the byte-identical reference environment on macOS arm64. Other operating systems or BLAS implementations are not verified byte-identical reproduction environments. The verifier stops on any TSV hash difference and has no implemented numeric-tolerance fallback. A pip installation alone does not establish equivalence to the recorded macOS arm64 Accelerate environment.

## One-command rerun

Extract this candidate with its three sibling directories intact: `Code_for_Review/`, `Source_Data/`, and `Supplementary_Tables/`. Run the following command from `Code_for_Review/`. The lower-case `code/` and `source_data/` layout belongs to the historical v1.0.0 release, not this candidate:

```bash
python reproduce_nargab.py
```

The full command, not rerun during Phase 2:

1. runs the available-case region analysis with seeds 30097890 and 30097860;
2. runs the complete-case population sensitivity;
3. checks the five TSV outputs from each mode byte-for-byte against the distributed frozen tables;
4. verifies the primary and region statistics used in the manuscript; and
5. rebuilds Figures 1, 2 and 4 from packaged source tables.

## Fast verification without recomputation

```bash
python verify_nargab_outputs.py
```

## Expected principal values

- Primary +75 nt: ER 0.089319 [−0.018090, 0.221255]; mitochondria 0.042555 [−0.168336, 0.261755]; cytosol 0.084269 [−0.010151, 0.171143].
- Available-case +90/+105/+120 interaction: 0.9325468 [0.7122142, 1.1328309]; replicate-label-matched estimates 1.0736264 and 0.7914672.
- Available-case +60/+75/+90/+105/+120 interaction: 0.6518202 [0.3573579, 1.0554072]; replicate-label-matched estimates 0.7840298 and 0.5196107.

These region summaries are outcome-informed post-reveal sensitivity analyses. The gene bootstrap addresses within-gene dependence but does not create biological replication.

## Output locations and limits

Run in a disposable extracted copy. Region TSVs are written to `Code_for_Review/reproduction_output/`. Figures 1, 2 and 4 overwrite the parent Figures directory; they are not written inside reproduction_output. For Figure 3, see `presentation_v5/README.md`: the included SVG is a historical input, and its separate helper applies the recorded terminology-only label change without recalculating data. A final PDF/PNG export is not bundled or produced by the main command. The approved submission Figure 3 remains separate. The main runner's historical completion message about frozen Figure 3 is not a claim that it creates or copies that artwork. Historical extension scripts have separate, incomplete upstream environment bindings, documented with Supplementary Table 11 and `provenance/VERIFIED_HISTORICAL_PROVENANCE.md`; they are not invoked by this command. Phase 2 regenerated only the eight ancillary fields from saved gene-deletion effects. The prior ten-file byte-identity test belongs to the frozen v1.0.0 run and does not certify the changed candidate. A full bootstrap/figure rerun was not performed in Phase 2.


## Phase 2 summary-only command

From code/ or Code_for_Review/: `python update_influence_summaries_phase2.py --output-dir ../phase2_summary_check`.

This uses only bundled gene-deletion effects and writes a separate output directory. The policy and exact tie rules are in nargab_region_composite/PHASE2_INFLUENCE_SUMMARY_AMENDMENT.md. The full runner calls the same helper. Primary/region values, populations, gene-deletion effects, maxima and CIs are not recomputed by this command. Cross-implementation summary tolerance is 1e-10; the historical full verifier still performs strict byte checks when given full-rerun directories and has no numeric fallback. Do not interpret a summary-only pass as a new full-rerun pass.

This tree is prepared as public release v1.0.1. The baseline v1.0.0 archive identifiers do not identify this changed tree; publication status must be checked on the platforms. Current conditional primary sensitivity evidence is in Supplementary Table S1; neither zero-excluding nor zero-crossing alternative is promoted to a replacement primary result.

## Added targeted sensitivity extension (scientific candidate 1.0.1-rc.2)

The post-reveal exploratory within-gene and composition sensitivities are a separate calculation path in `nargab_targeted_sensitivity/`. See that directory's README.md and Supplementary Table S13. The extension reconstructs unit densities from the packaged verified processed checkpoints; it does not replace frozen estimands, perform raw alignment, establish biological replication or validate biological confidence coverage. Its output directory is explicit and separate, and its command does not run the historical figure workflow. The public v1.0.0 repository/DOI baseline does not identify this v1.0.1 sensitivity-integrated tree.

## S17 targeted robustness addition

Delivery 1.0.1-rc.2+delivery.20261018.5 adds separate post-reveal annotation-overlap exclusion and matched-population operational UMI-like collapse calculations. Start with `nargab_targeted_robustness/README.md`. Its `--validate-only` mode checks delivered files and schemas without calculating effects; optional point/interval replay writes to a new output directory. The original main runner above does not invoke S17. Saved matrices, population flags and original results support the exact density-level estimands; the complete inventory is in sibling `Source_Data/Targeted_Robustness_S17/`, and the summary is in `Supplementary_Tables/Supplementary_Table_14/`. Regenerating annotation overlaps additionally requires the exact external Ensembl104 GTF. Molecular-UMI identity, source-author deduplication and a comparable read-length sensitivity remain unverified. The new restricted-population primary findings do not replace the frozen primary inference. No public v1.0.0 identifier is asserted to archive this local controlled revision.
