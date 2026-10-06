Current remediation scope: see REPRODUCE_NARGAB.md in Code_for_Review/code and the canonical manuscript availability statements. The region command is not universal raw-to-figure reproduction; it rebuilds Figures 1, 2 and 4 only. Author-owned code is MIT; third-party notices and source-data terms remain applicable. Public GitHub/Zenodo access is not verified and requires a separate release action.

# P3 NARGAB code and reproducibility package

This package supports the NARGAB Standard Paper revision. It preserves the historical analysis code and provides a deterministic, pinned rerun for the frozen region-composite analysis, checks of selected saved primary statistics, and Figures 1, 2 and 4 rendering.

Start with `REPRODUCE_NARGAB.md`.

For journal review, extract `Code_for_Review.zip` and `Source_Data.zip` into the same parent directory. The complete release archive uses the equivalent lowercase `code/` and `source_data/` directories; the verification and figure scripts recognize both layouts.

## What the one-command workflow does

Running `python reproduce_nargab.py` in the pinned Python 3.12 environment:

1. reruns the +90/+105/+120 available-case region summary;
2. reruns the +60/+75/+90/+105/+120 available-case sensitivity;
3. reruns both complete-case population sensitivities;
4. retains all records within each resampled gene in the whole-gene bootstrap;
5. checks ten core output tables byte-for-byte against the frozen audited outputs in the reference macOS arm64 Accelerate environment;
6. verifies the primary +75-nt and region-composite values reported in the manuscript; and
7. rebuilds main Figures 1, 2 and 4 from packaged source data. Figure 3 is the unchanged audited primary-result figure.

The workflow does not select windows, add subgroups or datasets, discover new results, reinterpret P values or treat gene resampling as biological replication.

## Inference boundaries

- The +75-nt center and [TMD+60, TMD+91) window are the prespecified primary analysis; the universal enrichment hypothesis was not supported.
- The TMD-order-by-compartment interaction is post-reveal and observational.
- The +90/+105/+120 region is outcome-informed and is used only as a region-level robustness summary.
- The +60/+75/+90/+105/+120 region is the broader sensitivity.
- Available-case is the principal region estimand; complete-case is only a population sensitivity.
- Whole-gene resampling addresses within-gene dependence but does not create additional biological replicates.
- Additional extensions tested distinct specificity/CHX questions, not the full compartment interaction; completed ordinal and TOM22 boundaries are disclosed separately.

## Reproducibility scope

The exact processed inputs required by the region-composite script are included in `frozen_intermediates/region_composite_project/`. Public accessions, source-data tables, seeds, contracts, environments and provenance are supplied. Historical raw-FASTQ-to-alignment workflows are retained under the existing workflow directories, but this package does not misrepresent them as a one-command raw-data rebuild.

The prepared release is intended for the repository with MIT licensing for original author-owned code and preserved third-party licenses `https://github.com/daming7977-byte/TMD-RiboSeq-Inference` and reserved DOI `10.5281/zenodo.22658221`. At package generation these records were not publicly accessible; activate and verify them before submission.
