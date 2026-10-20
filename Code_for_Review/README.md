# Code supporting scientific version 1.0.1-rc.2

Release **v1.0.1**, preserving scientific version **1.0.1-rc.2**, prepared from verified local candidate **1.0.1-rc.2+public-candidate.20261019.2**. Read the package-root README, component-rights statement, release metadata and identity manifest. No publication identifier is asserted by the preparation record.

- Frozen region-composite estimands, saved primary checks, and Figures 1, 2 and 4: `REPRODUCE_NARGAB.md` and `reproduce_nargab.py`.
- Post-reveal exploratory within-gene/composition sensitivities, Supplementary Methods S16 and Table S13: `nargab_targeted_sensitivity/README.md`.
- Post-reveal annotation-overlap exclusion and matched-population operational UMI-like collapse, S17 and Table S14: `nargab_targeted_robustness/README.md`.

- Figure 3 historical SVG and terminology-only presentation command: `presentation_v5/README.md`. The main reproduction command does not generate or copy Figure 3; this release package does not include a standalone final Figure 3 export.
- Verified Route A/B2 historical provenance and three-way digests: `provenance/VERIFIED_HISTORICAL_PROVENANCE.md`. The five Route A/B2 files and the five S17 privacy-normalized files are separate groups.

Keep sibling `Source_Data/` and `Supplementary_Tables/` directories. Original processed input matrices and scientific outputs are preserved. Five historical S17 sources contain only documented local-path-prefix normalization; their supported adapters use unchanged scientific functions or replace filesystem bindings before execution. Use the documented adapters, not historical scripts directly. Historical protocol/source hashes describe originals and are not current-copy checksum claims. `FILE_MANIFEST_RC2.tsv` and the package-root release manifest cover current bytes; older manifests retain historical scope.

The inherited coefficient-unit warning in `REPRODUCE_NARGAB.md` remains applicable and unchanged: a callable historical routine labels natural-log effects as log2, but did not produce reported P3 results. Do not use it to interpret reported log2 effects.

The public v1.0.0 baseline is at https://github.com/daming7977-byte/TMD-RiboSeq-Falsification/releases/tag/v1.0.0 and https://doi.org/10.5281/zenodo.22658221. Those identifiers do not archive this controlled revision. Original author-owned code is MIT; third-party licenses and source-data terms remain distinct. Raw-FASTQ-to-alignment historical workflows are retained for provenance and require external assets; complete raw-data reproduction is not claimed.
