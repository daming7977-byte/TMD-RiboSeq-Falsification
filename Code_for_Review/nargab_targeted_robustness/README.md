> Current release package: **v1.0.1**, preserving scientific version **1.0.1-rc.2**. Read the package-root README first. Historical delivery/version, original-byte, publication-state and path statements below describe the preserved reviewer history. The release identity manifest separates historical candidate hashes from current release hashes. Preparation did not publish files or rerun scientific analyses.

# Supplementary Methods S17: targeted robustness reconstruction

This directory supports Supplementary Methods S17 and Supplementary Table S14. It contains the saved inputs, completed calculation scripts with documented path-only normalization, and outputs for two **post-reveal exploratory sensitivities**: same-strand other-gene exon exclusion and matched-population operational UMI-like collapse. These analyses do not replace the frozen primary +75-nt inference or the original regional estimates, establish independent replication, or validate molecular-UMI deduplication.

## Contents and limits

- `inputs/PAIR_METADATA.tsv` and `PAIR_DENSITY_MATRICES.npz`: the exact saved center-specific regional population and paired real/pseudo density matrices from the completed verification. NPZ arrays `uncollapsed` and `collapsed` share metadata row order; `columns` names the 24 library/real-pseudo columns. Matched pseudo densities are arithmetic means over the original distinct pseudo geometries. The original geometry catalog, masks and uncollapsed geometry-density tables remain in `../frozen_intermediates/region_composite_project/`; they are not duplicated here.
- `inputs/pair_overlap_flags.tsv`: exclusion flags for the exact included primary and regional center-specific pairs. Complete real/pseudo pairs are removed if an unmasked base of either the real window or any matched pseudo window overlaps a same-strand exon of a gene outside the frozen owner-gene set. No partial masking or pseudo reweighting is performed.
- `inputs/PRIMARY_PLUS75_UNIT_DENSITIES_v1.tsv` and `PRIMARY_PLUS75_LOCATION_SUMMARY_v1.tsv`: original frozen primary paired densities and summary. `PRIMARY_PLUS75_UNIT_GEOMETRY_MAP_v1.tsv` supports optional annotation regeneration.
- `expected/`: byte-preserved completed primary baseline/exclusion, regional exclusion and matched-population collapse results. These contain the full-precision estimates, replicate-label effects, population sizes and conditional intervals supporting S14.
- `original_scripts/`: completed inventory, bridge and exclusion scripts with documented local-path-prefix normalization only, plus historical molecular-family collapse and geometry-counting source code. Original-source hashes remain for provenance; private absolute path prefixes are normalized in this public candidate. **Do not execute these originals directly on a different machine.** Use the portable adapters below. Historical BAM-processing source is included to document the operational definition, not as part of the density-level replay.
- `evidence/` and `historical_sources/`: original input hashes, density/denominator reconciliation, overlap summaries, coordinate evidence, analysis specification, seeds, environment and the exact historical collapse-family contract. Historical filenames or source text containing “UMI”, “duplicate” or “EXACT5UMI” do not constitute molecular-UMI validation; the current interpretation is **operational UMI-like collapse**.
- `PACKAGED_SOURCE_IDENTITY.tsv`: packaged paths, original paths, byte counts and verified source SHA256 values. `CHECKSUMS.tsv` covers this directory except itself.

The full per-window annotation inventory is provided once at `Source_Data/Targeted_Robustness_S17/P3_OVERLAP_INVENTORY.tsv` after reconstructing the reviewer bundle. The table is at `Supplementary_Tables/Supplementary_Table_14/SUPPLEMENTARY_TABLE_S14.tsv`. Each `source_output` in S14 resolves under this directory's `expected/`; the original A/B source location is retained in `PACKAGED_SOURCE_IDENTITY.tsv`.

## Environment and no-analysis verification

The completed calculations used Python 3.12.14, NumPy 2.3.5 and pandas 2.2.3. These are the recommended versions for exact conditional-interval replay. The scripts use only the Python standard library, NumPy and pandas. No BAM/FASTQ, reference download, alignment, counting or new molecular-family construction occurs in the density replay.

From this directory, check packaged hashes, matrix indexing, schemas and required flags without calculating effects:

```bash
python3 -B replay_density_sensitivities.py --validate-only
```

## Density-level replay

Write to a new output directory outside the package. The default replays the completed point estimates only. Adding `--bootstrap` also replays the completed conditional percentile intervals. Both modes compare their calculated results with the saved full-precision outputs, with absolute tolerance `1e-10`; any mismatch stops the run. No scientific input is overwritten.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python3 -B replay_density_sensitivities.py --output /path/to/new/S17_points

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python3 -B replay_density_sensitivities.py --bootstrap --output /path/to/new/S17_intervals
```

Regional replay uses the exact completed bridge functions. It preserves the frozen available-case, unique-gene MULTI_PASS FIRST/INTERNAL populations, the centerwise density-mean ratio estimands, replicate-label averaging and center averaging. The original baseline and collapsed analyses have identical populations and share sampled gene multiplicities. Overlap exclusion retains the original available-case convention on its explicitly restricted population. There are 10,000 multinomial gene-bootstrap draws in batches of 200, using sorted union-gene order and NumPy `default_rng` seeds 30097890 (plateau) and 30097860 (broader). The primary reconstruction resamples paired units, not genes, with 10,000 draws, batch 200 and seed 30097775. Percentile limits are 2.5 and 97.5. These conditional intervals do not represent biological-replicate uncertainty or establish cross-compartment culture/batch pairing.

Outputs are `REGIONAL_REPLAY.tsv`, `PRIMARY_REPLAY.tsv` and `REPLAY_RECORD.json`. A point-only run leaves interval fields blank; it does not report saved intervals as newly calculated. Primary restricted-population findings remain exploratory and do not replace the original frozen primary results.

## Optional annotation regeneration

The saved inventory and flags permit density-level replay without the GTF. Independent annotation regeneration additionally requires the exact external, uncompressed `Homo_sapiens.GRCh38.104.gtf`:

`SHA256 f3de2bccc02877709272a0ce66af59a2922d407a4076541bc5a65c9d3518581e`

The annotation adapter’s historical docstring uses “byte-preserved” to describe its source implementation. The packaged historical source has documented path-prefix normalization; the adapter also replaces filesystem bindings. That wording is not a claim of whole-file identity to the original script. Scientific coordinate and overlap logic are unchanged.

The large GTF is not duplicated in this package. The adapter neither downloads it nor accepts a different reference. It changes only filesystem bindings in the preserved inventory implementation; coordinate and overlap logic remain unchanged. It uses the existing reviewer frozen geometries, primary geometry map, annotations and masks, converts GTF 1-based inclusive coordinates to 0-based half-open coordinates, and intersects exact spliced bases. It does not use genomic bounding intervals or process reads.

```bash
python3 -B replay_annotation_inventory.py \
  --gtf /path/to/Homo_sapiens.GRCh38.104.gtf \
  --reviewer-root /path/to/Code_for_Review \
  --output /path/to/new/S17_annotation
```

Annotation overlap establishes potential transcript-origin ambiguity, not proven contamination. The examined frozen windows remain valid under their annotated transcript models, including shared readthrough-like structures. The historical five-base field is UMI-like; family collapse is not validated PCR-duplicate removal or reconstruction of the original data authors' procedure. Comparable footprint-length sensitivity remains unavailable from pooled density inputs. This package supports the stated density-level calculations and annotation review; it does not claim complete raw-data reproducibility or biological validation.
