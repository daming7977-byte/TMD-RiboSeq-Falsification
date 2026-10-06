# Reproduce the NARGAB revision

## Scope

This package reproduces the frozen post-reveal region-composite statistics, verifies the primary +75-nt statistics, and rebuilds main Figures 1, 2 and 4. Figure 3 is distributed as unchanged frozen artwork; this command neither regenerates nor copies it. The primary and benchmark statistics are read/checked as saved outputs, not recomputed by this command. The package does not perform discovery, choose windows or reinterpret P values.

The included 90-MB `frozen_intermediates/region_composite_project/` directory contains the exact processed inputs required by the region-composite script. Raw FASTQ-to-alignment workflows remain documented historical workflows and are not silently represented as a one-command raw-data reproduction.

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

For the journal review files, extract `Code_for_Review.zip` and `Source_Data.zip` into the same parent directory so that `Code_for_Review/` and `Source_Data/` are siblings. The complete GitHub/Zenodo release archive already contains the equivalent `code/` and `source_data/` layout. From `Code_for_Review/` in the journal layout, or `code/` in the release layout:

```bash
python reproduce_nargab.py
```

The command:

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
- Available-case +90/+105/+120 interaction: 0.9325468 [0.7122142, 1.1328309]; biological replicates 1.0736264 and 0.7914672.
- Available-case +60/+75/+90/+105/+120 interaction: 0.6518202 [0.3573579, 1.0554072]; biological replicates 0.7840298 and 0.5196107.

These region summaries are outcome-informed post-reveal sensitivity analyses. The gene bootstrap addresses within-gene dependence but does not create biological replication.

## Output locations and limits

Run in a disposable extracted copy. Region TSVs are written to code/reproduction_output (Code_for_Review/reproduction_output in the journal layout). Figures 1, 2 and 4 overwrite the parent Figures directory; they are not written inside reproduction_output. Figure 3 must be retained from the distribution. Historical extension scripts have separate, incomplete upstream environment bindings, documented with Supplementary Table 11; they are not invoked by this command. No scientific rerun was performed during forensic remediation. The prior recorded byte-identity test is historical evidence, not a fresh environment test.
