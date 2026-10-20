#!/usr/bin/env python3
"""Optional annotation replay: supply the exact external GTF; never downloads it.

Runs the byte-preserved completed inventory implementation after replacing only
its five filesystem bindings. No BAM/FASTQ, statistical resampling or read counts.
"""
from pathlib import Path
import argparse
import ast
import csv
import hashlib

HERE = Path(__file__).resolve().parent
GTF_SHA256 = 'f3de2bccc02877709272a0ce66af59a2922d407a4076541bc5a65c9d3518581e'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gtf', type=Path, required=True, help='Uncompressed Homo_sapiens.GRCh38.104.gtf')
    parser.add_argument('--reviewer-root', type=Path, default=HERE.parent, help='Extracted Code_for_Review directory')
    parser.add_argument('--output', type=Path, required=True, help='New directory outside this package')
    args = parser.parse_args()
    gtf = args.gtf.resolve()
    output = args.output.resolve()
    if output.exists() or output.is_relative_to(HERE):
        raise RuntimeError('Use a new output directory outside packaged inputs')
    with gtf.open('rb') as stream:
        if hashlib.file_digest(stream, 'sha256').hexdigest() != GTF_SHA256:
            raise RuntimeError('Exact Ensembl release 104 GTF hash mismatch')
    frozen = args.reviewer_root.resolve() / 'frozen_intermediates/region_composite_project'
    mapping = HERE / 'inputs/PRIMARY_PLUS75_UNIT_GEOMETRY_MAP_v1.tsv'
    bindings = {'OUT': output, 'F': frozen, 'H': HERE / 'inputs', 'GTF': gtf, 'PRIM': mapping}
    script = HERE / 'original_scripts/build_overlap_inventory.py'
    # Verify the preserved script and primary map before evaluating its logic.
    source_rows = list(csv.DictReader((HERE / 'PACKAGED_SOURCE_IDENTITY.tsv').open(), delimiter='\t'))
    for target in [script, mapping]:
        record = next(row for row in source_rows if HERE / row['packaged_path'] == target)
        with target.open('rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != record['sha256']:
                raise RuntimeError(f'Packaged source hash mismatch: {target}')
    # The corresponding frozen files are already present in the full reviewer
    # package, but must match the completed inventory's exact reference inputs.
    records = list(csv.DictReader((HERE / 'evidence/INVENTORY_ORIGINAL_INPUT_SHA256.tsv').open(), delimiter='\t'))
    for suffix in [
        'tmd_anchor_preregister/bridge/FORMAL_TMD_ANCHOR_CANDIDATES_v1p1.tsv',
        'crossmap_window_preregister/output/GSE300977_PAIRED_REAL_PSEUDO_CROSSMAP_MASK_AUDIT_v1.tsv',
        'tmd_spatial_secondary_scan_v1/output/SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv',
    ]:
        target = frozen / suffix
        record = next(row for row in records if row['path'].endswith('/' + suffix))
        with target.open('rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != record['sha256']:
                raise RuntimeError(f'Frozen source hash mismatch: {target}')
    tree = ast.parse(script.read_text(), filename=str(script))
    replaced = []
    for statement in tree.body:
        if isinstance(statement, ast.Assign) and len(statement.targets) == 1 and isinstance(statement.targets[0], ast.Name):
            name = statement.targets[0].id
            if name in bindings:
                statement.value = ast.parse(f'Path({str(bindings[name])!r})', mode='eval').body
                replaced.append(name)
    if set(replaced) != set(bindings) or len(replaced) != len(bindings):
        raise RuntimeError('Unexpected original script path bindings; stop rather than alter logic')
    ast.fix_missing_locations(tree)
    output.mkdir(parents=True)
    namespace = {'__name__': '__main__', '__file__': str(script)}
    exec(compile(tree, str(script), 'exec'), namespace)


if __name__ == '__main__':
    main()
