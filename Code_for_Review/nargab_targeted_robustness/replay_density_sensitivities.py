#!/usr/bin/env python3
"""Portable replay of the completed S17 density-level calculations; no read processing.

The imported bridge functions are the byte-preserved completed-verification code.
--validate-only checks files, schemas and indexing without calculating effects.
Replay writes only to a new output directory explicitly supplied by the reviewer.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import platform
import sys
import time
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TOL = 1e-10


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def verify_files():
    rows = list(csv.DictReader((HERE / 'CHECKSUMS.tsv').open(), delimiter='\t'))
    for row in rows:
        path = HERE / row['relative_path']
        require(path.is_file(), f'Missing packaged file: {path}')
        require(path.stat().st_size == int(row['bytes']), f'Size mismatch: {path}')
        require(digest(path) == row['sha256'], f'SHA256 mismatch: {path}')
    return len(rows)


def load_inputs():
    # Default pandas conversion intentionally matches the completed regional bridge.
    pairs = pd.read_csv(HERE / 'inputs/PAIR_METADATA.tsv', sep='\t', keep_default_na=False)
    with np.load(HERE / 'inputs/PAIR_DENSITY_MATRICES.npz', allow_pickle=False) as archive:
        matrices = {key: archive[key] for key in ['uncollapsed', 'collapsed']}
        columns = archive['columns']
    expected_columns = [f'B{i:02d}_{side}' for i in range(1, 13) for side in ['real', 'pseudo']]
    require(columns.tolist() == expected_columns, 'Density column order differs')
    require(np.array_equal(pairs.matrix_row, np.arange(len(pairs))), 'Matrix row index differs')
    require(not pairs.duplicated(['center', 'unit_key']).any(), 'Duplicate center/unit')
    require(not pairs.paired_real_unit_id.duplicated().any(), 'Duplicate pair ID')
    for key, matrix in matrices.items():
        require(matrix.shape == (len(pairs), 24), f'Matrix shape: {key}')
        require(np.isfinite(matrix).all() and (matrix >= 0).all(), f'Invalid densities: {key}')
    flags = pd.read_csv(HERE / 'inputs/pair_overlap_flags.tsv', sep='\t', keep_default_na=False)
    require(not flags.paired_real_unit_id.duplicated().any(), 'Duplicate annotation flag ID')
    flags = flags.set_index('paired_real_unit_id')
    selected = pairs.center.isin([60, 75, 90, 105, 120]).to_numpy()
    regional = pairs.loc[selected]
    require(set(regional.paired_real_unit_id) <= set(flags.index), 'Regional flags missing')
    require(np.array_equal(flags.loc[regional.paired_real_unit_id, 'center'], regional.center),
            'Regional flag center mismatch')
    # Round-trip conversion is the completed primary script's exact convention.
    primary = pd.read_csv(HERE / 'inputs/PRIMARY_PLUS75_UNIT_DENSITIES_v1.tsv',
                          sep='\t', float_precision='round_trip')
    require(len(primary) == 9098 and not primary.paired_real_unit_id.duplicated().any(),
            'Frozen primary population differs')
    require(set(primary.paired_real_unit_id) <= set(flags.index), 'Primary flags missing')
    return pairs, matrices, flags, primary


def load_bridge():
    spec = importlib.util.spec_from_file_location('completed_bridge', HERE / 'original_scripts/umi_bridge.py')
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)  # Functions only: original main() is guarded.
    return bridge


def close(actual, expected, label):
    require(abs(float(actual) - float(expected)) <= TOL,
            f'Completed-result mismatch for {label}: {actual} versus {expected}')


def regions(bridge, pairs, matrices, flags, with_bootstrap):
    expected = pd.read_csv(HERE / 'expected/P3_UMI_MATCHED_POPULATION_RESULTS.tsv', sep='\t', float_precision='round_trip')
    excluded_expected = pd.read_csv(HERE / 'expected/REGION_EXCLUSION_RESULTS.tsv', sep='\t', float_precision='round_trip')
    rows = []
    selected = pairs.center.isin([60, 75, 90, 105, 120]).to_numpy()
    reduced = pairs.loc[selected].reset_index(drop=True)
    removed = flags.loc[reduced.paired_real_unit_id, 'exclude_pair_unmasked_other_gene'].eq('YES').to_numpy()
    # Re-establish both frozen point gates before any collapsed estimate or
    # bootstrap replay, as in the completed verification.
    for region, centers in bridge.REGIONS.items():
        gate = bridge.prepare(pairs, matrices['uncollapsed'], centers)
        value = bridge.evaluate(gate, np.ones((1, len(gate['genes']))))['mean'][0]
        close(value, bridge.EXPECTED[region], f'{region}/frozen baseline gate')
    for region, centers in bridge.REGIONS.items():
        baseline = bridge.prepare(pairs, matrices['uncollapsed'], centers)
        collapsed = bridge.prepare(pairs, matrices['collapsed'], centers)
        keep = bridge.prepare(reduced.loc[~removed], matrices['uncollapsed'][selected][~removed], centers)
        require(baseline['genes'] == collapsed['genes'], 'Matched population changed')
        # Paired baseline/collapsed use the same multinomial draws. The exclusion
        # uses its explicitly reduced gene population and the same fixed seed.
        boot = [None, None, None]
        if with_bootstrap:
            boot[:2] = bridge.bootstrap_many([baseline, collapsed], bridge.SEEDS[region])
            boot[2] = bridge.bootstrap_many([keep], bridge.SEEDS[region])[0]
        for scenario, prepared, distribution in zip(
                ['uncollapsed', 'operationally_collapsed', 'overlap_exclusion'],
                [baseline, collapsed, keep], boot):
            obs = bridge.evaluate(prepared, np.ones((1, len(prepared['genes']))))
            if scenario == 'overlap_exclusion':
                ref = excluded_expected.loc[excluded_expected.region.eq(region)].iloc[0]
                r1, r2, lo, hi = ref.rep1, ref.rep2, ref.ci_low, ref.ci_high
            else:
                ref = expected.loc[expected.region.eq(region) & expected.density_track.eq(scenario)].iloc[0]
                r1, r2 = ref.replicate_1, ref.replicate_2
                lo, hi = ref.conditional_gene_bootstrap_95_lower, ref.conditional_gene_bootstrap_95_upper
            close(obs['mean'][0], ref.interaction, f'{region}/{scenario}/interaction')
            close(obs['rep1'][0], r1, f'{region}/{scenario}/replicate1')
            close(obs['rep2'][0], r2, f'{region}/{scenario}/replicate2')
            ci = [np.nan, np.nan] if distribution is None else np.percentile(distribution, [2.5, 97.5])
            if distribution is not None:
                close(ci[0], lo, f'{region}/{scenario}/lower')
                close(ci[1], hi, f'{region}/{scenario}/upper')
            rows.append({'analysis_status': 'POST-REVEAL EXPLORATORY SENSITIVITY' if scenario != 'uncollapsed' else 'frozen baseline reproduction',
                         'region': region, 'scenario': scenario, 'genes': len(prepared['genes']),
                         'union_units': prepared['meta'].unit_key.nunique(), 'center_unit_rows': len(prepared['meta']),
                         'replicate_1': obs['rep1'][0], 'replicate_2': obs['rep2'][0], 'interaction': obs['mean'][0],
                         'conditional_95_lower': ci[0], 'conditional_95_upper': ci[1],
                         'bootstrap_dependence_unit': 'gene', 'bootstrap_draws': 10000 if with_bootstrap else 0,
                         'bootstrap_seed': bridge.SEEDS[region], 'completed_result_tolerance': TOL})
    return rows


def primary_results(primary, flags, with_bootstrap):
    contrasts = {'ER': [(1, 2), (11, 12)], 'mitochondria': [(9, 10), (7, 8)], 'cytosol': [(5, 6), (3, 4)]}
    cols = [f'B{i:02d}_{side}_density' for i in range(1, 13) for side in ['real', 'pseudo']]
    matrix = primary[cols].to_numpy()
    require(np.isfinite(matrix).all() and (matrix >= 0).all(), 'Invalid primary densities')

    def evaluate(means):
        require((means > 0).all(), 'Nonpositive primary aggregate; no pseudocount rescue')
        return np.stack([np.stack([np.log2(means[:, 2*(a-1)] / means[:, 2*(b-1)]) -
                                   np.log2(means[:, 2*(a-1)+1] / means[:, 2*(b-1)+1])
                                   for a, b in reps], axis=1) for reps in contrasts.values()], axis=1)

    removed = flags.loc[primary.paired_real_unit_id, 'exclude_pair_unmasked_other_gene'].eq('YES').to_numpy()
    require(int(removed.sum()) == 348, 'Primary exclusion population differs')
    rows = []
    for scenario, values in [('baseline', matrix), ('exclusion', matrix[~removed])]:
        ref = pd.read_csv(HERE / f'expected/primary_{scenario}_results.tsv', sep='\t', float_precision='round_trip').set_index('location')
        n = len(values)
        obs = evaluate(values.mean(axis=0)[None, :])
        ci = np.full((2, 3), np.nan)
        if with_bootstrap:
            rng = np.random.default_rng(30097775)
            parts = []
            for start in range(0, 10000, 200):
                weights = rng.multinomial(n, np.full(n, 1/n), size=200)
                parts.append(evaluate((weights @ values) / n).mean(axis=2))
            ci = np.percentile(np.concatenate(parts), [2.5, 97.5], axis=0)
        for j, location in enumerate(contrasts):
            original = ref.loc[location]
            close(obs[0, j].mean(), original.estimate, f'primary/{scenario}/{location}/estimate')
            close(obs[0, j, 0], original.rep1, f'primary/{scenario}/{location}/replicate1')
            close(obs[0, j, 1], original.rep2, f'primary/{scenario}/{location}/replicate2')
            if with_bootstrap:
                close(ci[0, j], original.ci_low, f'primary/{scenario}/{location}/lower')
                close(ci[1, j], original.ci_high, f'primary/{scenario}/{location}/upper')
            rows.append({'analysis_status': 'frozen baseline reproduction' if scenario == 'baseline' else 'POST-REVEAL EXPLORATORY SENSITIVITY',
                         'analysis': 'PRIMARY_PLUS75', 'scenario': scenario, 'location': location,
                         'original_units': 9098, 'excluded_units': 0 if scenario == 'baseline' else 348,
                         'retained_units': n, 'replicate_1': obs[0, j, 0], 'replicate_2': obs[0, j, 1],
                         'estimate': obs[0, j].mean(), 'conditional_95_lower': ci[0, j], 'conditional_95_upper': ci[1, j],
                         'bootstrap_dependence_unit': 'paired unit', 'bootstrap_draws': 10000 if with_bootstrap else 0,
                         'bootstrap_seed': 30097775, 'completed_result_tolerance': TOL})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validate-only', action='store_true', help='Hash/schema/index checks only; no effect calculations')
    parser.add_argument('--output', type=Path, help='New, nonexisting directory outside this package')
    parser.add_argument('--bootstrap', action='store_true', help='Replay the completed 10,000-draw conditional intervals')
    args = parser.parse_args()
    count = verify_files()
    pairs, matrices, flags, primary = load_inputs()
    if args.validate_only:
        require(not args.bootstrap and args.output is None, '--validate-only cannot run or write analyses')
        print(json.dumps({'packaged_hash_checks': count, 'matrix_shape': list(matrices['uncollapsed'].shape),
                          'overlap_flag_rows': len(flags), 'primary_rows': len(primary),
                          'analysis_executed': False, 'validation': 'PASS'}, indent=2))
        return
    require(args.output is not None, '--output is required for replay')
    output = args.output.resolve()
    require(not output.exists(), 'Output must be a new directory')
    require(not output.is_relative_to(HERE), 'Write replay outputs outside the packaged inputs')
    started = time.monotonic()
    bridge = load_bridge()
    regional_rows = regions(bridge, pairs, matrices, flags, args.bootstrap)
    primary_rows = primary_results(primary, flags, args.bootstrap)
    output.mkdir(parents=True)
    pd.DataFrame(regional_rows).to_csv(output / 'REGIONAL_REPLAY.tsv', sep='\t', index=False, float_format='%.17g')
    pd.DataFrame(primary_rows).to_csv(output / 'PRIMARY_REPLAY.tsv', sep='\t', index=False, float_format='%.17g')
    record = {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__, 'platform': platform.platform(),
              'seconds': time.monotonic() - started, 'bootstrap': args.bootstrap, 'batch': 200,
              'packaged_hash_checks': count, 'tolerance': TOL, 'matched_completed_outputs': True,
              'input_boundaries': 'Saved paired densities and annotation flags; no BAM, FASTQ, or annotation reprocessing'}
    (output / 'REPLAY_RECORD.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
