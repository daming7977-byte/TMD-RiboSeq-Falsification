#!/usr/bin/env python3
import csv
import hashlib
from pathlib import Path

FREEZE = 'stage0c_deeptmhmm_v1_2026-08-24'
BATCH_SIZE = 100
REFERENCE = Path('REPRO_WORKSPACE/GSE133111_stage0C/reference/frozen/GSE133111_human_APPRIS19192_author_compatible_v1.fa')
REFERENCE_SHA256 = 'ecda9babdf126cffa6b91e5731c355db9e40db96e2c5c5f737125342961dbf56'
FEATURES = Path('REPRO_WORKSPACE/GSE133111_stage0C/reference/author_reference/human_features.tsv')
FEATURES_SHA256 = 'c54965e9bca3cca852f422402180b8f0811a9ed58e8c9fd7a18b92b0e3099f36'
ELIGIBILITY = Path('REPRO_WORKSPACE/GSE133111_stage0C/qc/tmd_coordinate_audit/deeptmhmm_pilot/deeptmhmm_input_eligibility.tsv')
TRANSLATION = Path('REPRO_WORKSPACE/GSE133111_stage0C/qc/tmd_coordinate_audit/transcript_translation_audit.tsv')
PROTEINS = Path('REPRO_WORKSPACE/GSE133111_stage0C/qc/tmd_coordinate_audit/transcript_derived_proteins_audit_only.fa')
OUT = Path('REPRO_WORKSPACE/GSE133111_stage0C/reference/deeptmhmm_raw')
QC = Path('REPRO_WORKSPACE/GSE133111_stage0C/qc/deeptmhmm_full')

def sha256_text(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()

def read_fasta(path):
    result = []
    rid = None
    seq = []
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip('\n\r')
            if line.startswith('>'):
                if rid is not None:
                    result.append((rid, ''.join(seq)))
                rid = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if rid is not None:
        result.append((rid, ''.join(seq)))
    return result

def read_tsv(path, key):
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 19192
    assert len({row[key] for row in rows}) == 19192
    return {row[key]: row for row in rows}

def write_tsv(path, rows, fields):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

OUT.mkdir(parents=True, exist_ok=True)
(OUT / 'inputs').mkdir(parents=True, exist_ok=True)
(OUT / 'jobs').mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

eligibility = read_tsv(ELIGIBILITY, 'transcript')
translation = read_tsv(TRANSLATION, 'transcript_id')
proteins = read_fasta(PROTEINS)
protein_map = dict(proteins)
assert len(proteins) == len(protein_map) == 18845
assert [rid for rid, _ in proteins] == sorted(protein_map)

eligible_ids = sorted(rid for rid, row in eligibility.items() if row['deeptmhmm_input_eligible'] == 'YES')
ineligible_ids = sorted(rid for rid, row in eligibility.items() if row['deeptmhmm_input_eligible'] == 'NO')
assert len(eligible_ids) == 18845
assert len(ineligible_ids) == 347
assert set(eligible_ids) == set(protein_map)
assert len(eligible_ids) + len(ineligible_ids) == 19192

reference_provenance = f'{REFERENCE};sha256={REFERENCE_SHA256};{FEATURES};sha256={FEATURES_SHA256}'
manifest_rows = []
for rid in sorted(eligibility):
    erow = eligibility[rid]
    trow = translation[rid]
    is_eligible = erow['deeptmhmm_input_eligible'] == 'YES'
    seq = protein_map.get(rid)
    if is_eligible:
        assert seq is not None
        assert int(erow['protein_length']) == len(seq)
        assert erow['protein_sha256'] == sha256_text(seq)
        assert trow['predictor_ready'] == 'TRUE'
        conversion_status = 'PREDICTOR_PROTEIN_CONSTRUCTED'
        exclusion_reason = ''
    else:
        assert seq is None
        assert trow['predictor_ready'] == 'FALSE'
        conversion_status = 'NO_UNAMBIGUOUS_PREDICTOR_PROTEIN'
        exclusion_reason = erow['technical_reason']
        assert exclusion_reason
    manifest_rows.append({
        'transcript': rid,
        'cds_protein_conversion_status': conversion_status,
        'deeptmhmm_input_eligible': erow['deeptmhmm_input_eligible'],
        'protein_length': erow['protein_length'] if is_eligible else '',
        'protein_sha256': erow['protein_sha256'] if is_eligible else '',
        'technical_exclusion_reason': exclusion_reason,
        'translation_terminal_stop_handling': trow['terminal_stop_symbol_handling'],
        'frozen_human_reference_provenance': reference_provenance,
        'deeptmhmm_freeze_version': FREEZE,
    })

manifest_fields = list(manifest_rows[0])
write_tsv(OUT / 'deeptmhmm_full_input_manifest.tsv', manifest_rows, manifest_fields)

batch_rows = []
for start in range(0, len(eligible_ids), BATCH_SIZE):
    ids = eligible_ids[start:start+BATCH_SIZE]
    batch_number = start // BATCH_SIZE + 1
    batch_id = f'batch_{batch_number:04d}'
    fasta_path = OUT / 'inputs' / f'{batch_id}.fasta'
    with fasta_path.open('w') as handle:
        for rid in ids:
            handle.write(f'>{rid}\n{protein_map[rid]}\n')
    fasta_hash = hashlib.sha256(fasta_path.read_bytes()).hexdigest()
    for offset, rid in enumerate(ids, start=1):
        batch_rows.append({
            'transcript': rid,
            'global_eligible_rank': start + offset,
            'batch_id': batch_id,
            'batch_number': batch_number,
            'batch_rank': offset,
            'batch_size_contract': BATCH_SIZE,
            'batch_input_records': len(ids),
            'batch_input_file': str(fasta_path),
            'batch_input_sha256': fasta_hash,
            'protein_length': len(protein_map[rid]),
            'protein_sha256': sha256_text(protein_map[rid]),
            'selection_rule': 'eligible transcripts sorted lexicographically; consecutive fixed blocks of 100; final remainder retained',
            'deeptmhmm_freeze_version': FREEZE,
        })

write_tsv(OUT / 'deeptmhmm_full_batch_manifest.tsv', batch_rows, list(batch_rows[0]))

batch_files = sorted((OUT / 'inputs').glob('batch_*.fasta'))
assert len(batch_files) == 189
assert sum(1 for row in batch_rows) == 18845

summary = [
    {'metric': 'frozen_transcripts', 'value': 19192},
    {'metric': 'deeptmhmm_input_eligible', 'value': 18845},
    {'metric': 'technically_ineligible', 'value': 347},
    {'metric': 'batch_size_contract', 'value': BATCH_SIZE},
    {'metric': 'deterministic_batches', 'value': 189},
    {'metric': 'full_batches_of_100', 'value': 188},
    {'metric': 'final_batch_records', 'value': 45},
]
write_tsv(QC / 'deeptmhmm_full_input_summary.tsv', summary, ['metric', 'value'])
print('frozen_transcripts=19192')
print('eligible=18845')
print('ineligible=347')
print('batches=189')
print('batch_size=100')
print('final_batch_records=45')
