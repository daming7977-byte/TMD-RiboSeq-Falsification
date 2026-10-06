#!/usr/bin/env python3
"""Audit-only validation/normalization for the frozen full DeepTMHMM raw run."""
import csv
import hashlib
import json
import os
import re
from collections import Counter
from pathlib import Path

APP = 'DTU/DeepTMHMM:1.0.24'
FREEZE = 'stage0c_deeptmhmm_v1_2026-08-24'
RAW = Path('REPRO_WORKSPACE/GSE133111_stage0C/reference/deeptmhmm_raw')
QC = Path('REPRO_WORKSPACE/GSE133111_stage0C/qc/deeptmhmm_full')
PROTEINS = Path('REPRO_WORKSPACE/GSE133111_stage0C/qc/tmd_coordinate_audit/transcript_derived_proteins_audit_only.fa')
TRANSCRIPTS = Path('REPRO_WORKSPACE/GSE133111_stage0C/reference/frozen/GSE133111_human_APPRIS19192_author_compatible_v1.fa')
FEATURES = Path('REPRO_WORKSPACE/GSE133111_stage0C/reference/author_reference/human_features.tsv')
ALLOWED_CLASSES = {'TM', 'SP+TM', 'SP', 'GLOB', 'BETA'}
ALLOWED_TOPOLOGY_CHARS = set('IOPSMB')
ALLOWED_GFF_LABELS = {'inside', 'outside', 'periplasm', 'signal', 'TMhelix', 'Beta sheet'}
LABEL_TO_CHAR = {'inside':'I','outside':'O','periplasm':'P','signal':'S','TMhelix':'M','Beta sheet':'B'}
GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def sha256_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

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

def read_tsv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write_tsv_atomic(path, rows, fields):
    tmp = path.with_suffix(path.suffix + '.tmp')
    with tmp.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)

def read_three_line(path):
    lines = path.read_text().splitlines()
    if len(lines) % 3:
        raise ValueError(f'malformed 3line output: {path}')
    records = []
    for i in range(0, len(lines), 3):
        match = re.fullmatch(r'>(\S+) \| (\S+)', lines[i])
        if not match:
            raise ValueError(f'malformed 3line header: {lines[i]}')
        records.append((match.group(1), match.group(2), lines[i+1], lines[i+2]))
    return records

def read_gff(path):
    records = {}
    lengths = {}
    declared_tmrs = {}
    current = None
    for line in path.read_text().splitlines():
        if line.startswith('##') or line == '//':
            continue
        if line.startswith('# '):
            match = re.fullmatch(r'# (\S+) Length: (\d+)', line)
            if match:
                current = match.group(1)
                lengths[current] = int(match.group(2))
                records.setdefault(current, [])
                continue
            match = re.fullmatch(r'# (\S+) Number of predicted TMRs: (\d+)', line)
            if match:
                declared_tmrs[match.group(1)] = int(match.group(2))
                continue
        fields = line.split('\t')
        if len(fields) < 4:
            raise ValueError(f'malformed GFF line: {line}')
        rid, label, start1, end1 = fields[:4]
        if rid != current:
            raise ValueError(f'GFF record/header mismatch: {rid} != {current}')
        records[rid].append((label, int(start1), int(end1)))
    return records, lengths, declared_tmrs

def translate(nt):
    if len(nt) % 3:
        raise ValueError('round-trip nucleotide interval is not triplet length')
    return ''.join(GENETIC_CODE[nt[i:i+3]] for i in range(0, len(nt), 3))

QC.mkdir(parents=True, exist_ok=True)
input_manifest = read_tsv(RAW / 'deeptmhmm_full_input_manifest.tsv')
batch_manifest = read_tsv(RAW / 'deeptmhmm_full_batch_manifest.tsv')
protein_records = read_fasta(PROTEINS)
protein_map = dict(protein_records)
tx_map = dict(read_fasta(TRANSCRIPTS))
with FEATURES.open(newline='') as handle:
    features = {row['Transcript']: int(row['UTR5_Length']) for row in csv.DictReader(handle, delimiter='\t')}

eligible_ids = [row['transcript'] for row in batch_manifest]
assert eligible_ids == sorted(eligible_ids)
assert len(eligible_ids) == len(set(eligible_ids)) == 18845
assert len(input_manifest) == 19192

batch_to_ids = {}
for row in batch_manifest:
    batch_to_ids.setdefault(row['batch_id'], []).append(row['transcript'])

normalized_rows = []
prediction_rows = []
accounting_rows = []
raw_inventory_rows = []
classes = Counter()
topology_chars = Counter()
gff_labels = Counter()
roundtrip_mismatches = 0
bounds_mismatches = 0
malformed_outputs = 0
successful_ids = set()
explicit_failure_ids = set()

for batch_number in range(1, 190):
    batch_id = f'batch_{batch_number:04d}'
    result_path = RAW / 'jobs' / batch_id / 'batch_result.json'
    expected_ids = batch_to_ids[batch_id]
    if not result_path.exists():
        explicit_failure_ids.update(expected_ids)
        for rid in expected_ids:
            accounting_rows.append({'transcript':rid,'batch_id':batch_id,'outcome':'EXPLICIT_TECHNICAL_FAILURE','source_job_id':'','detail':'missing batch_result.json'})
        continue
    result = json.loads(result_path.read_text())
    if result['final_status'] != 'SUCCESS':
        explicit_failure_ids.update(expected_ids)
        for rid in expected_ids:
            accounting_rows.append({'transcript':rid,'batch_id':batch_id,'outcome':'EXPLICIT_TECHNICAL_FAILURE','source_job_id':'','detail':result['final_status']})
        continue
    raw_dir = Path(result['raw_output_directory'])
    source_job_id = result['source_job_id']
    for path in sorted(p for p in raw_dir.rglob('*') if p.is_file()):
        raw_inventory_rows.append({
            'batch_id': batch_id,
            'source_job_id': source_job_id,
            'relative_path': str(path.relative_to(raw_dir)),
            'bytes': path.stat().st_size,
            'sha256': sha256_file(path),
            'official_app': APP,
            'freeze_version': FREEZE,
        })
    try:
        top_records = read_three_line(raw_dir / 'predicted_topologies.3line')
        gff_records, gff_lengths, declared_tmrs = read_gff(raw_dir / 'TMRs.gff3')
    except Exception:
        malformed_outputs += 1
        explicit_failure_ids.update(expected_ids)
        for rid in expected_ids:
            accounting_rows.append({'transcript':rid,'batch_id':batch_id,'outcome':'EXPLICIT_TECHNICAL_FAILURE','source_job_id':source_job_id,'detail':'malformed raw output'})
        continue
    output_ids = [record[0] for record in top_records]
    if output_ids != expected_ids or len(output_ids) != len(set(output_ids)) or set(gff_records) != set(expected_ids):
        malformed_outputs += 1
        explicit_failure_ids.update(expected_ids)
        for rid in expected_ids:
            accounting_rows.append({'transcript':rid,'batch_id':batch_id,'outcome':'EXPLICIT_TECHNICAL_FAILURE','source_job_id':source_job_id,'detail':'input/output ID accounting mismatch'})
        continue
    batch_valid = True
    batch_prediction_rows = []
    batch_normalized_rows = []
    for rid, protein_class, output_seq, topology in top_records:
        expected_seq = protein_map[rid]
        if protein_class not in ALLOWED_CLASSES or set(topology) - ALLOWED_TOPOLOGY_CHARS:
            batch_valid = False
        if output_seq != expected_seq or len(topology) != len(expected_seq):
            batch_valid = False
        if gff_lengths.get(rid) != len(expected_seq):
            batch_valid = False
        classes[protein_class] += 1
        topology_chars.update(topology)
        covered = [None] * len(expected_seq)
        segment_rows_for_rid = []
        for segment_index, (label, start1, end1) in enumerate(gff_records[rid], start=1):
            gff_labels[label] += 1
            if label not in ALLOWED_GFF_LABELS or not (1 <= start1 <= end1 <= len(expected_seq)):
                batch_valid = False
                bounds_mismatches += 1
                continue
            aa_start0 = start1 - 1
            aa_end0 = end1
            expected_topology_char = LABEL_TO_CHAR[label]
            if topology[aa_start0:aa_end0] != expected_topology_char * (aa_end0 - aa_start0):
                batch_valid = False
            for pos in range(aa_start0, aa_end0):
                if covered[pos] is not None:
                    batch_valid = False
                    bounds_mismatches += 1
                covered[pos] = label
            cds_nt_start0 = 3 * aa_start0
            cds_nt_end0 = 3 * aa_end0
            tx_nt_start0 = features[rid] + cds_nt_start0
            tx_nt_end0 = features[rid] + cds_nt_end0
            nt = tx_map[rid][tx_nt_start0:tx_nt_end0].upper().replace('U','T')
            translated = translate(nt)
            roundtrip = translated == expected_seq[aa_start0:aa_end0]
            if not roundtrip:
                batch_valid = False
                roundtrip_mismatches += 1
            segment_rows_for_rid.append({
                'transcript': rid,
                'segment_index_raw_order': segment_index,
                'raw_deeptmhmm_label': label,
                'source_residue_start': start1,
                'source_residue_end': end1,
                'normalized_aa_start0': aa_start0,
                'normalized_aa_end0_exclusive': aa_end0,
                'cds_nt_start0': cds_nt_start0,
                'cds_nt_end0_exclusive': cds_nt_end0,
                'transcript_nt_start0': tx_nt_start0,
                'transcript_nt_end0_exclusive': tx_nt_end0,
                'protein_length': len(expected_seq),
                'protein_sha256': hashlib.sha256(expected_seq.encode()).hexdigest(),
                'deeptmhmm_app_version': APP,
                'freeze_version': FREEZE,
                'source_job_id': source_job_id,
                'source_batch_id': batch_id,
                'roundtrip_match': 'YES' if roundtrip else 'NO',
            })
        if any(item is None for item in covered):
            batch_valid = False
            bounds_mismatches += 1
        parsed_tmrs = sum(1 for label, _, _ in gff_records[rid] if label in ('TMhelix','Beta sheet'))
        if declared_tmrs.get(rid) != parsed_tmrs:
            batch_valid = False
        batch_normalized_rows.extend(segment_rows_for_rid)
        batch_prediction_rows.append({
            'transcript': rid,
            'raw_predicted_protein_class': protein_class,
            'protein_length': len(expected_seq),
            'protein_sha256': hashlib.sha256(expected_seq.encode()).hexdigest(),
            'raw_topology_string_sha256': hashlib.sha256(topology.encode()).hexdigest(),
            'raw_topology_segment_count': len(gff_records[rid]),
            'raw_declared_tmr_count': declared_tmrs[rid],
            'source_job_id': source_job_id,
            'source_batch_id': batch_id,
            'deeptmhmm_app_version': APP,
            'freeze_version': FREEZE,
        })
    if not batch_valid:
        explicit_failure_ids.update(expected_ids)
        for rid in expected_ids:
            accounting_rows.append({'transcript':rid,'batch_id':batch_id,'outcome':'EXPLICIT_TECHNICAL_FAILURE','source_job_id':source_job_id,'detail':'post-output validation failure'})
        continue
    successful_ids.update(expected_ids)
    normalized_rows.extend(batch_normalized_rows)
    prediction_rows.extend(batch_prediction_rows)
    for rid in expected_ids:
        accounting_rows.append({'transcript':rid,'batch_id':batch_id,'outcome':'SUCCESS','source_job_id':source_job_id,'detail':''})

assert successful_ids.isdisjoint(explicit_failure_ids)
assert successful_ids | explicit_failure_ids == set(eligible_ids)

normalized_fields = ['transcript','segment_index_raw_order','raw_deeptmhmm_label','source_residue_start','source_residue_end','normalized_aa_start0','normalized_aa_end0_exclusive','cds_nt_start0','cds_nt_end0_exclusive','transcript_nt_start0','transcript_nt_end0_exclusive','protein_length','protein_sha256','deeptmhmm_app_version','freeze_version','source_job_id','source_batch_id','roundtrip_match']
prediction_fields = ['transcript','raw_predicted_protein_class','protein_length','protein_sha256','raw_topology_string_sha256','raw_topology_segment_count','raw_declared_tmr_count','source_job_id','source_batch_id','deeptmhmm_app_version','freeze_version']
accounting_fields = ['transcript','batch_id','outcome','source_job_id','detail']
inventory_fields = ['batch_id','source_job_id','relative_path','bytes','sha256','official_app','freeze_version']
write_tsv_atomic(RAW / 'deeptmhmm_full_normalized_topology_segments.tsv', normalized_rows, normalized_fields)
write_tsv_atomic(RAW / 'deeptmhmm_full_raw_prediction_classes.tsv', prediction_rows, prediction_fields)
write_tsv_atomic(QC / 'deeptmhmm_full_input_output_accounting.tsv', sorted(accounting_rows, key=lambda row: row['transcript']), accounting_fields)
write_tsv_atomic(QC / 'deeptmhmm_full_raw_output_inventory.tsv', raw_inventory_rows, inventory_fields)

summary = {
    'freeze_source_version': FREEZE,
    'official_app': APP,
    'frozen_transcripts': 19192,
    'deeptmhmm_input_eligible': 18845,
    'technically_ineligible': 347,
    'successful_outputs': len(successful_ids),
    'explicit_failures': len(explicit_failure_ids),
    'input_output_accounting_pass': len(successful_ids) + len(explicit_failure_ids) == 18845,
    'raw_topology_segments': len(normalized_rows),
    'coordinate_roundtrip_segments_tested': len(normalized_rows),
    'roundtrip_mismatches': roundtrip_mismatches,
    'bounds_mismatches': bounds_mismatches,
    'malformed_batch_outputs': malformed_outputs,
    'predicted_class_counts': dict(sorted(classes.items())),
    'topology_character_counts': dict(sorted(topology_chars.items())),
    'raw_gff_label_counts': dict(sorted(gff_labels.items())),
}
(QC / 'deeptmhmm_full_validation_summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
print(json.dumps(summary, indent=2, sort_keys=True))
