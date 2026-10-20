#!/usr/bin/env python3
"""Audit-only controller for the frozen GSE133111 DeepTMHMM raw run."""
import argparse
import csv
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import biolib
from biolib.biolib_api_client.api_client import BiolibApiClient

APP_URI = 'DTU/DeepTMHMM:1.0.24'
FREEZE_VERSION = 'stage0c_deeptmhmm_v1_2026-08-24'
RAW_ROOT = Path('REPRO_WORKSPACE/GSE133111_stage0C/reference/deeptmhmm_raw')
INPUT_ROOT = RAW_ROOT / 'inputs'
JOB_ROOT = RAW_ROOT / 'jobs'
QC_ROOT = Path('REPRO_WORKSPACE/GSE133111_stage0C/qc/deeptmhmm_full')
MAX_NEW_ATTEMPTS_PER_CONTROLLER_PASS = 2
POLL_SECONDS = 15
APP_OBJECT = None

def now():
    return datetime.now(timezone.utc).isoformat()

def sha256_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def fasta_count(path):
    with path.open() as handle:
        return sum(1 for line in handle if line.startswith('>'))

def write_json_atomic(path, data):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    os.replace(tmp, path)

def write_bytes(path, data):
    path.write_bytes(data if isinstance(data, bytes) else bytes(data))

def output_inventory(raw_dir):
    rows = []
    for path in sorted(p for p in raw_dir.rglob('*') if p.is_file()):
        rows.append({
            'relative_path': str(path.relative_to(raw_dir)),
            'bytes': path.stat().st_size,
            'sha256': sha256_file(path),
        })
    return rows

def run_attempt(batch_id, input_path, attempt_number):
    batch_dir = JOB_ROOT / batch_id
    attempt_dir = batch_dir / f'attempt_{attempt_number}'
    raw_dir = attempt_dir / 'raw'
    attempt_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    input_records = fasta_count(input_path)
    input_sha256 = sha256_file(input_path)
    command = f"biolib run '{APP_URI}' --fasta {input_path}"
    provenance_path = attempt_dir / 'job_provenance.json'
    provenance = {
        'batch_id': batch_id,
        'attempt_number': attempt_number,
        'freeze_version': FREEZE_VERSION,
        'official_app': APP_URI,
        'execution': 'BioLib Cloud',
        'command_equivalent': command,
        'input_file': str(input_path),
        'input_sha256': input_sha256,
        'input_records': input_records,
        'submission_time_utc': now(),
        'job_id': None,
        'status': 'SUBMITTING',
        'authentication': 'authenticated BioLib account; no token or account identifier persisted in project files',
        'retry_policy': 'up to two identical attempts per controller pass; every unresolved technical batch receives the same sequential/concurrency-4 recovery policy; app and input unchanged',
    }
    write_json_atomic(provenance_path, provenance)
    try:
        if APP_OBJECT is None:
            raise RuntimeError('controller app object not initialized')
        job = APP_OBJECT.cli(
            args=['--fasta', str(input_path)],
            blocking=False,
            result_prefix=f'gse133111-{FREEZE_VERSION}-{batch_id}-attempt-{attempt_number}',
        )
        provenance['job_id'] = job.id
        provenance['status'] = 'SUBMITTED'
        provenance['job_id_recorded_time_utc'] = now()
        write_json_atomic(provenance_path, provenance)
        while True:
            status = job.get_status()
            provenance['last_observed_status'] = status
            provenance['last_status_time_utc'] = now()
            write_json_atomic(provenance_path, provenance)
            if status not in ('in_progress', 'pending', 'created', 'queued'):
                break
            time.sleep(POLL_SECONDS)
        provenance['completion_time_utc'] = now()
        provenance['status'] = status
        try:
            provenance['exit_code'] = job.get_exit_code()
        except Exception as exc:
            provenance['exit_code_error'] = repr(exc)
            provenance['exit_code'] = None
        try:
            write_bytes(attempt_dir / 'stdout.log', job.get_stdout())
        except Exception as exc:
            provenance['stdout_capture_error'] = repr(exc)
        try:
            write_bytes(attempt_dir / 'stderr.log', job.get_stderr())
        except Exception as exc:
            provenance['stderr_capture_error'] = repr(exc)
        try:
            job.save_files(output_dir=str(raw_dir), overwrite=True)
            provenance['raw_output_inventory'] = output_inventory(raw_dir)
        except Exception as exc:
            provenance['output_save_error'] = repr(exc)
        required = {'TMRs.gff3', 'predicted_topologies.3line', 'deeptmhmm_results.md'}
        observed = {row['relative_path'] for row in provenance.get('raw_output_inventory', [])}
        provenance['required_outputs_present'] = required.issubset(observed)
        success = (
            status == 'completed'
            and provenance.get('exit_code') == 0
            and provenance['required_outputs_present']
        )
        provenance['technical_attempt_success'] = success
        write_json_atomic(provenance_path, provenance)
        return success, provenance
    except Exception as exc:
        provenance['completion_time_utc'] = now()
        provenance['status'] = 'CLIENT_EXCEPTION'
        provenance['client_exception'] = repr(exc)
        provenance['client_exception_detail'] = str(exc)
        provenance['technical_attempt_success'] = False
        write_json_atomic(provenance_path, provenance)
        return False, provenance

def run_batch(batch_id):
    input_path = INPUT_ROOT / f'{batch_id}.fasta'
    assert input_path.exists()
    batch_dir = JOB_ROOT / batch_id
    batch_dir.mkdir(parents=True, exist_ok=True)
    final_path = batch_dir / 'batch_result.json'
    if final_path.exists():
        existing = json.loads(final_path.read_text())
        if existing.get('final_status') == 'SUCCESS':
            return existing
    attempts = []
    existing_attempt_numbers = []
    for provenance_path in sorted(batch_dir.glob('attempt_*/job_provenance.json')):
        prior = json.loads(provenance_path.read_text())
        attempt_number = int(prior['attempt_number'])
        existing_attempt_numbers.append(attempt_number)
        attempts.append({
            'attempt_number': attempt_number,
            'job_id': prior.get('job_id'),
            'status': prior.get('status'),
            'technical_attempt_success': bool(prior.get('technical_attempt_success', False)),
        })
    next_attempt = max(existing_attempt_numbers, default=0) + 1
    for attempt_number in range(next_attempt, next_attempt + MAX_NEW_ATTEMPTS_PER_CONTROLLER_PASS):
        success, provenance = run_attempt(batch_id, input_path, attempt_number)
        attempts.append({
            'attempt_number': attempt_number,
            'job_id': provenance.get('job_id'),
            'status': provenance.get('status'),
            'technical_attempt_success': success,
        })
        if success:
            result = {
                'batch_id': batch_id,
                'input_file': str(input_path),
                'input_sha256': sha256_file(input_path),
                'input_records': fasta_count(input_path),
                'final_status': 'SUCCESS',
                'successful_attempt': attempt_number,
                'source_job_id': provenance.get('job_id'),
                'raw_output_directory': str(batch_dir / f'attempt_{attempt_number}' / 'raw'),
                'attempts': attempts,
            }
            write_json_atomic(final_path, result)
            return result
        if provenance.get('status') == 'CLIENT_EXCEPTION':
            time.sleep(60)
    result = {
        'batch_id': batch_id,
        'input_file': str(input_path),
        'input_sha256': sha256_file(input_path),
        'input_records': fasta_count(input_path),
        'final_status': 'EXPLICIT_FAILURE',
        'successful_attempt': None,
        'source_job_id': None,
        'raw_output_directory': None,
        'attempts': attempts,
    }
    write_json_atomic(final_path, result)
    return result

def rebuild_aggregate():
    rows = []
    for result_path in sorted(JOB_ROOT.glob('batch_*/batch_result.json')):
        result = json.loads(result_path.read_text())
        rows.append({
            'batch_id': result['batch_id'],
            'input_file': result['input_file'],
            'input_sha256': result['input_sha256'],
            'input_records': result['input_records'],
            'final_status': result['final_status'],
            'successful_attempt': result['successful_attempt'] if result['successful_attempt'] is not None else '',
            'source_job_id': result['source_job_id'] or '',
            'raw_output_directory': result['raw_output_directory'] or '',
            'attempt_count': len(result['attempts']),
        })
    aggregate = QC_ROOT / 'deeptmhmm_full_job_provenance.tsv'
    fields = ['batch_id','input_file','input_sha256','input_records','final_status','successful_attempt','source_job_id','raw_output_directory','attempt_count']
    tmp = aggregate.with_suffix('.tsv.tmp')
    with tmp.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, aggregate)
    return rows

def main():
    global APP_OBJECT
    parser = argparse.ArgumentParser()
    parser.add_argument('--start', type=int, default=1)
    parser.add_argument('--end', type=int, default=189)
    parser.add_argument('--concurrency', type=int, default=1)
    args = parser.parse_args()
    assert 1 <= args.start <= args.end <= 189
    assert 1 <= args.concurrency <= 8
    client = BiolibApiClient.get()
    if not client.is_signed_in:
        raise RuntimeError('BioLib authentication not loaded; refusing to submit anonymously')
    APP_OBJECT = biolib.load(APP_URI)
    batch_ids = [f'batch_{n:04d}' for n in range(args.start, args.end + 1)]
    with ThreadPoolExecutor(max_workers=args.concurrency) as executor:
        futures = {executor.submit(run_batch, batch_id): batch_id for batch_id in batch_ids}
        for future in as_completed(futures):
            batch_id = futures[future]
            try:
                result = future.result()
                print(json.dumps({'batch_id': batch_id, 'final_status': result['final_status'], 'job_id': result.get('source_job_id')}), flush=True)
            except Exception as exc:
                print(json.dumps({'batch_id': batch_id, 'controller_exception': repr(exc)}), file=sys.stderr, flush=True)
            rebuild_aggregate()
    rows = rebuild_aggregate()
    success = sum(int(r['input_records']) for r in rows if r['final_status'] == 'SUCCESS')
    failed = sum(int(r['input_records']) for r in rows if r['final_status'] == 'EXPLICIT_FAILURE')
    print(json.dumps({'recorded_batches': len(rows), 'successful_input_records': success, 'explicit_failure_records': failed}), flush=True)

if __name__ == '__main__':
    main()
