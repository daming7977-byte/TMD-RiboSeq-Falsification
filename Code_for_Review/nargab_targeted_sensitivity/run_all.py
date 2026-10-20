#!/usr/bin/env python3
"""Run only the packaged targeted sensitivity in an explicit new output directory."""
from pathlib import Path
import argparse,csv,hashlib,json,os,shutil,subprocess,sys,time
from datetime import datetime,timezone

P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',required=True,help='New, nonexistent output directory; existing paths are never overwritten.')
parser.add_argument('--verify-existing',action='store_true',help='Verify a completed existing calculation without re-executing it; preserves calculation outputs.')
args=parser.parse_args(); O=Path(args.output_dir).expanduser().resolve()
if O.exists() and not args.verify_existing:raise SystemExit('Refusing to overwrite existing output directory: '+str(O))
if args.verify_existing and not (O/'evidence/calculation_complete.json').is_file():raise SystemExit('Existing calculation is not complete.')
if O==P or P in O.parents:raise SystemExit('Output must be outside the immutable sensitivity extension directory.')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for row in csv.DictReader((P/'EXTENSION_FILE_SHA256.tsv').open(),delimiter='\t'):
 p=P/row['path']
 assert p.is_file() and p.stat().st_size==int(row['bytes']) and sha(p)==row['sha256'],('Extension integrity failure',row['path'])
for row in csv.DictReader((P/'INPUT_MANIFEST.tsv').open(),delimiter='\t'):
 p=P.parent/row['relative_to_code_package']
 assert p.is_file() and p.stat().st_size==int(row['bytes']) and sha(p)==row['sha256'],('Input integrity failure',row['relative_to_code_package'])
reference=json.loads((P/'FROZEN_REFERENCE_SHA256.json').read_text())
rp=P.parent/reference['relative_to_code_package']
assert rp.stat().st_size==reference['bytes'] and sha(rp)==reference['sha256'], 'Frozen reference checksum differs'
started=time.perf_counter()
if not args.verify_existing:
 O.mkdir(parents=True);(O/'evidence').mkdir();(O/'logs').mkdir()
 for rel in ['P3_WITHIN_GENE_PRESPECIFIED_PROTOCOL.md','evidence/PROTOCOL_FREEZE.json','evidence/PRE_RESULT_IMPLEMENTATION_CONVENTIONS.md']:shutil.copyfile(P/rel,O/rel)
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OMP_NUM_THREADS='1',P3_SENSITIVITY_OUTPUT=str(O))
 for script in ['prepare_inputs.py','run_sensitivity.py']:
  command=[sys.executable,str(P/'scripts'/script)]
  with (O/'logs'/(script+'.log')).open('w') as log:
   process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,env=env)
   for line in process.stdout:log.write(line);log.flush();print(line,end='',flush=True)
   assert process.wait()==0,('Execution failed',script)
# Comparisons do not alter results or tune the calculation.
import numpy as np
import pandas as pd
checks=[]
def compare_tsv(actual,expected):
 a=pd.read_csv(actual,sep='\t',dtype=str,keep_default_na=False);b=pd.read_csv(expected,sep='\t',dtype=str,keep_default_na=False)
 assert list(a.columns)==list(b.columns) and a.shape==b.shape,('Table layout mismatch',actual.name,a.shape,b.shape)
 max_error=0.;numeric_cells=0;exact_cells=0
 for col in a:
  try:x=pd.to_numeric(a[col]);y=pd.to_numeric(b[col])
  except (ValueError,TypeError):
   assert a[col].equals(b[col]),('Exact field differs',actual.name,col)
   exact_cells+=len(a);continue
  xv=x.to_numpy(dtype=float);yv=y.to_numpy(dtype=float)
  assert np.array_equal(np.isnan(xv),np.isnan(yv)) and np.array_equal(np.isposinf(xv),np.isposinf(yv)) and np.array_equal(np.isneginf(xv),np.isneginf(yv)),('Missing/nonfinite field differs',actual.name,col)
  finite=np.isfinite(xv);error=float(np.max(np.abs(xv[finite]-yv[finite]))) if finite.any() else 0.
  assert error<=1e-10,('Numeric field differs',actual.name,col,error)
  max_error=max(max_error,error);numeric_cells+=int(finite.sum())
 checks.append({'file':str(actual.relative_to(O)),'rows':len(a),'columns':len(a.columns),'numeric_cells_compared':numeric_cells,'exact_text_cells_compared':exact_cells,'max_absolute_error':max_error,'absolute_tolerance':1e-10,'status':'PASS'})
compare_tsv(O/'P3_WITHIN_GENE_INDEPENDENT_RESULTS.tsv',P/'expected/calculation_results.tsv')
compare_tsv(O/'P3_WITHIN_GENE_ELIGIBILITY_AND_BALANCE.tsv',P/'expected/eligibility_and_balance.tsv')
for name in ['balance_diagnostics.tsv','common_support_strata.tsv','balance_weights.tsv','gene_deletion_influence.tsv','effective_density_contributions.tsv','tmd_count_stratum_effects.tsv','frozen_reproduction_comparison.tsv','within_gene_contrasts_PLATEAU_90_120.tsv','within_gene_contrasts_BROAD_60_120.tsv','positivity_exclusions_PLATEAU_90_120.tsv','positivity_exclusions_BROAD_60_120.tsv']:
 compare_tsv(O/'evidence'/name,P/'expected'/name)
exact=[]
for name,expected in json.loads((P/'expected/EXACT_OUTPUT_SHA256.json').read_text()).items():
 actual=sha(O/'evidence'/name);exact.append({'file':name,'expected_sha256':expected,'actual_sha256':actual,'byte_identical':actual==expected})
# Metadata populations must be byte-identical. Bootstrap-array byte identity is an
# additional reference-platform result; numerical table comparisons above are mandatory.
assert all(c['byte_identical'] for c in exact if not c['file'].startswith('bootstrap_'))
result={'status':'PASS','completed_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-started,'verification_only':args.verify_existing,'source_protocol_sha256':sha(P/'P3_WITHIN_GENE_PRESPECIFIED_PROTOCOL.md'),'variants':30,'draws_per_variant':10000,'table_comparisons':checks,'exact_output_comparisons':exact,'all_reference_outputs_byte_identical':all(c['byte_identical'] for c in exact),'scope':'Re-execution of the previously independently implemented processed-input audit, with path/verification adapters only. Not new biological replication, raw reconstruction, empirical confidence coverage validation, or independently programmed verification of every bootstrap draw.'}
(O/'REVIEWER_REPRODUCTION_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
with (O/'OUTPUT_SHA256.tsv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t');writer.writeheader()
 for p in sorted(O.rglob('*')):
  if p.is_file() and p.name!='OUTPUT_SHA256.tsv':writer.writerow({'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':sha(p)})
print(json.dumps({'status':result['status'],'variants':30,'checked_tables':len(checks),'all_reference_outputs_byte_identical':result['all_reference_outputs_byte_identical'],'elapsed_seconds':result['elapsed_seconds'],'report':str(O/'REVIEWER_REPRODUCTION_CHECK.json')},indent=2),flush=True)
