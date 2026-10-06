"""Execute only frozen Task13 estimands in explicit freeze order."""
from pathlib import Path
import json,csv,hashlib,importlib.util,subprocess,sys
import numpy as np
R=Path(__file__).resolve().parent;D=R.parent/'ROUTE_A_BENCHMARK_DESIGN_20260910';T=R.parent/'A3_ALIGNMENT_CALIBRATION_20260910'
NODE='node'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def freeze(names,receipt):
 rows=[]
 for n in names:
  p=R/n;rows.append(sha(p)+'  '+n+'\n');p.chmod(0o444)
 (R/receipt).write_text(''.join(rows));(R/receipt).chmod(0o444)
def author(tables,stage):
 p=R/(stage+'_tables.json');p.write_text(json.dumps(tables)+'\n')
 subprocess.run([NODE,str(R/'author_stage.mjs'),str(p)],cwd=R,check=True)
def table(p):return list(csv.DictReader(p.open(),delimiter='\t'))
# Verify every frozen Task13 file before any numerical reconstruction.
inputs=[]
for line in (D/'SHA256SUMS.txt').read_text().splitlines():
 h,n=line.split('  ',1);p=(D/n).resolve();assert sha(p)==h,str(p);inputs.append([str(p),p.stat().st_size,h,'TASK13_FROZEN_INPUT'])
print('All Task13 frozen input hashes verified',flush=True)
mem=json.loads((D/'SUPPORT_MEMBERSHIP_FROZEN.json').read_text());PG=mem['primary_common_genes'];SG=mem['seven_dose_common_genes'];assert len(PG)==304 and len(SG)==280 and PG==sorted(PG) and SG==sorted(SG)
anchors={g:[] for g in PG}
for r in table(D/'ROUTE_A_REAL_ANCHOR_ELIGIBILITY_v0.tsv'):
 if r['eligibility']=='ELIGIBLE' and r['gene_locus'] in anchors:anchors[r['gene_locus']].append((r['anchor_id'],int(r['anchor_cds_codon_0based'])))
pseudo={}
for p in table(D/'ROUTE_A_PSEUDO_MATCHES_v0.tsv'):pseudo.setdefault(p['real_anchor_id'],[]).append(int(p['pseudo_cds_codon_0based']))
assert sum(map(len,anchors.values()))==946 and sum(len(anchors[g]) for g in SG)==879
assert all(len(pseudo[a])==5 for v in anchors.values() for a,_ in v)
counts_frozen=json.loads((D/'CDS_WIDE_TECHNICAL_SUPPORT.json').read_text())
old_manifest={n:h for h,n in (line.split('  ',1) for line in (T/'SHA256SUMS.txt').read_text().splitlines())}
assert sha(T/'b2_qc_definitions.py')==old_manifest['b2_qc_definitions.py']
spec=importlib.util.spec_from_file_location('base',T/'b2_qc_definitions.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base);base.ROOT=T;models,maps=base.arrays()
model={m['name']:m for m in models};K={g:model[g]['nt']//3 for g in PG}
RUNS=[f'SRR{1520311+i}' for i in range(8)];CONDITIONS=['noCHX','1/64x','1/16x','1/4x','1x','8x','100x_rep1','100x_rep2'];STATES=CONDITIONS[:6]+['100x']
END=[0,6,7];account=[]
def reconstruct(i,L):
 run=RUNS[i];o={27:14,28:15}[L];out={g:np.zeros(K[g],dtype=np.int64) for g in PG};group=[];old=None
 def finish(rows):
  if not rows or len(rows[0][9])!=L:return
  mapped=[a for a in rows if not int(a[1])&4]
  if len(mapped)!=1:return
  a=mapped[0];assert a[5]==f'{L}M' and 'NM:i:0' in a[11:]
  rev=bool(int(a[1])&16);p=int(a[3])-1+(L-1-o if rev else o);ids,qs,_,_=maps[a[2],rev]
  if not 0<=p<len(ids):return
  idx=int(ids[p]);q=int(qs[p])
  if idx>0 and q>=0 and q%3==0:
   g=models[idx-1]['name']
   if g in out:out[g][q//3]+=1
 with (T/'technical'/(run+'.genome.sam')).open() as f:
  for line in f:
   if line.startswith('@'):continue
   a=line.rstrip('\n').split('\t')
   if a[0]!=old:finish(group);group=[];old=a[0]
   group.append(a)
  finish(group)
 for g,v in out.items():assert int(v.sum())==counts_frozen[run][str(L)][g],(run,L,g)
 np.savez_compressed(R/'technical'/f'{run}.{L}nt.CDS_counts.npz',**out)
 account.append([run,CONDITIONS[i],L,o,len(out),sum(int(v.sum()) for v in out.values()),'EXACT_MATCH_TO_TASK13_CDS_TOTALS','NO_NEW_EXCLUSIONS'])
 print(f'{run} {L}nt reconstructed and whole-CDS totals verified',flush=True)
 return out
def gene_specificity(out,genes):
 vals=[]
 for g in genes:
  v=out[g];N=int(v.sum());assert N>0;X=v.astype(np.float64)*K[g]/N
  def d(a):return np.sum(X[a+30:a+61],dtype=np.float64)/31-np.sum(X[a-30:a+1],dtype=np.float64)/31
  vals.append(np.mean([d(a)-np.mean([d(q) for q in pseudo[aid]],dtype=np.float64) for aid,a in anchors[g]],dtype=np.float64))
 return np.asarray(vals,dtype=np.float64)
def bootstrap(v,seed):
 rng=np.random.Generator(np.random.PCG64(seed));n=len(v);assert n>=2
 draws=np.array([np.mean(v[rng.integers(0,n,size=n)],dtype=np.float64) for _ in range(10000)],dtype=np.float64)
 ci=np.quantile(draws,[.025,.975],method='linear');return float(np.mean(v)),float(ci[0]),float(ci[1]),draws
verdict=lambda lo,hi:'TECHNICAL_PERTURBATION_SUPPORT' if lo>0 or hi<0 else 'NO_DIRECTIONAL_CHANGE_SUPPORT'
results={};cache={};gs={}
header=['footprint_nt','A_offset_nt','genes','anchors','endpoint_statistic','CI95_lower','CI95_upper','bootstrap_draws','seed','frozen_verdict','null_interpretation','sign_vs_primary','sensitivity_assessment']
for L in [27,28]:
 if L==28:assert (R/'PRIMARY_OUTPUT_SHA256.txt').exists()
 cache[L]={};gs[L]={}
 for i in END:
  cache[L][i]=reconstruct(i,L);gs[L][i]=gene_specificity(cache[L][i],PG)
 endpoint=(gs[L][6]+gs[L][7])/2-gs[L][0];theta,lo,hi,draws=bootstrap(endpoint,5957314);v=verdict(lo,hi)
 concord='PRIMARY' if L==27 else 'CONCORDANT' if np.sign(theta)==np.sign(results[27]['statistic']) else 'DISCORDANT'
 assessment='PRIMARY' if L==27 else 'DISAGREES' if concord=='DISCORDANT' or v!=results[27]['verdict'] else 'SUPPORTS' if v=='TECHNICAL_PERTURBATION_SUPPORT' else 'INCONCLUSIVE_CONCORDANT_NULL'
 results[L]=dict(statistic=theta,lower=lo,upper=hi,verdict=v,sign=concord,assessment=assessment)
 np.save(R/'technical'/f'endpoint_{L}nt_bootstrap.npy',draws)
 name='ROUTE_A_PRIMARY_RESULT_v0.tsv' if L==27 else 'ROUTE_A_28NT_SENSITIVITY_RESULT_v0.tsv'
 tables={name:[header,[L,{27:14,28:15}[L],304,946,theta,lo,hi,10000,5957314,v,'NO_SUPPORTED_ENDPOINT_CHANGE' if v=='NO_DIRECTIONAL_CHANGE_SUPPORT' else 'NOT_APPLICABLE',concord,assessment]]}
 if L==27:
  geneheader=['gene_locus','eligible_anchors','noCHX_specificity','100x_rep1_specificity','100x_rep2_specificity','100x_equal_experiment_mean','endpoint_delta']
  tables['ROUTE_A_PRIMARY_GENE_LEVEL_v0.tsv']=[geneheader]+[[g,len(anchors[g]),float(gs[L][0][j]),float(gs[L][6][j]),float(gs[L][7][j]),float((gs[L][6][j]+gs[L][7][j])/2),float(endpoint[j])] for j,g in enumerate(PG)]
 else:
  (R/'technical/sensitivity_gene_values.json').write_text(json.dumps({'genes':PG,'noCHX':gs[L][0].tolist(),'rep1':gs[L][6].tolist(),'rep2':gs[L][7].tolist(),'delta':endpoint.tolist()})+'\n')
 author(tables,'primary' if L==27 else 'sensitivity')
 freeze(list(tables)+[f'technical/endpoint_{L}nt_bootstrap.npy'], 'PRIMARY_OUTPUT_SHA256.txt' if L==27 else 'SENSITIVITY_OUTPUT_SHA256.txt')
 print('FROZEN_ENDPOINT '+str(L)+' '+json.dumps(results[L]),flush=True)
# Explicit frozen disagreement stopping rule: do not continue secondary/visualization on discordance.
if results[28]['assessment']=='DISAGREES':
 results['stopped']='PRIMARY_SENSITIVITY_DISCORDANCE'
 author({'ROUTE_A_SEVEN_DOSE_SECONDARY_v0.tsv':[['status','reason'],['NOT_EXECUTED','Frozen disagreement stopping rule']], 'ROUTE_A_CONDITION_LEVEL_DESCRIPTIVE_v0.tsv':[['status','reason'],['NOT_EXECUTED','Frozen disagreement stopping rule']]},'stopped')
else:
 sec=[['footprint_nt','A_offset_nt','genes','anchors','mean_ordinal_slope','CI95_lower','CI95_upper','draws','seed','sign_vs_endpoint','CI_zero_exclusion_vs_endpoint']];desc=[['footprint_nt','dose_rank','dose_state','genes','anchors','equal_gene_specificity','purpose']];results['secondary']={}
 for L in [27,28]:
  for i in range(8):
   if i not in cache[L]:cache[L][i]=reconstruct(i,L)
  mat=np.stack([gene_specificity(cache[L][i],SG) for i in range(8)],axis=1);states=np.column_stack([mat[:,:6],(mat[:,6]+mat[:,7])/2]);slopes=np.sum(states*(np.arange(7,dtype=np.float64)-3),axis=1)/28
  theta,lo,hi,draws=bootstrap(slopes,5957315);np.save(R/'technical'/f'ordinal_{L}nt_bootstrap.npy',draws)
  sign='CONCORDANT' if np.sign(theta)==np.sign(results[L]['statistic']) else 'DISCORDANT';excl=(lo>0 or hi<0)==(results[L]['verdict']=='TECHNICAL_PERTURBATION_SUPPORT')
  sec.append([L,{27:14,28:15}[L],280,879,theta,lo,hi,10000,5957315,sign,'CONCORDANT' if excl else 'DISCORDANT'])
  desc.extend([[L,i,c,280,879,float(np.mean(states[:,i])), 'DESCRIPTIVE_ONLY_NO_PAIRWISE_TEST'] for i,c in enumerate(STATES)])
  results['secondary'][L]=dict(statistic=theta,lower=lo,upper=hi,sign_vs_endpoint=sign,CI_agreement=excl)
  (R/'technical'/f'ordinal_{L}nt_gene_values.json').write_text(json.dumps({'genes':SG,'states':states.tolist(),'slopes':slopes.tolist()})+'\n')
 author({'ROUTE_A_SEVEN_DOSE_SECONDARY_v0.tsv':sec,'ROUTE_A_CONDITION_LEVEL_DESCRIPTIVE_v0.tsv':desc},'secondary')
 freeze(['ROUTE_A_SEVEN_DOSE_SECONDARY_v0.tsv','ROUTE_A_CONDITION_LEVEL_DESCRIPTIVE_v0.tsv','technical/ordinal_27nt_bootstrap.npy','technical/ordinal_28nt_bootstrap.npy'],'SECONDARY_OUTPUT_SHA256.txt')
 print('FROZEN_SECONDARY '+json.dumps(results['secondary']),flush=True)
results['inputs_verified']=len(inputs)
(R/'NUMERICAL_RESULTS_FROZEN.json').write_text(json.dumps(results,indent=2)+'\n');freeze(['NUMERICAL_RESULTS_FROZEN.json'],'ALL_NUMERICAL_SHA256.txt')
author({'ROUTE_A_EXECUTION_ACCOUNTING_v0.tsv':[['run','condition','footprint_nt','A_offset_nt','primary_support_genes_reconstructed','accepted_CDS_reads_in_primary_gene_set','support_reconciliation','execution_exclusions']]+account,'ROUTE_A_EXECUTION_INPUT_MANIFEST_v0.tsv':[['path','bytes','SHA256','role']]+inputs+[[str(R/'run_benchmark.py'),(R/'run_benchmark.py').stat().st_size,sha(R/'run_benchmark.py'),'EXECUTION_IMPLEMENTATION']]},'accounting')
print('NUMERICAL_EXECUTION_COMPLETE',flush=True)
