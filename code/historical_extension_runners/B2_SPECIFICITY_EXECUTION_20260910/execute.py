"""Execute frozen B2 design. Numerical results first; visualization is a separate step."""
import csv,json,hashlib,datetime,collections
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent;T6=P.parent/'B2_ALIGNMENT_CALIBRATION_20260910';T7=P.parent/'B2_BENCHMARK_DESIGN_20260910'
assert np.__version__=='2.3.5'
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
inputs=[]
for root in (T6,T7):
 # Verify all upstream frozen files, retaining exact absolute identity manifest.
 for line in (root/'SHA256SUMS.txt').read_text().splitlines():
  sha,rel=line.split('  ',1);f=root/rel
  assert digest(f)==sha, str(f)
  inputs.append([str(f),sha,'UPSTREAM_FROZEN_FILE',root.name])
 inputs.append([str(root/'SHA256SUMS.txt'),digest(root/'SHA256SUMS.txt'),'UPSTREAM_MANIFEST',root.name])
print('All Task 6 and Task 7 hashes verified',flush=True)
models=json.loads((T6/'references/coding_models.json').read_text())
models=[m for m in models if m['chrom']!='chrM'];names=[m['name'] for m in models];assert len(set(names))==len(names)
lengths={}
for l in (T6/'references/genome.fa').read_text().splitlines():
 if l.startswith('>'):chrom=l[1:].split()[0];lengths[chrom]=0
 else:lengths[chrom]+=len(l)
maps={(c,s):(np.zeros(n,dtype=np.int32),np.full(n,-1,dtype=np.int32)) for c,n in lengths.items() for s in '+-'}
for i,m in enumerate(models):
 gid,qmap=maps[m['chrom'],m['strand']];off=0
 ex=sorted(m['exons'],reverse=m['strand']=='-')
 for a,b in ex:
  old=gid[a:b];free=old==0;old[~free]=-1;old[free]=i+1
  q=np.arange(off,off+b-a,dtype=np.int32);qmap[a:b]=q[::-1] if m['strand']=='-' else q;off+=b-a
 assert off==m['nt']
counts={(L,frac):[np.zeros(m['nt']//3-1,dtype=np.int64) for m in models] for L in (28,32) for frac in ('input','pulldown')}
account=[]
for run,frac in [('SRR1562895','input'),('SRR1562896','input'),('SRR1562897','pulldown'),('SRR1562898','pulldown')]:
 metrics={L:collections.Counter() for L in (28,32)};total=0
 def finish(g):
  global total
  if not g:return
  total+=1;L=len(g[0][9])
  if L not in metrics:return
  d=metrics[L];d['entering_genome_alignment']+=1
  aligned=[a for a in g if not int(a[1])&4]
  if not aligned:d['unmapped']+=1;return
  d['mapped']+=1
  if len(aligned)>1:d['multimapped_rejected']+=1;return
  a=aligned[0];flag=int(a[1]);assert not flag&256 and not flag&2048
  assert a[5]==f'{L}M' and 'NM:i:0' in a[11:]
  d['unique_exact']+=1;s='-' if flag&16 else '+'
  left=int(a[3])-1;site=left+(L-1 if s=='-' else 0)+(15 if L==28 else 16)*(1 if s=='+' else -1)
  assert 0<=site<lengths[a[2]]
  gid,qmap=maps[a[2],s];i=int(gid[site]);q=int(qmap[site])
  if i==0:d['outside_nuclear_sense_CDS']+=1;return
  if i<0:d['ambiguous_sense_CDS']+=1;return
  if q%3:d['off_frame']+=1;return
  j=q//3
  if j==len(counts[L,frac][i-1]):d['stop_codon']+=1;return
  assert 0<=j<len(counts[L,frac][i-1]);counts[L,frac][i-1][j]+=1;d['assigned_translated_CDS_N']+=1
 with (T6/'technical'/f'{run}.genome.sam').open() as f:
  group=[];old=None
  for line in f:
   if line.startswith('@'):continue
   a=line.rstrip('\n').split('\t')
   if a[0]!=old:finish(group);group=[];old=a[0]
   group.append(a)
  finish(group)
 prior=json.loads((T6/'technical'/f'{run}.qc.json').read_text())
 assert total==prior['totals']['entering']
 for L,d in metrics.items():
  assert d['unique_exact']==prior['lengths'][str(L)]['unique']
  assert d['mapped']==prior['lengths'][str(L)]['mapped']
  assert d['entering_genome_alignment']==d['unmapped']+d['mapped']
  assert d['mapped']==d['unique_exact']+d['multimapped_rejected']
  assert d['unique_exact']==sum(d[k] for k in ('outside_nuclear_sense_CDS','ambiguous_sense_CDS','off_frame','stop_codon','assigned_translated_CDS_N'))
  for k in ('entering_genome_alignment','mapped','unmapped','unique_exact','multimapped_rejected','outside_nuclear_sense_CDS','ambiguous_sense_CDS','off_frame','stop_codon','assigned_translated_CDS_N'):
   account.append([run,'GSM1495513' if frac=='input' else 'GSM1495514',frac,L,k,d[k]])
 print(run+' reconstructed; technical read accounting reconciled',flush=True)
anchors=[r for r in rows(T7/'B2_REAL_ANCHOR_ELIGIBILITY_v0.tsv') if r['eligibility']=='ELIGIBLE']
ps=collections.defaultdict(list)
for r in rows(T7/'B2_FINAL_PSEUDO_ANCHORS_v0.tsv'):ps[r['real_anchor_id']].append(int(r['pseudo_cds_codon_0based']))
ga=collections.defaultdict(list)
for a in anchors:ga[a['gene_locus']].append(a)
genes=sorted(ga);G=len(genes);assert G==338 and len(anchors)==1068 and sum(map(len,ps.values()))==5340
assert all(len(ps[a['anchor_id']])==5 for a in anchors)
idx={n:i for i,n in enumerate(names)}
draws=np.random.Generator(np.random.PCG64(6101208)).integers(0,G,size=(10000,G),dtype=np.int64)
tables={};results={};visual_input={};gene_data={};boots={}
result_header=['length_nt','A_site_offset_nt','genes','real_anchors','pseudo_definitions','N_input','N_pulldown','statistic','bootstrap_CI95_lower','bootstrap_CI95_upper','frozen_one_sided95_lower','frozen_directional_verdict','bootstrap_draws','bootstrap_seed','new_support_exclusions','biological_replication']
for L in (28,32):
 N={f:sum(int(c.sum()) for c in counts[L,f]) for f in ('input','pulldown')}
 assert min(N.values())>0
 for f in N:
  assert N[f]==sum(r[5] for r in account if r[2]==f and r[3]==L and r[4]=='assigned_translated_CDS_N')
 E={g:np.log2((1e6*(counts[L,'pulldown'][idx[g]]+.5)/N['pulldown'])/(1e6*(counts[L,'input'][idx[g]]+.5)/N['input'])) for g in genes}
 def z(v,a):return float(np.mean(v[a+30:a+61])-np.mean(v[a-30:a+1]))
 contributions=[];gr=[]
 for g in genes:
  real=[];pseudo=[]
  for a in ga[g]:
   q=int(a['anchor_cds_codon_0based']);v=E[g]
   assert q-60>=0 and q+100<len(v)
   real.append(z(v,q));pseudo.append(float(np.mean([z(v,p) for p in ps[a['anchor_id']]])))
  d=float(np.mean(np.array(real)-np.array(pseudo)));contributions.append(d)
  gr.append([g,len(real),float(np.mean(real)),float(np.mean(pseudo)),d])
 T=float(np.mean(contributions));boot=np.array(contributions)[draws].mean(axis=1)
 assert np.isfinite(boot).all()
 low,high=np.quantile(boot,[.025,.975],method='linear');lower=np.quantile(boot,.05,method='linear')
 verdict='DIRECTIONAL_SUPPORT' if lower>0 else 'NO_DIRECTIONAL_SUPPORT'
 row=[L,15 if L==28 else 16,G,len(anchors),5340,N['input'],N['pulldown'],T,float(low),float(high),float(lower),verdict,10000,6101208,0,'NO / NOT ESTABLISHED']
 name='B2_PRIMARY_SPECIFICITY_RESULT_v0.tsv' if L==28 else 'B2_32NT_SENSITIVITY_RESULT_v0.tsv'
 tables[name]=[result_header,row];results[str(L)]=dict(zip(result_header,row));gene_data[L]=gr;boots[L]=boot
 # Freeze numerical result before any aggregate visualization, including primary before sensitivity.
 frozen={'result':results[str(L)],'gene_statistics':gr,'bootstrap_sha256':hashlib.sha256(boot.tobytes()).hexdigest(),'gene_draw_indices_sha256':hashlib.sha256(draws.tobytes()).hexdigest(),'code_sha256':digest(P/'execute.py'),'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'central_CI_note':'Task8 requested central 95% interval added as reporting summary; frozen decision uses 5th percentile'}
 target=P/f'NUMERICAL_{L}NT_FREEZE.json';target.write_text(json.dumps(frozen,indent=2)+'\n');target.chmod(0o444)
 np.save(P/f'bootstrap_{L}nt.npy',boot)
 # Only vectors for already-eligible genes are retained for the separately authorized display.
 for g in genes:visual_input[f'{L}_{g}']=E[g]
 print(f'{L} nt numerical result frozen',flush=True)
tables['B2_PRIMARY_GENE_LEVEL_STATISTIC_v0.tsv']=[['gene_locus','eligible_real_anchors','mean_real_Z','mean_matched_pseudo_Z','paired_gene_specificity']]+gene_data[28]
tables['B2_PRIMARY_BOOTSTRAP_v0.tsv']=[['draw_1based','primary_28nt','sensitivity_32nt']]+[[i+1,float(boots[28][i]),float(boots[32][i])] for i in range(10000)]
tables['B2_EXECUTION_ACCOUNTING_v0.tsv']=[['run','GSM','fraction','length_nt','metric','read_count']]+account
tables['B2_EXECUTION_SOURCE_AND_INPUT_MANIFEST_v0.tsv']=[['input_path','sha256','role','source_stage']]+inputs
(P/'tables.json').write_text(json.dumps(tables));(P/'results.json').write_text(json.dumps(results,indent=2)+'\n')
np.savez_compressed(P/'visualization_inputs.npz',**visual_input)
print(json.dumps(results,indent=2),flush=True)
