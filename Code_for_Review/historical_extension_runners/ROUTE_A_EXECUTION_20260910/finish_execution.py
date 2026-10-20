"""Only frozen-window descriptive visualization after numerical freeze; final accounting."""
from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
R=Path(__file__).resolve().parent;D=R.parent/'ROUTE_A_BENCHMARK_DESIGN_20260910'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
for receipt in ['PRIMARY_OUTPUT_SHA256.txt','SENSITIVITY_OUTPUT_SHA256.txt','ALL_NUMERICAL_SHA256.txt']:
 for line in (R/receipt).read_text().splitlines():
  h,n=line.split('  ',1);assert sha(R/n)==h
res=json.loads((R/'NUMERICAL_RESULTS_FROZEN.json').read_text())
lines=['# Route A frozen execution verdict v0','', 'Role: KNOWN TECHNICAL-PERTURBATION / METHOD-EVALUATION. No established biological replication; whole-gene resampling is not experimental replication.','']
for L in ['27','28']:
 x=res[L];lines.append(f"{L}nt endpoint, 304 genes / 946 anchors: statistic {x['statistic']:.17g}; two-sided95%CI [{x['lower']:.17g}, {x['upper']:.17g}]; frozen label {x['verdict']}. Sign status {x['sign']}; assessment {x['assessment']}.")
 if x['verdict']=='NO_DIRECTIONAL_CHANGE_SUPPORT':lines.append('Interpretation: NO_SUPPORTED_ENDPOINT_CHANGE for this frozen specificity metric, not an assertion that CHX has no effects on ribosome profiling.')
 else:lines.append('Interpretation: supported change in this frozen real-TMD-versus-matched-pseudo specificity metric at100x relative to noCHX under this reconstruction. No mechanism or biological replication established.')
if 'stopped' in res:
 lines += ['', 'STOP: '+res['stopped']+'. The frozen design treats a different point-estimate sign or different support label as technical discordance. The primary is retained unchanged. Secondary and visualization were not executed because the explicit disagreement stopping rule applies. No rescue analysis.']
 (R/'ROUTE_A_FROZEN_VISUALIZATION_v0.md').write_text('NOT_GENERATED: primary/sensitivity discordance triggered the frozen stopping rule. No descriptive profiles calculated.\n')
else:
 for line in (R/'SECONDARY_OUTPUT_SHA256.txt').read_text().splitlines():
  h,n=line.split('  ',1);assert sha(R/n)==h
 for L,x in res['secondary'].items():lines += ['',f"{L}nt secondary, 280 genes / 879 anchors: ordinal mean slope {x['statistic']:.17g};95%CI [{x['lower']:.17g}, {x['upper']:.17g}]; sign versus endpoint {x['sign_vs_endpoint']}; interval-zero-exclusion agreement {x['CI_agreement']}. Prespecified secondary only; never replaces primary."]
 lines += ['', 'Condition-level values use the fixed280-gene/879-anchor seven-dose population for both footprints; primary endpoint uses304genes/946anchors. These populations are not silently interchanged. No additional dose pairwise tests were performed.']
 # All numerical results are already immutable before any profile calculation.
 import subprocess
 genes=json.loads((D/'SUPPORT_MEMBERSHIP_FROZEN.json').read_text())['seven_dose_common_genes'];anchors={g:[] for g in genes};pseudo={}
 for r in csv.DictReader((D/'ROUTE_A_REAL_ANCHOR_ELIGIBILITY_v0.tsv').open(),delimiter='\t'):
  if r['eligibility']=='ELIGIBLE' and r['gene_locus'] in anchors:anchors[r['gene_locus']].append((r['anchor_id'],int(r['anchor_cds_codon_0based'])))
 for r in csv.DictReader((D/'ROUTE_A_PSEUDO_MATCHES_v0.tsv').open(),delimiter='\t'):pseudo.setdefault(r['real_anchor_id'],[]).append(int(r['pseudo_cds_codon_0based']))
 x=np.arange(-60,101);profiles={};labels=['noCHX','1/64x CHX','1/16x CHX','1/4x CHX','1x CHX','8x CHX','100x CHX (equal experiment mean)']
 for L in [27,28]:
  real=[];ps=[]
  for i in range(8):
   data=np.load(R/'technical'/f'SRR{1520311+i}.{L}nt.CDS_counts.npz');rg=[];pg=[]
   for g in genes:
    v=data[g].astype(np.float64);v=v*len(v)/np.sum(v);rg.append(np.mean([v[a+x] for aid,a in anchors[g]],axis=0));pg.append(np.mean([np.mean([v[q+x] for q in pseudo[aid]],axis=0) for aid,a in anchors[g]],axis=0))
   real.append(np.asarray(rg));ps.append(np.asarray(pg))
  real=real[:6]+[(real[6]+real[7])/2];ps=ps[:6]+[(ps[6]+ps[7])/2]
  profiles[L]={'real':[np.mean(v,axis=0).tolist() for v in real],'pseudo':[np.mean(v,axis=0).tolist() for v in ps]}
 (R/'technical/frozen_visualization_values.json').write_text(json.dumps({'codon_positions':x.tolist(),'conditions':labels,'genes':280,'anchors':879,'profiles':profiles})+'\n')
 plotdata=np.column_stack([x]+[np.asarray(profiles[L][kind][row]) for L in [27,28] for row in range(7) for kind in ['real','pseudo']])
 np.savetxt(R/'technical/frozen_profiles.plotdata.dat',plotdata,fmt='%.17g',delimiter='\t')
 subprocess.run(['/usr/local/bin/Rscript',str(R/'draw_profiles.R')],cwd=R,check=True)
 lines += ['', 'Visualization is restricted to[-60,+100] with the frozen real/pseudo geometry and gene normalization. Baseline[-30,0] and downstream[30,60] bands are marked for orientation; no alternative windows, centers, smoothing or stratification.']
lines += ['', 'Technical accounting: reconstructed CDS totals exactly matched Task13 for every reconstructed run/length/gene. No new support exclusion or parameter change. Original Task12/13 input hashes verified before execution. Primary outputs were hash-frozen before28nt calculation; secondary, if permitted, was frozen before visualization. All original contracts, B2 and P3 materials remain unchanged.', '', 'Runtime: Python '+sys.version.split()[0]+'; NumPy '+np.__version__+'. Endpoint bootstrap seed5957314; ordinal seed5957315;10000whole-gene draws, PCG64, linear percentile intervals. Exact execution implementation, bootstrap draw arrays and frozen input identifiers are retained.', '', 'Stopped after authorized frozen execution. No post-result exploratory analysis.']
(R/'ROUTE_A_EXECUTION_VERDICT_v0.md').write_text('\n'.join(lines)+'\n')
print('FINAL_REPORT_AND_AUTHORIZED_VISUALIZATION_COMPLETE')
