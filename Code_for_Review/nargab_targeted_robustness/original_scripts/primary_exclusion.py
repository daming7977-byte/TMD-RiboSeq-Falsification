import pathlib,hashlib,sys,json,time
import numpy as np,pandas as pd
W=pathlib.Path('historical_location/Documents/P3_FINAL_TARGETED_ROBUSTNESS_VERIFICATION_20261018_v1'); H=pathlib.Path('historical_location/paper3_GSE300977_raw')
inputs=[H/'tmd_spatial_primary_reveal_v1/output/PRIMARY_PLUS75_UNIT_DENSITIES_v1.tsv',H/'tmd_spatial_primary_reveal_v1/output/PRIMARY_PLUS75_LOCATION_SUMMARY_v1.tsv',H/'tmd_spatial_analysis_preregister/GSE300977_PRIMARY_SPATIAL_ESTIMAND_CONTRACT_v1.txt']
mode=sys.argv[1]; start=time.time(); p=pd.read_csv(inputs[0],sep='\t',float_precision='round_trip');assert len(p)==9098 and not p.paired_real_unit_id.duplicated().any()
cols=[f'B{i:02d}_{k}_density' for i in range(1,13) for k in ('real','pseudo')]; x=p[cols].to_numpy();assert np.isfinite(x).all() and (x>=0).all()
contrasts={'ER':[(1,2),(11,12)],'mitochondria':[(9,10),(7,8)],'cytosol':[(5,6),(3,4)]}
def evaluate(means):
 assert (means>0).all()
 return np.stack([np.stack([np.log2(means[:,2*(a-1)]/means[:,2*(b-1)])-np.log2(means[:,2*(a-1)+1]/means[:,2*(b-1)+1]) for a,b in reps],axis=1) for reps in contrasts.values()],axis=1)
saved=pd.read_csv(inputs[1],sep='\t',float_precision='round_trip');base=evaluate(x.mean(axis=0)[None,:]);assert np.max(np.abs(base.mean(axis=2)[0]-saved.D_location_mean.to_numpy()))<1e-10
if mode=='exclusion':
 flags=W/'A/pair_overlap_flags.tsv';inputs.append(flags); f=pd.read_csv(flags,sep='\t',dtype=str).set_index('paired_real_unit_id'); assert f.index.is_unique
 s=f.loc[p.paired_real_unit_id,'exclude_pair_unmasked_other_gene'].eq('YES').to_numpy();x=x[~s];p=p.loc[~s]
else:s=np.zeros(len(x),bool)
n=len(x); observed=evaluate(x.mean(axis=0)[None,:]);rng=np.random.default_rng(30097775);boot=[]
for i in range(0,10000,200):
 weights=rng.multinomial(n,np.full(n,1/n),size=200)
 boot.append(evaluate((weights@x)/n).mean(axis=2))
 if i%2000==0:print(f'{mode} {i+200}/10000',flush=True)
b=np.concatenate(boot); ci=np.percentile(b,[2.5,97.5],axis=0)
rows=[]
for j,loc in enumerate(contrasts):
 rows.append({'analysis':'PRIMARY_PLUS75','scenario':mode,'location':loc,'n_original_units':9098,'n_excluded_units':int(s.sum()),'n_retained_units':n,'rep1':observed[0,j,0],'rep2':observed[0,j,1],'estimate':observed[0,j].mean(),'ci_low':ci[0,j],'ci_high':ci[1,j],'interval_type':'conditional unit-bootstrap percentile','bootstrap_n':10000,'seed':30097775,'strong_support_rule_satisfied':bool(np.sign(observed[0,j,0])==np.sign(observed[0,j,1]) and (ci[0,j]>0 or ci[1,j]<0))})
result=pd.DataFrame(rows);result.to_csv(W/'A'/f'primary_{mode}_results.tsv',sep='\t',index=False)
if mode=='baseline':assert np.max(np.abs(ci.T-saved[['bootstrap_ci_2.5','bootstrap_ci_97.5']].to_numpy()))<1e-10
pd.DataFrame({'paired_real_unit_id':p.paired_real_unit_id}).to_csv(W/'A'/f'primary_{mode}_retained_ids.tsv',sep='\t',index=False)
pd.DataFrame([{'path':str(z),'bytes':z.stat().st_size,'sha256':hashlib.file_digest(z.open('rb'),'sha256').hexdigest()} for z in inputs]).to_csv(W/'A'/f'primary_{mode}_inputs.tsv',sep='\t',index=False)
print(result.to_string(index=False));print(json.dumps({'elapsed_seconds':time.time()-start,'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'baseline_gate':'PASS'}))
