from pathlib import Path
import importlib.util,sys,json,time,hashlib
import numpy as np,pandas as pd
W=Path('historical_location/Documents/P3_FINAL_TARGETED_ROBUSTNESS_VERIFICATION_20261018_v1')
spec=importlib.util.spec_from_file_location('bridge',W/'B/umi_bridge.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
p=pd.read_csv(W/'B/PAIR_METADATA.tsv',sep='\t',keep_default_na=False); a=np.load(W/'B/PAIR_DENSITY_MATRICES.npz');x=a['uncollapsed']; flags=pd.read_csv(W/'A/pair_overlap_flags.tsv',sep='\t',keep_default_na=False).set_index('paired_real_unit_id');assert flags.index.is_unique
select=p.center.isin([60,75,90,105,120]).to_numpy();p=p.loc[select].reset_index(drop=True);x=x[select]
f=flags.loc[p.paired_real_unit_id];assert np.array_equal(f.center.to_numpy(),p.center.to_numpy());exclude=f.exclude_pair_unmasked_other_gene.eq('YES').to_numpy()
rows=[];support=[];t=time.monotonic()
for region,centers in b.REGIONS.items():
 base=b.prepare(p,x,centers);orig=b.evaluate(base,np.ones((1,len(base['genes']))));assert abs(orig['mean'][0]-b.EXPECTED[region])<1e-10
 q=p.loc[~exclude].copy();z=x[~exclude];pr=b.prepare(q,z,centers)
 obs=b.evaluate(pr,np.ones((1,len(pr['genes']))));boot=b.bootstrap_many([pr],b.SEEDS[region])[0];ci=np.percentile(boot,[2.5,97.5])
 np.savez_compressed(W/'A'/f'{region}_EXCLUSION_BOOTSTRAP.npz',interaction=boot)
 rows.append({'analysis':'POST-REVEAL EXPLORATORY SENSITIVITY','region':region,'scenario':'exclude pairs with any unmasked same-strand other-gene exon overlap','original_genes':len(base['genes']),'retained_genes':len(pr['genes']),'original_union_units':base['meta'].unit_key.nunique(),'retained_union_units':pr['meta'].unit_key.nunique(),'original_center_unit_rows':len(base['meta']),'retained_center_unit_rows':len(pr['meta']),'frozen_interaction':orig['mean'][0],'rep1':obs['rep1'][0],'rep2':obs['rep2'][0],'interaction':obs['mean'][0],'ci_low':ci[0],'ci_high':ci[1],'interval_type':'conditional gene-bootstrap percentile','bootstrap_draws':10000,'bootstrap_seed':b.SEEDS[region],'interval_includes_zero':bool(ci[0]<=0<=ci[1])})
 for c in centers:
  for k in ['FIRST','INTERNAL']:
   idx=((p.center==c)&p.position_class_consensus.eq(k)).to_numpy();part=p.loc[idx];affected=part.loc[exclude[idx]]
   support.append({'region':region,'center':c,'class':k,'original_units':int(idx.sum()),'removed_units':int((exclude&idx).sum()),'remaining_units':int((~exclude&idx).sum()),'original_genes':part.gene_id.nunique(),'affected_genes':affected.gene_id.nunique(),'retained_genes':part.loc[~exclude[idx]].gene_id.nunique()})
 print(region,rows[-1],flush=True)
pd.DataFrame(rows).to_csv(W/'A/REGION_EXCLUSION_RESULTS.tsv',sep='\t',index=False,float_format='%.17g');pd.DataFrame(support).to_csv(W/'A/REGION_EXCLUSION_SUPPORT.tsv',sep='\t',index=False)
# Frozen per-gene deletion results are only indexed, not recomputed.
fp=Path('historical_location/Documents/P3_FINAL_COVER_LETTER_MICRO_EDIT_20261018_v1/check/reviewer_parts/Code_for_Review/reproduction_output/available_case/REGION_COMPOSITE_GENE_INFLUENCE_v1.tsv')
loo=pd.read_csv(fp,sep='\t');found=[]
for region,centers in b.REGIONS.items():
 gs=set(p.loc[exclude&p.center.isin(centers).to_numpy(),'gene_id']);sub=loo.loc[loo.region.eq(region)&loo.gene_id.isin(gs)].copy();sub['evidence']='existing frozen whole-gene deletion, not the new pair exclusion';found.append(sub.sort_values('absolute_change',ascending=False))
pd.concat(found).to_csv(W/'A/OVERLAP_GENES_EXISTING_INFLUENCE.tsv',sep='\t',index=False)
inputs=[W/'B/PAIR_METADATA.tsv',W/'B/PAIR_DENSITY_MATRICES.npz',W/'B/umi_bridge.py',W/'A/pair_overlap_flags.tsv',fp]
pd.DataFrame([{'path':str(z),'bytes':z.stat().st_size,'sha256':hashlib.file_digest(z.open('rb'),'sha256').hexdigest()} for z in inputs]).to_csv(W/'A/region_exclusion_inputs.tsv',sep='\t',index=False)
print(json.dumps({'elapsed_seconds':time.monotonic()-t,'numpy':np.__version__,'pandas':pd.__version__,'python':sys.version,'threads':1,'same_estimand_changed_exclusion_population':True}))
