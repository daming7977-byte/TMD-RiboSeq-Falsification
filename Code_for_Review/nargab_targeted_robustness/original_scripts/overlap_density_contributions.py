from pathlib import Path
import numpy as np,pandas as pd,hashlib,json
W=Path('historical_location/Documents/P3_FINAL_TARGETED_ROBUSTNESS_VERIFICATION_20261018_v1');ip=W/'A/P3_OVERLAP_INVENTORY.tsv';gp=W/'B/GEOMETRY_DENSITIES.npz'
p=pd.read_csv(ip,sep='\t',keep_default_na=False,low_memory=False);geo=np.load(gp);lookup={g:i for i,g in enumerate(geo['geometry_ids'])};ids=np.array([lookup[g] for g in p.geometry_id]);density=geo['uncollapsed'][ids]*p.pseudo_aggregation_weight.to_numpy()[:,None]
# Density contributions are whole-window, pseudo-component-weighted contributions, NOT counts attributable to overlap bases.
p['flag']=p.unmasked_other_gene_overlap_nt>0
libmap={1:('ER',1,'pulldown'),2:('ER',1,'input'),11:('ER',2,'pulldown'),12:('ER',2,'input'),9:('mitochondria',1,'pulldown'),10:('mitochondria',1,'input'),7:('mitochondria',2,'pulldown'),8:('mitochondria',2,'input'),5:('cytosol',1,'pulldown'),6:('cytosol',1,'input'),3:('cytosol',2,'pulldown'),4:('cytosol',2,'input')}
rows=[];gene_rows=[]
for key,idx in p.groupby(['analysis','center','topology_class','position_class','side']).groups.items():
 ii=np.array(idx);sub=p.loc[ii];flag=sub.flag.to_numpy();wt=sub.pseudo_aggregation_weight.to_numpy();wtshare=wt[flag].sum()/wt.sum()
 for libnum,(comp,rep,role) in libmap.items():
  if key[0]!='PRIMARY_PLUS75' and comp=='cytosol':continue
  vals=density[ii,libnum-1];total=vals.sum();marked=vals[flag].sum()
  rows.append(dict(zip(['analysis','center','topology_class','position_class','side'],key))|{'compartment':comp,'replicate_label':rep,'role':role,'library':f'B{libnum:02d}','pair_equivalent_windows':wt.sum(),'affected_pair_equivalent_windows':wt[flag].sum(),'affected_window_weight_fraction':wtshare,'total_weighted_window_density':total,'affected_whole_window_density':marked,'affected_whole_window_density_fraction':marked/total if total>0 else np.nan,'density_share_over_window_share':marked/total/wtshare if total>0 and wtshare>0 else np.nan})
  for gene,gidx in sub.loc[flag].groupby('owner_gene_ids').groups.items():
   locs=p.index.get_indexer(gidx);value=density[locs,libnum-1].sum()
   gene_rows.append(dict(zip(['analysis','center','topology_class','position_class','side'],key))|{'owner_gene_ids':gene,'owner_gene_names':p.loc[gidx[0],'owner_gene_names'],'compartment':comp,'replicate_label':rep,'role':role,'library':f'B{libnum:02d}','affected_whole_window_density':value,'fraction_of_full_factor':value/total if total>0 else np.nan})
r=pd.DataFrame(rows);r.to_csv(W/'A/DENSITY_CONTRIBUTIONS_BY_CLASS_CENTER_LIBRARY_SIDE.tsv',sep='\t',index=False)
g=pd.DataFrame(gene_rows);g.to_csv(W/'A/FLAGGED_GENE_DENSITY_CONTRIBUTIONS.tsv',sep='\t',index=False)
# Primary actual estimand pools all topology classes; this separate summary uses that exact pooling.
pr=p.loc[p.analysis.eq('PRIMARY_PLUS75')];actual=[]
for side,idx in pr.groupby('side').groups.items():
 ii=np.array(idx);f=p.loc[ii,'flag'].to_numpy();wt=p.loc[ii,'pseudo_aggregation_weight'].to_numpy()
 for j,(comp,rep,role) in libmap.items():
  v=density[ii,j-1];actual.append({'analysis':'PRIMARY_PLUS75','side':side,'compartment':comp,'replicate_label':rep,'role':role,'library':f'B{j:02d}','affected_window_weight_fraction':wt[f].sum()/wt.sum(),'affected_whole_window_density_fraction':v[f].sum()/v.sum(),'total_weighted_window_density':v.sum(),'affected_whole_window_density':v[f].sum()})
pd.DataFrame(actual).to_csv(W/'A/PRIMARY_POOLED_DENSITY_CONTRIBUTIONS.tsv',sep='\t',index=False)
# Strongest contributing overlap genes for each analysis, indexed existing factors only.
g.sort_values('fraction_of_full_factor',ascending=False).groupby(['analysis','owner_gene_ids'],sort=False).head(1).sort_values(['analysis','fraction_of_full_factor'],ascending=[True,False]).to_csv(W/'A/OVERLAP_GENE_MAX_FACTOR_SHARE.tsv',sep='\t',index=False)
pd.DataFrame([{'path':str(z),'sha256':hashlib.file_digest(z.open('rb'),'sha256').hexdigest()} for z in [ip,gp]]).to_csv(W/'A/density_contribution_inputs.tsv',sep='\t',index=False)
print('Windowdensityshares; not overlap-base specific counts. No new recount.')
print(r.groupby('analysis').affected_whole_window_density_fraction.agg(['min','max']).to_string());print(pd.DataFrame(actual).to_string(index=False))
