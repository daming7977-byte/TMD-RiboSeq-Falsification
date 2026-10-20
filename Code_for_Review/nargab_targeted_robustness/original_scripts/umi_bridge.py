#!/usr/bin/env python3
"""Independent bounded point/cluster-bootstrap bridge; never processes BAM/FASTQ."""
from pathlib import Path
import hashlib,json,sys,platform,time,re
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
C=Path('historical_location/Documents/P3_FINAL_COVER_LETTER_MICRO_EDIT_20261018_v1/check/reviewer_parts/Code_for_Review')
F=C/'frozen_intermediates/region_composite_project'
H=Path('historical_location/paper3_GSE300977_raw')
U=H/'P3_ELIFE_UPGRADE'
REGIONS={'PLATEAU_90_120':[90,105,120],'BROAD_60_120':[60,75,90,105,120]}
SEEDS={'PLATEAU_90_120':30097890,'BROAD_60_120':30097860}
EXPECTED={'PLATEAU_90_120':0.9325467977159678,'BROAD_60_120':0.6518202467809296}
LIBS=[f'B{i:02d}' for i in range(1,13)]
COLS=[f'{l}_{s}' for l in LIBS for s in ['real','pseudo']]
CONTRASTS=[('ER',1,'B01','B02'),('ER',2,'B11','B12'),('mitochondria',1,'B09','B10'),('mitochondria',2,'B07','B08')]
INPUTS=[]
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(1048576),b''):h.update(x)
 return h.hexdigest()
def register(p):
 p=Path(p);INPUTS.append({'path':str(p),'bytes':p.stat().st_size,'sha256':digest(p)});return p
def read(p,**kw):return pd.read_csv(register(p),sep='\t',**kw)
def require(c,msg):
 if not c:raise RuntimeError(msg)
def make_pair_matrix(pairs,catalog,values):
 """values in catalog row order, one density per library. Pseudo means preserve cardinality."""
 lookup={(r.strand,r.geometry,r.mask_bits):i for i,r in enumerate(catalog.itertuples(index=False))}
 out=np.empty((len(pairs),24));realidx=[];pseudolist=[]
 for row in pairs.itertuples(index=False):
  r=lookup[(row.strand,row.real_window_geometry,row.real_crossmap_mask_bits)]
  gs=row.pseudo_window_geometries.split(';');ms=row.pseudo_crossmap_mask_bits.split(';')
  require(len(gs)==len(ms) and len(gs)>0,'pseudo mask cardinality')
  ps=[lookup[(row.strand,g,m)] for g,m in zip(gs,ms)]
  require(len(ps)==int(row.n_distinct_pseudo_geometries),'pseudo count differs')
  realidx.append(r);pseudolist.append(ps)
 for i,(r,ps) in enumerate(zip(realidx,pseudolist)):
  out[i,0::2]=values[r,:];out[i,1::2]=np.mean(values[ps,:],axis=0)
 require(np.isfinite(out).all() and (out>=0).all(),'invalid pair density')
 return out

def prepare(pairs,matrix,centers):
 """Return sorted-gene sums/counts for EXACT available-case retained rows."""
 selection=pairs.center.isin(centers).to_numpy();meta=pairs.loc[selection].copy();vals=matrix[selection]
 require(not meta.duplicated(['center','unit_key']).any(),'duplicate center-unit')
 genes=sorted(meta.gene_id.unique());gi={g:i for i,g in enumerate(genes)}
 groups={}
 for c in centers:
  for k in ['FIRST','INTERNAL']:
   idx=((meta.center==c)&(meta.position_class_consensus==k)).to_numpy()
   require(idx.any(),f'empty group {c} {k}')
   groupmeta=meta.loc[idx];v=vals[idx]
   order=np.argsort(groupmeta.unit_key.to_numpy());v=v[order];grows=groupmeta.iloc[order].gene_id
   sums=np.zeros((len(genes),24));counts=np.zeros(len(genes))
   for g,row in zip(grows,v):sums[gi[g]]+=row;counts[gi[g]]+=1
   groups[(c,k)]=(sums,counts)
 return {'genes':genes,'groups':groups,'centers':centers,'meta':meta}

def evaluate(prepared,weights):
 effects={};fi={};groups=prepared['groups'];centers=prepared['centers'];ci={x:i for i,x in enumerate(COLS)}
 for (c,k),(s,n) in groups.items():
  denom=weights@n;require((denom>0).all(),'empty weighted class')
  means=(weights@s)/denom[:,None]
  for loc,rep,pdlib,inlib in CONTRASTS:
   rpd=means[:,ci[f'{pdlib}_real']];rin=means[:,ci[f'{inlib}_real']];ppd=means[:,ci[f'{pdlib}_pseudo']];pin=means[:,ci[f'{inlib}_pseudo']]
   require(all((x>0).all() for x in (rpd,rin,ppd,pin)),'zero aggregate; no rescue')
   effects[(c,k,loc,rep)]=np.log2(rpd/rin)-np.log2(ppd/pin)
 for loc in ['ER','mitochondria']:
  for rep in [1,2]:fi[(loc,rep)]=np.mean(np.vstack([effects[(c,'FIRST',loc,rep)]-effects[(c,'INTERNAL',loc,rep)] for c in centers]),axis=0)
 r1=fi[('mitochondria',1)]-fi[('ER',1)];r2=fi[('mitochondria',2)]-fi[('ER',2)]
 return {'rep1':r1,'rep2':r2,'mean':(r1+r2)/2,'ER':(fi[('ER',1)]+fi[('ER',2)])/2,'mitochondria':(fi[('mitochondria',1)]+fi[('mitochondria',2)])/2,'effects':effects}

def bootstrap_many(prepared_list,seed,n=10000,batch=200):
 genes=prepared_list[0]['genes'];require(all(p['genes']==genes for p in prepared_list),'unmatched gene bootstrap populations')
 rng=np.random.default_rng(seed);prob=np.full(len(genes),1/len(genes));parts=[[] for _ in prepared_list]
 for start in range(0,n,batch):
  weights=rng.multinomial(len(genes),prob,size=min(batch,n-start))
  for i,p in enumerate(prepared_list):parts[i].append(evaluate(p,weights)['mean'])
 return [np.concatenate(x) for x in parts]

def main():
 t=time.monotonic()
 log=register(U/'U2B2_v1/U2B2_RUN_LOG_v1.txt').read_text();historical={Path(p).name:h for h,p in re.findall(r'INPUT_SHA256 ([0-9a-f]{64}) (.+)',log)}
 formalpath=F/'tmd_anchor_preregister/bridge/FORMAL_TMD_ANCHOR_CANDIDATES_v1p1.tsv'
 pairpath=F/'crossmap_window_preregister/output/GSE300977_PAIRED_REAL_PSEUDO_CROSSMAP_MASK_AUDIT_v1.tsv'
 catpath=F/'tmd_spatial_secondary_scan_v1/output/SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv'
 for p in [formalpath,pairpath,catpath]:require(digest(p)==historical[p.name],f'historical identity {p.name}')
 formal=read(formalpath,dtype=str);pairs=read(pairpath,dtype=str,keep_default_na=False);cat=read(catpath,dtype=str)
 require(not cat.geometry_id.duplicated().any(),'catalog duplicate IDs')
 require(not cat.duplicated(['strand','geometry','mask_bits']).any(),'catalog duplicate key')
 cat.unmasked_nt=cat.unmasked_nt.astype(int)
 require(all(m.count('0')==n for m,n in zip(cat.mask_bits,cat.unmasked_nt)),'mask denominator mismatch')
 n=len(cat);base=np.empty((n,12));collapsed=np.empty_like(base);compat=[]
 for j,lib in enumerate(LIBS):
  bp=F/f'tmd_spatial_secondary_scan_v1/geometry_counts/{lib}.SECONDARY_GEOMETRY_DENSITIES_v1.tsv'
  cp=U/f'U2B2A_v1/geometry_counts/{lib}.U2B2_EXACT5UMI_COLLAPSED_GEOMETRY_DENSITIES_v1.tsv'
  require(digest(cp)==historical[cp.name],f'collapse historical hash {lib}')
  b=read(bp);d=read(cp)
  for frame,label in [(b,'baseline'),(d,'collapsed')]:
   require(not frame.geometry_id.duplicated().any(),f'duplicate {lib} {label}')
   require(set(frame.geometry_id)==set(cat.geometry_id),f'geometry key coverage {lib} {label}')
  b=b.set_index('geometry_id').loc[cat.geometry_id];d=d.set_index('geometry_id').loc[cat.geometry_id]
  for frame,label in [(b,'baseline'),(d,'collapsed')]:
   require(np.array_equal(frame.unmasked_nt.to_numpy(),cat.unmasked_nt.to_numpy()),f'denominator {lib} {label}')
   require(np.isfinite(frame.density).all() and (frame.density>=0).all(),f'density {lib} {label}')
  base[:,j]=b.density;collapsed[:,j]=d.density
  require((collapsed[:,j]<=base[:,j]+1e-8).all(),f'collapsed density increased {lib}')
  compat.append({'library':lib,'catalog_geometries':n,'uncollapsed_rows':len(b),'collapsed_rows':len(d),'missing_extra_duplicate_keys':0,'denominator_mismatches':0,'historical_collapsed_hash':'PASS','uncollapsed_zero_geometries':int((base[:,j]==0).sum()),'collapsed_zero_geometries':int((collapsed[:,j]==0).sum()),'positive_to_zero_geometries':int(((base[:,j]>0)&(collapsed[:,j]==0)).sum()),'density_increase_geometries':int((collapsed[:,j]>base[:,j]+1e-8).sum()),'uncollapsed_density_sum':float(base[:,j].sum()),'collapsed_density_sum':float(collapsed[:,j].sum())})
 pd.DataFrame(compat).to_csv(HERE/'INPUT_COMPATIBILITY.tsv',sep='\t',index=False)
 np.savez_compressed(HERE/'GEOMETRY_DENSITIES.npz',geometry_ids=cat.geometry_id.to_numpy(dtype=str),unmasked_nt=cat.unmasked_nt.to_numpy(),uncollapsed=base,collapsed=collapsed,libraries=np.array(LIBS))
 catalogout=cat.copy();catalogout.to_csv(HERE/'GEOMETRY_CATALOG.tsv',sep='\t',index=False)
 pairs.center=pairs.center.astype(int)
 for col in ['real_crossmap_unmasked_nt','pseudo_min_unmasked_nt','n_gene_ids']:pairs[col]=pairs[col].astype(int)
 tx= formal[['target_transcript','pass_class']].dropna().drop_duplicates();require((tx.groupby('target_transcript').pass_class.nunique()<=1).all(),'pass class conflict')
 txmap=dict(zip(tx.target_transcript,tx.pass_class))
 pairs['topology_class']=[('MULTI_PASS' if s=={'MULTI_PASS'} else 'SINGLE_PASS' if s=={'SINGLE_PASS'} else 'MIXED_TOPOLOGY') for s in [{txmap[t] for t in ts.split(',')} for ts in pairs.target_transcripts]]
 pairs=pairs[(pairs.real_crossmap_unmasked_nt>=1)&(pairs.pseudo_min_unmasked_nt>=1)&(pairs.topology_class=='MULTI_PASS')&pairs.position_class_consensus.isin(['FIRST','INTERNAL'])&(pairs.unique_gene_assignment=='YES')&(pairs.n_gene_ids==1)].copy()
 pairs['gene_id']=pairs.gene_ids
 keycols=['chrom','strand','real_anchor_genomic_start0','target_transcripts','gene_ids','position_class_consensus','tmd_indices']
 pairs['unit_key']=pairs[keycols].astype(str).agg('|'.join,axis=1)
 require(not pairs.duplicated(['center','unit_key']).any(),'center key duplicate')
 pairs=pairs.reset_index(drop=True);pairs.insert(0,'matrix_row',np.arange(len(pairs)))
 matrix=make_pair_matrix(pairs,cat,base);cmatrix=make_pair_matrix(pairs,cat,collapsed)
 pairs.to_csv(HERE/'PAIR_METADATA.tsv',sep='\t',index=False)
 np.savez_compressed(HERE/'PAIR_DENSITY_MATRICES.npz',uncollapsed=matrix,collapsed=cmatrix,columns=np.array(COLS))
 print('MATRICES_READY',len(pairs),flush=True)
 # Gate is completed for BOTH baseline regions before any collapsed effect calculation.
 preps={};gate=[]
 frozen=read(C/'reproduction_output/available_case/REGION_COMPOSITE_SUMMARY_v1.tsv')
 for region,centers in REGIONS.items():
  p=prepare(pairs,matrix,centers);preps[region]=p;obs=evaluate(p,np.ones((1,len(p['genes']))))
  value=float(obs['mean'][0]);error=abs(value-EXPECTED[region]);require(error<=1e-10,f'baseline failed {region} {error}')
  gate.append({'region':region,'reconstructed':value,'frozen':EXPECTED[region],'absolute_error':error,'gate':'PASS'})
 pd.DataFrame(gate).to_csv(HERE/'BASELINE_POINT_GATES.tsv',sep='\t',index=False)
 print('BOTH_BASELINE_POINT_GATES_PASS',gate,flush=True)
 results=[];centerrows=[];support=[]
 for region,centers in REGIONS.items():
  p=preps[region];q=prepare(pairs,cmatrix,centers);require(p['meta'].paired_real_unit_id.tolist()==q['meta'].paired_real_unit_id.tolist(),'matched population mismatch')
  boot=bootstrap_many([p,q],SEEDS[region]);np.savez_compressed(HERE/f'{region}_BOOTSTRAP.npz',uncollapsed=boot[0],collapsed=boot[1],difference=boot[1]-boot[0])
  for label,prepared,dist in [('uncollapsed',p,boot[0]),('operationally_collapsed',q,boot[1])]:
   obs=evaluate(prepared,np.ones((1,len(prepared['genes']))));ci=np.percentile(dist,[2.5,97.5]);meta=prepared['meta']
   results.append({'analysis':'POST-REVEAL EXPLORATORY SENSITIVITY','region':region,'density_track':label,'population':'frozen AVAILABLE_CASE unique-gene MULTI_PASS FIRST/INTERNAL','N_genes':len(prepared['genes']),'N_units_union':meta.unit_key.nunique(),'N_center_unit_rows':len(meta),'replicate_1':float(obs['rep1'][0]),'replicate_2':float(obs['rep2'][0]),'interaction':float(obs['mean'][0]),'conditional_gene_bootstrap_95_lower':float(ci[0]),'conditional_gene_bootstrap_95_upper':float(ci[1]),'CI_includes_zero':bool(ci[0]<=0<=ci[1]),'bootstrap_seed':SEEDS[region],'bootstrap_draws':10000,'bootstrap_batch':200})
   for (center,cls,loc,rep),vals in obs['effects'].items():centerrows.append({'region':region,'density_track':label,'center':center,'class':cls,'compartment':loc,'replicate':rep,'D':float(vals[0])})
  fr=frozen.loc[frozen.region==region].iloc[0];uci=np.percentile(boot[0],[2.5,97.5]);require(np.all(np.abs(uci-[fr['interaction_CI_2.5'],fr['interaction_CI_97.5']])<=1e-10),'baseline interval mismatch')
  for center in centers:
   for cls in ['FIRST','INTERNAL']:
    ids=((pairs.center==center)&(pairs.position_class_consensus==cls)).to_numpy()
    for j,col in enumerate(COLS):
     bv=matrix[ids,j];cv=cmatrix[ids,j]
     support.append({'region':region,'center':center,'class':cls,'library_status':col,'eligible_pairs':int(ids.sum()),'uncollapsed_positive_pairs':int((bv>0).sum()),'collapsed_positive_pairs':int((cv>0).sum()),'positive_to_zero_pairs':int(((bv>0)&(cv==0)).sum()),'uncollapsed_sum_density':float(bv.sum()),'collapsed_sum_density':float(cv.sum()),'ratio_sum_density':float(cv.sum()/bv.sum())})
  print('REGION_COMPLETE',region,results[-2:],flush=True)
 pd.DataFrame(results).to_csv(HERE/'P3_UMI_MATCHED_POPULATION_RESULTS.tsv',sep='\t',index=False,float_format='%.17g')
 pd.DataFrame(centerrows).to_csv(HERE/'CENTER_SPECIFIC_ESTIMATES.tsv',sep='\t',index=False,float_format='%.17g')
 pd.DataFrame(support).to_csv(HERE/'MATCHED_POPULATION_SUPPORT.tsv',sep='\t',index=False,float_format='%.17g')
 pd.DataFrame(INPUTS).drop_duplicates('path').to_csv(HERE/'INPUT_SHA256.tsv',sep='\t',index=False)
 record={'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'platform':platform.platform(),'elapsed_seconds':time.monotonic()-t,'script_sha256':digest(__file__),'inputs':len({x['path'] for x in INPUTS}),'bootstrap_n':10000,'batch':200,'seeds':SEEDS,'threads':'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1','baseline_point_tolerance':1e-10,'baseline_interval_tolerance':1e-10,'point_gates':'PASS','baseline_CI_gates':'PASS','no_BAM_processing':True}
 (HERE/'RUN_RECORD.json').write_text(json.dumps(record,indent=2)+'\n')
 print(json.dumps(record,indent=2),flush=True)
if __name__=='__main__':main()
