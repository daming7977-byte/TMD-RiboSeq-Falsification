#!/usr/bin/env python3
"""Independent prespecified gene/composition sensitivity. No project imports."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,math
import numpy as np
import pandas as pd

import os
P=Path(__file__).resolve().parents[1]
W=Path(os.environ['P3_SENSITIVITY_OUTPUT']).resolve(); E=W/'evidence'; R=P.parent
freeze=json.loads((E/'PROTOCOL_FREEZE.json').read_text())
assert hashlib.sha256((W/'P3_WITHIN_GENE_PRESPECIFIED_PROTOCOL.md').read_bytes()).hexdigest()==freeze['sha256']
u=pd.read_csv(E/'verified_units.tsv',sep='\t',dtype={'chrom':str,'gene_ids':str})
X=np.load(E/'verified_unit_densities.npy')
ann=pd.read_csv(E/'gene_annotation.tsv',sep='\t').set_index('gene_id')
manifest=pd.read_csv(E/'sample_manifest.tsv',sep='\t',dtype=str)
coefficient=np.zeros((24,4))
for j,(loc,rep) in enumerate([('ER','1'),('ER','2'),('mitochondria','1'),('mitochondria','2')]):
    for fraction,sgn in [('pulldown',1),('input',-1)]:
        bid=manifest[(manifest.location==loc)&(manifest.replicate==rep)&(manifest.fraction==fraction)].iloc[0].blind_id
        ix=(int(bid[1:])-1)*2
        coefficient[ix,j]=sgn;coefficient[ix+1,j]=-sgn
active=np.where((coefficient!=0).any(axis=1))[0];C=coefficient[active]
inputcols=np.array([i for b in [2,10,12,8] for i in [(b-1)*2,(b-1)*2+1]])
metrics=['ER_FI_rep1','ER_FI_rep2','MITO_FI_rep1','MITO_FI_rep2','interaction_rep1','interaction_rep2','interaction_mean']
names=['FROZEN_UNIT_ALL','EQUAL_GENE_ALL','UNIT_SHARED','EQUAL_GENE_SHARED','UNIT_POSITIVE_SHARED','EQUAL_GENE_POSITIVE_SHARED','WITHIN_GENE_MEAN','WITHIN_GENE_MEDIAN','WITHIN_GENE_COMMON_TRANSCRIPT','BALANCE_CORE_RESTRICTED_UNIT','BALANCE_CORE_RESTRICTED_EQUAL_GENE','BALANCE_CORE_WEIGHTED','BALANCE_EXTENDED_RESTRICTED_UNIT','BALANCE_EXTENDED_RESTRICTED_EQUAL_GENE','BALANCE_EXTENDED_WEIGHTED']
results=[]; eligibility=[]; influence=[]; diagnostics=[]; balance_records=[]; support_rows=[]; comparisons=[]; strata_results=[]; unit_contrib=[]

def effects(t):
    """t (..., center, class, 24), result (...,7). Class normalizers cancel."""
    a=t[...,active]
    valid=np.all(a>0,axis=(-3,-2,-1))
    with np.errstate(divide='ignore',invalid='ignore'):
        d=np.log2(a)@C
        fi=(d[...,0,:]-d[...,1,:]).mean(axis=-2)
        it=fi[...,2:4]-fi[...,:2]
        out=np.concatenate([fi,it,it.mean(axis=-1,keepdims=True)],axis=-1)
    return np.where(np.asarray(valid)[...,None],out,np.nan)

def gene_effects(m):
    return effects(m) # g, center, class, 24 -> g,7

def binom_quantile(n,p,a):
    # Stable recurrence for p=.025; symmetry supplies .975. Exact finite-binomial CDF.
    if p>.5: return n-binom_quantile(n,1-p,1-a)
    prob=math.exp(n*math.log1p(-p)); cdf=prob
    for k in range(n+1):
        if cdf>=a: return k
        prob*=((n-k)/(k+1))*p/(1-p);cdf+=prob
    return n

def quantile_mc_bounds(draws,p):
    z=np.sort(draws);n=len(z)
    lo=max(0,binom_quantile(n,p,.005)-1);hi=min(n-1,binom_quantile(n,p,.995))
    return float(z[lo]),float(z[hi])

def weighted_medians(v, multiplicities):
    n=int(multiplicities[0].sum());out=[]
    for j in range(v.shape[1]):
        order=np.argsort(v[:,j],kind='stable');cs=np.cumsum(multiplicities[:,order],axis=1)
        il=np.argmax(cs>((n-1)//2),axis=1);ih=np.argmax(cs>(n//2),axis=1)
        zz=v[order,j];out.append((zz[il]+zz[ih])/2)
    return np.stack(out,axis=1)

def summary_cov(region,design,center,variable,values,weights):
    vals=[]
    for cl in range(2):
        v=values[:,cl];w=weights[:,cl];ok=(w>0)&np.isfinite(v);v=v[ok];w=w[ok]
        if not len(v): return
        w=w/w.sum();mu=float(w@v);variance=float(w@((v-mu)**2))
        vals.append((mu,variance,len(v),float(1/(w@w)),float(w.max())))
    den=math.sqrt((vals[0][1]+vals[1][1])/2);delta=vals[0][0]-vals[1][0]
    smd=delta/den if den>0 else (0.0 if delta==0 else math.copysign(float('inf'),delta))
    diagnostics.append(dict(region=region,design=design,center=center,variable=variable,FIRST_mean=vals[0][0],INTERNAL_mean=vals[1][0],SMD=smd,abs_SMD_gt_0_1=abs(smd)>.1,FIRST_genes=vals[0][2],INTERNAL_genes=vals[1][2],FIRST_Kish_ESS=vals[0][3],INTERNAL_Kish_ESS=vals[1][3],FIRST_max_weight=vals[0][4],INTERNAL_max_weight=vals[1][4]))

def evaluate(region,rindex,name,genes,T=None,V=None,method='linear',supported=True,reason='',frozen_point=None):
    G=len(genes);vid=names.index(name);seed=30097890 if (vid==0 and rindex==1) else 30097860 if vid==0 else 202610160000+rindex*100+vid
    meta=dict(region=region,estimand=name,n_genes=G,seed=seed,bootstrap_draws=10000,support_status='SUPPORTED' if supported and G>=30 else 'SPARSE_OR_UNSUPPORTED',support_reason=reason,protocol_sha256=freeze['sha256'])
    if G==0:
        meta.update(status='NOT_ESTIMABLE_EMPTY_POPULATION');results.append(meta);print(region,name,'EMPTY',flush=True);return
    if method=='linear':
        point=effects(T.sum(axis=0));loo=effects(T.sum(axis=0)[None]-T)
    elif method=='mean':
        point=V.mean(axis=0);loo=(V.sum(axis=0)[None]-V)/(G-1) if G>1 else np.full_like(V,np.nan)
    else:
        point=np.median(V,axis=0);loo=np.stack([np.median(np.delete(V,g,axis=0),axis=0) for g in range(G)]) if G>1 else np.full_like(V,np.nan)
    if not np.isfinite(point).all():
        meta.update(status='NOT_ESTIMABLE_NONPOSITIVE_AGGREGATE');results.append(meta);return
    rng=np.random.default_rng(seed);draws=[]
    for start in range(0,10000,200):
        multiplicity=rng.multinomial(G,np.full(G,1/G),size=min(200,10000-start))
        if method=='linear': d=effects((multiplicity@T.reshape(G,-1)).reshape((-1,)+T.shape[1:]))
        elif method=='mean': d=(multiplicity@V)/G
        else:d=weighted_medians(V,multiplicity)
        draws.append(d)
    draws=np.concatenate(draws);valid=np.isfinite(draws).all(axis=1)
    np.save(E/f'bootstrap_{region}_{name}.npy',draws)
    meta.update(status='COMPLETE' if valid.all() else 'INTERVAL_UNAVAILABLE_INVALID_DRAWS',invalid_draws=int((~valid).sum()),interval_scope='conditional_gene_resampling_fixed_population_and_design')
    for j,key in enumerate(metrics):meta[key]=float(point[j])
    if valid.all():
        q=np.quantile(draws,[.025,.975],axis=0,method='linear')
        for j,key in enumerate(metrics):
            meta[key+'_percentile_lo']=float(q[0,j]);meta[key+'_percentile_hi']=float(q[1,j])
            meta[key+'_basic_lo']=float(2*point[j]-q[1,j]);meta[key+'_basic_hi']=float(2*point[j]-q[0,j])
        mc1=quantile_mc_bounds(draws[:,-1],.025);mc2=quantile_mc_bounds(draws[:,-1],.975)
        meta.update(percentile_lo_MC99_lo=mc1[0],percentile_lo_MC99_hi=mc1[1],percentile_hi_MC99_lo=mc2[0],percentile_hi_MC99_hi=mc2[1],percentile_crosses_zero=bool(q[0,-1]<=0<=q[1,-1]),basic_crosses_zero=bool(2*point[-1]-q[1,-1]<=0<=2*point[-1]-q[0,-1]))
    delta=loo[:,-1]-point[-1]
    for g,lv,ch in zip(genes,loo[:,-1],delta):
        influence.append(dict(region=region,estimand=name,gene_id=g,full=point[-1],leave_one_out=lv,signed_change=ch,absolute_change=abs(ch),relative_change=abs(ch)/abs(point[-1]) if point[-1]!=0 else np.nan,sign_reversal=bool(np.sign(lv)!=np.sign(point[-1])) if np.isfinite(lv) else None))
    finite=np.isfinite(delta)
    if finite.any():
        best=sorted(np.where(finite)[0],key=lambda i:(-round(abs(float(delta[i])),12),genes[i]))[0]
        meta.update(max_LOO_abs_change=float(abs(delta[best])),max_LOO_gene=genes[best],max_LOO_relative_change=float(abs(delta[best])/abs(point[-1])) if point[-1]!=0 else np.nan,any_LOO_sign_reversal=bool((np.sign(loo[finite,-1])!=np.sign(point[-1])).any()),invalid_LOO=int((~finite).sum()))
    meta.update(positive_direction=bool(point[-1]>0),replicates_both_positive=bool((point[4:6]>0).all()))
    if frozen_point is not None:
        meta.update(change_from_frozen=float(point[-1]-frozen_point),relative_magnitude_change=float((point[-1]-frozen_point)/abs(frozen_point)),material_magnitude_flag=bool(abs(point[-1]-frozen_point)>=.10 or abs(point[-1]-frozen_point)/abs(frozen_point)>=.25))
    results.append(meta)
    print(region,name,'N',G,'POINT',round(point[-1],9),'CI',[round(meta.get('interaction_mean_percentile_'+a,float('nan')),9) for a in ['lo','hi']],'support',meta['support_status'],flush=True)

def analyze(region,cs,rindex):
    z=u[u.center.isin(cs)].copy();genes=sorted(z.gene_ids.unique());G=len(genes);gi={g:i for i,g in enumerate(genes)}
    S=np.zeros((G,len(cs),2,24));N=np.zeros((G,len(cs),2),dtype=int)
    pos=np.zeros((G,len(cs),2));mask=pos.copy();ordinal=pos.copy();absolute=pos.copy()
    txsets=[[] for _ in genes]
    for ix,r in z.sort_values('unit_key').iterrows():
        g=gi[r.gene_ids];c=cs.index(r.center);k=0 if r.position_class_consensus=='FIRST' else 1
        S[g,c,k]+=X[ix];N[g,c,k]+=1;pos[g,c,k]+=r.relative_cds_position;mask[g,c,k]+=r.mask_fraction;ordinal[g,c,k]+=r.ordinal_mean;absolute[g,c,k]+=r.absolute_cds_start
        txsets[g].append(set(r.target_transcripts.split(',')))
    M=np.divide(S,N[...,None],out=np.zeros_like(S),where=N[...,None]>0)
    for arr in [pos,mask,ordinal,absolute]:np.divide(arr,N,out=arr,where=N>0)
    shared=(N>0).all(axis=(1,2));positive=shared&(M[...,active]>0).all(axis=(1,2,3))
    coherent=np.array([bool(set.intersection(*ts)) for ts in txsets])
    J=gene_effects(M[positive]);positivegenes=np.array(genes)[positive]
    np.savez_compressed(E/f'gene_tensor_{region}.npz',genes=np.array(genes),centers=cs,sums=S,counts=N,means=M,shared=shared,positive=positive,coherent=coherent,within_gene_effects=J)
    gann=ann.loc[genes];ct=gann.annotated_tmd_count.to_numpy();length=gann.cds_length.to_numpy()
    coverage=M[...,inputcols].mean(axis=-1)
    for g,gid in enumerate(genes):
        unitset=z[z.gene_ids==gid].drop_duplicates('unit_key')
        eligibility.append(dict(region=region,gene_id=gid,gene_name=gann.iloc[g].gene_name,n_FIRST_union=int((unitset.position_class_consensus=='FIRST').sum()),n_INTERNAL_union=int((unitset.position_class_consensus=='INTERNAL').sum()),FIRST_centers=int((N[g,:,0]>0).sum()),INTERNAL_centers=int((N[g,:,1]>0).sum()),n_region_centers=len(cs),shared_every_center=bool(shared[g]),positive_every_required_density=bool(positive[g]),common_transcript=bool(coherent[g]),common_transcript_ids=','.join(sorted(set.intersection(*txsets[g]))),exclusion_reason='ELIGIBLE_POSITIVE_SHARED' if positive[g] else 'MISSING_CLASS_AT_ONE_OR_MORE_CENTERS' if not shared[g] else 'ZERO_CLASS_DENSITY_IN_REQUIRED_LIBRARY',zero_required_densities=int((M[g,...,active]==0).sum()) if shared[g] else np.nan,annotated_tmd_count=ct[g],cds_length=length[g],n_annotation_transcripts=int(gann.iloc[g].n_annotation_transcripts),tmd_count_min=gann.iloc[g].tmd_count_min,tmd_count_max=gann.iloc[g].tmd_count_max,cds_length_min=gann.iloc[g].cds_length_min,cds_length_max=gann.iloc[g].cds_length_max))
    pd.DataFrame(dict(gene_id=positivegenes,**{m:J[:,j] for j,m in enumerate(metrics)})).to_csv(E/f'within_gene_contrasts_{region}.tsv',sep='\t',index=False)
    zero=[]
    for g in np.where(shared & ~positive)[0]:
        for c,k,col in zip(*np.where(M[g][...,active]==0)):
            orig=active[col];zero.append(dict(region=region,gene_id=genes[g],center=cs[c],topology_class=['FIRST','INTERNAL'][k],blind_id=f'B{orig//2+1:02d}',density_type=['real','pseudo'][orig%2]))
    pd.DataFrame(zero).to_csv(E/f'positivity_exclusions_{region}.tsv',sep='\t',index=False)
    frozen=effects(S.sum(axis=0))[-1]
    evaluate(region,rindex,names[0],genes,T=S,frozen_point=frozen)
    ref=pd.read_csv(R/'reproduction_output/available_case/REGION_COMPOSITE_SUMMARY_v1.tsv',sep='\t').set_index('region').loc[region]
    got=results[-1]
    for key,rkey in [('interaction_mean','MITO_minus_ER_mean'),('interaction_rep1','MITO_minus_ER_rep1'),('interaction_rep2','MITO_minus_ER_rep2'),('interaction_mean_percentile_lo','interaction_CI_2.5'),('interaction_mean_percentile_hi','interaction_CI_97.5')]:
        error=abs(got[key]-float(ref[rkey]));comparisons.append(dict(region=region,metric=key,independent=got[key],frozen=float(ref[rkey]),absolute_error=error,tolerance=1e-10,status='PASS' if error<=1e-10 else 'FAIL'))
        assert error<=1e-10,(region,key,error)
    evaluate(region,rindex,names[1],genes,T=M,frozen_point=frozen)
    for name,array,sel in [(names[2],S,shared),(names[3],M,shared),(names[4],S,positive),(names[5],M,positive)]:evaluate(region,rindex,name,list(np.array(genes)[sel]),T=array[sel],frozen_point=frozen)
    evaluate(region,rindex,names[6],list(positivegenes),V=J,method='mean',frozen_point=frozen)
    evaluate(region,rindex,names[7],list(positivegenes),V=J,method='median',frozen_point=frozen)
    sel=coherent[positive]
    evaluate(region,rindex,names[8],list(positivegenes[sel]),V=J[sel],method='mean',frozen_point=frozen)
    # Continuous baseline/composition diagnostics for all, shared and positive populations.
    variables={'annotated_tmd_count':np.broadcast_to(ct[:,None,None],N.shape),'cds_length':np.broadcast_to(length[:,None,None],N.shape),'relative_cds_position':pos,'absolute_cds_start':absolute,'ordinal_index':ordinal,'input_coverage':coverage,'mask_fraction':mask,'retained_unit_count':N.astype(float)}
    for label,sel0,w0 in [('ALL_UNIT',np.ones(G,bool),N),('ALL_EQUAL_GENE',np.ones(G,bool),(N>0).astype(float)),('SHARED_EQUAL_GENE',shared,(N>0).astype(float)),('POSITIVE_SHARED_EQUAL_GENE',positive,(N>0).astype(float))]:
        for ci,c in enumerate(cs):
            weight=w0[:,ci]*sel0[:,None]
            for key,value in variables.items():summary_cov(region,label,c,key,value[:,ci],weight)
    # Save genuine density contribution concentration separately from nominal gene weights.
    for label,T,selection in [('FROZEN',S,np.ones(G,bool)),('EQUAL_GENE',M,np.ones(G,bool)),('SHARED_EQUAL_GENE',M,shared),('POSITIVE_EQUAL_GENE',M,positive)]:
        for ci,c in enumerate(cs):
            for k,cl in enumerate(['FIRST','INTERNAL']):
                for col in active:
                    v=T[:,ci,k,col]*selection;total=v.sum();v=v/total if total>0 else v
                    unit_contrib.append(dict(region=region,design=label,center=c,topology_class=cl,blind_id=f'B{col//2+1:02d}',density_type=['real','pseudo'][col%2],effective_density_contributing_genes=float(1/(v@v)) if v@v>0 else 0,max_gene_density_share=float(v.max()),max_contributing_gene=genes[int(v.argmax())]))
    for design in ['CORE','EXTENDED']:
        kept=np.zeros(N.shape,dtype=bool);weights=np.zeros(N.shape);supported=True;reasons=[]
        for ci,c in enumerate(cs):
            positive_cov=coverage[:,ci][(N[:,ci]>0)&(coverage[:,ci]>0)];cut=float(np.median(positive_cov)) if len(positive_cov) else np.nan
            strata={}
            for g in range(G):
                for k in range(2):
                    if N[g,ci,k]==0:continue
                    if not (np.isfinite(ct[g]) and ct[g]>=2 and np.isfinite(length[g])):continue
                    tb=int(np.searchsorted([3,5,9],ct[g],side='right'));lb=int(np.searchsorted([1000,2000],length[g],side='right'))
                    key=(tb,lb)
                    if design=='EXTENDED':
                        if not (np.isfinite(pos[g,ci,k]) and 0<=pos[g,ci,k]<=1):continue
                        pb=int(np.searchsorted([.25,.5,.75],pos[g,ci,k],side='right'));cv=coverage[g,ci,k];cb=0 if cv==0 else 1 if cv<=cut else 2;mb=int(mask[g,ci,k]>=.5)
                        key=key+(pb,cb,mb)
                    strata.setdefault(key,[[],[]])[k].append(g)
            eligible_strata={key:members for key,members in strata.items() if min(map(len,members))>=5}
            mass=sum(min(map(len,v)) for v in eligible_strata.values())
            for key,members in strata.items():
                retained=key in eligible_strata;q=min(map(len,members))/mass if retained else 0
                support_rows.append(dict(region=region,design=design,center=c,stratum=';'.join(map(str,key)),n_FIRST=len(members[0]),n_INTERNAL=len(members[1]),retained=retained,target_mass=q,positive_input_median_cut=cut))
                for k,ids in enumerate(members):
                    for g in ids:
                        w=q/len(ids) if retained else 0;weights[g,ci,k]=w;kept[g,ci,k]=retained
                        balance_records.append(dict(region=region,design=design,center=c,gene_id=genes[g],topology_class=['FIRST','INTERNAL'][k],stratum=';'.join(map(str,key)),retained=retained,weight=w,n_units=int(N[g,ci,k]),relative_cds_position=pos[g,ci,k],annotated_tmd_count=ct[g],cds_length=length[g],input_coverage=coverage[g,ci,k],mask_fraction=mask[g,ci,k]))
            for k in range(2):
                w=weights[:,ci,k];n=int((w>0).sum());ess=float(w.sum()**2/(w@w)) if w@w>0 else 0;mx=float(w.max()/w.sum()) if w.sum()>0 else 1
                if n<30 or ess<30 or mx>.10:supported=False;reasons.append(f'center{c}/class{k}:n={n},ESS={ess:.3f},max={mx:.3f}')
            for label,w0 in [('RESTRICTED_UNIT',N[:,ci]*kept[:,ci]),('RESTRICTED_EQUAL_GENE',kept[:,ci].astype(float)),('WEIGHTED',weights[:,ci])]:
                for key,value in variables.items():summary_cov(region,design+'_'+label,c,key,value[:,ci],w0)
        gids=kept.any(axis=(1,2));glist=list(np.array(genes)[gids])
        for suffix,T in [('RESTRICTED_UNIT',S*kept[...,None]),('RESTRICTED_EQUAL_GENE',M*kept[...,None]),('WEIGHTED',M*weights[...,None])]:
            name='BALANCE_'+design+'_'+suffix
            if not supported:
                # Report points only within observed retained strata; no unsupported bootstrap claim.
                p=effects(T.sum(axis=0)) if len(glist) else np.full(7,np.nan)
                results.append(dict(region=region,estimand=name,n_genes=len(glist),status='UNSUPPORTED_COMMON_SUPPORT',support_status='UNSUPPORTED',support_reason=';'.join(reasons),**{m:float(p[i]) for i,m in enumerate(metrics)},protocol_sha256=freeze['sha256']))
                print(region,name,'UNSUPPORTED',reasons,flush=True)
            else:evaluate(region,rindex,name,glist,T=T[gids],supported=supported,reason=';'.join(reasons),frozen_point=frozen)
    # Prespecified descriptive TMD count strata, no new intervals or searches.
    cb=np.searchsorted([3,5,9],ct,side='right')
    for s in range(4):
        for label,select in [('SHARED_EQUAL_GENE',shared&(cb==s)&np.isfinite(ct)),('WITHIN_GENE_MEAN',positive&(cb==s)&np.isfinite(ct))]:
            n=int(select.sum());status='COMPLETE' if n>=30 else 'SPARSE_NOT_REPORTED';pt=np.nan
            if n>=30:
                pt=float(effects(M[select].sum(axis=0))[-1]) if label=='SHARED_EQUAL_GENE' else float(gene_effects(M[select])[:,-1].mean())
                if not np.isfinite(pt):status='NOT_ESTIMABLE_NONPOSITIVE_AGGREGATE'
            strata_results.append(dict(region=region,estimand=label,tmd_count_bin=['[2,3)','[3,5)','[5,9)','[9,infinity)'][s],n_genes=n,status=status,interaction=pt))
    # Save at each completed region so interruptions preserve results.
    write_outputs()

def write_outputs():
    for rows,path in [(results,W/'P3_WITHIN_GENE_INDEPENDENT_RESULTS.tsv'),(eligibility,W/'P3_WITHIN_GENE_ELIGIBILITY_AND_BALANCE.tsv'),(influence,E/'gene_deletion_influence.tsv'),(diagnostics,E/'balance_diagnostics.tsv'),(balance_records,E/'balance_weights.tsv'),(support_rows,E/'common_support_strata.tsv'),(comparisons,E/'frozen_reproduction_comparison.tsv'),(strata_results,E/'tmd_count_stratum_effects.tsv'),(unit_contrib,E/'effective_density_contributions.tsv')]:pd.DataFrame(rows).to_csv(path,sep='\t',index=False)

(E/'calculation_start.json').write_text(json.dumps({'utc':datetime.now(timezone.utc).isoformat(),'protocol_sha256':freeze['sha256'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'conventions_sha256':hashlib.sha256((E/'PRE_RESULT_IMPLEMENTATION_CONVENTIONS.md').read_bytes()).hexdigest()},indent=2)+'\n')
for args in [('PLATEAU_90_120',[90,105,120],1),('BROAD_60_120',[60,75,90,105,120],2)]:analyze(*args)
write_outputs()
(E/'calculation_complete.json').write_text(json.dumps({'utc':datetime.now(timezone.utc).isoformat(),'variants':len(results),'protocol_sha256':freeze['sha256'],'frozen_comparisons_passed':all(r['status']=='PASS' for r in comparisons)},indent=2)+'\n')
print('COMPLETED_ALL_PRESPECIFIED_BRANCHES',flush=True)
