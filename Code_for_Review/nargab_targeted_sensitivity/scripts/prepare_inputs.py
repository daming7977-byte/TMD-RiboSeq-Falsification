#!/usr/bin/env python3
"""Independent source-density reconstruction and exact annotation joins; no effects."""
from pathlib import Path
import hashlib, json, sys, platform
from datetime import datetime, timezone
import numpy as np
import pandas as pd

import os
P=Path(__file__).resolve().parents[1]
W=Path(os.environ['P3_SENSITIVITY_OUTPUT']).resolve()
R=P.parent; F=R/'frozen_intermediates/region_composite_project'; E=W/'evidence'
hashes=[]
release=pd.read_csv(P/'INPUT_MANIFEST.tsv',sep='\t').set_index('relative_to_code_package')
def load(path, expected=None):
    b=path.read_bytes(); sha=hashlib.sha256(b).hexdigest(); rel=str(path.relative_to(R))
    assert rel in release.index, rel
    assert release.loc[rel,'sha256']==sha, (rel,'audit input checksum')
    assert int(release.loc[rel,'bytes'])==len(b)
    if expected: assert sha==expected
    hashes.append(dict(path=rel,bytes=len(b),sha256=sha,packaged_audit_manifest_verified=True,explicit_baseline_verified=bool(expected)))
    return pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)

freeze=json.loads((E/'PROTOCOL_FREEZE.json').read_text())
assert hashlib.sha256((W/'P3_WITHIN_GENE_PRESPECIFIED_PROTOCOL.md').read_bytes()).hexdigest()==freeze['sha256']
manifest=load(R/'example_manifests/GSE300977/biological_reveal/GSE300977_BIOLOGICAL_REVEAL_MANIFEST_v1r1.tsv')
assert set(manifest.blind_id)=={f'B{i:02d}' for i in range(1,13)}
assert not manifest.duplicated(['location','fraction','replicate']).any()
expected={'ER':{'1':('B01','B02'),'2':('B11','B12')},'mitochondria':{'1':('B09','B10'),'2':('B07','B08')},'cytosol':{'1':('B05','B06'),'2':('B03','B04')}}
for loc, rr in expected.items():
    for rep,(pdid,inid) in rr.items():
        for frac,bid in [('pulldown',pdid),('input',inid)]:
            row=manifest[(manifest.location==loc)&(manifest.replicate==rep)&(manifest.fraction==frac)]
            assert len(row)==1 and row.iloc[0].blind_id==bid
manifest.to_csv(E/'sample_manifest.tsv',sep='\t',index=False)
catalog=load(F/'tmd_spatial_secondary_scan_v1/output/SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv')
assert not catalog.geometry_id.duplicated().any()
key=list(zip(catalog.strand,catalog.geometry,catalog.mask_bits))
assert len(set(key))==len(key)
for row in catalog.itertuples(index=False):
    assert len(row.mask_bits)==31 and set(row.mask_bits)<={'0','1'}
    assert row.mask_bits.count('0')==int(row.unmasked_nt)>0
    assert sum(int(p.split(':')[1].split('-')[1])-int(p.split(':')[1].split('-')[0]) for p in row.geometry.split('^'))==31
density=[]
for i in range(1,13):
    d=load(F/f'tmd_spatial_secondary_scan_v1/geometry_counts/B{i:02d}.SECONDARY_GEOMETRY_DENSITIES_v1.tsv')
    assert not d.geometry_id.duplicated().any() and set(d.geometry_id)==set(catalog.geometry_id)
    d=d.set_index('geometry_id').loc[catalog.geometry_id]
    assert np.array_equal(d.unmasked_nt.astype(int),catalog.unmasked_nt.astype(int))
    v=d.density.astype(float).to_numpy();assert np.isfinite(v).all() and (v>=0).all()
    density.append(v)
density=np.stack(density,axis=1)
geoid={k:i for i,k in enumerate(key)}
formal=load(F/'tmd_anchor_preregister/bridge/FORMAL_TMD_ANCHOR_CANDIDATES_v1p1.tsv')
fidx=formal.set_index(['target_transcript','TMD_index']);assert not fidx.index.duplicated().any()
tx=formal[['target_transcript','gene_id','pass_class','n_primary_TMDs','target_cds_length']].drop_duplicates().set_index('target_transcript')
assert not tx.index.duplicated().any()
genomic=load(P/'inputs/GSE300977_FORMAL_GENOMIC_TMD_ANCHORS_v1.tsv',expected='6f065fad3c851cbd6bcc07e2871d3459e1b395a684b7ec57ad19ea440f0db86f')
for row in genomic.itertuples(index=False):
    f=fidx.loc[(row.target_transcript,row.TMD_index)]
    assert f.gene_id==row.gene_id and float(f.target_cds_length)==float(row.target_cds_length)
    assert int(f.target_tmd_cds_nt_start0)==int(row.tmd_cds_nt_start0)
    assert f.TMD_position_class==row.TMD_position_class
gidx=genomic.set_index(['target_transcript','chrom','strand','anchor_genomic_start0'])
assert not gidx.index.duplicated().any()
pair=load(F/'crossmap_window_preregister/output/GSE300977_PAIRED_REAL_PSEUDO_CROSSMAP_MASK_AUDIT_v1.tsv')
assert not pair.paired_real_unit_id.duplicated().any()
pair['center']=pair.center.astype(int)
pair=pair[pair.center.isin([60,75,90,105,120])].copy()
row_records=[]; mats=[]; missing=[]
for r in pair.itertuples(index=False):
    targets=r.target_transcripts.split(',')
    linked=[gidx.loc[(t,r.chrom,r.strand,r.real_anchor_genomic_start0)] for t in targets]
    assert {x.gene_id for x in linked}==set(r.gene_ids.split(','))
    assert {str(x.TMD_index) for x in linked}==set(r.tmd_indices.split(','))
    classes=set(tx.loc[targets,'pass_class']); topo=next(iter(classes)) if len(classes)==1 else 'MIXED_TOPOLOGY'
    structural=(topo=='MULTI_PASS' and r.position_class_consensus in ['FIRST','INTERNAL'] and r.unique_gene_assignment=='YES' and r.n_gene_ids=='1')
    maskok=int(r.real_crossmap_unmasked_nt)>=1 and int(r.pseudo_min_unmasked_nt)>=1
    missing.append({'paired_real_unit_id':r.paired_real_unit_id,'center':r.center,'gene_ids':r.gene_ids,'topology':topo,'class':r.position_class_consensus,'structural_eligible':structural,'mask_eligible':maskok,'real_unmasked':int(r.real_crossmap_unmasked_nt),'pseudo_min_unmasked':int(r.pseudo_min_unmasked_nt)})
    if not (structural and maskok): continue
    real=geoid[(r.strand,r.real_window_geometry,r.real_crossmap_mask_bits)]
    pg=r.pseudo_window_geometries.split(';'); pm=r.pseudo_crossmap_mask_bits.split(';')
    assert len(pg)==len(pm)==int(r.n_distinct_pseudo_geometries) and len(set(zip(pg,pm)))==len(pg)
    pseudo=[geoid[(r.strand,g,m)] for g,m in zip(pg,pm)]
    assert int(catalog.iloc[real].unmasked_nt)==int(r.real_crossmap_unmasked_nt)
    assert min(int(catalog.iloc[i].unmasked_nt) for i in pseudo)==int(r.pseudo_min_unmasked_nt)
    # Library x real/pseudo: even entries real, odd pseudo.
    mats.append(np.column_stack((density[real],density[pseudo].mean(axis=0))).ravel())
    rr=r._asdict();rr.update(topology_class=topo,
        relative_cds_position=np.mean([int(x.tmd_cds_nt_start0)/float(x.target_cds_length) for x in linked]),
        absolute_cds_start=np.mean([int(x.tmd_cds_nt_start0) for x in linked]),
        ordinal_mean=np.mean([int(x.TMD_index) for x in linked]),
        mask_fraction=min(int(r.real_crossmap_unmasked_nt),int(r.pseudo_min_unmasked_nt))/31,
        unit_key='|'.join(str(getattr(r,k)) for k in ['chrom','strand','real_anchor_genomic_start0','target_transcripts','gene_ids','position_class_consensus','tmd_indices']))
    row_records.append(rr)
u=pd.DataFrame(row_records); X=np.asarray(mats)
assert not u.duplicated(['center','unit_key']).any()
assert np.isfinite(X).all() and (X>=0).all()
# Exact reconstruction comparison to the independently verified audit density identity.
identity=json.loads((P/'expected/density_identity.json').read_text())
assert list(X.shape)==identity['shape'] and str(X.dtype)==identity['dtype']
assert hashlib.sha256(np.ascontiguousarray(X,dtype=np.float64).tobytes()).hexdigest()==identity['contiguous_float64_sha256']
maxerr=0.0

gene_rows=[]
for gene,grp in u.groupby('gene_ids',sort=True):
    tids=sorted(set(t for s in grp.target_transcripts for t in s.split(',')))
    t=tx.loc[tids];assert set(t.gene_id)=={gene}
    n=t.n_primary_TMDs.astype(float);l=t.target_cds_length.astype(float)
    gene_rows.append(dict(gene_id=gene,gene_name=','.join(sorted(set(grp.gene_names))),annotation_transcripts=','.join(tids),n_annotation_transcripts=len(tids),annotated_tmd_count=float(n.median()),tmd_count_min=float(n.min()),tmd_count_max=float(n.max()),cds_length=float(l.median()),cds_length_min=float(l.min()),cds_length_max=float(l.max())))
u.to_csv(E/'verified_units.tsv',sep='\t',index=False);np.save(E/'verified_unit_densities.npy',X)
pd.DataFrame(gene_rows).to_csv(E/'gene_annotation.tsv',sep='\t',index=False)
pd.DataFrame(missing).to_csv(E/'pre_filter_eligibility.tsv',sep='\t',index=False)
pd.DataFrame(hashes).to_csv(E/'INPUT_SHA256.tsv',sep='\t',index=False)
env={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':'not installed; not required; binomial ranks use standard-library recurrence','started_after_protocol_freeze':freeze['frozen_utc'],'completed_utc':datetime.now(timezone.utc).isoformat()}
(E/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
checks={'verified_input_files':len(hashes),'geometry_count':len(catalog),'five_center_prefilter_rows':len(pair),'eligible_region_center_units':len(u),'unique_genes':u.gene_ids.nunique(),'density_shape':list(X.shape),'previous_independent_density_max_abs_error':maxerr,'protocol_sha256':freeze['sha256'],'compartment_mapping':'PASS','exact_transcript_anchor_linkage':'PASS','finite_nonnegative_density':'PASS'}
(E/'input_verification.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2),flush=True)
