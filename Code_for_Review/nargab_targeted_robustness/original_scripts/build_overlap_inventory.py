#!/usr/bin/env python3
"""Read-only annotation inventory on frozen 31-nt spliced window components."""
from pathlib import Path
import csv, re, hashlib, json, sys, platform, time, bisect
from collections import defaultdict, Counter
from functools import lru_cache

OUT=Path(__file__).resolve().parent
F=Path('historical_location/Documents/P3_FINAL_COVER_LETTER_MICRO_EDIT_20261018_v1/check/reviewer_parts/Code_for_Review/frozen_intermediates/region_composite_project')
H=Path('historical_location/paper3_GSE300977_raw')
GTF=H/'mapping_reference_blinded/Homo_sapiens.GRCh38.104.gtf'
PAIR=F/'crossmap_window_preregister/output/GSE300977_PAIRED_REAL_PSEUDO_CROSSMAP_MASK_AUDIT_v1.tsv'
FORMAL=F/'tmd_anchor_preregister/bridge/FORMAL_TMD_ANCHOR_CANDIDATES_v1p1.tsv'
CAT=F/'tmd_spatial_secondary_scan_v1/output/SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv'
PRIM=H/'tmd_spatial_primary_reveal_v1/output/PRIMARY_PLUS75_UNIT_GEOMETRY_MAP_v1.tsv'
START=time.time(); INPUTS={}; counts=Counter()
def log(s):print(s,flush=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 INPUTS[str(p)]=sha(p)
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields=None):
 tmp=p.with_name(p.name+'.tmp')
 with tmp.open('w') as f:
  w=csv.DictWriter(f,delimiter='\t',fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
 tmp.replace(p)
def yn(x):return 'YES' if x else 'NO'
def tokens(s):return [t for t in s.split(',') if t]
def blocks(s):
 out=[]
 for token in s.split('^'):
  chrom,coords=token.split(':');a,b=map(int,coords.split('-'));assert b>a;out.append((chrom,a,b))
 assert sum(b-a for _,a,b in out)==31,s
 return out
@lru_cache(None)
def positions(geom,strand):
 return tuple((c,p) for c,a,b in blocks(geom) for p in (range(a,b) if strand=='+' else range(b-1,a-1,-1)))
formal=read(FORMAL); classes=defaultdict(set)
for r in formal:classes[r['target_transcript']].add(r['pass_class'])
assert all(len(v)==1 for v in classes.values())
prim_map={r['paired_real_unit_id']:r for r in read(PRIM)}
cat={ (r['strand'],r['geometry'],r['mask_bits']):r['geometry_id'] for r in read(CAT)}
pairs=[];txs=set()
for r in read(PAIR):
 c=int(r['center'])
 if c not in [60,75,90,105,120]:continue
 if int(r['real_crossmap_unmasked_nt'])<1 or int(r['pseudo_min_unmasked_nt'])<1:continue
 cs=set().union(*(classes[t] for t in tokens(r['target_transcripts'])))
 r['topology_class']=next(iter(cs)) if len(cs)==1 else 'MIXED'
 r['primary_included']=yn(r['paired_real_unit_id'] in prim_map)
 region=r['topology_class']=='MULTI_PASS' and r['position_class_consensus'] in ['FIRST','INTERNAL'] and r['unique_gene_assignment']=='YES' and int(r['n_gene_ids'])==1
 r['plateau_included']=yn(region and c in [90,105,120]);r['broad_included']=yn(region)
 if r['primary_included']==r['plateau_included']==r['broad_included']=='NO':continue
 r['unit_key']='|'.join(r[k] for k in ['chrom','strand','real_anchor_genomic_start0','target_transcripts','gene_ids','position_class_consensus','tmd_indices'])
 pairs.append(r);txs.update(tokens(r['target_transcripts']))
assert sum(r['primary_included']=='YES' for r in pairs)==len(prim_map)==9098
assert len({r['paired_real_unit_id'] for r in pairs})==len(pairs)
log(f'Loaded {len(pairs)} selected center-specific pairs, {len(txs)} target transcripts')

# Pass 1: exact target paths and metadata. Hash GTF while streaming it.
ATTR=re.compile(r'(\w+) "([^"]*)"')
def attrs(s):
 result=defaultdict(list)
 for k,v in ATTR.findall(s):result[k].append(v)
 return result
features=defaultdict(lambda:defaultdict(list));txmeta={};genemeta={};h=hashlib.sha256()
with GTF.open('rb') as f:
 for ln,raw in enumerate(f,1):
  h.update(raw)
  if raw.startswith(b'#'):continue
  s=raw.decode().rstrip('\n').split('\t'); typ=s[2]
  if typ not in ['gene','transcript','exon','CDS','stop_codon']:continue
  a=attrs(s[8]);gid=a['gene_id'][0];tx=a.get('transcript_id',[''])[0]
  if typ=='gene':genemeta[gid]={'name':a.get('gene_name',[''])[0],'version':a.get('gene_version',[''])[0],'biotype':a.get('gene_biotype',[''])[0],'line':ln,'chrom':s[0],'strand':s[6]}
  if tx not in txs:continue
  if typ=='transcript':txmeta[tx]={'gene_id':gid,'version':a.get('transcript_version',[''])[0],'tags':a.get('tag',[]),'support':a.get('transcript_support_level',[''])[0],'line':ln,'chrom':s[0],'strand':s[6]}
  elif typ in ['exon','CDS','stop_codon']:features[tx][typ].append((s[0],int(s[3])-1,int(s[4]),s[6],ln))
INPUTS[str(GTF)]=h.hexdigest();assert h.hexdigest()=='f3de2bccc02877709272a0ce66af59a2922d407a4076541bc5a65c9d3518581e'
assert txs<=txmeta.keys(),txs-txmeta.keys()
log(f'GTF pass1 done; hash verified; {len(genemeta)} annotated genes')

# Bins are a retrieval accelerator ONLY. All overlap tests below use exact exon/base coordinates.
BINSIZE=1000;querybins=defaultdict(set)
for r in pairs:
 for g in [r['real_window_geometry']]+r['pseudo_window_geometries'].split(';'):
  for c,a,b in blocks(g):querybins[(c,r['strand'])].update(range(a//BINSIZE,(b-1)//BINSIZE+1))
# Full target exons are also queried to identify structural shared-exon/readthrough-like candidates.
for tx in txs:
 for c,a,b,st,ln in features[tx]['exon']:querybins[(c,st)].update(range(a//BINSIZE,(b-1)//BINSIZE+1))
querybins={k:sorted(v) for k,v in querybins.items()};idx=defaultdict(list);exons=[]
with GTF.open() as f:
 for ln,line in enumerate(f,1):
  if line.startswith('#'):continue
  s=line.rstrip('\n').split('\t')
  if s[2]!='exon':continue
  key=(s[0],s[6]);q=querybins.get(key)
  if not q:continue
  start,end=int(s[3])-1,int(s[4]);i=bisect.bisect_left(q,start//BINSIZE);j=bisect.bisect_right(q,(end-1)//BINSIZE)
  if i==j:continue
  a=attrs(s[8]);eid=len(exons);exons.append((start,end,a['gene_id'][0],a['transcript_id'][0],ln,a.get('exon_id',[''])[0]))
  for bn in q[i:j]:idx[(s[0],s[6],bn)].append(eid)
log(f'GTF pass2 done; {len(exons)} relevant exon records indexed')

def query(c,st,a,b):
 ids=set()
 for bn in range(a//BINSIZE,(b-1)//BINSIZE+1):ids.update(idx.get((c,st,bn),[]))
 return (exons[i] for i in ids if exons[i][0]<b and exons[i][1]>a)
# Distinguish explicit tags from an operational structural candidate: target transcript
# shares complete exons with two other protein-coding genes whose shared exon
# ranges are spatially separated. Shared constant exons alone are not readthrough.
struct={};struct_rows=[]
for tx in sorted(txs):
 own=txmeta[tx]['gene_id'];gene_exons=defaultdict(set);ev=defaultdict(set)
 for c,a,b,st,ln in features[tx]['exon']:
  for x,y,g,t,l,e in query(c,st,a,b):
   if g!=own and x==a and y==b:gene_exons[g].add((c,a,b));ev[g].add(l)
 allshared=set().union(*gene_exons.values()) if gene_exons else set()
 coding=[g for g in gene_exons if genemeta[g]['biotype']=='protein_coding']
 candidate=False
 for ga in coding:
  for gb in coding:
   if ga==gb:continue
   aa=gene_exons[ga];bb=gene_exons[gb]
   if {x[0] for x in aa}=={x[0] for x in bb} and max(x[2] for x in aa)<=min(x[1] for x in bb):candidate=True
 explicit=any('readthrough' in z.lower() for z in txmeta[tx]['tags'])
 struct[tx]=(candidate,explicit,tuple(sorted(gene_exons)))
 if gene_exons:
  struct_rows.append({'target_transcript':tx,'target_gene':own,'transcript_version':txmeta[tx]['version'],'gtf_transcript_line':txmeta[tx]['line'],'tags':','.join(txmeta[tx]['tags']),'explicit_readthrough_tag':yn(explicit),'structural_readthrough_like_candidate':yn(candidate),'other_gene_ids_exact_shared_exons':','.join(sorted(gene_exons)),'other_gene_names':','.join(genemeta[g]['name'] or g for g in sorted(gene_exons)),'n_exact_shared_exons':len(allshared),'shared_exon_gtf_lines':';'.join(f'{g}:'+','.join(map(str,sorted(ev[g]))) for g in sorted(ev))})
write(OUT/'target_transcript_shared_exon_structure.tsv',struct_rows)

# Confirm a frozen 31-nt vector is a continuous path through annotated target CDS(+stop).
paths={}
for tx in txs:
 fs=features[tx]['CDS']+features[tx]['stop_codon'];st=txmeta[tx]['strand']
 fs=sorted(fs,key=lambda x:x[1],reverse=(st=='-'));cum=0;p=[]
 for c,a,b,st,ln in fs:p.append((c,a,b,st,cum));cum+=b-a
 paths[tx]=p
@lru_cache(None)
def valid_path(tx,geom,strand):
 pp=positions(geom,strand);mapped=[]
 for c,pos in pp:
  choices=[cum+(pos-a if st=='+' else b-1-pos) for cc,a,b,st,cum in paths[tx] if cc==c and a<=pos<b and st==strand]
  if len(choices)!=1:return False
  mapped.append(choices[0])
 return all(b==a+1 for a,b in zip(mapped,mapped[1:]))

@lru_cache(None)
def inspect(geom,strand,bits,owners_s,txs_s):
 pp=positions(geom,strand);assert len(pp)==len(bits)==31 and set(bits)<=set('01')
 owners=set(tokens(owners_s));targets=set(tokens(txs_s));unmasked={i for i,b in enumerate(bits) if b=='0'}
 by_gene=defaultdict(set);by_gene_unmask=defaultdict(set);same_other_tx=set();own_target=set();evidence=defaultdict(set);other_tx=defaultdict(set)
 for c,a,b in blocks(geom):
  for x,y,g,t,l,e in query(c,strand,a,b):
   inds={i for i,(cc,p) in enumerate(pp) if cc==c and x<=p<y}
   if not inds:continue
   if g in owners:
    if t not in targets:same_other_tx.update(inds)
    else:own_target.update(inds)
   else:
    by_gene[g].update(inds);by_gene_unmask[g].update(inds&unmasked);evidence[g].add(l);other_tx[g].add(t)
 full=set().union(*by_gene.values()) if by_gene else set();um=set().union(*by_gene_unmask.values()) if by_gene_unmask else set()
 compatible=[t for t in sorted(targets) if valid_path(t,geom,strand)]
 candidates=[t for t in compatible if struct[t][0]];explicit=[t for t in compatible if struct[t][1]]
 return {'full_other_gene_overlap_nt':len(full),'unmasked_other_gene_overlap_nt':len(um),'other_gene_ids_full':','.join(sorted(by_gene)),'other_gene_ids_unmasked':','.join(g for g in sorted(by_gene_unmask) if by_gene_unmask[g]),'other_gene_names_full':','.join(genemeta[g]['name'] or g for g in sorted(by_gene)),'other_gene_nt_full':';'.join(f'{g}:{len(by_gene[g])}' for g in sorted(by_gene)),'other_gene_nt_unmasked':';'.join(f'{g}:{len(by_gene_unmask[g])}' for g in sorted(by_gene_unmask) if by_gene_unmask[g]),'other_transcript_ids':';'.join(g+':'+','.join(sorted(other_tx[g])) for g in sorted(other_tx)),'gtf_other_exon_lines':';'.join(g+':'+','.join(map(str,sorted(evidence[g]))) for g in sorted(evidence)),'same_gene_other_transcript_overlap_nt':len(same_other_tx),'same_gene_other_transcript_unmasked_nt':len(same_other_tx&unmasked),'target_transcript_exon_coverage_nt':len(own_target),'continuous_target_cds_stop_path_valid':yn(bool(compatible)),'compatible_target_transcripts':','.join(compatible),'n_compatible_target_transcripts':len(compatible),'all_named_targets_path_compatible':yn(len(compatible)==len(targets)),'structural_readthrough_like_compatible_transcripts':','.join(candidates),'explicit_readthrough_tag_compatible_transcripts':','.join(explicit),'read_origin_status':'UNRESOLVED_DIFFERENT_GENE_SHARED_EXON' if um else ('OTHER_GENE_OVERLAP_MASKED' if full else 'NO_SAME_STRAND_OTHER_GENE_EXON_OVERLAP_AT_COUNTED_BASES')}

inventory=[];flags=[]
for no,r in enumerate(pairs):
 analyses=[]
 if r['primary_included']=='YES':analyses.append('PRIMARY_PLUS75')
 if r['plateau_included']=='YES':analyses.append('PLATEAU_90_120')
 if r['broad_included']=='YES':analyses.append('BROAD_60_120')
 pgs=r['pseudo_window_geometries'].split(';');pms=r['pseudo_crossmap_mask_bits'].split(';');assert len(pgs)==len(pms)==int(r['n_distinct_pseudo_geometries'])
 windows=[('real',0,r['real_window_geometry'],r['real_crossmap_mask_bits'])]+[('pseudo',i,g,m) for i,(g,m) in enumerate(zip(pgs,pms),1)]
 checks=[]
 for side,pi,g,m in windows:
  d=inspect(g,r['strand'],m,r['gene_ids'],r['target_transcripts']);checks.append(d)
  row={'paired_real_unit_id':r['paired_real_unit_id'],'center':r['center'],'unit_key':r['unit_key'],'chrom':r['chrom'],'strand':r['strand'],'topology_class':r['topology_class'],'position_class':r['position_class_consensus'],'owner_gene_ids':r['gene_ids'],'owner_gene_names':r['gene_names'],'target_transcripts':r['target_transcripts'],'side':side,'pseudo_index':pi,'n_pseudos':len(pgs),'pseudo_aggregation_weight':1/len(pgs) if side=='pseudo' else 1,'geometry_id':cat[(r['strand'],g,m)],'geometry':g,'mask_bits':m,'unmasked_nt':m.count('0'),**d}
  for analysis in analyses:inventory.append({'analysis':analysis,**row})
 real=checks[0];ps=checks[1:]
 flags.append({'paired_real_unit_id':r['paired_real_unit_id'],'center':r['center'],'unit_key':r['unit_key'],'owner_gene_ids':r['gene_ids'],'owner_gene_names':r['gene_names'],'target_transcripts':r['target_transcripts'],'topology_class':r['topology_class'],'position_class':r['position_class_consensus'],'primary_included':r['primary_included'],'plateau_included':r['plateau_included'],'broad_included':r['broad_included'],'real_unmasked_other_gene':yn(real['unmasked_other_gene_overlap_nt']>0),'pseudo_any_unmasked_other_gene':yn(any(x['unmasked_other_gene_overlap_nt']>0 for x in ps)),'exclude_pair_unmasked_other_gene':yn(any(x['unmasked_other_gene_overlap_nt']>0 for x in checks)),'real_full_other_gene':yn(real['full_other_gene_overlap_nt']>0),'pseudo_any_full_other_gene':yn(any(x['full_other_gene_overlap_nt']>0 for x in ps)),'any_invalid_target_cds_stop_path':yn(any(x['continuous_target_cds_stop_path_valid']=='NO' for x in checks)),'other_gene_ids_unmasked':','.join(sorted(set(g for x in checks for g in tokens(x['other_gene_ids_unmasked']))))})
 if no%5000==0:log(f'Processed {no}/{len(pairs)} paired units')
write(OUT/'pair_overlap_flags.tsv',flags);write(OUT/'P3_OVERLAP_INVENTORY.tsv',inventory)
log(f'FLAGS_READY {OUT / "pair_overlap_flags.tsv"}')

summary=[]
for analysis in ['PRIMARY_PLUS75','PLATEAU_90_120','BROAD_60_120']:
 ar=[r for r in inventory if r['analysis']==analysis]
 for cl in sorted(set(r['position_class'] for r in ar)):
  for top in sorted(set(r['topology_class'] for r in ar if r['position_class']==cl)):
   for side in ['real','pseudo']:
    sr=[r for r in ar if r['position_class']==cl and r['topology_class']==top and r['side']==side]
    if not sr:continue
    full=[r for r in sr if r['full_other_gene_overlap_nt']>0];u=[r for r in sr if r['unmasked_other_gene_overlap_nt']>0]
    genes=lambda rr:set(g for r in rr for g in tokens(r['owner_gene_ids']))
    for compartment in (['ER','mitochondria','cytosol'] if analysis=='PRIMARY_PLUS75' else ['ER','mitochondria']):
     summary.append({'analysis':analysis,'compartment':compartment,'topology_class':top,'position_class':cl,'side':side,'window_component_records':len(sr),'center_specific_pairs':len({r['paired_real_unit_id'] for r in sr}),'unique_TMD_units':len({r['unit_key'] for r in sr}),'owner_genes':len(genes(sr)),'full_overlap_window_components':len(full),'counted_overlap_window_components':len(u),'affected_center_specific_pairs':len({r['paired_real_unit_id'] for r in u}),'affected_unique_TMD_units':len({r['unit_key'] for r in u}),'affected_owner_genes':len(genes(u)),'counted_overlap_base_instances':sum(r['unmasked_other_gene_overlap_nt'] for r in sr),'unmasked_base_instances':sum(r['unmasked_nt'] for r in sr),'continuous_target_path_failures':sum(r['continuous_target_cds_stop_path_valid']=='NO' for r in sr)})
write(OUT/'overlap_summary_by_analysis_compartment_class_side.tsv',summary)
genesummary=[]
for analysis in ['PRIMARY_PLUS75','PLATEAU_90_120','BROAD_60_120']:
 gr=defaultdict(list)
 for r in inventory:
  if r['analysis']==analysis:
   for g in tokens(r['owner_gene_ids']):gr[g].append(r)
 for g,rows in sorted(gr.items()):
  affected=[r for r in rows if r['unmasked_other_gene_overlap_nt']>0]
  if affected:
   genesummary.append({'analysis':analysis,'owner_gene_id':g,'owner_gene_name':genemeta[g]['name'],'center_specific_pairs':len({r['paired_real_unit_id'] for r in rows}),'affected_center_specific_pairs':len({r['paired_real_unit_id'] for r in affected}),'affected_unique_TMD_units':len({r['unit_key'] for r in affected}),'affected_real_components':sum(r['side']=='real' for r in affected),'affected_pseudo_components':sum(r['side']=='pseudo' for r in affected),'overlapping_other_gene_ids':','.join(sorted(set(z for r in affected for z in tokens(r['other_gene_ids_unmasked'])))),'readthrough_like_target_transcripts':','.join(sorted(set(z for r in affected for z in tokens(r['structural_readthrough_like_compatible_transcripts']))))})
write(OUT/'affected_gene_summary.tsv',genesummary)
write(OUT/'input_sha256.tsv',[{'path':p,'sha256':v} for p,v in INPUTS.items()])
meta={'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.time()-START,'selected_pair_rows':len(pairs),'inventory_rows':len(inventory),'target_transcripts':len(txs),'invalid_target_path_window_records':sum(r['continuous_target_cds_stop_path_valid']=='NO' for r in inventory),'gtf_sha256':INPUTS[str(GTF)],'distinct_inspection_calls':inspect.cache_info()._asdict(),'pair_summary':{}}
for a,key in [('PRIMARY_PLUS75','primary_included'),('PLATEAU_90_120','plateau_included'),('BROAD_60_120','broad_included')]:
 rr=[r for r in flags if r[key]=='YES'];aff=[r for r in rr if r['exclude_pair_unmasked_other_gene']=='YES']
 meta['pair_summary'][a]={'pairs':len(rr),'units':len({r['unit_key'] for r in rr}),'genes':len(set(g for r in rr for g in tokens(r['owner_gene_ids']))),'affected_pairs':len(aff),'affected_units':len({r['unit_key'] for r in aff}),'affected_genes':len(set(g for r in aff for g in tokens(r['owner_gene_ids'])))}
(OUT/'inventory_metadata.json').write_text(json.dumps(meta,indent=2)+'\n');log(json.dumps(meta,indent=2))
