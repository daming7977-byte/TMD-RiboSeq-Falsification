"""Presentation-only export of frozen TSV results. No inference or counting.
Run with --input-root (frozen supplement directory) and --output-root (new directory).
"""
import argparse,csv,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from PIL import Image

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True);args=p.parse_args()
IN=args.input_root;OUT=args.output_root;OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':10,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','savefig.transparent':False,'figure.facecolor':'white'})
COL=['#0072B2','#B45F06','#16856B'];LOCS=['ER','mitochondria','cytosol'];NAMES=['ER','Mitochondria','Cytosol'];CENTERS=[-30,-15,0,15,30,45,60,75,90,105,120];mapping=[];files=[]
oldmanifest=list(csv.DictReader((IN/'SUPPLEMENTARY_FILE_MANIFEST.tsv').open(),delimiter='\t'));orig={r['packaged_file']:r for r in oldmanifest}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name):
 path=IN/name
 if name in orig:assert digest(path)==orig[name]['packaged_sha256'],name
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def map_panel(n,panel,name,cols,selector,tier):
 r=orig[name];mapping.append([n,panel,name,r['source_frozen_file'],cols,selector,tier,'EXISTING_FROZEN_VALUES_NO_RECALCULATION',r['source_sha256'],digest(IN/name)])
def label(ax,letter,title):ax.set_title(title,loc='left',pad=8);ax.figure.text(ax.get_position().x0-.040,ax.get_position().y1+.012,letter,fontweight='bold',fontsize=12)
def axis(ax,ylab):
 ax.axhline(0,color='#777777',lw=.6);ax.axvline(75,color='#555555',lw=.75,ls=':');ax.set_xticks([-30,0,30,60,90,120]);ax.set_xlabel('Window center relative to anchor (nt)');ax.set_ylabel(ylab);ax.grid(axis='y',alpha=.15);ax.set_xlim(-34,124)
def plot(ax,rows,col,mean='D_mean',r1='D_rep1',r2='D_rep2',ylab='Real-minus-pseudo D'):
 rows=sorted(rows,key=lambda r:int(r['center']));assert [int(r['center']) for r in rows]==CENTERS
 x=[int(r['center']) for r in rows]
 ax.fill_between(x,[float(r['CI_2.5']) for r in rows],[float(r['CI_97.5']) for r in rows],color=col,alpha=.16,lw=0)
 ax.plot(x,[float(r[mean]) for r in rows],'-',color=col,lw=1.45)
 ax.plot(x,[float(r[r1]) for r in rows],'o',color=col,ms=3.5,lw=.6)
 ax.plot(x,[float(r[r2]) for r in rows],'s',mfc='white',mec=col,ms=3.5,mew=.8)
 axis(ax,ylab)
def start(n,title,tier,size):
 f=plt.figure(figsize=size);f.suptitle(f'Supplementary Figure {n} | {title}',x=.055,y=.978,ha='left',fontsize=13,fontweight='bold');f.text(.055,.936,tier,fontsize=9,color='#444444');return f
def legend(f):
 f.legend([Line2D([],[],color='#333333',lw=1.5),Line2D([],[],color='#333333',marker='o',ls='none',ms=4),Line2D([],[],color='#333333',marker='s',mfc='white',ls='none',ms=4)],['Mean estimate','Biological replicate 1','Biological replicate 2'],loc='lower center',bbox_to_anchor=(.53,.062),ncol=3,frameon=False,fontsize=8)
def save(f,n,point=True):
 if point:f.text(.055,.021,'Shading: pointwise 95% bootstrap intervals, not multiplicity-corrected across positions.\nDotted line: +75 nt, the sole primary center. Unit-bootstrap uncertainty does not replace biological replication.',fontsize=7.3,color='#444444')
 else:f.text(.055,.023,'Recorded technical QC and eligibility fields only. No new statistical tests, counts from alignments or exclusions.',fontsize=7.5,color='#444444')
 stem=OUT/f'Supplementary_Figure_{n}_v1'
 for ext in ['pdf','svg','png']:
  path=stem.with_suffix('.'+ext);f.savefig(path,dpi=320,facecolor='white',transparent=False);files.append([n,path.name,digest(path)])
 im=Image.open(stem.with_suffix('.png')).convert('RGBA');assert im.getchannel('A').getextrema()==(255,255);assert im.getpixel((0,0))==(255,255,255,255)
 plt.close(f)
# Fig S1: complete fixed-common-anchor replicate profiles and explicitly recorded common N.
name='Table_3/COMMON_ANCHOR_11CENTER_PROFILE_v1.tsv';d=load(name);f=start(1,'Fixed-common-anchor sensitivity','Post-reveal sensitivity | Same fixed anchor population across all 11 centers',(10,6.4));g=f.add_gridspec(2,3,left=.075,right=.97,bottom=.18,top=.85,hspace=.7,height_ratios=[2.2,1])
for j,loc in enumerate(LOCS):
 a=f.add_subplot(g[0,j]);plot(a,[r for r in d if r['location']==loc],COL[j]);a.set_ylim(-.7,.75);label(a,'ABC'[j],NAMES[j]);map_panel(1,'ABC'[j],name,'center; D_mean; D_rep1; D_rep2; CI_2.5; CI_97.5','location='+loc,'POST_REVEAL_SENSITIVITY')
a=f.add_subplot(g[1,:]);r=[r for r in d if r['location']=='ER'];a.plot([int(x['center']) for x in r],[int(x['N_common_anchors']) for x in r],'o-',color='#555555',ms=3);a.set_ylim(0,4000);a.set_yticks([0,3510]);axis(a,'Common anchors');label(a,'D','Fixed common-anchor population (N = 3,510)');map_panel(1,'D',name,'center; N_common_anchors','location=ER (same recorded N across compartments)','POST_REVEAL_SENSITIVITY');legend(f);save(f,1)
# S2/S3: complete class profiles, replicate-specific contrasts and recorded eligibility.
for n,multi in [(2,False),(3,True)]:
 folder='Table_5/' if multi else 'Table_4/';prof=folder+('MULTIPASS_FIRST_INTERNAL_11CENTER_PROFILE_v1.tsv' if multi else 'FIRST_INTERNAL_11CENTER_PROFILE_v1.tsv');con=folder+('MULTIPASS_FIRST_MINUS_INTERNAL_11CENTER_v1.tsv' if multi else 'FIRST_MINUS_INTERNAL_11CENTER_CONTRAST_v1.tsv');acct=folder+('MULTIPASS_TOPOLOGY_ACCOUNTING_v1.tsv' if multi else 'FIRST_INTERNAL_ELIGIBILITY_ACCOUNTING_v1.tsv');d=load(prof);c=load(con);acc=load(acct);tier='POST_REVEAL_SENSITIVITY' if multi else 'PREREGISTERED_SECONDARY'
 f=start(n,'Multi-pass-only TMD-order profiles' if multi else 'Complete FIRST versus INTERNAL profiles','Post-reveal multi-pass-only sensitivity' if multi else 'Preregistered secondary stratification | Original FIRST and INTERNAL populations',(10,11));g=f.add_gridspec(4,3,left=.085,right=.97,bottom=.155,top=.865,hspace=.76,height_ratios=[1.5,1.5,1.5,1.05])
 for row,cl in enumerate(['FIRST','INTERNAL']):
  for j,loc in enumerate(LOCS):
   a=f.add_subplot(g[row,j]);r=[x for x in d if x['location']==loc and x['TMD_class']==cl];plot(a,r,COL[j]);a.set_ylim(-1.1,1.1);letter='ABCDEFGHIJ'[row*3+j];label(a,letter,NAMES[j]+' | '+('First TMD' if cl=='FIRST' else 'Internal TMD'));map_panel(n,letter,prof,'center; D_mean; D_rep1; D_rep2; CI_2.5; CI_97.5',f'location={loc}; TMD_class={cl}',tier)
 for j,loc in enumerate(LOCS):
  a=f.add_subplot(g[2,j]);plot(a,[r for r in c if r['location']==loc],COL[j],mean='FIRST_minus_INTERNAL_mean',r1='FIRST_minus_INTERNAL_rep1',r2='FIRST_minus_INTERNAL_rep2',ylab='First - internal D');a.set_ylim(-1.35,1.2);letter='GHI'[j];label(a,letter,NAMES[j]+' | First minus internal');map_panel(n,letter,con,'center; FIRST_minus_INTERNAL_mean; FIRST_minus_INTERNAL_rep1; FIRST_minus_INTERNAL_rep2; CI_2.5; CI_97.5','location='+loc,tier)
 a=f.add_subplot(g[3,:]);fields=['MULTI_FIRST','MULTI_INTERNAL','MULTI_MIXED_POSITION'] if multi else ['FIRST','INTERNAL','MIXED']
 for col,field,lab in zip(['#0072B2','#D55E00','#777777'],fields,['First','Internal','Mixed (excluded from class contrast)']):a.plot([int(r['center']) for r in acc],[int(r[field]) for r in acc],'-o',color=col,ms=3,label=lab)
 axis(a,'Eligible units');a.set_ylim(bottom=0);label(a,'J','Recorded class accounting');a.legend(loc='upper center',bbox_to_anchor=(.59,1.55),ncol=3,frameon=False,fontsize=7);map_panel(n,'J',acct,'center; '+'; '.join(fields),'all centers; saved class counts',tier);legend(f);save(f,n)
# S4: exactly the two non-focal direct interactions already specified by manifest.
name='Table_6/MULTIPASS_TOPOLOGY_COMPARTMENT_INTERACTIONS_v1.tsv';d=load(name);acct='Table_6/TOPOLOGY_COMPARTMENT_INTERACTION_ACCOUNTING_v1.tsv';acc=load(acct);f=start(4,'Additional topology-by-compartment interactions','Post-reveal direct follow-up | Multi-pass proteins',(10,6.7));g=f.add_gridspec(2,2,left=.085,right=.97,bottom=.20,top=.84,hspace=.68,height_ratios=[2,1])
for j,(key,title) in enumerate([('CYTOSOL_MINUS_ER','Cytosolic versus ER'),('MITO_MINUS_CYTOSOL','Mitochondrial versus cytosolic')]):
 a=f.add_subplot(g[0,j]);plot(a,[r for r in d if r['interaction']==key],COL[j],mean='mean',r1='rep1',r2='rep2',ylab='Difference in first-minus-internal effects');a.set_ylim(-1.4,2);label(a,'AB'[j],title);map_panel(4,'AB'[j],name,'center; rep1; rep2; mean; CI_2.5; CI_97.5','interaction='+key,'POST_REVEAL_FOLLOW_UP')
a=f.add_subplot(g[1,:])
for cl,col in [('FIRST',COL[0]),('INTERNAL',COL[1])]:
 r=[x for x in acc if x['TMD_class']==cl];a.plot([int(x['center']) for x in r],[int(x['N_units']) for x in r],'-o',color=col,ms=3,label='First TMD' if cl=='FIRST' else 'Internal TMD')
axis(a,'Units per class');a.set_ylim(bottom=0);a.legend(frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.72,1.45));label(a,'C','Recorded multi-pass class populations');map_panel(4,'C',acct,'center; TMD_class; N_units','all centers and classes','POST_REVEAL_FOLLOW_UP');legend(f);save(f,4)
# S5: QC flow and recorded counts, without recounting or aggregation.
f=start(5,'Engineering and technical QC','Blinded technical processing and frozen window eligibility',(10,8.2));g=f.add_gridspec(3,2,left=.085,right=.97,bottom=.10,top=.86,hspace=.70,height_ratios=[.65,1.4,1.4])
a=f.add_subplot(g[0,:]);a.set_axis_off();a.set_xlim(-.025,1.02);a.text(-.06,1.06,'A',transform=a.transAxes,fontsize=12,fontweight='bold')
steps=['Blinded libraries\nB01-B12','Terminal read\npreprocessing','Unique alignment\nand length QC','Frozen geometry\nand paired masks']
for j,t in enumerate(steps):
 x=j*.255;a.add_patch(FancyBboxPatch((x,.1),.215,.70,boxstyle='round,pad=.008',fc='#EDF2F6',ec='#82919F',lw=.7));a.text(x+.1075,.45,t,ha='center',va='center',fontsize=9)
 if j<3:a.annotate('',xy=(x+.247,.45),xytext=(x+.22,.45),arrowprops={'arrowstyle':'->','color':'#65717D'})
map_panel(5,'A','Table_8/BLINDED_RUN_MANIFEST.tsv','blind_id; run_accession','all rows; library IDs only','ENGINEERING_QC');map_panel(5,'A','Table_8/BLIND_MATRIX_QC_FINAL_AUDIT_v2p1r1.tsv','recorded_field; recorded_value','recorded QC PASS stages; qualitative flow','ENGINEERING_QC')
pre='Table_8/FORMAL_PREPROCESSING_SUMMARY_v3.tsv';leng='Table_8/BLIND_READ_LENGTH_SUMMARY.tsv';d=load(pre);l=load(leng);ids=[r['blind_id'] for r in d];lookup={r['blind_id']:r for r in l};x=list(range(len(ids)))
a=f.add_subplot(g[1,:]);
for delta,col,key,title in [(-.25,'#AAB5BF','raw_reads','Raw reads'),(0,'#457B9D','passing_reads','After preprocessing'),(.25,'#16856B','unique_reads','Unique alignments')]:
 vals=[int(lookup[r['blind_id']][key] if key=='unique_reads' else r[key]) for r in d];a.bar([i+delta for i in x],vals,width=.23,color=col,label=title)
a.set_xticks(x,ids);a.set_ylabel('Recorded reads / alignments');a.ticklabel_format(axis='y',style='sci',scilimits=(6,6));a.grid(axis='y',alpha=.15);label(a,'B','Recorded per-library processing counts');a.legend(frameon=False,ncol=3,loc='upper right');map_panel(5,'B',pre,'blind_id; raw_reads; passing_reads','all 12 libraries','ENGINEERING_QC');map_panel(5,'B',leng,'blind_id; unique_reads','all 12 libraries; matched by blind_id','ENGINEERING_QC')
a=f.add_subplot(g[2,0]);a.bar(x,[float(lookup[k]['fraction_27_36']) for k in ids],color='#16856B');a.set_xticks(x,ids,rotation=45);a.set_ylim(0,1);a.set_ylabel('Recorded fraction of unique reads');label(a,'C','Read lengths 27-36 nt');map_panel(5,'C',leng,'blind_id; fraction_27_36','all libraries','ENGINEERING_QC')
a=f.add_subplot(g[2,1]);geo='Table_8/GSE300977_TMD_GEOMETRY_RESOLVED_UNIT_SUMMARY_v2.tsv';pair='Table_8/GSE300977_REAL_PSEUDO_PAIRED_GEOMETRY_SUMMARY_v1.tsv'
for file,field,col,title in [(geo,'geometry_resolved_units','#457B9D','Structural real units'),(pair,'unique_real_geometry_units','#B45F06','Structurally paired real units')]:
 rr=load(file);a.plot([int(r['center']) for r in rr],[int(r[field]) for r in rr],'-o',color=col,ms=3,label=title);map_panel(5,'D',file,'center; '+field,'all centers; structural counts before final Crossmap eligibility','ENGINEERING_QC')
axis(a,'Recorded units');a.set_ylim(bottom=0);a.legend(frameon=False,fontsize=7,loc='lower right');label(a,'D','Structural eligibility before mask filtering');save(f,5,False)
with (OUT/'panel_sources.tsv').open('w') as h:
 w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['figure','panel','packaged_input','frozen_source','columns','selection','tier','value_mode','frozen_sha256','input_sha256']);w.writerows(mapping)
with (OUT/'artifact_hashes.tsv').open('w') as h:
 w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['figure','file','sha256']);w.writerows(files)
print(json.dumps({'figures':5,'artifacts':len(files),'panel_source_rows':len(mapping),'new_statistics':0}))
