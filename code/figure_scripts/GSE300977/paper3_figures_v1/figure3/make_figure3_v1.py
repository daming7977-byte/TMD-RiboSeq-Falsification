"""Figure 3: render frozen tables only; no estimation or biological analysis."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SECONDARY=Path('tmd_spatial_secondary_scan_v1/output/SECONDARY_11CENTER_LOCATION_PROFILE_v1.tsv')
COMMON=Path('tmd_spatial_common_anchor_sensitivity_v1/output/COMMON_ANCHOR_11CENTER_PROFILE_v1.tsv')
OUT=Path('paper3_figures_v1/figure3'); OUT.mkdir(parents=True,exist_ok=True)
s=pd.read_csv(SECONDARY,sep='\t'); c=pd.read_csv(COMMON,sep='\t')
CENTERS=[-30,-15,0,15,30,45,60,75,90,105,120]
LOCATIONS=[('ER','ER','#009E73','o'),('mitochondria','Mitochondria','#D55E00','s'),('cytosol','Cytosol','#0072B2','^')]
for df,lo,hi,ncol in [(s,'pointwise_bootstrap_ci_2.5','pointwise_bootstrap_ci_97.5','N_units'),(c,'CI_2.5','CI_97.5','N_common_anchors')]:
    assert {'center','location','D_mean','D_rep1','D_rep2',lo,hi,ncol} <= set(df.columns)
    assert len(df)==33 and not df.duplicated(['center','location']).any()
    assert set(df.location)=={x[0] for x in LOCATIONS}
    for location,_,_,_ in LOCATIONS:
        z=df[df.location.eq(location)].sort_values('center')
        assert z.center.tolist()==CENTERS
        assert np.isfinite(z[['D_mean','D_rep1','D_rep2',lo,hi]].to_numpy()).all()
        assert (z[lo]<=z.D_mean).all() and (z.D_mean<=z[hi]).all()
assert c.N_common_anchors.eq(3510).all(), 'COMMON_ANCHOR_UNIVERSE_CONFLICT'
print('FIGURE3_SCHEMA_ACCOUNTING_AND_N3510 = PASS')
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':10,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','savefig.transparent':False})
fig,axs=plt.subplots(2,2,figsize=(7.3,6.7))
fig.subplots_adjust(left=.095,right=.975,bottom=.20,top=.865,hspace=.54,wspace=.29)
fig.text(.095,.97,'GSE300977 | Positional structure across spatial states',fontsize=12,weight='bold',va='top')
fig.text(.095,.932,'A-C  Preregistered secondary scan     D  Post-reveal common-anchor sensitivity',fontsize=8,va='top')
for ax,letter in zip(axs.flat,'ABCD'):
    ax.text(-.16,1.14,letter,transform=ax.transAxes,fontsize=12,weight='bold')
    ax.spines[['top','right']].set_visible(False)
    ax.axhline(0,color='.63',lw=.8,zorder=0)
    ax.axvline(0,color='.86',lw=.7,zorder=0)
    ax.axvline(75,color='.53',ls=':',lw=.8,zorder=0)
    ax.set_xlim(-36,126); ax.set_ylim(-.62,.57)
    ax.set_xticks(CENTERS,labels=[f'{x:+d}' if x else '0' for x in CENTERS],rotation=45,ha='right')
    ax.set_xlabel('Position relative to TMD start (nt)')
    ax.set_ylabel(r'Mean $D$ (log$_2$)')
for ax,(location,label,color,marker) in zip(axs.flat,LOCATIONS):
    z=s[s.location.eq(location)].sort_values('center'); x=z.center.to_numpy(); y=z.D_mean.to_numpy()
    ax.fill_between(x,z['pointwise_bootstrap_ci_2.5'].to_numpy(),z['pointwise_bootstrap_ci_97.5'].to_numpy(),color=color,alpha=.15,lw=0)
    ax.plot(x,y,color=color,marker=marker,ms=3.5,lw=1.5)
    ax.scatter(x,z.D_rep1,s=10,color=color,alpha=.4,zorder=3)
    ax.scatter(x,z.D_rep2,s=10,facecolors='white',edgecolors=color,linewidths=.65,zorder=3)
    ax.set_title(label,loc='left',color=color,pad=10,weight='bold')
ax=axs[1,1]
for location,label,color,marker in LOCATIONS:
    z=c[c.location.eq(location)].sort_values('center'); x=z.center.to_numpy(); y=z.D_mean.to_numpy()
    ax.fill_between(x,z['CI_2.5'].to_numpy(),z['CI_97.5'].to_numpy(),color=color,alpha=.09,lw=0)
    ax.plot(x,y,color=color,marker=marker,ms=3,lw=1.3,label=label)
ax.set_title('Fixed common anchors | N = 3,510',loc='left',pad=10,weight='bold')
ax.legend(frameon=False,loc='lower right',fontsize=7)
handles=[Line2D([],[],color='.3',marker='o',ms=3,label='Mean'),Line2D([],[],ls='',color='.5',marker='o',ms=3,label='Replicate 1 (A-C)'),Line2D([],[],ls='',color='.4',marker='o',mfc='white',ms=3,label='Replicate 2 (A-C)')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.073),ncol=3,frameon=False,fontsize=7)
fig.text(.5,.055,'Shading: pointwise 95% bootstrap intervals; not multiplicity-corrected.',ha='center',fontsize=7)
fig.text(.5,.033,'Dotted line: frozen +75-nt primary position. No secondary center is promoted to primary.',ha='center',fontsize=7)
fig.text(.5,.011,'A-C use center-specific eligible units; D uses the same 3,510 anchors at all 11 centers.',ha='center',fontsize=7,color='.3')
for ext in ['pdf','svg','png']:
    p=OUT/f'FIGURE3_POSITIONAL_COMMON_ANCHOR_v1.{ext}'
    fig.savefig(p,dpi=450,facecolor='white',edgecolor='white',transparent=False)
    print('SHA256',hashlib.sha256(p.read_bytes()).hexdigest(),p)
plt.close(fig)
print('FIGURE3_GENERATION_GATE = PASS')
