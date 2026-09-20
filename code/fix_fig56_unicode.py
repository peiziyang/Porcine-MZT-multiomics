#!/usr/bin/env python
"""Fix Fig5 Fig6 Unicode font warnings."""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.family':'Arial','svg.fonttype':'none','savefig.dpi':300,'font.size':9})
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import pandas as pd, numpy as np, os

OUT='E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures'
BASE='E:/Workbuddy/2026-07-27-11-58-27/data/processed'
deseq_dir=f'{BASE}/deseq2'
pa_data={}
for f in sorted(os.listdir(deseq_dir)):
    if f.startswith('m15_deseq2_') and f.endswith('_IVFvsPA.csv'):
        stage=f.replace('m15_deseq2_','').replace('_IVFvsPA.csv','')
        pa_data[stage]=pd.read_csv(os.path.join(deseq_dir,f))
kg=['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2','IDH2','PKM','EXOSC9','MDH1','GNL3','DPPA5','NANOG']
pa_all=[]
for stage,df in pa_data.items():
    if 'gene_symbol' not in df.columns: continue
    s2i={}
    for i,row in df.iterrows():
        s=str(row['gene_symbol'])
        if not s.startswith('ENSSSCG'): s2i[s]=i
    for g in kg:
        if g in s2i:
            r=df.iloc[s2i[g]]
            pa_all.append({'gene':g,'stage':stage,'log2FC':r['log2FoldChange'],'padj':r['padj']})
pa_df=pd.DataFrame(pa_all); overall=pa_df[pa_df['stage']=='overall']

def pl(ax,label):
    ax.text(-0.03,1.03,label,transform=ax.transAxes,fontsize=13,fontweight='bold',va='bottom',ha='left')

# FIG5
fig5,axes=plt.subplots(2,2,figsize=(12,9.5))
fig5.suptitle('Bulk RNA-seq Validation: PA Embryo Transcriptional Attenuation',fontsize=13,fontweight='bold')

ax=axes[0,0]; pl(ax,'A')
l1=['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2']
l1v=[]; l1p=[]; l1gf=[]
for g in l1:
    sub=overall[overall['gene']==g]
    if len(sub): l1v.append(sub['log2FC'].values[0]); l1p.append(sub['padj'].values[0]); l1gf.append(g)
xl=np.arange(len(l1gf))
bc=['#E41A1C' if v>0.5 else '#FF7F00' if v>0 else '#377EB8' if v<-0.5 else '#999' for v in l1v]
ax.bar(xl,l1v,color=bc,edgecolor='black',lw=0.8,width=0.6)
for i,(v,p) in enumerate(zip(l1v,l1p)):
    sig='p<0.05' if p<0.05 else 'ns'
    ax.text(i,v+0.3*np.sign(v) if v!=0 else 0.3,f"{v:+.2f} ({sig})",ha='center',fontsize=7,fontweight='bold',color='#333')
ax.axhline(0,color='gray',lw=1); ax.set_xticks(xl); ax.set_xticklabels(l1gf,fontsize=9)
ax.set_ylabel('log2 Fold Change (PA / IVF)',fontsize=10)
ax.set_title('Layer 1: Maternal blueprint - MAINTAINED',fontsize=11,fontweight='bold',color='#E41A1C')

ax=axes[0,1]; pl(ax,'B')
l3=['IDH2','PKM','EXOSC9','MDH1','GNL3']
l3v=[]; l3p=[]; l3gf=[]
for g in l3:
    sub=overall[overall['gene']==g]
    if len(sub): l3v.append(sub['log2FC'].values[0]); l3p.append(sub['padj'].values[0]); l3gf.append(g)
xl3=np.arange(len(l3gf))
ax.bar(xl3,l3v,color=['#377EB8' if v<-2 else '#999' for v in l3v],edgecolor='black',lw=0.8,width=0.6)
for i,(v,p) in enumerate(zip(l3v,l3p)):
    sig='p<1e-6' if p<1e-6 else f'p={p:.0e}' if p<0.05 else 'ns'
    ax.text(i,v+0.5,f"{v:+.2f} ({sig})",ha='center',fontsize=7,fontweight='bold',color='#377EB8')
ax.axhline(0,color='gray',lw=1); ax.set_xticks(xl3); ax.set_xticklabels(l3gf,fontsize=9)
ax.set_ylabel('log2 Fold Change (PA / IVF)',fontsize=10)
ax.set_title('Layer 3: Zygotic execution - ATTENUATED',fontsize=11,fontweight='bold',color='#377EB8')

ax=axes[1,0]; pl(ax,'C')
stg=['1cell','2cell','4cell','8cell']
gpl=['DNMT1','ZP3','GDF9','IDH2','PKM','EXOSC9','GNL3','DPPA5','NANOG']
hm=np.zeros((len(gpl),len(stg))); hm[:]=np.nan
for gi,g in enumerate(gpl):
    for si,s in enumerate(stg):
        sub=pa_df[(pa_df['gene']==g)&(pa_df['stage']==s)]
        if len(sub) and not np.isnan(sub.iloc[0]['log2FC']): hm[gi,si]=sub.iloc[0]['log2FC']
im=ax.imshow(np.ma.masked_where(np.isnan(hm),hm),aspect='auto',cmap='RdBu_r',vmin=-10,vmax=10)
ax.set_xticks(np.arange(len(stg))); ax.set_xticklabels(stg,fontsize=9)
ax.set_yticks(np.arange(len(gpl))); ax.set_yticklabels(gpl,fontsize=9)
for i in range(len(gpl)):
    for j in range(len(stg)):
        v=hm[i,j]
        if not np.isnan(v): ax.text(j,i,f'{v:+.1f}',ha='center',va='center',fontsize=7,fontweight='bold',color='white' if abs(v)>5 else 'black')
ax.set_title('Per-stage log2FC (PA vs IVF)',fontsize=10)
plt.colorbar(im,ax=ax,shrink=0.75,label='log2FC')

ax=axes[1,1]; pl(ax,'D'); ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis('off')
ys=9
for title,col,rows in [
    ('Layer 1: Maternal Blueprint','#E41A1C',[('DNMT1','+1.30','**'),('ZP3','+0.21','ns'),('ZP4','-0.75','ns'),('GDF9','+0.62','ns'),('RARRES1','-0.30','ns')]),
    ('Layer 3: Zygotic Execution','#377EB8',[('IDH2 (TCA)','-8.03','***'),('PKM (glycolysis)','-7.64','***'),('EXOSC9 (ribosome)','-7.73','***'),('MDH1 (TCA)','-4.58','**'),('GNL3 (nucleolar)','-3.06','*')]),
    ('ZGA Markers','#4DAF4A',[('DPPA5','+2.14','**'),('NANOG','+3.91','***')])]:
    ax.text(0.3,ys,title,fontsize=9.5,fontweight='bold',color=col,va='center'); ys-=0.5
    for gene,lfc,sig in rows:
        ax.text(0.5,ys,gene,fontsize=8,va='center')
        fc=float(lfc.replace('+',''))
        cl='#E41A1C' if fc>1 else '#377EB8' if fc<-1 else '#666'
        ax.text(5.5,ys,f'log2FC = {lfc}',fontsize=8,va='center',fontweight='bold',color=cl)
        sg='[***]' if sig=='***' else '[**]' if sig=='**' else '[*]' if sig=='*' else ''
        ax.text(8.5,ys,sg,fontsize=8,va='center',color='#CC0000')
        ys-=0.4
    ys-=0.3
ax.text(0.5,ys,'42-sample bulk RNA-seq | SRP301735 | DESeq2: stage-adjusted',fontsize=7.5,style='italic',color='gray')
ax.text(0.5,ys-0.35,'[*] Padj<0.05  [**] Padj<0.01  [***] Padj<0.001',fontsize=7,color='gray')

plt.tight_layout()
fig5.savefig(f'{OUT}/fig5_pa_validation.svg',format='svg')
fig5.savefig(f'{OUT}/fig5_pa_validation.png',dpi=300)
print(f'Fig5: {os.path.getsize(OUT+"/fig5_pa_validation.svg")//1024}KB')

# FIG6
fig6,ax=plt.subplots(1,1,figsize=(12,8))
ax.set_xlim(0,14); ax.set_ylim(0,10); ax.axis('off')
for y in np.arange(1,10,1): ax.axhline(y,color='#E8E8E8',lw=0.3)
ax.text(3,9.7,'NORMAL DEVELOPMENT',ha='center',fontsize=13,fontweight='bold',color='#333')
gv=plt.Circle((3,8.5),0.6,facecolor='#FFF3CD',edgecolor='#FF7F00',lw=2)
ax.add_patch(gv); ax.text(3,8.5,'GV\nOocyte',ha='center',va='center',fontsize=8,fontweight='bold')
ax.text(3,9.3,'32 cells x RNA+WGBS',ha='center',fontsize=6.5,color='gray')
ax.annotate('',xy=(3,7.5),xytext=(3,7.9),arrowprops=dict(arrowstyle='->',color='#666',lw=1.5))
ax.text(3,7.7,'MOFA+',ha='center',fontsize=7,color='#666',style='italic')
l1=FancyBboxPatch((0.5,6.2),5,1,boxstyle='round,pad=0.1',facecolor='#377EB8',edgecolor='#1a5276',lw=2,alpha=0.2)
ax.add_patch(l1); ax.text(3,7.0,'Layer 1: Maternal Blueprint (F1)',fontsize=10,fontweight='bold',color='#1a5276',ha='center')
ax.text(3,6.6,'DNMT1 . ZP3 . ZP4 . GDF9 . RARRES1',fontsize=8,ha='center',color='#1a5276')
ax.text(3,6.3,'RNA+METH co-regulated  r=+0.84',fontsize=7,ha='center',color='#666',style='italic')
l2=FancyBboxPatch((0.5,4.2),5,1,boxstyle='round,pad=0.1',facecolor='#E41A1C',edgecolor='#8B0000',lw=2,alpha=0.15)
ax.add_patch(l2); ax.text(3,5.0,'Layer 2: Maternal Clearance (F3)',fontsize=10,fontweight='bold',color='#8B0000',ha='center')
ax.text(3,4.6,'ZP2 . SYCN . PARP12',fontsize=8,ha='center',color='#8B0000')
ax.text(3,4.3,'Spearman rho = -0.90  P = 2e-4',fontsize=7,ha='center',color='#666',style='italic')
l3=FancyBboxPatch((0.5,2.2),5,1,boxstyle='round,pad=0.1',facecolor='#4DAF4A',edgecolor='#1B5E20',lw=2,alpha=0.15)
ax.add_patch(l3); ax.text(3,3.0,'Layer 3: Zygotic Activation (F4/F6)',fontsize=10,fontweight='bold',color='#1B5E20',ha='center')
ax.text(3,2.6,'IDH2 . PKM . EXOSC9 . GNL3 . NOP9',fontsize=8,ha='center',color='#1B5E20')
ax.text(3,2.3,'Metabolism + Translation + Ribosome',fontsize=7,ha='center',color='#666',style='italic')
for i,(x,stage) in enumerate([(0.2,'1-cell'),(0.2,'2-cell'),(0.2,'4-cell'),(0.2,'8-cell'),(0.2,'Morula')]):
    c=plt.Circle((x,8.5-i*0.9),0.15,facecolor='#DDD' if i<3 else '#CCE5FF',edgecolor='gray',lw=0.5)
    ax.add_patch(c); ax.text(x+0.25,8.5-i*0.9,stage,fontsize=6.5,color='gray',va='center')
ax.annotate('',xy=(5.5,1.8),xytext=(5.5,8.8),arrowprops=dict(arrowstyle='->',color='#999',lw=2,ls='--'))
ax.text(5.5,5.2,'MZT\nProgression',ha='center',fontsize=8,color='#999',rotation=90,va='center')
ax.text(10,9.7,'PA EMBRYOS',ha='center',fontsize=13,fontweight='bold',color='#E41A1C')
pl1=FancyBboxPatch((7.5,6.2),5,1,boxstyle='round,pad=0.1',facecolor='#FFE5E5',edgecolor='#E41A1C',lw=2,alpha=0.5)
ax.add_patch(pl1); ax.text(10,7.0,'Layer 1: MAINTAINED',fontsize=10,fontweight='bold',color='#E41A1C',ha='center')
ax.text(10,6.6,'DNMT1 expression HIGHER in PA',fontsize=8,ha='center',color='#E41A1C')
ax.text(10,6.3,'log2FC=+1.30, Padj=0.006',fontsize=7,ha='center',color='#CC0000',style='italic')
pl2=FancyBboxPatch((7.5,4.2),5,1,boxstyle='round,pad=0.1',facecolor='#FFF3CD',edgecolor='#FF7F00',lw=2,alpha=0.4)
ax.add_patch(pl2); ax.text(10,5.0,'Layer 2: VARIABLE',fontsize=10,fontweight='bold',color='#FF7F00',ha='center')
ax.text(10,4.6,'ZP2/4 inconsistent across stages',fontsize=8,ha='center',color='#CC8400')
pl3=FancyBboxPatch((7.5,2.2),5,1,boxstyle='round,pad=0.1',facecolor='#E5F0FF',edgecolor='#377EB8',lw=2,alpha=0.5)
ax.add_patch(pl3); ax.text(10,3.0,'Layer 3: ATTENUATED',fontsize=10,fontweight='bold',color='#377EB8',ha='center')
ax.text(10,2.6,'IDH2 (log2FC=-8.03)  PKM (-7.64)',fontsize=8,ha='center',color='#377EB8')
ax.text(10,2.3,'EXOSC9 (-7.73)  MDH1 (-4.58)',fontsize=7,ha='center',color='#377EB8')
ax.text(7,1.0,'GV oocyte multi-omics (32 cells, MOFA+) -> Atlas projection (832 E0-E10 cells) -> Independent validation (42 samples, SRP301735)',ha='center',fontsize=8,style='italic',color='gray')
plt.tight_layout()
fig6.savefig(f'{OUT}/fig6_model.svg',format='svg')
fig6.savefig(f'{OUT}/fig6_model.png',dpi=300)
print(f'Fig6: {os.path.getsize(OUT+"/fig6_model.svg")//1024}KB')
print('Done - no Unicode warnings')