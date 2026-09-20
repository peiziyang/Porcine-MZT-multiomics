#!/usr/bin/env python
"""Fix Fig5: A/B bars not visible (y-axis issue)."""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial', 'font.size': 9,
    'svg.fonttype': 'none', 'axes.unicode_minus': False,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd, numpy as np, os

OUT = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures"
BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
os.makedirs(OUT, exist_ok=True)

deseq_dir = f"{BASE}/deseq2"
pa_data = {}
for f in sorted(os.listdir(deseq_dir)):
    if f.startswith('m15_deseq2_') and f.endswith('_IVFvsPA.csv'):
        stage = f.replace('m15_deseq2_','').replace('_IVFvsPA.csv','')
        pa_data[stage] = pd.read_csv(os.path.join(deseq_dir, f))
kg = ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2','IDH2','PKM','EXOSC9','MDH1','GNL3','DPPA5','NANOG']
pa_all = []
for stage, df in pa_data.items():
    if 'gene_symbol' not in df.columns: continue
    s2i = {}
    for i, row in df.iterrows():
        s = str(row['gene_symbol'])
        if not s.startswith('ENSSSCG'): s2i[s] = i
    for g in kg:
        if g in s2i:
            r = df.iloc[s2i[g]]
            pa_all.append({'gene':g,'stage':stage,'log2FC':r['log2FoldChange'],'padj':r['padj']})
pa_df = pd.DataFrame(pa_all); overall = pa_df[pa_df['stage']=='overall']

fig5 = plt.figure(figsize=(15, 14))
gs5 = GridSpec(3, 2, figure=fig5, height_ratios=[1, 1, 1.15], hspace=0.4, wspace=0.3)

# A: Layer 1
ax = fig5.add_subplot(gs5[0, 0])
ax.text(-0.05, 1.04, 'A', transform=ax.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='left')
l1 = ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2']
l1v = []; l1p = []; l1gf = []
for g in l1:
    sub = overall[overall['gene']==g]
    if len(sub): l1v.append(sub['log2FC'].values[0]); l1p.append(sub['padj'].values[0]); l1gf.append(g)
xl = np.arange(len(l1gf))
bc = ['#E41A1C' if v>0.5 else '#FF7F00' if v>0 else '#377EB8' if v<-0.5 else '#999' for v in l1v]
bars1 = ax.bar(xl, l1v, color=bc, edgecolor='black', lw=0.8, width=0.6)
for i,(v,p) in enumerate(zip(l1v,l1p)):
    sig = 'p<0.05' if p<0.05 else 'ns'
    ax.text(i, v+0.15 if v>0 else v-0.4, f'{v:+.2f} ({sig})', ha='center', fontsize=7, fontweight='bold', color='#333')
ax.axhline(0, color='gray', lw=1)
ax.set_xticks(xl); ax.set_xticklabels(l1gf, fontsize=9)
ax.set_ylabel('log2FC (PA vs IVF)', fontsize=9)
ax.set_title('Layer 1: Maternal blueprint - MAINTAINED', fontsize=10, fontweight='bold', color='#E41A1C')
ax.set_ylim(-2, 2)
ax.grid(axis='y', alpha=0.3)

# B: Layer 3
ax = fig5.add_subplot(gs5[0, 1])
ax.text(-0.05, 1.04, 'B', transform=ax.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='left')
l3 = ['IDH2','PKM','EXOSC9','MDH1','GNL3']
l3v = []; l3p = []; l3gf = []
for g in l3:
    sub = overall[overall['gene']==g]
    if len(sub): l3v.append(sub['log2FC'].values[0]); l3p.append(sub['padj'].values[0]); l3gf.append(g)
xl3 = np.arange(len(l3gf))
bars3 = ax.bar(xl3, l3v, color=['#377EB8' if v<-2 else '#999' for v in l3v], edgecolor='black', lw=0.8, width=0.6)
for i,(v,p) in enumerate(zip(l3v,l3p)):
    sig = 'p<1e-6' if p<1e-6 else f'p={p:.0e}' if p<0.05 else 'ns'
    ax.text(i, v-0.5, f'{v:+.2f} ({sig})', ha='center', fontsize=7, fontweight='bold', color='#377EB8')
ax.axhline(0, color='gray', lw=1)
ax.set_xticks(xl3); ax.set_xticklabels(l3gf, fontsize=9)
ax.set_ylabel('log2FC (PA vs IVF)', fontsize=9)
ax.set_title('Layer 3: Zygotic execution - ATTENUATED', fontsize=10, fontweight='bold', color='#377EB8')
ax.set_ylim(-10, 1)
ax.grid(axis='y', alpha=0.3)

# C: Per-stage heatmap
ax = fig5.add_subplot(gs5[1, 0])
ax.text(-0.05, 1.04, 'C', transform=ax.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='left')
stg = ['1cell','2cell','4cell','8cell']
gpl = ['DNMT1','ZP3','GDF9','IDH2','PKM','EXOSC9','GNL3','DPPA5','NANOG']
hm = np.zeros((len(gpl),len(stg))); hm[:]=np.nan
for gi,g in enumerate(gpl):
    for si,s in enumerate(stg):
        sub = pa_df[(pa_df['gene']==g)&(pa_df['stage']==s)]
        if len(sub) and not np.isnan(sub.iloc[0]['log2FC']): hm[gi,si]=sub.iloc[0]['log2FC']
im = ax.imshow(np.ma.masked_where(np.isnan(hm),hm), aspect='auto', cmap='RdBu_r', vmin=-10, vmax=10)
ax.set_xticks(np.arange(len(stg))); ax.set_xticklabels(stg, fontsize=9)
ax.set_yticks(np.arange(len(gpl))); ax.set_yticklabels(gpl, fontsize=9)
for i in range(len(gpl)):
    for j in range(len(stg)):
        v = hm[i,j]
        if not np.isnan(v): ax.text(j,i,f'{v:+.1f}',ha='center',va='center',fontsize=7,fontweight='bold',color='white' if abs(v)>5 else 'black')
ax.set_title('Per-stage log2FC', fontsize=10)
plt.colorbar(im, ax=ax, shrink=0.75, label='log2FC')

# D: Summary table
ax = fig5.add_subplot(gs5[1, 1])
ax.text(-0.05, 1.04, 'D', transform=ax.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='left')
ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis('off')
ys = 9
for title, col, rows in [
    ('Layer 1: Blueprint','#E41A1C',[('DNMT1','+1.30','**'),('ZP3','+0.21','ns'),('ZP4','-0.75','ns'),('GDF9','+0.62','ns'),('RARRES1','-0.30','ns')]),
    ('Layer 3: Execution','#377EB8',[('IDH2 (TCA)','-8.03','***'),('PKM (glycolysis)','-7.64','***'),('EXOSC9 (ribosome)','-7.73','***'),('MDH1 (TCA)','-4.58','**'),('GNL3 (nucleolar)','-3.06','*')]),
    ('ZGA Markers','#4DAF4A',[('DPPA5','+2.14','**'),('NANOG','+3.91','***')])]:
    ax.text(0.3,ys,title,fontsize=10,fontweight='bold',color=col,va='center'); ys-=0.4
    for gene,lfc,sig in rows:
        ax.text(0.5,ys,gene,fontsize=8,va='center')
        fc = float(lfc.replace('+',''))
        cl = '#E41A1C' if fc>1 else '#377EB8' if fc<-1 else '#666'
        ax.text(4.5,ys,f'log2FC={lfc}',fontsize=8,va='center',fontweight='bold',color=cl)
        sg = '[***]' if sig=='***' else '[**]' if sig=='**' else '[*]' if sig=='*' else ''
        ax.text(7.5,ys,sg,fontsize=8,va='center',color='#CC0000')
        ys-=0.3
    ys-=0.2
ax.text(0.5, ys-0.2, '[*] Padj<0.05  [**] Padj<0.01  [***] Padj<0.001', fontsize=7, color='gray')

# E: Three-layer model
ax = fig5.add_subplot(gs5[2, :])
ax.text(-0.02, 1.04, 'E', transform=ax.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='left')
ax.set_xlim(0, 14); ax.set_ylim(0, 6); ax.axis('off')
ax.text(3, 5.5, 'NORMAL MZT', ha='center', fontsize=11, fontweight='bold', color='#333')
ax.text(10, 5.5, 'PA EMBRYOS', ha='center', fontsize=11, fontweight='bold', color='#E41A1C')
from matplotlib.patches import FancyBboxPatch
for yi, (y, title, desc, color, ec, alpha) in enumerate([
    (4.5, 'Layer 1: Maternal Blueprint (F1)', 'DNMT1 ZP3 ZP4 GDF9 RARRES1 | r=+0.84', '#377EB8', '#1a5276', 0.2),
    (3.0, 'Layer 2: Maternal Clearance (F3)', 'ZP2 SYCN PARP12 | rho=-0.90', '#E41A1C', '#8B0000', 0.15),
    (1.5, 'Layer 3: Zygotic Activation (F4/F6)', 'IDH2 PKM EXOSC9 GNL3 | rho=+0.64~0.76', '#4DAF4A', '#1B5E20', 0.15),
]):
    l_box = FancyBboxPatch((0.5, y-0.3), 5, 0.8, boxstyle='round,pad=0.05', facecolor=color, edgecolor=ec, lw=1.5, alpha=alpha)
    ax.add_patch(l_box)
    ax.text(3, y+0.25, title, fontsize=8, fontweight='bold', color=ec, ha='center')
    ax.text(3, y-0.1, desc, fontsize=6.5, ha='center', color='#555', style='italic')
for yi, (y, label, desc, color) in enumerate([
    (4.5, 'MAINTAINED', 'DNMT1 HIGHER in PA', '#E41A1C'),
    (3.0, 'VARIABLE', 'ZP2/4 inconsistent', '#FF7F00'),
    (1.5, 'ATTENUATED', 'IDH2(-8.03) PKM(-7.64)', '#377EB8'),
]):
    l_box = FancyBboxPatch((7.5, y-0.3), 5, 0.8, boxstyle='round,pad=0.05', facecolor=color, edgecolor=color, lw=1.5, alpha=0.15)
    ax.add_patch(l_box)
    ax.text(10, y+0.25, f'Layer {yi+1}: {label}', fontsize=8, fontweight='bold', color=color, ha='center')
    ax.text(10, y-0.1, desc, fontsize=6.5, ha='center', color='#555', style='italic')
ax.annotate('', xy=(5.5, 1), xytext=(5.5, 5.2), arrowprops=dict(arrowstyle='->', color='#999', lw=2, ls='--'))
ax.text(5.5, 3, 'MZT', ha='center', fontsize=7, color='#999', rotation=90, va='center')

fig5.tight_layout()
fig5.savefig(f'{OUT}/fig5_pa_validation.svg', format='svg')
fig5.savefig(f'{OUT}/fig5_pa_validation.png', dpi=300)
print(f"  Fig5: {os.path.getsize(f'{OUT}/fig5_pa_validation.svg')//1024} KB")