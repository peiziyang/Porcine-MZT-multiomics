#!/usr/bin/env python
"""Fix Fig5 (merge Fig6 model), Fig7C (split metabolism into sub-panels), Fig8A (split divergence).

Changes:
  Fig5: Add panel E as the 3-layer model schematic (from Fig6), keep A-D same
  Fig7C: Reconstruct metabolism as 3 sub-panels from data (not PNG paste)
  Fig8A: Reconstruct trajectory divergence sub-panels from M5 data
  Fig6: Retire — content now in Fig5E
"""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial', 'font.size': 8,
    'svg.fonttype': 'none', 'axes.unicode_minus': False,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
})
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch
import pandas as pd, numpy as np, os
from scipy.stats import spearmanr

OUT = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures"
BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
os.makedirs(OUT, exist_ok=True)

def pl(ax, label):
    ax.text(-0.03, 1.03, label, transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')

# ================================================================
# FIGURE 5: PA validation + 3-layer model (merge Fig6 into panel E)
# ================================================================
print("Rebuilding Fig5 (add 3-layer model as panel E) ...")
from matplotlib.gridspec import GridSpec

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
gs5 = GridSpec(3, 2, figure=fig5, height_ratios=[1, 1, 1.15], hspace=0.35, wspace=0.3)

# A: Layer 1
ax = fig5.add_subplot(gs5[0, 0]); pl(ax, 'A')
l1 = ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2']
l1v = []; l1p = []; l1gf = []
for g in l1:
    sub = overall[overall['gene']==g]
    if len(sub): l1v.append(sub['log2FC'].values[0]); l1p.append(sub['padj'].values[0]); l1gf.append(g)
xl = np.arange(len(l1gf))
bc = ['#E41A1C' if v>0.5 else '#FF7F00' if v>0 else '#377EB8' if v<-0.5 else '#999' for v in l1v]
ax.bar(xl, l1v, color=bc, edgecolor='black', lw=0.8, width=0.6)
for i,(v,p) in enumerate(zip(l1v,l1p)):
    sig = 'p<0.05' if p<0.05 else 'ns'
    ax.text(i, v+0.3*np.sign(v) if v!=0 else 0.3, f'{v:+.2f} ({sig})', ha='center', fontsize=7, fontweight='bold', color='#333')
ax.axhline(0, color='gray', lw=1); ax.set_xticks(xl); ax.set_xticklabels(l1gf, fontsize=8)
ax.set_ylabel('log2FC (PA vs IVF)', fontsize=9)
ax.set_title('Layer 1: Maternal blueprint - MAINTAINED', fontsize=10, fontweight='bold', color='#E41A1C')

# B: Layer 3
ax = fig5.add_subplot(gs5[0, 1]); pl(ax, 'B')
l3 = ['IDH2','PKM','EXOSC9','MDH1','GNL3']
l3v = []; l3p = []; l3gf = []
for g in l3:
    sub = overall[overall['gene']==g]
    if len(sub): l3v.append(sub['log2FC'].values[0]); l3p.append(sub['padj'].values[0]); l3gf.append(g)
xl3 = np.arange(len(l3gf))
ax.bar(xl3, l3v, color=['#377EB8' if v<-2 else '#999' for v in l3v], edgecolor='black', lw=0.8, width=0.6)
for i,(v,p) in enumerate(zip(l3v,l3p)):
    sig = 'p<1e-6' if p<1e-6 else f'p={p:.0e}' if p<0.05 else 'ns'
    ax.text(i, v+0.5, f'{v:+.2f} ({sig})', ha='center', fontsize=7, fontweight='bold', color='#377EB8')
ax.axhline(0, color='gray', lw=1); ax.set_xticks(xl3); ax.set_xticklabels(l3gf, fontsize=8)
ax.set_ylabel('log2FC (PA vs IVF)', fontsize=9)
ax.set_title('Layer 3: Zygotic execution - ATTENUATED', fontsize=10, fontweight='bold', color='#377EB8')

# C: Per-stage heatmap
ax = fig5.add_subplot(gs5[1, 0]); pl(ax, 'C')
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
        if not np.isnan(v): ax.text(j,i,f'{v:+.1f}',ha='center',va='center',fontsize=6,fontweight='bold',color='white' if abs(v)>5 else 'black')
ax.set_title('Per-stage log2FC', fontsize=10)
plt.colorbar(im, ax=ax, shrink=0.75, label='log2FC')

# D: Summary table
ax = fig5.add_subplot(gs5[1, 1]); pl(ax, 'D')
ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis('off')
ys = 9
for title, col, rows in [
    ('Layer 1: Blueprint','#E41A1C',[('DNMT1','+1.30','**'),('ZP3','+0.21','ns'),('ZP4','-0.75','ns'),('GDF9','+0.62','ns'),('RARRES1','-0.30','ns')]),
    ('Layer 3: Execution','#377EB8',[('IDH2 (TCA)','-8.03','***'),('PKM (glycolysis)','-7.64','***'),('EXOSC9 (ribosome)','-7.73','***'),('MDH1 (TCA)','-4.58','**'),('GNL3 (nucleolar)','-3.06','*')]),
    ('ZGA Markers','#4DAF4A',[('DPPA5','+2.14','**'),('NANOG','+3.91','***')])]:
    ax.text(0.3,ys,title,fontsize=9,fontweight='bold',color=col,va='center'); ys-=0.4
    for gene,lfc,sig in rows:
        ax.text(0.5,ys,gene,fontsize=7.5,va='center')
        fc = float(lfc.replace('+',''))
        cl = '#E41A1C' if fc>1 else '#377EB8' if fc<-1 else '#666'
        ax.text(4.5,ys,f'log2FC={lfc}',fontsize=7.5,va='center',fontweight='bold',color=cl)
        sg = '[***]' if sig=='***' else '[**]' if sig=='**' else '[*]' if sig=='*' else ''
        ax.text(7.5,ys,sg,fontsize=7.5,va='center',color='#CC0000')
        ys-=0.3
    ys-=0.2
ax.text(0.5, ys-0.2, '[*] Padj<0.05  [**] Padj<0.01  [***] Padj<0.001', fontsize=6.5, color='gray')

# E: Three-layer model (from Fig6, compact)
ax = fig5.add_subplot(gs5[2, :]); pl(ax, 'E')
ax.set_xlim(0, 14); ax.set_ylim(0, 6); ax.axis('off')
ax.text(3, 5.5, 'NORMAL MZT', ha='center', fontsize=11, fontweight='bold', color='#333')
ax.text(10, 5.5, 'PA EMBRYOS', ha='center', fontsize=11, fontweight='bold', color='#E41A1C')

# Three layers compact
for yi, (y, title, desc, color, ec, alpha) in enumerate([
    (4.5, 'Layer 1: Maternal Blueprint (F1)', 'DNMT1 ZP3 ZP4 GDF9 RARRES1 | r=+0.84', '#377EB8', '#1a5276', 0.2),
    (3.0, 'Layer 2: Maternal Clearance (F3)', 'ZP2 SYCN PARP12 | rho=-0.90', '#E41A1C', '#8B0000', 0.15),
    (1.5, 'Layer 3: Zygotic Activation (F4/F6)', 'IDH2 PKM EXOSC9 GNL3 | rho=+0.64~0.76', '#4DAF4A', '#1B5E20', 0.15),
]):
    l_box = FancyBboxPatch((0.5, y-0.3), 5, 0.8, boxstyle='round,pad=0.05', facecolor=color, edgecolor=ec, lw=1.5, alpha=alpha)
    ax.add_patch(l_box)
    ax.text(3, y+0.25, title, fontsize=8, fontweight='bold', color=ec, ha='center')
    ax.text(3, y-0.1, desc, fontsize=6.5, ha='center', color='#555', style='italic')

# PA overlay (right side)
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

# ================================================================
# FIGURE 7C: Split metabolism into sub-panels (NOT PNG paste)
# ================================================================
print("Rebuilding Fig7C (metabolism sub-panels from data) ...")

# Load metabolism data from M6
metab_path = f"{BASE}/m6_metabolism_ratio.csv"
reg_heat = f"{BASE}/pyscenic_out/fig5b_regulon_heatmap_stage.png"
reg_stage = f"{BASE}/pyscenic_out/fig5b_regulon_by_stage.png"

fig7 = plt.figure(figsize=(14, 14))
gs7 = fig7.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

# A: Regulon heatmap
ax = fig7.add_subplot(gs7[0, 0])
if os.path.exists(reg_heat):
    ax.imshow(mpimg.imread(reg_heat))
ax.set_title('A: Stage-conserved regulon activity (AUCell)', fontsize=10, fontweight='bold')
ax.axis('off')

# B: Regulon by stage
ax = fig7.add_subplot(gs7[0, 1])
if os.path.exists(reg_stage):
    ax.imshow(mpimg.imread(reg_stage))
ax.set_title('B: Key regulon temporal dynamics', fontsize=10, fontweight='bold')
ax.axis('off')

# C: Metabolism — 3 sub-panels from data
ax_c = fig7.add_subplot(gs7[1, :])
ax_c.set_title('C: Transcriptome-based metabolic scoring', fontsize=11, fontweight='bold', loc='center', pad=10)

# If metabolism data exists, plot from scratch
if os.path.exists(metab_path):
    met = pd.read_csv(metab_path)
    # Create 3 inline subplots within this axis using inset_axes
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
    # Sub-panel i: Glycolysis/OxPhos ratio barplot
    ax_c1 = inset_axes(ax_c, width='30%', height='70%', loc='lower left', borderpad=0)
    ax_c1.set_title('(i) Glycolysis/OxPhos ratio', fontsize=8)
    if 'stage' in met.columns and 'ratio' in met.columns:
        stages_met = met['stage'].values[:6]
        ratios = met['ratio'].values[:6]
        ax_c1.bar(range(len(stages_met)), ratios, color='#377EB8', edgecolor='black', width=0.6)
        ax_c1.set_xticks(range(len(stages_met))); ax_c1.set_xticklabels(stages_met, fontsize=6, rotation=45)
        ax_c1.set_ylabel('Ratio', fontsize=7); ax_c1.axhline(1, color='red', ls='--', lw=0.5)
    else:
        ax_c1.text(0.5, 0.5, 'Metabolism data\nnot available', ha='center', va='center', fontsize=7)
    ax_c1.axis('off') if 'stage' not in met.columns else None
    
    # Sub-panel ii: Glycolysis score per stage
    ax_c2 = inset_axes(ax_c, width='30%', height='70%', loc='lower center', borderpad=0)
    ax_c2.set_title('(ii) Glycolysis module score', fontsize=8)
    if 'glycolysis' in met.columns:
        gl = met['glycolysis'].values[:6]
        ax_c2.bar(range(len(stages_met)), gl, color='#E41A1C', edgecolor='black', width=0.6)
        ax_c2.set_xticks(range(len(stages_met))); ax_c2.set_xticklabels(stages_met, fontsize=6, rotation=45)
        ax_c2.set_ylabel('Score', fontsize=7)
    else:
        ax_c2.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=7)
    ax_c2.axis('off') if 'glycolysis' not in met.columns else None
    
    # Sub-panel iii: OxPhos score per stage
    ax_c3 = inset_axes(ax_c, width='30%', height='70%', loc='lower right', borderpad=0)
    ax_c3.set_title('(iii) OxPhos module score', fontsize=8)
    if 'oxphos' in met.columns:
        ox = met['oxphos'].values[:6]
        ax_c3.bar(range(len(stages_met)), ox, color='#4DAF4A', edgecolor='black', width=0.6)
        ax_c3.set_xticks(range(len(stages_met))); ax_c3.set_xticklabels(stages_met, fontsize=6, rotation=45)
        ax_c3.set_ylabel('Score', fontsize=7)
    else:
        ax_c3.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=7)
    ax_c3.axis('off') if 'oxphos' not in met.columns else None
else:
    ax_c.text(0.5, 0.5, 'Metabolism data file not found', ha='center', va='center', fontsize=9)

ax_c.axis('off')

fig7.tight_layout()
fig7.savefig(f'{OUT}/fig7_regulon_metabolism.svg', format='svg')
fig7.savefig(f'{OUT}/fig7_regulon_metabolism.png', dpi=300)
print(f"  Fig7: {os.path.getsize(f'{OUT}/fig7_regulon_metabolism.svg')//1024} KB")

# ================================================================
# FIGURE 8A: Split trajectory divergence into sub-panels
# ================================================================
print("Rebuilding Fig8A (trajectory sub-panels from M5 data) ...")

div_path = f"{BASE}/m5_gene_divergence_v2.csv"
volc_path = f"{BASE}/m3_ivf_vs_pa_volcano.png"
traj_path = f"{BASE}/m5_trajectory_comparison.png"

fig8 = plt.figure(figsize=(15, 9))
gs8 = fig8.add_gridspec(1, 2, wspace=0.3)

# A: Trajectory divergence — 3 sub-panels
ax_a = fig8.add_subplot(gs8[0, 0])
ax_a.set_title('A: Pseudotime-based cross-condition divergence', fontsize=11, fontweight='bold', loc='center', pad=8)

if os.path.exists(div_path):
    div_df = pd.read_csv(div_path)
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
    # (i) Top 50 divergent genes: mean divergence per comparison
    ax_a1 = inset_axes(ax_a, width='28%', height='45%', loc='upper left', borderpad=0.5)
    ax_a1.set_title('(i) Mean divergence', fontsize=7)
    if 'PA_invivo' in div_df.columns:
        cats = ['IVF vs PA', 'IVF vs InVivo', 'PA vs InVivo']
        vals = [div_df['IVF_PA'].mean() if 'IVF_PA' in div_df.columns else 0,
                div_df['IVF_invivo'].mean() if 'IVF_invivo' in div_df.columns else 0,
                div_df['PA_invivo'].mean() if 'PA_invivo' in div_df.columns else 0]
        if sum(vals) > 0:
            ax_a1.bar(cats, vals, color=['#377EB8','#4DAF4A','#E41A1C'], edgecolor='black', width=0.5)
            ax_a1.set_ylabel('Divergence', fontsize=6)
            ax_a1.tick_params(labelsize=5)
        else:
            ax_a1.text(0.5,0.5,'N/A',ha='center',va='center',fontsize=7)
    else:
        ax_a1.text(0.5,0.5,'N/A',ha='center',va='center',fontsize=7)
    
    # (ii) PA_max_genes: in how many of top 50 is PA the max deviate?
    ax_a2 = inset_axes(ax_a, width='28%', height='45%', loc='upper center', borderpad=0.5)
    ax_a2.set_title('(ii) Max deviation source', fontsize=7)
    if 'PA_invivo' in div_df.columns and 'IVF_invivo' in div_df.columns:
        top50 = div_df.sort_values('PA_invivo', ascending=False).head(50)
        pa_max = (top50['PA_invivo'] > top50['IVF_invivo']).sum()
        ivf_max = (top50['IVF_invivo'] > top50['PA_invivo']).sum()
        ax_a2.bar(['PA', 'IVF'], [pa_max, ivf_max], color=['#E41A1C','#377EB8'], edgecolor='black', width=0.4)
        ax_a2.tick_params(labelsize=5)
    else:
        ax_a2.text(0.5,0.5,'N/A',ha='center',va='center',fontsize=7)
    
    # (iii) Top divergent genes barplot
    ax_a3 = inset_axes(ax_a, width='38%', height='45%', loc='upper right', borderpad=0.5)
    ax_a3.set_title('(iii) Top genes: PA divergence', fontsize=7)
    if 'gene' in div_df.columns and 'PA_invivo' in div_df.columns:
        top10 = div_df.sort_values('PA_invivo', ascending=False).head(10)
        ax_a3.barh(range(10), top10['PA_invivo'].values, color='#E41A1C', edgecolor='black', height=0.6)
        ax_a3.set_yticks(range(10)); ax_a3.set_yticklabels(top10['gene'].values, fontsize=5)
        ax_a3.invert_yaxis(); ax_a3.set_xlabel('Divergence', fontsize=5)
        ax_a3.tick_params(labelsize=5)
    else:
        ax_a3.text(0.5,0.5,'N/A',ha='center',va='center',fontsize=7)
    
    # (iv) If trajectory comparison PNG exists, show it in bottom half
    ax_a4 = inset_axes(ax_a, width='90%', height='40%', loc='lower center', borderpad=0.5)
    if os.path.exists(traj_path):
        ax_a4.imshow(mpimg.imread(traj_path))
        ax_a4.set_title('(iv) Cross-condition trajectory comparison', fontsize=7)
    else:
        ax_a4.text(0.5,0.5,'Trajectory image unavailable',ha='center',va='center',fontsize=7)
    ax_a4.axis('off')

ax_a.axis('off')

# B: IVF vs PA volcano
ax_b = fig8.add_subplot(gs8[0, 1])
if os.path.exists(volc_path):
    ax_b.imshow(mpimg.imread(volc_path))
ax_b.set_title('B: IVF vs PA differential expression (public scRNA)', fontsize=11, fontweight='bold')
ax_b.axis('off')

fig8.tight_layout()
fig8.savefig(f'{OUT}/fig8_divergence_volcano.svg', format='svg')
fig8.savefig(f'{OUT}/fig8_divergence_volcano.png', dpi=300)
print(f"  Fig8: {os.path.getsize(f'{OUT}/fig8_divergence_volcano.svg')//1024} KB")

# ================================================================
# Update manuscript: remove Fig6 references, renumber Fig7→Fig6, Fig8→Fig7
# ================================================================
print("\nDone. Fig6 merged into Fig5E. Now Fig7=regulon+metab, Fig8→renumber.")
print("Remember: manuscript needs Fig6→Fig5E, Fig7→Fig6, Fig8→Fig7 in text+legends.")