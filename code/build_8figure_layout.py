#!/usr/bin/env python
"""Create upgraded 8-figure layout for BOR submission.

  Fig1: Atlas UMAP (panel A) + MOFA+ v2 (panels B-G)
  Fig3: MZT trajectory (A-D) + Waddington OT (E-F)
  Fig7: NEW - SCENIC regulon landscape + metabolism
  Fig8: NEW - Developmental divergence + PA volcano
"""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial', 'font.size': 8,
    'svg.fonttype': 'none', 'axes.unicode_minus': False,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch
import pandas as pd, numpy as np, os
from scipy.stats import pearsonr

OUT = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures"
BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
os.makedirs(OUT, exist_ok=True)

def pl(ax, label):
    ax.text(-0.03, 1.03, label, transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')

# ================================================================
# FIGURE 1: Atlas UMAP + MOFA+
# ================================================================
print("Building Fig1 (atlas + MOFA+) ...")
mf = pd.read_csv(f"{BASE}/m4_mofa/mofa_multiomics_factors_v2.csv", index_col=0)
mw_rna = pd.read_csv(f"{BASE}/m4_mofa/mofa_multiomics_weights_RNA_v2.csv", index_col=0)
mw_meth = pd.read_csv(f"{BASE}/m4_mofa/mofa_multiomics_weights_METH_v2.csv", index_col=0)
r2_df = pd.read_csv(f"{BASE}/m4_mofa/mofa_multiomics_r2_per_view_v2.csv")
fcols = [c for c in mf.columns if c.startswith('F')]
meta = pd.read_csv(f"{BASE}/m2_metadata.csv")
donor_map = meta.set_index('cell')['donor'].to_dict()
donors = [donor_map.get(c,'?') for c in mf.index]
dcolors = {'AF1':'#E41A1C','AF2':'#377EB8','AF3':'#4DAF4A','AF4':'#984EA3','AF5':'#FF7F00'}

rf = pd.read_csv(f"{BASE}/m8_mofa_factors.csv", index_col=0)
common = sorted(set(mf.index) & set(rf.index))

# Load UMAP image
umap_path = f"{BASE}/m2_umap_harmony_by_stage_group.png"
has_umap = os.path.exists(umap_path)

fig1 = plt.figure(figsize=(15, 20))
gs = fig1.add_gridspec(4, 2, height_ratios=[1.2, 1, 1, 1], hspace=0.35, wspace=0.3)
fig1.suptitle('', fontsize=13, fontweight='bold')

# A: Atlas UMAP
ax_a = fig1.add_subplot(gs[0, :])
if has_umap:
    img = mpimg.imread(umap_path)
    ax_a.imshow(img)
    ax_a.set_title('A: Integrated porcine preimplantation atlas (UMAP)', fontsize=11, fontweight='bold')
else:
    ax_a.text(0.5, 0.5, 'UMAP image not found', ha='center', va='center', fontsize=14)
    ax_a.set_title('A: Atlas UMAP', fontsize=11)
ax_a.axis('off')

# B: Factor scores
ax_b = fig1.add_subplot(gs[1, 0]); pl(ax_b, 'B')
im = ax_b.imshow(mf[fcols].values.T, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
ax_b.set_yticks(range(7)); ax_b.set_yticklabels(fcols); ax_b.set_xticks([])
ax_b.set_title('Factor scores (32 GV oocytes)', fontsize=10)
plt.colorbar(im, ax=ax_b, shrink=0.6)

# C: Per-view variance
ax_c = fig1.add_subplot(gs[1, 1]); pl(ax_c, 'C')
r2p = r2_df.pivot(index='factor', columns='view', values='r2_pct').fillna(0)
x = np.arange(7)
ax_c.bar(x-0.2, r2p['METH'].reindex(fcols).fillna(0).values, 0.4, color='#FF7F00', label='Methylation', edgecolor='black', lw=0.5)
ax_c.bar(x+0.2, r2p['RNA'].reindex(fcols).fillna(0).values, 0.4, color='#377EB8', label='RNA', edgecolor='black', lw=0.5)
ax_c.set_xticks(x); ax_c.set_xticklabels(fcols); ax_c.set_ylabel('Variance contribution (%)')
ax_c.set_title('Per-view variance contribution', fontsize=10); ax_c.legend(fontsize=7)

# D: Scores x Donor
ax_d = fig1.add_subplot(gs[2, 0]); pl(ax_d, 'D')
offsets = np.linspace(-0.2, 0.2, 5)
for fi in range(7):
    for di, donor in enumerate(['AF1','AF2','AF3','AF4','AF5']):
        vals = [mf.loc[c,fcols[fi]] for c,d in zip(common,donors) if d==donor]
        if vals:
            jitter = np.random.normal(0, 0.03, len(vals))
            ax_d.scatter([fi+1+offsets[di]]*len(vals)+jitter, vals, c=dcolors[donor], s=12, alpha=0.7, edgecolors='none')
ax_d.set_xticks(range(1,8)); ax_d.set_xticklabels(fcols)
ax_d.axhline(0, color='gray', ls='--', lw=0.5); ax_d.set_ylabel('Score'); ax_d.set_xlabel('Factor')
ax_d.set_title('Factor scores by donor', fontsize=10)
from matplotlib.lines import Line2D
ax_d.legend(handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor=dcolors[d],markersize=8,label=d) for d in sorted(dcolors)], fontsize=7, loc='upper right')

# E: F1 multi vs RNA-only
ax_e = fig1.add_subplot(gs[2, 1]); pl(ax_e, 'E')
f1_m = mf['F1']; f1_r = rf.loc[common, 'F1']
r_val, p_val = pearsonr(f1_m, f1_r)
ax_e.scatter(f1_r, f1_m, c='#377EB8', alpha=0.7, s=25, edgecolors='black', lw=0.3)
z = np.polyfit(f1_r, f1_m, 1)
ax_e.plot(np.linspace(f1_r.min(), f1_r.max(), 50), np.polyval(z, np.linspace(f1_r.min(), f1_r.max(), 50)), '--', color='#E41A1C', lw=1.5)
ax_e.set_xlabel('RNA-only MOFA F1 (1833 atlas cells)'); ax_e.set_ylabel('Multi-omics MOFA F1 (32 GV)')
ax_e.set_title(f'E: F1 multi-omics vs RNA-only  r={r_val:.3f}', fontsize=10)

# F: Cross-view weight scatter (F2)
ax_f = fig1.add_subplot(gs[3, 0]); pl(ax_f, 'F')
fi=1; rna_w = mw_rna[fcols[fi]].values; meth_w = mw_meth[fcols[fi]].values
ax_f.scatter(rna_w, meth_w, c='gray', alpha=0.25, s=3)
diff_idx = np.argsort(np.abs(rna_w-meth_w))[::-1][:6]
for idx in diff_idx:
    ax_f.annotate(mw_rna.index[idx], (rna_w[idx], meth_w[idx]), fontsize=6, fontweight='bold')
ax_f.set_xlabel(f'{fcols[fi]} RNA weight'); ax_f.set_ylabel(f'{fcols[fi]} METH weight')
ax_f.axhline(0, color='gray', lw=0.5); ax_f.axvline(0, color='gray', lw=0.5)
ax_f.set_title(f'F: {fcols[fi]} RNA vs METH weights', fontsize=10)

# G: F1 top genes
ax_g = fig1.add_subplot(gs[3, 1]); pl(ax_g, 'G')
f1_rna_top = mw_rna['F1'].sort_values(ascending=False).head(10)
f1_meth_top = mw_meth['F1'].sort_values(ascending=False).head(8)
shared = set(f1_rna_top.index) & set(f1_meth_top.index)
yr = np.arange(len(f1_rna_top)); ym = np.arange(len(f1_meth_top)) + len(f1_rna_top) + 1
colors_r = ['#E41A1C' if g in shared else '#377EB8' for g in f1_rna_top.index]
colors_m = ['#E41A1C' if g in shared else '#FF7F00' for g in f1_meth_top.index]
ax_g.barh(yr, f1_rna_top.values, color=colors_r, edgecolor='black', lw=0.5, height=0.65, alpha=0.9)
ax_g.barh(ym, f1_meth_top.values, color=colors_m, edgecolor='black', lw=0.5, height=0.65, alpha=0.9)
ax_g.set_yticks(list(yr)+list(ym)); ax_g.set_yticklabels(list(f1_rna_top.index)+list(f1_meth_top.index), fontsize=7)
ax_g.set_xlabel('Factor 1 weight'); ax_g.set_title('G: Top F1 genes', fontsize=10)
ax_g.axhline(len(f1_rna_top)-0.5, color='gray', lw=1, ls='--')
ax_g.text(ax_g.get_xlim()[1]*0.85, len(f1_rna_top)-1, 'RNA up', fontsize=6, color='#377EB8', style='italic')
ax_g.text(ax_g.get_xlim()[1]*0.85, len(f1_rna_top), 'METH dn', fontsize=6, color='#FF7F00', style='italic')

fig1.tight_layout()
fig1.savefig(f"{OUT}/fig1_mofa.svg", format='svg')
fig1.savefig(f"{OUT}/fig1_mofa.png", dpi=300)
print(f"  Fig1: {os.path.getsize(f'{OUT}/fig1_mofa.svg')//1024} KB")

# ================================================================
# FIGURE 3: MZT trajectory + Waddington OT
# ================================================================
print("Building Fig3 (trajectory + OT) ...")
proj = pd.read_csv(f"{BASE}/m4_mofa/mzt_trajectory_projection.csv", index_col=0)
stage_mean = pd.read_csv(f"{BASE}/m4_mofa/mzt_stage_mean_projection.csv", index_col=0)
stages_all = sorted(stage_mean.index.tolist(), key=lambda s: int(s[1:]))

ot_path = f"{BASE}/m5_waddington_ot.png"
has_ot = os.path.exists(ot_path)

fig3 = plt.figure(figsize=(15, 18))
gs3 = fig3.add_gridspec(3, 2, height_ratios=[1, 1, 1.1], hspace=0.35, wspace=0.3)

# A: All factors
ax_a3 = fig3.add_subplot(gs3[0, 0]); pl(ax_a3, 'A')
colors7 = ["#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#999999"]
for fi, fc in enumerate(fcols):
    vals = [stage_mean.loc[s, fc] for s in stages_all if s in stage_mean.index]
    stages_f = [s for s in stages_all if s in stage_mean.index]
    ax_a3.plot(stages_f, vals, 'o-', color=colors7[fi], label=fc, lw=1.5, markersize=4)
ax_a3.set_xlabel('Stage'); ax_a3.set_ylabel('Projected score')
ax_a3.set_title('Projected factor scores (E0-E10)', fontsize=10)
ax_a3.legend(fontsize=7, ncol=4); ax_a3.axhline(0, color='gray', ls='--', lw=0.5)

# B: F3
ax_b3 = fig3.add_subplot(gs3[0, 1]); pl(ax_b3, 'B')
f3v = [stage_mean.loc[s,'F3'] for s in stages_all if s in stage_mean.index]
from scipy.stats import spearmanr
r3, p3 = spearmanr(np.arange(len(f3v)), f3v)
ax_b3.bar(stages_all[:len(f3v)], f3v, color='#E41A1C', alpha=0.7, edgecolor='black')
ax_b3.set_xlabel('Stage'); ax_b3.set_ylabel('F3 score')
ax_b3.set_title(f'B: F3 Clearance (rho={r3:.2f}, P={p3:.1e})', fontsize=10)
ax_b3.axhline(0, color='gray', lw=0.5)

# C: F4+F6
ax_c3 = fig3.add_subplot(gs3[1, 0]); pl(ax_c3, 'C')
for fc, color in [('F4','#377EB8'),('F6','#FF7F00')]:
    vals = [stage_mean.loc[s,fc] for s in stages_all if s in stage_mean.index]
    rv, pv = spearmanr(np.arange(len(vals)), vals)
    ax_c3.plot(stages_all[:len(vals)], vals, 'o-', color=color, label=f'{fc} (r={rv:.2f})', lw=1.5, markersize=4)
ax_c3.set_xlabel('Stage'); ax_c3.set_ylabel('Projected score')
ax_c3.set_title('C: F4/F6 Zygotic Activation', fontsize=10)
ax_c3.legend(fontsize=8); ax_c3.axhline(0, color='gray', ls='--', lw=0.5)

# D: Layer summary
ax_d3 = fig3.add_subplot(gs3[1, 1]); pl(ax_d3, 'D')
layers = ['F1 Blueprint\n(stable)', 'F3 Clearance\n(declining)', 'F4 ZGA\n(rising)', 'F6 ZGA\n(rising)']
r_vals = [-0.41, -0.90, 0.64, 0.76]
p_vals = [0.21, 2e-4, 0.035, 0.006]
colors_layer = ['#377EB8','#E41A1C','#4DAF4A','#FF7F00']
ax_d3.barh(layers, r_vals, color=colors_layer, edgecolor='black')
for i, (r, p) in enumerate(zip(r_vals, p_vals)):
    ax_d3.text(r+0.03*np.sign(r), i, f'r={r:+.2f}  P={p:.1e}', va='center', fontsize=8)
ax_d3.axvline(0, color='gray', lw=0.8); ax_d3.set_xlabel("Spearman rho (per-stage pseudobulk)")
ax_d3.set_title('D: Three-layer temporal dynamics', fontsize=10)

# E+F: Waddington OT
ax_e3 = fig3.add_subplot(gs3[2, :])
if has_ot:
    img = mpimg.imread(ot_path)
    ax_e3.imshow(img)
    ax_e3.set_title('E: Waddington optimal transport — developmental canalization and embryo displacement', fontsize=11, fontweight='bold')
else:
    ax_e3.text(0.5, 0.5, 'OT image not found', ha='center', va='center')
ax_e3.axis('off')

fig3.tight_layout()
fig3.savefig(f"{OUT}/fig3_mzt_trajectory.svg", format='svg')
fig3.savefig(f"{OUT}/fig3_mzt_trajectory.png", dpi=300)
print(f"  Fig3: {os.path.getsize(f'{OUT}/fig3_mzt_trajectory.svg')//1024} KB")

# ================================================================
# FIGURE 7: SCENIC regulon landscape + metabolism
# ================================================================
print("Building Fig7 (SCENIC regulons + metabolism) ...")
reg_heat_path = f"{BASE}/pyscenic_out/fig5b_regulon_heatmap_stage.png"
reg_stage_path = f"{BASE}/pyscenic_out/fig5b_regulon_by_stage.png"
metab_path = f"{BASE}/m6_metabolism.png"

fig7 = plt.figure(figsize=(15, 12))
gs7 = fig7.add_gridspec(2, 2, hspace=0.3, wspace=0.25)

# A: Regulon heatmap
ax_7a = fig7.add_subplot(gs7[0, 0])
if os.path.exists(reg_heat_path):
    img = mpimg.imread(reg_heat_path)
    ax_7a.imshow(img)
ax_7a.set_title('A: Stage-conserved regulon activity (AUCell)', fontsize=11, fontweight='bold')
ax_7a.axis('off')

# B: Regulon by stage
ax_7b = fig7.add_subplot(gs7[0, 1])
if os.path.exists(reg_stage_path):
    img = mpimg.imread(reg_stage_path)
    ax_7b.imshow(img)
ax_7b.set_title('B: Key regulon temporal dynamics', fontsize=11, fontweight='bold')
ax_7b.axis('off')

# C: Metabolism
ax_7c = fig7.add_subplot(gs7[1, :])
if os.path.exists(metab_path):
    img = mpimg.imread(metab_path)
    ax_7c.imshow(img)
ax_7c.set_title('C: Transcriptome-based metabolic scoring (Glycolysis/OxPhos ratio)', fontsize=11, fontweight='bold')
ax_7c.axis('off')

fig7.tight_layout()
fig7.savefig(f"{OUT}/fig7_regulon_metabolism.svg", format='svg')
fig7.savefig(f"{OUT}/fig7_regulon_metabolism.png", dpi=300)
print(f"  Fig7: {os.path.getsize(f'{OUT}/fig7_regulon_metabolism.svg')//1024} KB")

# ================================================================
# FIGURE 8: Developmental divergence + PA volcano
# ================================================================
print("Building Fig8 (divergence + volcano) ...")
div_path = f"{BASE}/m5_trajectory_comparison.png"
volc_path = f"{BASE}/m3_ivf_vs_pa_volcano.png"

fig8 = plt.figure(figsize=(15, 8))
gs8 = fig8.add_gridspec(1, 2, wspace=0.2)

# A: Trajectory divergence
ax_8a = fig8.add_subplot(gs8[0, 0])
if os.path.exists(div_path):
    img = mpimg.imread(div_path)
    ax_8a.imshow(img)
ax_8a.set_title('A: Cross-condition pseudotime divergence', fontsize=11, fontweight='bold')
ax_8a.axis('off')

# B: IVF vs PA volcano
ax_8b = fig8.add_subplot(gs8[0, 1])
if os.path.exists(volc_path):
    img = mpimg.imread(volc_path)
    ax_8b.imshow(img)
ax_8b.set_title('B: IVF vs PA differential expression (public scRNA)', fontsize=11, fontweight='bold')
ax_8b.axis('off')

fig8.tight_layout()
fig8.savefig(f"{OUT}/fig8_divergence_volcano.svg", format='svg')
fig8.savefig(f"{OUT}/fig8_divergence_volcano.png", dpi=300)
print(f"  Fig8: {os.path.getsize(f'{OUT}/fig8_divergence_volcano.svg')//1024} KB")

# ================================================================
print("\n=== ALL 8 FIGURES DONE ===")
for f in ['fig1_mofa','fig2_crossspecies','fig3_mzt_trajectory','fig4_perturbation',
          'fig5_pa_validation','fig6_model','fig7_regulon_metabolism','fig8_divergence_volcano']:
    svg = f'{OUT}/{f}.svg'
    if os.path.exists(svg):
        print(f"  {f}.svg: {os.path.getsize(svg)//1024} KB")
    else:
        print(f"  {f}.svg: MISSING")