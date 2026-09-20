#!/usr/bin/env python
"""Regenerate Fig6 and Fig7 with correct data column names."""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial', 'font.size': 9,
    'svg.fonttype': 'none', 'axes.unicode_minus': False,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd, numpy as np, os

OUT = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures"
BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"

# ================================================================
# FIGURE 6: SCENIC regulon + metabolism (3 sub-panels)
# ================================================================
print("Rebuilding Fig6 with correct metabolism data ...")
reg_heat_path = f"{BASE}/pyscenic_out/fig5b_regulon_heatmap_stage.png"
reg_stage_path = f"{BASE}/pyscenic_out/fig5b_regulon_by_stage.png"
met = pd.read_csv(f"{BASE}/m6_metabolism_ratio.csv")

fig6 = plt.figure(figsize=(14, 14))
gs6 = fig6.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

# A: Regulon heatmap
ax = fig6.add_subplot(gs6[0, 0])
ax.text(-0.03, 1.03, 'A', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')
if os.path.exists(reg_heat_path):
    ax.imshow(mpimg.imread(reg_heat_path))
ax.set_title('Stage-conserved regulon activity (AUCell)', fontsize=10)
ax.axis('off')

# B: Regulon by stage
ax = fig6.add_subplot(gs6[0, 1])
ax.text(-0.03, 1.03, 'B', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')
if os.path.exists(reg_stage_path):
    ax.imshow(mpimg.imread(reg_stage_path))
ax.set_title('Key regulon temporal dynamics', fontsize=10)
ax.axis('off')

# C: 3 sub-panels from metabolism data — use a fresh subgridspec
gs_c_inner = fig6.add_gridspec(2, 3, left=0.04, bottom=0.04, right=0.99, top=0.92, hspace=0.5, wspace=0.35)

stages = met['stage'].values
x = np.arange(len(stages))

# (i) Ratio per stage
ax_i = fig6.add_subplot(gs_c_inner[0, 0])
ax_i.bar(x, met['ratio'].values, color='#377EB8', edgecolor='black', width=0.7)
ax_i.axhline(1, color='red', ls='--', lw=1, label='equal')
ax_i.set_xticks(x[::2]); ax_i.set_xticklabels(stages[::2], fontsize=6, rotation=45, ha='right')
ax_i.set_ylabel('Glycolysis/OxPhos ratio', fontsize=7)
ax_i.set_title('(i) Ratio (Glyc/OxPhos)', fontsize=9, fontweight='bold')
ax_i.tick_params(labelsize=6)
ax_i.set_ylim(0, 0.7)

# (ii) Glycolysis per stage
ax_ii = fig6.add_subplot(gs_c_inner[0, 1])
ax_ii.bar(x, met['gly_mean'].values, color='#E41A1C', edgecolor='black', width=0.7)
ax_ii.set_xticks(x[::2]); ax_ii.set_xticklabels(stages[::2], fontsize=6, rotation=45, ha='right')
ax_ii.set_ylabel('Module score', fontsize=7)
ax_ii.set_title('(ii) Glycolysis module score', fontsize=9, fontweight='bold')
ax_ii.tick_params(labelsize=6)

# (iii) OxPhos per stage
ax_iii = fig6.add_subplot(gs_c_inner[0, 2])
ax_iii.bar(x, met['ox_mean'].values, color='#4DAF4A', edgecolor='black', width=0.7)
ax_iii.set_xticks(x[::2]); ax_iii.set_xticklabels(stages[::2], fontsize=6, rotation=45, ha='right')
ax_iii.set_ylabel('Module score', fontsize=7)
ax_iii.set_title('(iii) OxPhos module score', fontsize=9, fontweight='bold')
ax_iii.tick_params(labelsize=6)

# Add a 4th panel: summary table or visualization
ax_summary = fig6.add_subplot(gs_c_inner[1, :])
ax_summary.axis('off')
summary_text = ('Metabolic scoring summary:\n'
                f'  Mean Glyc/OxPhos ratio across E0-E14: {met["ratio"].mean():.3f}\n'
                f'  Lowest ratio: {met["ratio"].min():.3f} (E{met.loc[met["ratio"].idxmin(),"stage"][1:]})\n'
                f'  Highest ratio: {met["ratio"].max():.3f} (E{met.loc[met["ratio"].idxmax(),"stage"][1:]})\n'
                '  Interpretation: preimplantation embryos maintain OxPhos-skewed\n'
                '  transcriptional profile throughout E0-E14 (ratio < 1 in all stages)')
ax_summary.text(0.0, 0.95, summary_text, transform=ax_summary.transAxes, fontsize=8,
                family='monospace', va='top',
                bbox=dict(boxstyle='round', facecolor='#F0F8FF', edgecolor='#1a5276', lw=1))

fig6.tight_layout()
fig6.savefig(f'{OUT}/fig6_regulon_metabolism.svg', format='svg')
fig6.savefig(f'{OUT}/fig6_regulon_metabolism.png', dpi=300)
print(f"  Fig6: {os.path.getsize(f'{OUT}/fig6_regulon_metabolism.svg')//1024} KB")

# ================================================================
# FIGURE 7: Trajectory divergence + IVF/PA volcano
# ================================================================
print("Rebuilding Fig7 with correct m5 column names ...")
volc_path = f"{BASE}/m3_ivf_vs_pa_volcano.png"

# Load divergence with CORRECT columns
div_df = pd.read_csv(f"{BASE}/m5_gene_divergence_v2.csv")
print(f"  Total genes: {len(div_df)}")
print(f"  IVF_vs_PA mean: {div_df['IVF_vs_PA'].mean():.4f}")
print(f"  IVF_vs_vivo mean: {div_df['IVF_vs_vivo'].mean():.4f}")
print(f"  PA_vs_vivo mean: {div_df['PA_vs_vivo'].mean():.4f}")

fig7 = plt.figure(figsize=(15, 12))
gs7 = fig7.add_gridspec(2, 2, height_ratios=[1, 1.1], hspace=0.4, wspace=0.3)

# A: 4 sub-panels from m5 data — use nested gridspec for proper sizing
ax = fig7.add_subplot(gs7[0, 0])
ax.text(-0.05, 1.02, 'A', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')
ax.set_title('Pseudotime-based cross-condition divergence (m5, n=8,276 genes)', fontsize=10, fontweight='bold')
ax.axis('off')

# Use proper subgridspec
gs_a = ax.inset_axes([0, 0, 1, 1]).get_subplotspec()
# Subdivide with proper axes
ax_i = fig7.add_subplot(gs7[0, 0])
# Actually use subgridspec for 2x2 within
gs7a = fig7.add_gridspec(2, 2, left=0.04, bottom=0.04, right=0.48, top=0.92, hspace=0.5, wspace=0.35)

ax_i = fig7.add_subplot(gs7a[0, 0])
ivf_pa = div_df['IVF_vs_PA'].mean()
ivf_vivo = div_df['IVF_vs_vivo'].mean()
pa_vivo = div_df['PA_vs_vivo'].mean()
ax_i.bar(['IVF-PA', 'IVF-InVivo', 'PA-InVivo'], [ivf_pa, ivf_vivo, pa_vivo],
          color=['#377EB8','#4DAF4A','#E41A1C'], edgecolor='black', width=0.5)
ax_i.set_ylabel('Mean divergence', fontsize=7)
ax_i.set_title('(i) Mean divergence per comparison', fontsize=8, fontweight='bold')
ax_i.tick_params(labelsize=6)
for i, v in enumerate([ivf_pa, ivf_vivo, pa_vivo]):
    ax_i.text(i, v+0.005, f'{v:.3f}', ha='center', fontsize=6, fontweight='bold')

ax_ii = fig7.add_subplot(gs7a[0, 1])
top50 = div_df.sort_values('PA_vs_vivo', ascending=False).head(50)
pa_max = (top50['PA_vs_vivo'] > top50['IVF_vs_vivo']).sum()
ivf_max = (top50['IVF_vs_vivo'] > top50['PA_vs_vivo']).sum()
ax_ii.bar(['PA highest', 'IVF highest'], [pa_max, ivf_max],
          color=['#E41A1C','#377EB8'], edgecolor='black', width=0.4)
ax_ii.set_ylabel('Number of genes (top 50)', fontsize=7)
ax_ii.set_title('(ii) Top-50 max deviation source', fontsize=8, fontweight='bold')
ax_ii.tick_params(labelsize=6)
for i, v in enumerate([pa_max, ivf_max]):
    ax_ii.text(i, v+1, str(v), ha='center', fontsize=7, fontweight='bold')

ax_iii = fig7.add_subplot(gs7a[1, 0])
top15 = div_df.sort_values('PA_vs_vivo', ascending=False).head(15)
colors_g = ['#E41A1C' if v > 0.5 else '#FF7F00' for v in top15['PA_vs_vivo'].values]
ax_iii.barh(range(15), top15['PA_vs_vivo'].values, color=colors_g, edgecolor='black', height=0.7)
ax_iii.set_yticks(range(15))
ax_iii.set_yticklabels(top15['gene'].values, fontsize=6)
ax_iii.invert_yaxis()
ax_iii.set_xlabel('PA vs in vivo divergence', fontsize=7)
ax_iii.set_title('(iii) Top 15 PA-divergent genes', fontsize=8, fontweight='bold')
ax_iii.tick_params(labelsize=6)

ax_iv = fig7.add_subplot(gs7a[1, 1])
top10 = div_df.sort_values('PA_vs_vivo', ascending=False).head(10)
x_pos = np.arange(3)
width = 0.08
for gi, (_, row) in enumerate(top10.iterrows()):
    vals = [row['IVF_vs_PA'], row['IVF_vs_vivo'], row['PA_vs_vivo']]
    ax_iv.bar(x_pos + gi*width, vals, width, color=plt.cm.tab20(gi/10), edgecolor='black', linewidth=0.2)
ax_iv.set_xticks(x_pos + 4.5*width)
ax_iv.set_xticklabels(['IVF-PA', 'IVF-InVivo', 'PA-InVivo'], fontsize=6)
ax_iv.set_ylabel('Divergence', fontsize=7)
ax_iv.set_title('(iv) Top 10 across 3 comparisons', fontsize=8, fontweight='bold')
ax_iv.tick_params(labelsize=6)
# Add legend
ax_iv.legend(top10['gene'].values, fontsize=5, loc='upper right', ncol=2)

# Remove the placeholder ax
ax.remove()

# B: Volcano
ax = fig7.add_subplot(gs7[0, 1])
ax.text(-0.05, 1.02, 'B', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')
if os.path.exists(volc_path):
    ax.imshow(mpimg.imread(volc_path))
ax.set_title('IVF vs PA differential expression (public scRNA, GSE164812)', fontsize=10, fontweight='bold')
ax.axis('off')

# Add a 5th panel: IVF-PA and PA-vivo side-by-side comparison for documentation
ax_extra = fig7.add_subplot(gs7[1, :])
ax_extra.text(-0.05, 1.02, 'C', transform=ax_extra.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')
# Use a 2-panel comparison (use proper subgridspec)
gs_c_sub = fig7.add_gridspec(1, 2, left=0.04, bottom=0.04, right=0.98, top=0.92, wspace=0.3)
ax_1 = fig7.add_subplot(gs_c_sub[0, 0])
ax_1.scatter(div_df['IVF_vs_vivo'], div_df['PA_vs_vivo'], alpha=0.3, s=3, c='gray')
ax_1.set_xlabel('IVF vs in vivo', fontsize=7)
ax_1.set_ylabel('PA vs in vivo', fontsize=7)
ax_1.set_title('(i) Gene-wise scatter', fontsize=8, fontweight='bold')
ax_1.axhline(0, color='gray', lw=0.3); ax_1.axvline(0, color='gray', lw=0.3)
ax_1.tick_params(labelsize=6)
mx = max(div_df['IVF_vs_vivo'].max(), div_df['PA_vs_vivo'].max())
ax_1.plot([0, mx], [0, mx], 'r--', lw=0.5, alpha=0.5)

ax_2 = fig7.add_subplot(gs_c_sub[0, 1])
ax_2.hist(div_df['PA_vs_vivo'], bins=80, color='#E41A1C', alpha=0.6, label='PA vs in vivo', edgecolor='black', linewidth=0.3)
ax_2.hist(div_df['IVF_vs_vivo'], bins=80, color='#377EB8', alpha=0.6, label='IVF vs in vivo', edgecolor='black', linewidth=0.3)
ax_2.set_xlabel('Divergence', fontsize=7)
ax_2.set_ylabel('Gene count', fontsize=7)
ax_2.set_title('(ii) Distribution comparison', fontsize=8, fontweight='bold')
ax_2.legend(fontsize=6)
ax_2.tick_params(labelsize=6)

ax_extra.set_title('Gene-wise divergence: PA vs in vivo is greater than IVF vs in vivo', fontsize=10, fontweight='bold')
ax_extra.axis('off')

fig7.tight_layout()
fig7.savefig(f'{OUT}/fig7_divergence_volcano.svg', format='svg')
fig7.savefig(f'{OUT}/fig7_divergence_volcano.png', dpi=300)
print(f"  Fig7: {os.path.getsize(f'{OUT}/fig7_divergence_volcano.svg')//1024} KB")

print("\n=== DONE ===")
print(f"Fig6 and Fig7 regenerated with correct data.")