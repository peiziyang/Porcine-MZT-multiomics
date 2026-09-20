#!/usr/bin/env python
"""
Reviewer-eye fix pass for 7 main figures.
Specific issues to fix:
1. Fig1A: enlarge UMAP markers; use actual `stage` column
2. Fig2B/C: separate overlapping gene labels
3. Fig4D: replace with per-factor magnitude breakdown (informative)
4. Fig5A/B: fix y-axis to show bars
5. Fig6C: remove parent panel title (subpanels have their own)
6. Fig7B: convert Ensembl IDs to gene symbols
"""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial',
    'svg.fonttype': 'none',
    'savefig.dpi': 300,
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
})
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import FancyBboxPatch, Circle
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
import os, warnings
warnings.filterwarnings('ignore')

BASE  = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
MF    = os.path.join(BASE, 'm4_mofa')
PERT  = os.path.join(BASE, 'perturbation')
DESEQ = os.path.join(BASE, 'deseq2')
SCEN  = os.path.join(BASE, 'pyscenic_out')
OUT   = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures'

# ── Load data ──
meta    = pd.read_csv(os.path.join(BASE, 'm2_metadata.csv'))
umap_df = pd.read_csv(os.path.join(BASE, 'm2_harmony_umap.csv'))
factor  = pd.read_csv(os.path.join(MF, 'mofa_multiomics_factors_v2.csv'), index_col=0)
r2_df   = pd.read_csv(os.path.join(MF, 'mofa_multiomics_r2_per_view_v2.csv'))
cross_v = pd.read_csv(os.path.join(MF, 'factor_cross_view_v2.csv'))
wrna    = pd.read_csv(os.path.join(MF, 'mofa_multiomics_weights_RNA_v2.csv'), index_col=0)
wmeth   = pd.read_csv(os.path.join(MF, 'mofa_multiomics_weights_METH_v2.csv'), index_col=0)
ortho   = pd.read_csv(os.path.join(MF, 'cross_species_orthologs.csv'))
stg_m   = pd.read_csv(os.path.join(MF, 'mzt_stage_mean_projection.csv'), index_col=0)
tproj   = pd.read_csv(os.path.join(MF, 'mzt_trajectory_projection.csv'), index_col=0)
pert_v  = pd.read_csv(os.path.join(PERT, 'perturbation_vectors.csv'))
pert_pa = pd.read_csv(os.path.join(PERT, 'perturbation_vs_pa.csv'))
div_df  = pd.read_csv(os.path.join(BASE, 'm5_gene_divergence_v2.csv'))
met     = pd.read_csv(os.path.join(BASE, 'm6_metabolism_ratio.csv'))
m3_de   = pd.read_csv(os.path.join(BASE, 'm3_ivf_vs_pa_results.csv'))
# For symbol mapping
gene_id_map = pd.read_csv(os.path.join(DESEQ, 'm15_deseq2_overall_IVFvsPA_stageAdj.csv'))
geneid2sym = dict(zip(gene_id_map['gene'], gene_id_map['gene_symbol']))

# Merge UMAP into metadata
meta_umap = meta.copy()
meta_umap['UMAP1'] = umap_df['UMAP1_harmony'].values
meta_umap['UMAP2'] = umap_df['UMAP2_harmony'].values

def panel_label(ax, label, x=-0.04, y=1.04, fontsize=13):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize,
            fontweight='bold', va='bottom', ha='left')

FCOLS = ['F1','F2','F3','F4','F5','F6','F7']

# ═══════════════════════════════════════════════════════════════════
# FIGURE 1: Fixed UMAP with proper stage coloring
# ═══════════════════════════════════════════════════════════════════
print("=== Fig1: Fixed ===")
fig1 = plt.figure(figsize=(14, 18))
gs1 = fig1.add_gridspec(3, 2, hspace=0.35, wspace=0.3,
                        top=0.96, bottom=0.04, left=0.06, right=0.98)

# ── A: Atlas UMAP (fixed: bigger dots, proper stages) ──
ax = fig1.add_subplot(gs1[0, 0])
panel_label(ax, 'A')
# Use stage column; build proper sorting
stage_list = sorted(meta_umap['stage'].dropna().unique(),
                   key=lambda s: ({'GV': -2, 'MII': -1}.get(s, 99*10)
                                  if not s.startswith('E') and not s.endswith('C')
                                  else int(s.replace('E','').replace('C','')) if s.startswith('E')
                                  else (int(s[0]) * 10 + 1) if s.endswith('C')
                                  else 999))
# Just show all unique stages with unique colors
unique_stages = sorted(meta_umap['stage'].dropna().unique())
cmap_stages = plt.cm.tab20(np.linspace(0, 1, len(unique_stages)))
for i, s in enumerate(unique_stages):
    sub = meta_umap[meta_umap['stage'] == s]
    ax.scatter(sub['UMAP1'], sub['UMAP2'], s=8, c=[cmap_stages[i]],
              alpha=0.7, label=f'{s} (n={len(sub)})', rasterized=True, edgecolors='none')
ax.set_title('Porcine preimplantation atlas (1,955 cells)', fontsize=11, fontweight='bold')
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=5.5, ncol=1, markerscale=1.5, frameon=False)
ax.tick_params(labelsize=7)

# ── B: Factor scores heatmap ──
ax = fig1.add_subplot(gs1[0, 1])
panel_label(ax, 'B')
fmat = factor[FCOLS].values.T
im = ax.imshow(fmat, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
ax.set_yticks(range(7)); ax.set_yticklabels(FCOLS, fontsize=8)
ax.set_xticks([]); ax.set_xlabel('32 GV oocytes', fontsize=8)
ax.set_title('Factor scores (z-scaled)', fontsize=10, fontweight='bold')
fig1.colorbar(im, ax=ax, shrink=0.8, aspect=20).ax.tick_params(labelsize=7)

# ── C: Per-view variance ──
ax = fig1.add_subplot(gs1[1, 0])
panel_label(ax, 'C')
x = np.arange(7); w = 0.35
r2p = r2_df.pivot(index='factor', columns='view', values='r2_pct').fillna(0)
rna_vals = [r2p.loc[f, 'RNA'] if f in r2p.index else 0 for f in FCOLS]
meth_vals = [r2p.loc[f, 'METH'] if f in r2p.index else 0 for f in FCOLS]
ax.bar(x - w/2, meth_vals, w, color='#FF7F00', edgecolor='black', lw=0.5, label='METH')
ax.bar(x + w/2, rna_vals, w, color='#377EB8', edgecolor='black', lw=0.5, label='RNA')
for i in range(7):
    if meth_vals[i] > 0.5:
        ax.text(i - w/2, meth_vals[i] + 1.5, f'{meth_vals[i]:.1f}%', ha='center', fontsize=6)
    if rna_vals[i] > 0.5:
        ax.text(i + w/2, rna_vals[i] + 1.5, f'{rna_vals[i]:.1f}%', ha='center', fontsize=6)
ax.set_xticks(x); ax.set_xticklabels(FCOLS, fontsize=8)
ax.set_ylabel('Variance contribution (%)', fontsize=9)
ax.set_title('Per-view variance contribution', fontsize=10, fontweight='bold')
ax.legend(fontsize=8, loc='upper right')
ax.set_ylim(0, max(max(rna_vals), max(meth_vals)) * 1.3)

# ── D: Cross-view correlation ──
ax = fig1.add_subplot(gs1[1, 1])
panel_label(ax, 'D')
r_vals = cross_v['rna_meth_corr_r'].values
colors_d = ['#E41A1C' if r > 0.5 else '#377EB8' if r < -0.3 else '#999999' for r in r_vals]
ax.bar(FCOLS, r_vals, color=colors_d, edgecolor='black', lw=0.5, width=0.6)
ax.axhline(0, color='gray', lw=0.8)
for i, r in enumerate(r_vals):
    ax.text(i, r + 0.04 if r >= 0 else r - 0.06,
            f'r={r:+.2f}', ha='center', fontsize=7, fontweight='bold')
ax.set_ylabel('RNA-METH weight Pearson r', fontsize=9)
ax.set_title('Cross-view weight correlation', fontsize=10, fontweight='bold')
ax.set_ylim(-0.15, 1.05)

# ── E: Factor scores by donor ──
ax = fig1.add_subplot(gs1[2, 0])
panel_label(ax, 'E')
donor_map = meta.set_index('cell')['donor'].to_dict()
donor_colors = {'AF1': '#E41A1C', 'AF2': '#377EB8', 'AF3': '#4DAF4A',
                'AF4': '#984EA3', 'AF5': '#FF7F00'}
offsets = np.linspace(-0.2, 0.2, 5)
for fi in range(7):
    for di, donor in enumerate(['AF1','AF2','AF3','AF4','AF5']):
        vals = []
        for c in factor.index:
            d = donor_map.get(c, 'unknown')
            if d == donor:
                vals.append(factor.loc[c, FCOLS[fi]])
        if vals:
            jitter = np.random.RandomState(42+fi*5+di).normal(0, 0.02, len(vals))
            ax.scatter(np.full(len(vals), fi+1+offsets[di]) + jitter, vals,
                      c=donor_colors[donor], s=10, alpha=0.6, edgecolors='none',
                      label=donor if fi == 0 else '')
ax.set_xticks(range(1, 8)); ax.set_xticklabels(FCOLS, fontsize=8)
ax.axhline(0, color='gray', ls='--', lw=0.5)
ax.set_xlabel('Factor', fontsize=9); ax.set_ylabel('Score', fontsize=9)
ax.set_title('Factor scores by donor', fontsize=10, fontweight='bold')
ax.legend(fontsize=6, loc='upper right', ncol=5, markerscale=0.8)
ax.tick_params(labelsize=7)

# ── F: F1 RNA vs METH weights ──
ax = fig1.add_subplot(gs1[2, 1])
panel_label(ax, 'F')
r1, p1 = pearsonr(wrna['F1'].values, wmeth['F1'].values)
ax.scatter(wrna['F1'].values, wmeth['F1'].values, c='gray', alpha=0.2, s=3)
shared_idx = np.argsort(np.abs(wrna['F1'].values) + np.abs(wmeth['F1'].values))[-8:]
for idx in shared_idx:
    ax.annotate(wrna.index[idx], (wrna['F1'].iloc[idx], wmeth['F1'].iloc[idx]),
               fontsize=6, fontweight='bold', color='#E41A1C')
ax.set_xlabel('F1 RNA weight', fontsize=9); ax.set_ylabel('F1 METH weight', fontsize=9)
ax.set_title(f'F1: RNA vs METH weights (r={r1:.3f})', fontsize=10, fontweight='bold')
ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
ax.tick_params(labelsize=7)

fig1.savefig(os.path.join(OUT, 'fig1_mofa.svg'), format='svg', bbox_inches='tight')
fig1.savefig(os.path.join(OUT, 'fig1_mofa.png'), dpi=300, bbox_inches='tight')
plt.close(fig1)
print("  -> fig1 OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 2: Fixed label overlap
# ═══════════════════════════════════════════════════════════════════
print("=== Fig2: Fixed ===")
fig2, axes2 = plt.subplots(2, 2, figsize=(12, 10))
axes2 = axes2.flatten()

# A
ax = axes2[0]; panel_label(ax, 'A')
pig_only = ortho['pig_gene'].nunique()
in_human = ortho[ortho['in_human']].shape[0]
in_mouse = ortho[ortho['in_mouse']].shape[0]
in_both  = ortho[ortho['in_human'] & ortho['in_mouse']].shape[0]
human_only = in_human - in_both
mouse_only = in_mouse - in_both
pig_excl = pig_only - in_human - in_mouse + in_both
cats = ['Pig+Human\n+Mouse', 'Pig+Human\nonly', 'Pig+Mouse\nonly', 'Pig only']
vals = [in_both, human_only, mouse_only, pig_excl]
cc = ['#377EB8', '#4DAF4A', '#FF7F00', '#999999']
ax.bar(cats, vals, color=cc, edgecolor='black', lw=0.5, width=0.6)
for i, v in enumerate(vals):
    ax.text(i, v + 0.5, str(v), ha='center', fontsize=9, fontweight='bold')
ax.set_ylabel('Number of F1 genes', fontsize=9)
ax.set_title(f'F1 maternal blueprint genes (n={pig_only})', fontsize=10, fontweight='bold')

# B: Pig vs Human with separated labels
ax = axes2[1]; panel_label(ax, 'B')
hu_genes = ortho[ortho['in_human']].dropna(subset=['pig_GV_mean','human_mean'])
r_hu, p_hu = spearmanr(hu_genes['pig_GV_mean'], hu_genes['human_mean'])
ax.scatter(hu_genes['pig_GV_mean'], hu_genes['human_mean'],
           c='#377EB8', alpha=0.6, s=25, edgecolors='black', lw=0.3)
# Label key genes with offsets to avoid overlap
label_offsets_b = {'DNMT1': (5, 5), 'ZP3': (-30, 5), 'ZP4': (-30, -10), 'GDF9': (5, 5), 'RARRES1': (5, -10)}
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    sub = hu_genes[hu_genes['pig_gene'] == g]
    if len(sub):
        off = label_offsets_b.get(g, (5, 5))
        ax.annotate(g, (sub['pig_GV_mean'].values[0], sub['human_mean'].values[0]),
                   fontsize=8, fontweight='bold',
                   xytext=(sub['pig_GV_mean'].values[0]+off[0]*0.005*hu_genes['pig_GV_mean'].max(),
                          sub['human_mean'].values[0]+off[1]*0.005*hu_genes['human_mean'].max()),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
ax.set_xlabel('Pig GV mean expression', fontsize=9)
ax.set_ylabel('Human preimplantation mean expr.', fontsize=9)
ax.set_title(f'Pig vs Human (rho={r_hu:.3f}, p={p_hu:.2e})', fontsize=10, fontweight='bold')
ax.tick_params(labelsize=7)

# C: Pig vs Mouse with separated labels
ax = axes2[2]; panel_label(ax, 'C')
mo_genes = ortho[ortho['in_mouse']].dropna(subset=['pig_GV_mean','mouse_mean'])
r_mo, p_mo = spearmanr(mo_genes['pig_GV_mean'], mo_genes['mouse_mean'])
ax.scatter(mo_genes['pig_GV_mean'], mo_genes['mouse_mean'],
           c='#FF7F00', alpha=0.6, s=25, edgecolors='black', lw=0.3)
label_offsets_c = {'DNMT1': (5, 5), 'ZP3': (5, 5), 'ZP4': (-30, -10), 'GDF9': (5, 5), 'RARRES1': (-30, 5)}
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    sub = mo_genes[mo_genes['pig_gene'] == g]
    if len(sub):
        off = label_offsets_c.get(g, (5, 5))
        ax.annotate(g, (sub['pig_GV_mean'].values[0], sub['mouse_mean'].values[0]),
                   fontsize=8, fontweight='bold',
                   xytext=(sub['pig_GV_mean'].values[0]+off[0]*0.005*mo_genes['pig_GV_mean'].max(),
                          sub['mouse_mean'].values[0]+off[1]*0.005*mo_genes['mouse_mean'].max()),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
ax.set_xlabel('Pig GV mean expression', fontsize=9)
ax.set_ylabel('Mouse preimplantation mean expr.', fontsize=9)
ax.set_title(f'Pig vs Mouse (rho={r_mo:.3f}, p={p_mo:.2e})', fontsize=10, fontweight='bold')
ax.tick_params(labelsize=7)

# D
ax = axes2[3]; panel_label(ax, 'D'); ax.axis('off')
summary = (
    f'Cross-species MZT conservation\n\n'
    f'F1 maternal blueprint genes: {pig_only}\n'
    f'  Conserved in human: {in_human}/{pig_only} ({100*in_human//pig_only}%)\n'
    f'  Conserved in mouse: {in_mouse}/{pig_only} ({100*in_mouse//pig_only}%)\n'
    f'  Conserved in both:  {in_both}/{pig_only} ({100*in_both//pig_only}%)\n\n'
    f'Pig-Human Spearman rho = {r_hu:.3f}, p = {p_hu:.2e}\n'
    f'Pig-Mouse Spearman rho = {r_mo:.3f}, p = {p_mo:.2e}\n\n'
    f'Core conserved across 3 species:\n'
    f'  DNMT1, ZP3, ZP4, GDF9, RARRES1\n\n'
    f'Data: GSE44183 (human+mouse preimplantation)\n'
    f'Matching: orthologous gene symbols'
)
ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=9, va='top',
        family='monospace', bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='gray'))

fig2.tight_layout()
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.svg'), format='svg', bbox_inches='tight')
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.png'), dpi=300, bbox_inches='tight')
plt.close(fig2)
print("  -> fig2 OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 3: keep same (looked good)
# ═══════════════════════════════════════════════════════════════════
print("=== Fig3: rebuild ===")
fig3, axes3 = plt.subplots(3, 2, figsize=(14, 16))
axes3 = axes3.flatten()
stage_order = ['E0','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10']
stage_num = np.arange(len(stage_order))

# A
ax = axes3[0]; panel_label(ax, 'A')
f3_vals = [stg_m.loc[s, 'F3'] if s in stg_m.index else np.nan for s in stage_order]
r_f3, p_f3 = spearmanr(stage_num[~np.isnan(f3_vals)], np.array(f3_vals)[~np.isnan(f3_vals)])
ax.plot(stage_num, f3_vals, '-o', color='#E41A1C', lw=2, markersize=6)
ax.fill_between(stage_num, np.array(f3_vals)-np.std(f3_vals)*0.5,
                np.array(f3_vals)+np.std(f3_vals)*0.5, alpha=0.15, color='#E41A1C')
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=7)
ax.set_ylabel('Projected score', fontsize=9)
ax.set_title(f'F3: Maternal clearance (rho={r_f3:+.3f}, p={p_f3:.1e})', fontsize=10, fontweight='bold', color='#E41A1C')
ax.tick_params(labelsize=7)

# B
ax = axes3[1]; panel_label(ax, 'B')
f4_vals = [stg_m.loc[s, 'F4'] if s in stg_m.index else np.nan for s in stage_order]
r_f4, p_f4 = spearmanr(stage_num[~np.isnan(f4_vals)], np.array(f4_vals)[~np.isnan(f4_vals)])
ax.plot(stage_num, f4_vals, '-o', color='#4DAF4A', lw=2, markersize=6)
ax.fill_between(stage_num, np.array(f4_vals)-np.std(f4_vals)*0.5,
                np.array(f4_vals)+np.std(f4_vals)*0.5, alpha=0.15, color='#4DAF4A')
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=7)
ax.set_ylabel('Projected score', fontsize=9)
ax.set_title(f'F4: Zygotic metabolic ({r_f4:+.3f}, p={p_f4:.2e})', fontsize=10, fontweight='bold', color='#4DAF4A')
ax.tick_params(labelsize=7)

# C
ax = axes3[2]; panel_label(ax, 'C')
f6_vals = [stg_m.loc[s, 'F6'] if s in stg_m.index else np.nan for s in stage_order]
r_f6, p_f6 = spearmanr(stage_num[~np.isnan(f6_vals)], np.array(f6_vals)[~np.isnan(f6_vals)])
ax.plot(stage_num, f6_vals, '-o', color='#FF7F00', lw=2, markersize=6)
ax.fill_between(stage_num, np.array(f6_vals)-np.std(f6_vals)*0.5,
                np.array(f6_vals)+np.std(f6_vals)*0.5, alpha=0.15, color='#FF7F00')
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=7)
ax.set_ylabel('Projected score', fontsize=9)
ax.set_title(f'F6: Zygotic translational ({r_f6:+.3f}, p={p_f6:.2e})', fontsize=10, fontweight='bold', color='#FF7F00')
ax.tick_params(labelsize=7)

# D
ax = axes3[3]; panel_label(ax, 'D')
f1_vals = [stg_m.loc[s, 'F1'] if s in stg_m.index else np.nan for s in stage_order]
r_f1, p_f1 = spearmanr(stage_num[~np.isnan(f1_vals)], np.array(f1_vals)[~np.isnan(f1_vals)])
ax.plot(stage_num, f1_vals, '-o', color='#377EB8', lw=2, markersize=6)
ax.fill_between(stage_num, np.array(f1_vals)-np.std(f1_vals)*0.5,
                np.array(f1_vals)+np.std(f1_vals)*0.5, alpha=0.15, color='#377EB8')
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=7)
ax.set_ylabel('Projected score', fontsize=9)
ax.set_title(f'F1: Maternal blueprint (STABLE) (rho={r_f1:+.3f}, p={p_f1:.2f})', fontsize=10, fontweight='bold', color='#377EB8')
ax.tick_params(labelsize=7)

# E
ax = axes3[4]; panel_label(ax, 'E')
for fc, col in [('F1','#377EB8'),('F3','#E41A1C'),('F4','#4DAF4A'),('F6','#FF7F00')]:
    vals = [stg_m.loc[s, fc] if s in stg_m.index else np.nan for s in stage_order]
    vals_arr = np.array(vals)
    if vals_arr.std() > 0:
        vals_norm = (vals_arr - vals_arr.mean()) / vals_arr.std()
    else:
        vals_norm = vals_arr
    ax.plot(stage_num, vals_norm, '-o', lw=1.5, markersize=4, color=col, label=fc, alpha=0.8)
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=7)
ax.set_ylabel('Normalized score', fontsize=9)
ax.set_title('All factor trajectories (z-normalized)', fontsize=10, fontweight='bold')
ax.axhline(0, color='gray', lw=0.5); ax.legend(fontsize=7); ax.tick_params(labelsize=7)

# F
ax = axes3[5]; panel_label(ax, 'F')
stages_in_data = sorted(tproj['stage'].unique(), key=lambda s: int(s[1:]))
stage_idx_map = {s: i for i, s in enumerate(stages_in_data)}
tproj_sorted = tproj.copy()
for s in stages_in_data:
    sub = tproj_sorted[tproj_sorted['stage'] == s]
    if len(sub) > 0:
        ax.violinplot(sub['F1'].values, positions=[stage_idx_map[s]],
                     showmeans=True, showmedians=False, widths=0.7)
ax.set_xticks(range(len(stages_in_data)))
ax.set_xticklabels(stages_in_data, fontsize=6, rotation=45, ha='right')
ax.set_ylabel('F1 projected score', fontsize=9)
ax.set_title('Per-cell F1 projection (n=832 E0-E10 cells)', fontsize=10, fontweight='bold')
ax.tick_params(labelsize=7)

fig3.tight_layout()
fig3.savefig(os.path.join(OUT, 'fig3_mzt_trajectory.svg'), format='svg', bbox_inches='tight')
fig3.savefig(os.path.join(OUT, 'fig3_mzt_trajectory.png'), dpi=300, bbox_inches='tight')
plt.close(fig3)
print("  -> fig3 OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 4: replace Panel D with per-factor magnitude breakdown
# ═══════════════════════════════════════════════════════════════════
print("=== Fig4: Fix Panel D ===")
fig4, axes4 = plt.subplots(2, 2, figsize=(12, 10))
axes4 = axes4.flatten()

# A
ax = axes4[0]; panel_label(ax, 'A')
dnmt1 = pert_v[pert_v['tf'] == 'DNMT1']
for label, color, marker in [('100%_KD','#E41A1C','o'),('50%_KD','#FF7F00','s'),
                              ('50%_OE','#377EB8','^'),('100%_OE','#4DAF4A','D')]:
    sub = dnmt1[dnmt1['perturbation'] == label]
    ax.scatter(sub['shift_F1'], sub['shift_F2'], c=color, marker=marker, s=20, alpha=0.7,
              label=label, edgecolors='black', lw=0.3)
ax.axhline(0, color='gray', ls='--', lw=0.5); ax.axvline(0, color='gray', ls='--', lw=0.5)
ax.set_xlabel('Shift in F1', fontsize=9); ax.set_ylabel('Shift in F2', fontsize=9)
ax.set_title('DNMT1 perturbation in MOFA+ space', fontsize=10, fontweight='bold')
ax.legend(fontsize=7); ax.tick_params(labelsize=7)

# B
ax = axes4[1]; panel_label(ax, 'B')
kd_cos = pert_pa[pert_pa['perturbation'].str.contains('KD')]
kd_avg = kd_cos.groupby('tf')['cos_sim'].mean().sort_values()
colors_b = ['#E41A1C' if v > 0.3 else '#FF7F00' if v > 0 else '#999999' for v in kd_avg.values]
ax.barh(range(len(kd_avg)), kd_avg.values, color=colors_b, edgecolor='black', lw=0.5, height=0.6)
ax.set_yticks(range(len(kd_avg))); ax.set_yticklabels(kd_avg.index, fontsize=8)
ax.axvline(0, color='gray', lw=0.8)
for i, v in enumerate(kd_avg.values):
    ax.text(v + 0.02, i, f'{v:+.3f}', fontsize=7, va='center', fontweight='bold')
ax.set_xlabel('Cosine similarity with PA deviation', fontsize=9)
ax.set_title('TF KD alignment with PA deviation', fontsize=10, fontweight='bold')
ax.tick_params(labelsize=7)

# C
ax = axes4[2]; panel_label(ax, 'C')
kds = pert_v[pert_v['perturbation'] == '100%_KD']
tfs = sorted(kds['tf'].unique())
hm = np.zeros((len(tfs), 4))
for ti, tf in enumerate(tfs):
    sub = kds[kds['tf'] == tf]
    for fi in range(4):
        hm[ti, fi] = sub[f'shift_{FCOLS[fi]}'].mean()
im = ax.imshow(hm.T, aspect='auto', cmap='RdBu_r', vmin=-0.025, vmax=0.025)
ax.set_xticks(range(len(tfs))); ax.set_xticklabels(tfs, fontsize=7, rotation=45, ha='right')
ax.set_yticks(range(4)); ax.set_yticklabels(FCOLS[:4], fontsize=8)
for i in range(4):
    for j in range(len(tfs)):
        ax.text(j, i, f'{hm[j,i]:.4f}', ha='center', va='center', fontsize=5.5,
               color='white' if abs(hm[j,i]) > 0.015 else 'black')
ax.set_title('100% KD: per-factor mean shift (n=32 cells)', fontsize=10, fontweight='bold')
fig4.colorbar(im, ax=ax, shrink=0.7, aspect=20).ax.tick_params(labelsize=7)

# D: Per-FACTOR magnitude breakdown (KD) - shows which factor each TF affects most
ax = axes4[3]; panel_label(ax, 'D')
tfs_show = ['DNMT1','ATF3','POU5F1','ESRRA','TFAP2C']
# Compute mean absolute shift per factor per TF
abs_shift = np.zeros((len(tfs_show), 4))
for ti, tf in enumerate(tfs_show):
    sub_kd = pert_v[(pert_v['tf']==tf) & (pert_v['perturbation']=='100%_KD')]
    for fi in range(4):
        abs_shift[ti, fi] = sub_kd[f'shift_{FCOLS[fi]}'].abs().mean()
x_pos = np.arange(len(tfs_show))
bottom = np.zeros(len(tfs_show))
factor_colors = {'F1':'#E41A1C', 'F2':'#377EB8', 'F3':'#4DAF4A', 'F4':'#FF7F00'}
for fi in range(4):
    ax.bar(x_pos, abs_shift[:, fi], 0.6, bottom=bottom, color=factor_colors[FCOLS[fi]],
          edgecolor='black', lw=0.4, label=FCOLS[fi])
    bottom += abs_shift[:, fi]
ax.set_xticks(x_pos); ax.set_xticklabels(tfs_show, fontsize=8)
ax.set_ylabel('Sum of |shift| per factor', fontsize=9)
ax.set_title('Per-factor |shift| contribution (100% KD)', fontsize=10, fontweight='bold')
ax.legend(fontsize=7, loc='upper right', ncol=4)
ax.tick_params(labelsize=7)

fig4.tight_layout()
fig4.savefig(os.path.join(OUT, 'fig4_perturbation.svg'), format='svg', bbox_inches='tight')
fig4.savefig(os.path.join(OUT, 'fig4_perturbation.png'), dpi=300, bbox_inches='tight')
plt.close(fig4)
print("  -> fig4 OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 5: Critical fix — set ylim properly to show bars
# ═══════════════════════════════════════════════════════════════════
print("=== Fig5: Fix ylim ===")
pa_all = []
# Load both per-stage and stage-adjusted files
for f in sorted(os.listdir(DESEQ)):
    if not f.startswith('m15_deseq2_'):
        continue
    if not (f.endswith('_IVFvsPA.csv') or f.endswith('_IVFvsPA_stageAdj.csv')):
        continue
    # Extract stage name
    name = f.replace('m15_deseq2_', '').replace('_IVFvsPA.csv', '').replace('_stageAdj.csv', '')
    df = pd.read_csv(os.path.join(DESEQ, f))
    if 'gene_symbol' not in df.columns:
        continue
    s2i = {}
    for i, row in df.iterrows():
        s = str(row['gene_symbol'])
        if not s.startswith('ENSSSCG'):
            s2i[s] = i
    for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2',
              'IDH2','PKM','EXOSC9','MDH1','GNL3','DPPA5','NANOG']:
        if g in s2i:
            r = df.iloc[s2i[g]]
            pa_all.append({'gene': g, 'stage': name,
                          'log2FC': r['log2FoldChange'],
                          'padj': r['padj']})
print(f'Loaded {len(pa_all)} PA records from DESeq2 files')
print(f'Stages found: {sorted(set(r["stage"] for r in pa_all))}')
pa_df = pd.DataFrame(pa_all)
overall = pa_df[pa_df['stage'].str.startswith('overall')]

fig5 = plt.figure(figsize=(14, 16))
gs5 = fig5.add_gridspec(3, 2, hspace=0.4, wspace=0.3)

# A: Layer 1
ax = fig5.add_subplot(gs5[0, 0]); panel_label(ax, 'A')
l1_genes = ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2']
l1_val = []; l1_padj = []; l1_names = []
for g in l1_genes:
    sub = overall[overall['gene'] == g]
    if len(sub) and not np.isnan(sub.iloc[0]['log2FC']):
        l1_val.append(sub.iloc[0]['log2FC']); l1_padj.append(sub.iloc[0]['padj']); l1_names.append(g)
x1 = np.arange(len(l1_names))
c1 = ['#E41A1C' if v > 0.5 else '#FF7F00' if v > -0.5 else '#377EB8' for v in l1_val]
bars = ax.bar(x1, l1_val, color=c1, edgecolor='black', lw=0.5, width=0.6)
ax.axhline(0, color='gray', lw=0.8)
for i, (v, p) in enumerate(zip(l1_val, l1_padj)):
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    yoff = 0.15 if v >= 0 else -0.35
    ax.text(i, v + yoff, f'{v:+.2f} {sig}', ha='center', fontsize=7.5, fontweight='bold')
ax.set_xticks(x1); ax.set_xticklabels(l1_names, fontsize=8)
ax.set_ylabel('log2 Fold Change (PA / IVF)', fontsize=9)
ax.set_title('Layer 1: Maternal blueprint -- MAINTAINED', fontsize=10, fontweight='bold', color='#E41A1C')
# KEY FIX: explicit ylim
ax.set_ylim(-1.5, 2.5)

# B: Layer 3
ax = fig5.add_subplot(gs5[0, 1]); panel_label(ax, 'B')
l3_genes = ['IDH2','PKM','EXOSC9','MDH1','GNL3']
l3_val = []; l3_padj = []; l3_names = []
for g in l3_genes:
    sub = overall[overall['gene'] == g]
    if len(sub) and not np.isnan(sub.iloc[0]['log2FC']):
        l3_val.append(sub.iloc[0]['log2FC']); l3_padj.append(sub.iloc[0]['padj']); l3_names.append(g)
x3 = np.arange(len(l3_names))
c3 = ['#377EB8' if v < -2 else '#999999' for v in l3_val]
ax.bar(x3, l3_val, color=c3, edgecolor='black', lw=0.5, width=0.6)
ax.axhline(0, color='gray', lw=0.8)
for i, (v, p) in enumerate(zip(l3_val, l3_padj)):
    sig = '***' if p < 1e-6 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    ax.text(i, 0.4, f'{v:+.2f} {sig}', ha='center', fontsize=7.5, fontweight='bold',
           color='#377EB8')
ax.set_xticks(x3); ax.set_xticklabels(l3_names, fontsize=8)
ax.set_ylabel('log2 Fold Change (PA / IVF)', fontsize=9)
ax.set_title('Layer 3: Zygotic execution -- ATTENUATED', fontsize=10, fontweight='bold', color='#377EB8')
ax.set_ylim(-10, 1)

# C: Per-stage heatmap (rebuilt)
ax = fig5.add_subplot(gs5[1, 0]); panel_label(ax, 'C')
stg = ['1cell','2cell','4cell','8cell']
gpl = ['DNMT1','ZP3','GDF9','IDH2','PKM','EXOSC9','GNL3','DPPA5','NANOG']
hm = np.full((len(gpl), len(stg)), np.nan)
for gi, g in enumerate(gpl):
    for si, s in enumerate(stg):
        sub = pa_df[(pa_df['gene'] == g) & (pa_df['stage'] == s)]
        if len(sub) and not np.isnan(sub.iloc[0]['log2FC']):
            hm[gi, si] = sub.iloc[0]['log2FC']
hm_masked = np.ma.masked_where(np.isnan(hm), hm)
im = ax.imshow(hm_masked, aspect='auto', cmap='RdBu_r', vmin=-10, vmax=10)
ax.set_xticks(np.arange(len(stg))); ax.set_xticklabels(stg, fontsize=8)
ax.set_yticks(np.arange(len(gpl))); ax.set_yticklabels(gpl, fontsize=8)
for i in range(len(gpl)):
    for j in range(len(stg)):
        v = hm[i, j]
        if not np.isnan(v):
            ax.text(j, i, f'{v:+.1f}', ha='center', va='center', fontsize=6.5,
                   fontweight='bold', color='white' if abs(v) > 5 else 'black')
ax.set_title('Per-stage log2FC (PA vs IVF)', fontsize=10, fontweight='bold')
fig5.colorbar(im, ax=ax, shrink=0.75, aspect=20).ax.tick_params(labelsize=7)

# D: Summary
ax = fig5.add_subplot(gs5[1, 1]); panel_label(ax, 'D')
ax.set_xlim(0, 12); ax.set_ylim(0, 12); ax.axis('off')
ys = 11.5
for title, col, rows in [
    ('Layer 1 -- Maternal Blueprint (MAINTAINED)', '#E41A1C',
     [('DNMT1','+1.30','0.006','**'),
      ('ZP3','+0.21','0.52','ns'),
      ('ZP4','-0.75','0.08','ns'),
      ('GDF9','+0.62','0.11','ns'),
      ('RARRES1','-0.30','0.45','ns'),
      ('EIF4G2','+0.15','0.72','ns')]),
    ('Layer 3 -- Zygotic Execution (ATTENUATED)', '#377EB8',
     [('IDH2 (TCA)','-8.03','3.5e-9','***'),
      ('PKM (glycolysis)','-7.64','5.3e-12','***'),
      ('EXOSC9 (ribosome)','-7.73','3.9e-13','***'),
      ('MDH1 (TCA)','-4.58','2.4e-4','**'),
      ('GNL3 (nucleolar)','-3.06','0.014','*')]),
    ('ZGA Markers (ABERRANTLY ACTIVATED)', '#4DAF4A',
     [('DPPA5','+2.14','0.008','**'),
      ('NANOG','+3.91','1.4e-4','***')])]:
    ax.text(0.2, ys, title, fontsize=9, fontweight='bold', color=col, va='center')
    ys -= 0.5
    for gene, lfc, pval, sig in rows:
        ax.text(0.4, ys, gene, fontsize=7.5, va='center')
        cl = '#E41A1C' if lfc.startswith('+') and float(lfc) > 1 else '#377EB8' if lfc.startswith('-') and float(lfc) < -1 else '#666'
        ax.text(4.5, ys, f'log2FC = {lfc}', fontsize=7.5, va='center', fontweight='bold', color=cl)
        ax.text(7.5, ys, f'padj = {pval}', fontsize=7.5, va='center', color='#666')
        ax.text(10, ys, f'[{sig}]', fontsize=7.5, va='center', color='#CC0000', fontweight='bold')
        ys -= 0.38
    ys -= 0.3
ax.text(0.3, ys, '42-sample bulk RNA-seq (SRP301735) | DESeq2: stage-adjusted design', fontsize=7, style='italic', color='gray')

# E: 3-layer model (fixed - no overlap with arrow)
ax = fig5.add_subplot(gs5[2, :]); panel_label(ax, 'E')
ax.set_xlim(0, 18); ax.set_ylim(0, 10); ax.axis('off')

ax.text(4, 9.7, 'NORMAL MZT', ha='center', fontsize=12, fontweight='bold', color='#333')
gv_c = Circle((4, 8.5), 0.4, facecolor='#FFF3CD', edgecolor='#FF7F00', lw=2)
ax.add_patch(gv_c)
ax.text(4, 8.5, 'GV', ha='center', va='center', fontsize=8, fontweight='bold')

for y, label, fc, ec, genes, note in [
    (7.0, 'Layer 1: Maternal Blueprint', '#E5F0FF', '#377EB8',
     'DNMT1 . ZP3 . ZP4 . GDF9 . RARRES1', 'RNA+METH co-regulated r=+0.84'),
    (5.3, 'Layer 2: Maternal Clearance', '#FFE5E5', '#E41A1C',
     'ZP2 . SYCN . PARP12', 'rho=-0.90  p=2e-4'),
    (3.6, 'Layer 3: Zygotic Activation', '#E5FFE5', '#4DAF4A',
     'IDH2 . PKM . EXOSC9 . NOP9 . GNL3', 'Metabolism + Translation'),
]:
    fb = FancyBboxPatch((0.5, y-0.5), 7, 0.95, boxstyle='round,pad=0.08',
                        facecolor=fc, edgecolor=ec, lw=1.5, alpha=0.7)
    ax.add_patch(fb)
    ax.text(4, y+0.2, label, fontsize=9, fontweight='bold', color=ec, ha='center')
    ax.text(4, y-0.05, genes, fontsize=7.5, ha='center', color=ec)
    ax.text(4, y-0.32, note, fontsize=6.5, ha='center', color='#666', style='italic')

ax.annotate('', xy=(4, 2.8), xytext=(4, 9.0),
           arrowprops=dict(arrowstyle='->', color='#666', lw=1.5))

ax.text(13.5, 9.7, 'PA EMBRYOS', ha='center', fontsize=12, fontweight='bold', color='#E41A1C')
pa_data_e = [
    (7.0, 'Layer 1: MAINTAINED', '#E5F0FF', '#377EB8', 'DNMT1 HIGHER in PA\nlog2FC=+1.30, padj=0.006'),
    (5.3, 'Layer 2: VARIABLE', '#FFF3CD', '#FF7F00', 'ZP2/4 inconsistent\nacross stages'),
    (3.6, 'Layer 3: ATTENUATED', '#FFE5E5', '#E41A1C', 'IDH2 log2FC=-8.03 ***\nPKM log2FC=-7.64 ***\nEXOSC9 log2FC=-7.73 ***'),
]
for y, label, fc, ec, detail in pa_data_e:
    fb2 = FancyBboxPatch((10, y-0.5), 7, 0.95, boxstyle='round,pad=0.08',
                        facecolor=fc, edgecolor=ec, lw=1.5, alpha=0.7)
    ax.add_patch(fb2)
    ax.text(13.5, y+0.2, label, fontsize=9, fontweight='bold', color=ec, ha='center')
    ax.text(13.5, y-0.2, detail, fontsize=6.5, ha='center', color='#333', va='center')

ax.text(9, 1.5, 'GV oocyte multi-omics (MOFA+) -> Atlas projection -> 42-sample independent validation (SRP301735)',
        ha='center', fontsize=8, style='italic', color='gray')

fig5.savefig(os.path.join(OUT, 'fig5_pa_validation.svg'), format='svg', bbox_inches='tight')
fig5.savefig(os.path.join(OUT, 'fig5_pa_validation.png'), dpi=300, bbox_inches='tight')
plt.close(fig5)
print("  -> fig5 OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 6: Fix Panel C overlap
# ═══════════════════════════════════════════════════════════════════
print("=== Fig6: Fix Panel C ===")
regulon_activity = pd.read_csv(os.path.join(BASE, 'm7_regulon_activity.csv'))
reg_names = regulon_activity['Unnamed: 0'].str.extract(r'^(\w+)')[0].fillna('Unknown')
ra_num_cols = [c for c in regulon_activity.columns
               if c.startswith('E') and c[1:].isdigit()]
regulon_stage = regulon_activity[ra_num_cols].copy()
regulon_stage.index = reg_names.values
ra_stages = sorted(ra_num_cols, key=lambda s: int(s[1:]))

fig6 = plt.figure(figsize=(14, 14))
gs6 = fig6.add_gridspec(3, 1, height_ratios=[1.2, 1, 1.2], hspace=0.5)

# A
ax = fig6.add_subplot(gs6[0]); panel_label(ax, 'A')
top_20 = regulon_stage[ra_stages].mean(axis=1).sort_values(ascending=False).head(20).index.tolist()
hm_data = regulon_stage.loc[top_20, ra_stages].values
im = ax.imshow(hm_data, aspect='auto', cmap='YlOrRd')
ax.set_xticks(range(len(ra_stages)))
ax.set_xticklabels(ra_stages, fontsize=7, rotation=45, ha='right')
ax.set_yticks(range(len(top_20)))
ax.set_yticklabels([str(t)[:20] for t in top_20], fontsize=6.5)
ax.set_title('SCENIC regulon activity: top 20 across E0-E8 stages', fontsize=10, fontweight='bold')
fig6.colorbar(im, ax=ax, shrink=0.8, aspect=20).ax.tick_params(labelsize=7)

# B
ax = fig6.add_subplot(gs6[1]); panel_label(ax, 'B')
key_tfs = ['DNMT1','POU5F1','ESRRA','TFAP2C','GATA4','ATF3','NFE2L3','XBP1']
found_tfs = [t for t in key_tfs if t in regulon_stage.index]
if found_tfs:
    colors_tf = plt.cm.tab10(np.linspace(0, 1, len(found_tfs)))
    for ti, tf in enumerate(found_tfs):
        vals = regulon_stage.loc[tf, ra_stages].values
        ax.plot(range(len(ra_stages)), vals, '-o', lw=1.5, markersize=4,
               color=colors_tf[ti], label=tf)
    ax.set_xticks(range(len(ra_stages)))
    ax.set_xticklabels(ra_stages, fontsize=7, rotation=45, ha='right')
    ax.set_ylabel('Regulon activity', fontsize=9)
    ax.set_title('Key regulon temporal dynamics', fontsize=10, fontweight='bold')
    ax.legend(fontsize=7, ncol=4)
    ax.tick_params(labelsize=7)
else:
    ax.text(0.5, 0.5, 'No key regulons found in data', ha='center', va='center', fontsize=12)
    ax.set_title('Key regulon dynamics', fontsize=10, fontweight='bold')

# C: Metabolism — use proper subplot_spec, no parent title
gs_c = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs6[2], wspace=0.3)
stages_m = met['stage'].values
x_met = np.arange(len(stages_m))

ax1 = fig6.add_subplot(gs_c[0])
panel_label(ax1, 'C1')
ax1.bar(x_met, met['ratio'].values, color='#377EB8', edgecolor='black', lw=0.3, width=0.7)
ax1.axhline(1, color='red', ls='--', lw=0.8, alpha=0.5)
ax1.set_xticks(x_met[::2]); ax1.set_xticklabels(stages_m[::2], fontsize=6, rotation=45, ha='right')
ax1.set_ylabel('Glycolysis / OxPhos', fontsize=8)
ax1.set_title('Glycolysis/OxPhos ratio', fontsize=9, fontweight='bold')
ax1.tick_params(labelsize=6)
ax1.set_ylim(0, max(met['ratio'].values) * 1.3)

ax2 = fig6.add_subplot(gs_c[1])
panel_label(ax2, 'C2')
ax2.bar(x_met, met['gly_mean'].values, color='#E41A1C', edgecolor='black', lw=0.3, width=0.7)
ax2.set_xticks(x_met[::2]); ax2.set_xticklabels(stages_m[::2], fontsize=6, rotation=45, ha='right')
ax2.set_ylabel('Module score', fontsize=8)
ax2.set_title('Glycolysis module score', fontsize=9, fontweight='bold')
ax2.tick_params(labelsize=6)

ax3 = fig6.add_subplot(gs_c[2])
panel_label(ax3, 'C3')
ax3.bar(x_met, met['ox_mean'].values, color='#4DAF4A', edgecolor='black', lw=0.3, width=0.7)
ax3.set_xticks(x_met[::2]); ax3.set_xticklabels(stages_m[::2], fontsize=6, rotation=45, ha='right')
ax3.set_ylabel('Module score', fontsize=8)
ax3.set_title('OxPhos module score', fontsize=9, fontweight='bold')
ax3.tick_params(labelsize=6)

fig6.savefig(os.path.join(OUT, 'fig6_regulon_metabolism.svg'), format='svg', bbox_inches='tight')
fig6.savefig(os.path.join(OUT, 'fig6_regulon_metabolism.png'), dpi=300, bbox_inches='tight')
plt.close(fig6)
print("  -> fig6 OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 7: Map Ensembl IDs to symbols
# ═══════════════════════════════════════════════════════════════════
print("=== Fig7: Map IDs to symbols ===")
fig7, axes7 = plt.subplots(2, 2, figsize=(14, 12))
axes7 = axes7.flatten()

# A: Mean divergence
ax = axes7[0]; panel_label(ax, 'A')
ivf_pa = div_df['IVF_vs_PA'].mean()
ivf_vivo = div_df['IVF_vs_vivo'].mean()
pa_vivo = div_df['PA_vs_vivo'].mean()
ax.bar(['IVF-PA','IVF-InVivo','PA-InVivo'], [ivf_pa, ivf_vivo, pa_vivo],
       color=['#377EB8','#4DAF4A','#E41A1C'], edgecolor='black', lw=0.5, width=0.5)
for i, v in enumerate([ivf_pa, ivf_vivo, pa_vivo]):
    ax.text(i, v + 0.005, f'{v:.3f}', ha='center', fontsize=8, fontweight='bold')
ax.set_ylabel('Mean divergence (8,276 genes)', fontsize=9)
ax.set_title('Cross-condition divergence', fontsize=10, fontweight='bold')

# B: Top 15 PA-divergent genes (with SYMBOL mapping)
ax = axes7[1]; panel_label(ax, 'B')
top15 = div_df.sort_values('PA_vs_vivo', ascending=False).head(15)
# Map Ensembl to symbol
top15_sym = top15.copy()
top15_sym['gene'] = top15_sym['gene'].map(lambda g: geneid2sym.get(g, g))
# If multiple same symbol, add counter
sym_counts = {}
def uniquify(name):
    if name in sym_counts:
        sym_counts[name] += 1
        return f'{name}.{sym_counts[name]}'
    sym_counts[name] = 1
    return name
top15_sym['gene'] = top15_sym['gene'].apply(uniquify)
colors_g = ['#E41A1C' if v > 1.0 else '#FF7F00' if v > 0.5 else '#999' for v in top15['PA_vs_vivo'].values]
ax.barh(range(15), top15['PA_vs_vivo'].values, color=colors_g, edgecolor='black', lw=0.3, height=0.7)
ax.set_yticks(range(15))
ax.set_yticklabels(top15_sym['gene'].values, fontsize=6.5)
ax.invert_yaxis()
ax.set_xlabel('PA vs in vivo divergence', fontsize=9)
ax.set_title('Top 15 PA-divergent genes (mapped to symbols)', fontsize=10, fontweight='bold')
ax.tick_params(labelsize=6.5)

# C: Gene-wise scatter
ax = axes7[2]; panel_label(ax, 'C')
ax.scatter(div_df['IVF_vs_vivo'].values, div_df['PA_vs_vivo'].values,
           alpha=0.15, s=2, c='gray', rasterized=True)
mx = max(div_df['IVF_vs_vivo'].max(), div_df['PA_vs_vivo'].max()) * 1.05
ax.plot([0, mx], [0, mx], 'r--', lw=0.8, alpha=0.5, label='equality')
ax.set_xlabel('IVF vs in vivo divergence', fontsize=9)
ax.set_ylabel('PA vs in vivo divergence', fontsize=9)
ax.set_title('Gene-wise divergence comparison', fontsize=10, fontweight='bold')
ax.legend(fontsize=8); ax.tick_params(labelsize=7)

# D: Volcano
ax = axes7[3]; panel_label(ax, 'D')
m3_de = pd.read_csv(os.path.join(BASE, 'm3_ivf_vs_pa_results.csv'))
ns = m3_de[(m3_de['padj'] >= 0.05) | (m3_de['significant'] != True)]
ax.scatter(ns['log2FC_PA_vs_IVF'].values, -np.log10(ns['pvalue'].clip(1e-300)),
          c='gray', alpha=0.2, s=2, rasterized=True, label='NS')
sig_data = m3_de[(m3_de['padj'] < 0.05) | (m3_de['significant'] == True)]
if len(sig_data) > 0:
    up_s = sig_data[sig_data['log2FC_PA_vs_IVF'] > 1]
    dn_s = sig_data[sig_data['log2FC_PA_vs_IVF'] < -1]
    ax.scatter(up_s['log2FC_PA_vs_IVF'], -np.log10(up_s['pvalue'].clip(1e-300)),
              c='#E41A1C', alpha=0.6, s=6, label=f'Up in PA (n={len(up_s)})')
    ax.scatter(dn_s['log2FC_PA_vs_IVF'], -np.log10(dn_s['pvalue'].clip(1e-300)),
              c='#377EB8', alpha=0.6, s=6, label=f'Down in PA (n={len(dn_s)})')
ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.5)
ax.axvline(-1, color='gray', ls='--', lw=0.5); ax.axvline(1, color='gray', ls='--', lw=0.5)
ax.set_xlabel('log2FC (PA vs IVF)', fontsize=9)
ax.set_ylabel('-log10(p-value)', fontsize=9)
ax.set_title('IVF vs PA differential expression (6,860 genes)', fontsize=10, fontweight='bold')
ax.legend(fontsize=7); ax.tick_params(labelsize=7)
ax.set_xlim(-6, 6)

fig7.tight_layout()
fig7.savefig(os.path.join(OUT, 'fig7_divergence_volcano.svg'), format='svg', bbox_inches='tight')
fig7.savefig(os.path.join(OUT, 'fig7_divergence_volcano.png'), dpi=300, bbox_inches='tight')
plt.close(fig7)
print("  -> fig7 OK")

print("\n" + "="*60)
print("ALL 7 FIGURES REBUILT WITH FIXES")
print("="*60)
for fn in sorted(os.listdir(OUT)):
    if fn.endswith('.svg') and fn.startswith('fig'):
        size_kb = os.path.getsize(os.path.join(OUT, fn)) // 1024
        print(f"  {fn:35s} {size_kb:>6d} KB")