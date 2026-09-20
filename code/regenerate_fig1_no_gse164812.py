#!/usr/bin/env python
"""重新生成 Fig 1（移除 GSE164812 后，1,914 细胞图谱）。
只改 Fig 1A：标题 1,955 -> 1,914，数据用更新后的 m2_metadata/m2_harmony_umap。
其余 panel B-G 与原来一致。"""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial', 'svg.fonttype': 'none', 'savefig.dpi': 300,
    'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7,
})
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
from adjustText import adjust_text
import os

BASE  = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
MF    = os.path.join(BASE, 'm4_mofa')
OUT   = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/submission/figures'

meta    = pd.read_csv(os.path.join(BASE, 'm2_metadata.csv'))
umap_df = pd.read_csv(os.path.join(BASE, 'm2_harmony_umap.csv'))
factor  = pd.read_csv(os.path.join(MF, 'mofa_multiomics_factors_v2.csv'), index_col=0)
r2_df   = pd.read_csv(os.path.join(MF, 'mofa_multiomics_r2_per_view_v2.csv'))
cross_v = pd.read_csv(os.path.join(MF, 'factor_cross_view_v2.csv'))
wrna    = pd.read_csv(os.path.join(MF, 'mofa_multiomics_weights_RNA_v2.csv'), index_col=0)
wmeth   = pd.read_csv(os.path.join(MF, 'mofa_multiomics_weights_METH_v2.csv'), index_col=0)

FCOLS = ['F1','F2','F3','F4','F5','F6','F7']

def panel_label(ax, label, x=0.02, y=0.97, fontsize=6):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize,
            fontweight='bold', va='top', ha='left')

print("=== Fig1 ===")
fig1 = plt.figure(figsize=(5.8, 9.0))
gs1 = fig1.add_gridspec(4, 2, hspace=0.35, wspace=0.3,
                        top=0.97, bottom=0.03, left=0.06, right=0.98)

# ── A: Atlas UMAP ──
ax = fig1.add_subplot(gs1[0, 0])
panel_label(ax, 'A')
meta_umap = meta.copy()
meta_umap['UMAP1'] = umap_df['UMAP1_harmony'].values
meta_umap['UMAP2'] = umap_df['UMAP2_harmony'].values
stage_order = sorted(meta_umap['stage_group'].dropna().unique(), key=str)
cmap = plt.cm.tab10
for i, st in enumerate(stage_order):
    sub = meta_umap[meta_umap['stage_group'] == st]
    if len(sub):
        ax.scatter(sub['UMAP1'], sub['UMAP2'], s=8, alpha=0.6,
                   c=[cmap(i % 10)], label=st, rasterized=True)
ax.set_title('Porcine preimplantation atlas (1,910 cells)', fontsize=8, fontweight='bold')
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.legend(loc='lower right', fontsize=6, ncol=2, markerscale=2, title='Stage', framealpha=0.85)
ax.tick_params(labelsize=7)

# ── B: Factor scores heatmap ──
ax = fig1.add_subplot(gs1[0, 1])
panel_label(ax, 'B')
fmat = factor[FCOLS].values.T
im = ax.imshow(fmat, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
ax.set_yticks(range(7)); ax.set_yticklabels(FCOLS, fontsize=6)
ax.set_xticks([]); ax.set_xlabel('32 GV oocytes', fontsize=7)
ax.set_title('Factor scores (z-scaled)', fontsize=8, fontweight='bold')
cbar = fig1.colorbar(im, ax=ax, shrink=0.8, aspect=20)
cbar.ax.tick_params(labelsize=7)

# ── C: Per-view variance explained ──
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
ax.set_xticks(x); ax.set_xticklabels(FCOLS, fontsize=6)
ax.set_ylabel('Variance contribution (%)', fontsize=7)
ax.set_title('Per-view variance contribution', fontsize=8, fontweight='bold')
ax.legend(fontsize=6, loc='upper right')
ax.set_ylim(0, max(max(rna_vals), max(meth_vals)) * 1.3)

# ── D: Cross-view weight correlation ──
ax = fig1.add_subplot(gs1[1, 1])
panel_label(ax, 'D')
r_vals = cross_v['rna_meth_corr_r'].values
colors_d = ['#E41A1C' if r > 0.5 else '#377EB8' if r < -0.3 else '#999999' for r in r_vals]
ax.bar(FCOLS, r_vals, color=colors_d, edgecolor='black', lw=0.5, width=0.6)
ax.axhline(0, color='gray', lw=0.8)
for i, r in enumerate(r_vals):
    ax.text(i, r + 0.04 if r >= 0 else r - 0.06, f'r={r:+.2f}', ha='center', fontsize=6, fontweight='bold')
ax.set_ylabel('RNA-METH weight Pearson r', fontsize=7)
ax.set_title('Cross-view weight correlation', fontsize=8, fontweight='bold')
ax.set_ylim(-0.15, 1.05)

# ── E: Factor scores by donor ──
ax = fig1.add_subplot(gs1[2, 0])
panel_label(ax, 'E')
donor_map = meta.set_index('cell')['donor'].to_dict()
donor_colors = {'AF1': '#E41A1C', 'AF2': '#377EB8', 'AF3': '#4DAF4A', 'AF4': '#984EA3', 'AF5': '#FF7F00'}
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
ax.set_xticks(range(1, 8)); ax.set_xticklabels(FCOLS, fontsize=6)
ax.axhline(0, color='gray', ls='--', lw=0.5)
ax.set_xlabel('Factor', fontsize=7); ax.set_ylabel('Score', fontsize=7)
ax.set_title('Factor scores by donor', fontsize=8, fontweight='bold')
ax.legend(fontsize=6, loc='upper right', ncol=5, markerscale=0.8)
ax.tick_params(labelsize=7)

# ── F: F1 RNA vs METH weight scatter ──
ax = fig1.add_subplot(gs1[2, 1])
panel_label(ax, 'F')
r1, p1 = pearsonr(wrna['F1'].values, wmeth['F1'].values)
ax.scatter(wrna['F1'].values, wmeth['F1'].values, c='gray', alpha=0.2, s=3)
shared_idx = np.argsort(np.abs(wrna['F1'].values) + np.abs(wmeth['F1'].values))[-6:]
texts = []
for idx in shared_idx:
    texts.append(ax.annotate(wrna.index[idx], (wrna['F1'].iloc[idx], wmeth['F1'].iloc[idx]),
               fontsize=7, fontweight='bold', color='#E41A1C'))
adjust_text(texts, ax=ax, force_text=0.5, force_points=0.3, force_static=0.2,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.4))
ax.set_xlabel('F1 RNA weight', fontsize=7); ax.set_ylabel('F1 METH weight', fontsize=7)
ax.set_title(f'F1: RNA vs METH weights (r={r1:.3f}, p={p1:.1e})', fontsize=8, fontweight='bold')
ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
ax.tick_params(labelsize=7)

# ── G: GO enrichment ──
ax = fig1.add_subplot(gs1[3, :])
panel_label(ax, 'G')
go = pd.read_csv('E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa/factor_go_enrichment_v2.csv')
go['log10p'] = -np.log10(go['p_raw'].clip(lower=1e-300))
top = go.groupby('factor').apply(lambda x: x.nsmallest(3, 'p_raw')).reset_index(drop=True)
colors = {'F1':'#377EB8','F2':'#FF7F00','F3':'#E41A1C','F4':'#4DAF4A','F5':'#984EA3','F6':'#A65628','F7':'#F781BF'}
for f in sorted(top['factor'].unique()):
    sub = top[top['factor']==f]
    ax.barh(sub.index[::-1], sub['log10p'][::-1], color=colors.get(f,'#999'), label=f)
terms = [f"{r['factor']} | {str(r['GO_term'])[:48]}" for _, r in top.iterrows()]
ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms, fontsize=6)
ax.set_xlabel('-log10(raw p)', fontsize=7); ax.invert_yaxis()
ax.set_title('GO enrichment (top 3 terms per factor)', fontsize=8, fontweight='bold')
ax.legend(fontsize=6, ncol=7, loc='lower center', bbox_to_anchor=(0.5, -0.18),
           frameon=True, framealpha=0.9, edgecolor='gray')
ax.tick_params(labelsize=6)

fig1.savefig(os.path.join(OUT, 'fig1_mofa.svg'), format='svg', bbox_inches='tight')
fig1.savefig(os.path.join(OUT, 'fig1_mofa.png'), dpi=300, bbox_inches='tight')
fig1.savefig(os.path.join(OUT, 'fig1_mofa.pdf'), format='pdf', bbox_inches='tight')
plt.close(fig1)
print("-> fig1_mofa.{svg,png,pdf} OK (1,910 cells)")
