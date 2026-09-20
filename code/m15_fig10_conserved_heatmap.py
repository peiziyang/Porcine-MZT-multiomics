"""
m15 Fig 10: conserved DEG heatmap across cleavage stages (IVF vs PA).
Conserved DEG = gene significantly & directionally consistent across >=3 of 4 stages.
"""
import pandas as pd, numpy as np, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

BASE = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed/deseq2'
OUT  = os.path.join(BASE, 'm15_fig10_conserved_heatmap.png')
stages = ['1cell', '2cell', '4cell', '8cell']

# ---- load per-stage, dedup gene_symbol (keep most significant) ----
lfc = {}
padj = {}
for s in stages:
    df = pd.read_csv(f'{BASE}/m15_deseq2_{s}_IVFvsPA.csv')
    df = df.sort_values('padj').drop_duplicates('gene_symbol', keep='first')
    lfc[s]  = df.set_index('gene_symbol')['log2FoldChange']
    padj[s] = df.set_index('gene_symbol')['padj']

lfc  = pd.DataFrame(lfc)
padj = pd.DataFrame(padj)
lfc  = lfc.loc[lfc.index.intersection(padj.index)]
padj = padj.loc[lfc.index]

# ---- define conserved DEG ----
TH_FC, TH_P = 1.0, 0.05
sig   = (padj < TH_P) & (lfc.abs() > TH_FC)
n_sig = sig.sum(axis=1)
pos   = (lfc > 0).sum(axis=1)
neg   = (lfc < 0).sum(axis=1)
consistent = (pos >= 3) | (neg >= 3)        # >=3 of 4 stages same direction
mask = (n_sig >= 2) & consistent            # robust in >=2 stages AND direction-consistent

info = pd.DataFrame({
    'n_sig': n_sig,
    'pos': pos, 'neg': neg,
    'mean_abs': lfc.abs().mean(axis=1),
    'mean_lfc': lfc.mean(axis=1),
})
sel = info[mask].sort_values(['n_sig', 'mean_abs'], ascending=False)
top = sel.index[:30]
mat = lfc.loc[top, stages]

# ---- row order by hierarchical clustering on LFC profile ----
Z = linkage(pdist(mat.values), method='average')
order = leaves_list(Z)
mat_o = mat.iloc[order]

# ---- plot ----
vmax = np.nanmax(np.abs(mat_o.values))
fig, ax = plt.subplots(figsize=(6.5, 11))
im = ax.imshow(mat_o.values, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
ax.set_xticks(range(len(stages)))
ax.set_xticklabels([f'{s}\n(PA vs IVF)' for s in stages], fontsize=11)
ax.set_yticks(range(len(mat_o)))
ax.set_yticklabels(mat_o.index, fontsize=8)
ax.set_title('Figure 10. Top 30 conserved DEGs across cleavage stages\n(log2FC, PA vs IVF; red=PA-up, blue=PA-down)',
             fontsize=11, pad=10)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cb.set_label('log2 fold change (PA / IVF)', fontsize=9)
# direction annotation per gene
for i, g in enumerate(mat_o.index):
    direction = 'PA-up' if info.loc[g, 'mean_lfc'] > 0 else 'PA-down'
    ax.text(len(stages) + 0.35, i, direction, va='center', fontsize=6.5, color='#444')
fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches='tight')
print('Saved', OUT)
print(f'Conserved DEG candidates (n_sig>=2 & consistent): {int(mask.sum())}')
print(f'Top-30 selected. n_sig distribution among them: {sorted(sel.loc[top,"n_sig"].tolist(), reverse=True)}')
print('\nTop 30 conserved DEGs (gene | meanLFC | n_sig | dominant direction):')
for g in top:
    d = 'PA-up' if info.loc[g, 'mean_lfc'] > 0 else 'PA-down'
    print(f'  {g:12s} meanLFC={info.loc[g,"mean_lfc"]:+.2f}  n_sig={int(info.loc[g,"n_sig"])}  {d}')
