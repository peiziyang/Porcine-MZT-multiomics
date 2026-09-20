#!/usr/bin/env python
"""Regenerate supp_fig_cpg_informative with bar chart for panel D."""
import matplotlib; matplotlib.use('Agg')
matplotlib.rcParams.update({'font.family':'Arial','svg.fonttype':'none','savefig.dpi':300})
import matplotlib.pyplot as plt
import pandas as pd, numpy as np, os
from matplotlib.colors import LinearSegmentedColormap

OUT = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/submission/figures'

gs = pd.read_csv(r'E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa/per_gene_cpg_stats.csv')
qc = pd.read_csv(r'E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa/cpg_qc_report.csv')

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('CpG-level methylation QC (promoter regions, TSS +/- 2 kb)', fontsize=13, fontweight='bold')

# A: Histogram
ax = axes[0,0]; ax.set_title('A: Informative CpG sites per gene')
ax.hist(np.clip(gs['mean_promoter_ncpg'], 0, 200), bins=80, color='#377EB8', edgecolor='black', lw=0.3)
ax.axvline(gs['mean_promoter_ncpg'].median(), color='#E41A1C', lw=2, ls='--', label=f'median={gs["mean_promoter_ncpg"].median():.1f}')
ax.set_xlabel('Informative CpG sites per promoter'); ax.set_ylabel('Gene count'); ax.legend(fontsize=8)

# B: Samples with data
ax = axes[0,1]; ax.set_title('B: Samples with promoter methylation data per gene')
ax.hist(gs['samples_with_data'], bins=45, color='#4DAF4A', edgecolor='black', lw=0.3)
ax.axvline(gs['samples_with_data'].median(), color='#E41A1C', lw=2, ls='--', label=f'median={gs["samples_with_data"].median():.0f}')
ax.set_xlabel('Number of samples with data'); ax.set_ylabel('Gene count'); ax.legend(fontsize=8)

# C: Boxplot
ax = axes[1,0]; ax.set_title('C: F1 genes CpG count (n=41)')
f1_w = pd.read_csv(r'E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa/mofa_multiomics_weights_RNA_v2.csv', index_col=0)
f1_genes = list(f1_w['F1'].abs().sort_values(ascending=False).head(44).index)
f1_stats = gs[gs['gene_name'].isin(f1_genes)]
all_stats = gs
ax.boxplot([all_stats['mean_promoter_ncpg'].values, f1_stats['mean_promoter_ncpg'].values],
           tick_labels=['All genes\n(n=34,151)', 'F1 genes\n(n=41)'], patch_artist=True,
           boxprops=dict(facecolor='#E5F0FF'), medianprops=dict(color='#E41A1C', lw=2))
ax.set_ylabel('Informative CpG sites per promoter')

# D: Bar chart for core 5 genes
ax = axes[1,1]; ax.set_title('D: Core five F1 gene CpG statistics')
core = []
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    m = gs[gs['gene_name'] == g]
    if len(m): core.append((g, m.iloc[0]['mean_promoter_ncpg']))
core_genes = [c[0] for c in core]
core_vals = [c[1] for c in core]
colors = ['#E41A1C' if v < 30 else '#377EB8' for v in core_vals]
bars = ax.barh(core_genes, core_vals, color=colors, edgecolor='black', lw=0.5)
# Annotate values
for bar, v in zip(bars, core_vals):
    ax.text(v + 5, bar.get_y() + bar.get_height()/2, f'{v:.0f}',
            va='center', fontsize=9, fontweight='bold')
# Reference line: median = 31.5
ax.axvline(31.5, color='gray', ls='--', lw=1.5, alpha=0.7, label='All-genes median (31.5)')
ax.set_xlabel('Informative CpG sites per promoter')
ax.set_xlim(0, 180)
ax.legend(loc='lower right', fontsize=8)
ax.invert_yaxis()  # top-to-bottom matching paper order
# Add annotation box for context
ax.text(0.97, 0.05, 'n = 41 of 44 F1 genes\nmapped in CGmap data',
        transform=ax.transAxes, fontsize=8, ha='right',
        bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='gray'))

plt.tight_layout()
fig.savefig(os.path.join(OUT, 'Fig8_cpg_informative.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'Fig8_cpg_informative.svg'), format='svg', bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'Fig8_cpg_informative.pdf'), format='pdf', bbox_inches='tight')
print(f'Saved: {os.path.getsize(os.path.join(OUT, "Fig8_cpg_informative.png"))//1024} KB')