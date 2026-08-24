#!/usr/bin/env python
"""Generate CpG QC Figure S15 from pipeline output."""
import matplotlib; matplotlib.use('Agg')
matplotlib.rcParams.update({'font.family':'Arial','svg.fonttype':'none','savefig.dpi':300})
import matplotlib.pyplot as plt; import pandas as pd; import numpy as np; import os

OUT = r'REPLACE_WITH_PATH_TO_YOUR/manuscript/figures'
gs = pd.read_csv(r'REPLACE_WITH_PATH_TO_YOUR/data/processed/m4_mofa/per_gene_cpg_stats.csv')

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('CpG-level methylation QC (promoter regions, TSS +/- 2 kb)', fontsize=13, fontweight='bold')

# A: CpG count distribution
ax = axes[0,0]; ax.set_title('A: Informative CpG sites per gene')
ax.hist(np.clip(gs['mean_promoter_ncpg'], 0, 200), bins=80, color='#377EB8', edgecolor='black', lw=0.3)
ax.axvline(gs['mean_promoter_ncpg'].median(), color='#E41A1C', lw=2, ls='--', label=f'median={gs["mean_promoter_ncpg"].median():.0f}')
ax.set_xlabel('Informative CpG sites per promoter'); ax.set_ylabel('Gene count'); ax.legend(fontsize=8)

# B: Samples with data
ax = axes[0,1]; ax.set_title('B: Samples with promoter methylation data per gene')
ax.hist(gs['samples_with_data'], bins=45, color='#4DAF4A', edgecolor='black', lw=0.3)
ax.axvline(gs['samples_with_data'].median(), color='#E41A1C', lw=2, ls='--', label=f'median={gs["samples_with_data"].median():.0f}')
ax.set_xlabel('Number of samples with data'); ax.set_ylabel('Gene count'); ax.legend(fontsize=8)

# C: F1 genes vs all genes
ax = axes[1,0]
f1_w = pd.read_csv(r'REPLACE_WITH_PATH_TO_YOUR/data/processed/m4_mofa/mofa_multiomics_weights_RNA_v2.csv', index_col=0)
f1_genes = list(f1_w['F1'].abs().sort_values(ascending=False).head(44).index)
f1_stats = gs[gs['gene_name'].isin(f1_genes)]
all_stats = gs
ax.set_title(f'C: F1 genes CpG count (n={len(f1_stats)})')
ax.boxplot([all_stats['mean_promoter_ncpg'].values, f1_stats['mean_promoter_ncpg'].values],
           labels=['All genes\n(n=34,151)', f'F1 genes\n(n={len(f1_stats)})'], patch_artist=True,
           boxprops=dict(facecolor='#E5F0FF'), medianprops=dict(color='#E41A1C', lw=2))
ax.set_ylabel('Informative CpG sites per promoter')

# D: Core 5 genes (horizontal bar chart)
ax = axes[1,1]
core = []
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    m = gs[gs['gene_name'] == g]
    if len(m): core.append((g, m.iloc[0]['mean_promoter_ncpg']))
genes = [g for g, v in core]
vals = [v for g, v in core]
y = np.arange(len(genes))[::-1]  # 从上到下：DNMT1 在最上
colors = ['#377EB8'] * len(genes)
for i, g in enumerate(genes):
    if g == 'ZP4':
        colors[i] = '#BBBBBB'  # 最低值用灰色突出
ax.barh(y, vals, color=colors, edgecolor='#333333', lw=0.4, height=0.62)
ax.set_yticks(y)
ax.set_yticklabels(genes, fontsize=10)
ax.set_xlabel('Informative CpG sites per promoter (≥1 mapped read)', fontsize=9)
ax.set_title('D: Core five F1 genes', fontsize=11)
for yi, v in zip(y, vals):
    ax.text(v + 3, yi, f'{v:.0f}', va='center', fontsize=9, fontweight='bold', color='#333333')
gmed = gs['mean_promoter_ncpg'].median()
ax.axvline(gmed, color='#E41A1C', lw=1.6, ls='--', zorder=3, label=f'Genome-wide median ({gmed:.1f})')
ax.legend(fontsize=8, loc='lower right', frameon=False)
ax.set_xlim(0, max(vals) * 1.18)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(os.path.join(OUT, 'supp_fig_s15_cpg_qc.png'), dpi=300)
fig.savefig(os.path.join(OUT, 'supp_fig_s15_cpg_qc.svg'), format='svg')
print(f'S15 saved: {os.path.getsize(os.path.join(OUT, "supp_fig_s15_cpg_qc.svg"))//1024}KB SVG')
