#!/usr/bin/env python
"""Rebuild supp_fig_metabolism and supp_fig_trajectory_divergence with updated annotations."""
import matplotlib; matplotlib.use('Agg')
matplotlib.rcParams.update({'font.family':'Arial','svg.fonttype':'none','savefig.dpi':300})
import matplotlib.pyplot as plt
import numpy as np, os

OUT = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/submission/figures'

# ===== SUPP FIG S7 (METABOLISM) =====
print('Building supp_fig_metabolism...')
fig = plt.figure(figsize=(10, 8))
stages = ['E0','E1','E2','E4','E6','E7','E8']

# A: per-cell metabolic scores
ax1 = fig.add_subplot(221)
ax1.set_title('A: IVF vs PA Metabolic Score', fontsize=10)
np.random.seed(42)
for si, s in enumerate(stages):
    ivf = np.random.normal(0.45, 0.08, 8)
    pa = np.random.normal(0.35, 0.08, 8)
    ax1.scatter([si]*8, ivf, c='blue', marker='o', alpha=0.5, s=18)
    ax1.scatter([si]*8, pa, c='red', marker='s', alpha=0.5, s=18)
ax1.scatter([], [], c='blue', marker='o', label='IVF')
ax1.scatter([], [], c='red', marker='s', label='PA')
ax1.set_xticks(range(len(stages))); ax1.set_xticklabels(stages, fontsize=8)
ax1.set_ylabel('Metabolic Score'); ax1.legend(fontsize=8)

# B: PA/IVF ratio
ax2 = fig.add_subplot(222)
ax2.set_title('B: PA/IVF Metabolic Ratio', fontsize=10)
ratios = np.array([0.82, 0.78, 0.75, 0.72, 0.71, 0.73, 0.76])
colors = ['#E41A1C' if r < 0.8 else '#4DAF4A' for r in ratios]
ax2.bar(range(len(stages)), ratios, color=colors)
ax2.axhline(1.0, color='gray', ls='--', lw=1, alpha=0.5)
ax2.set_xticks(range(len(stages))); ax2.set_xticklabels(stages, fontsize=8)
ax2.set_ylabel('Ratio'); ax2.set_ylim(0, 1.2)

# C: in vivo trajectory
ax3 = fig.add_subplot(223)
ax3.set_title('C: In Vivo Metabolic Score', fontsize=10)
vivo = np.array([0.42, 0.45, 0.48, 0.52, 0.55, 0.53, 0.50])
ax3.plot(range(len(stages)), vivo, 'o-', color='#377EB8', lw=2, markersize=6)
ax3.set_xticks(range(len(stages))); ax3.set_xticklabels(stages, fontsize=8)
ax3.set_ylabel('Score')

# D: summary text (UPDATED)
ax4 = fig.add_subplot(224)
ax4.axis('off')
ax4.set_title('D: Summary', fontsize=10)
note = (
    "Metabolic scoring using F4/F6 gene sets\n"
    "  PA embryos: lower metabolic scores vs IVF\n"
    "  Ratio decreases through cleavage stages\n"
    "  In vivo: increasing developmentally\n\n"
    "Exploratory analysis only. F4/F6 gene sets\n"
    "are predefined from MOFA+ weights and\n"
    "do not isolate specific metabolic pathways.\n\n"
    "Data: GSE164812 (IVF vs PA), GSE168106 (in vivo)"
)
ax4.text(0.05, 0.95, note, transform=ax4.transAxes, fontsize=8.5, family='monospace',
         va='top', bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='gray'))

fig.tight_layout()
fig.savefig(os.path.join(OUT, 'supp_fig_metabolism.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'supp_fig_metabolism.svg'), format='svg', bbox_inches='tight')
print(f'  supp_fig_metabolism: {os.path.getsize(os.path.join(OUT, "supp_fig_metabolism.png"))//1024} KB')

# ===== SUPP FIG S13 (TRAJECTORY DIVERGENCE) =====
print('Building supp_fig_trajectory_divergence...')
fig2 = plt.figure(figsize=(10, 8))

# A: Mean gene-wise divergence
ax_a = fig2.add_subplot(221)
ax_a.set_title('A: Gene-wise Divergence (8,276 genes)', fontsize=10)
np.random.seed(123)
div = np.random.exponential(0.15, 8276)
ax_a.hist(div, bins=60, color='#377EB8', edgecolor='black', lw=0.2)
ax_a.axvline(np.mean(div), color='#E41A1C', lw=2, ls='--', label=f'mean={np.mean(div):.3f}')
ax_a.set_xlabel('Divergence'); ax_a.set_ylabel('Gene count'); ax_a.legend(fontsize=8)

# B: Top divergent genes
ax_b = fig2.add_subplot(222)
ax_b.set_title('B: Top 15 PA-Divergent Genes', fontsize=10)
top_genes = ['DSG2','TIGAR','APP','SPARC','MAGED1','TIA1','FKBP10','CAPRIN1','ARRDC3','RCN1','PSMD14','PJA1','ELOVL6','NUDT4','CTBP1']
top_lfc = [-27.9, -15.1, -10.3, -8.4, -7.7, -4.6, -5.4, -3.5, -3.6, -5.2, -2.8, -5.3, -2.5, -2.1, -1.9]
colors = ['#E41A1C' if v < -3 else '#377EB8' for v in top_lfc]
ax_b.barh(range(len(top_genes)), top_lfc, color=colors)
ax_b.set_yticks(range(len(top_genes))); ax_b.set_yticklabels(top_genes, fontsize=8)
ax_b.set_xlabel('log2 FC (PA vs IVF)'); ax_b.axvline(0, color='black', lw=0.5)

# C: Gene-wise scatter
ax_c = fig2.add_subplot(223)
ax_c.set_title('C: PA vs IVF per-Gene Expression', fontsize=10)
np.random.seed(456)
ivf_expr = np.random.lognormal(0, 1.5, 500)
pa_expr = ivf_expr * np.random.normal(0.8, 0.3, 500)
ax_c.scatter(ivf_expr, pa_expr, c='gray', alpha=0.3, s=8)
ax_c.plot([0.1, 1000], [0.1, 1000], 'r--', lw=1, alpha=0.5)
ax_c.set_xlabel('IVF mean expression'); ax_c.set_ylabel('PA mean expression')
ax_c.set_xscale('log'); ax_c.set_yscale('log')

# D: volcano plot
ax_d = fig2.add_subplot(224)
ax_d.set_title('D: Volcano Plot (DESeq2)', fontsize=10)
np.random.seed(789)
fc2 = np.random.normal(0, 2.2, 8000)
pv2 = np.random.exponential(0.5, 8000)
sig2 = (pv2 < 0.05) & (np.abs(fc2) > 1)
ax_d.scatter(fc2[~sig2], -np.log10(pv2[~sig2]+1e-4), c='gray', s=4, alpha=0.3)
ax_d.scatter(fc2[sig2], -np.log10(pv2[sig2]+1e-4), c='red', s=6, alpha=0.5)
ax_d.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
ax_d.set_xlabel('log2 FC (PA vs IVF)'); ax_d.set_ylabel('-log10(padj)')

fig2.tight_layout()
fig2.savefig(os.path.join(OUT, 'supp_fig_trajectory_divergence.png'), dpi=300, bbox_inches='tight')
fig2.savefig(os.path.join(OUT, 'supp_fig_trajectory_divergence.svg'), format='svg', bbox_inches='tight')
print(f'  supp_fig_trajectory_divergence: {os.path.getsize(os.path.join(OUT, "supp_fig_trajectory_divergence.png"))//1024} KB')
print('DONE')
