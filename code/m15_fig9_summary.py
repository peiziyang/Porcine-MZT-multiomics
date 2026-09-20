#!/usr/bin/env python
"""Generate Fig 9 summary figure for the 42-sample bulk RNA-seq IVF-vs-PA validation cohort.
Panel A: grouped bar of per-stage DEG counts (PA-down / PA-up) at padj<0.05 & |log2FC|>1.
Panel B: the most robust stage (4cell, n=9 per group) volcano.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/deseq2"
OUT = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/deseq2/m15_fig9_summary.png"

stages = ['1cell', '2cell', '4cell', '8cell']
up, down = [], []
for s in stages:
    df = pd.read_csv(os.path.join(BASE, f"m15_deseq2_{s}_IVFvsPA.csv"))
    sig = (df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 1)
    up.append(int((sig & (df['log2FoldChange'] > 0)).sum()))
    down.append(int((sig & (df['log2FoldChange'] < 0)).sum()))

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: grouped bars
x = np.arange(len(stages))
w = 0.38
axes[0].bar(x - w/2, down, w, label='PA-down (PA<IVF)', color='#d62728')
axes[0].bar(x + w/2, up, w, label='PA-up (PA>IVF)', color='#1f77b4')
axes[0].set_xticks(x); axes[0].set_xticklabels([f"{s}\n(n=3/3)" if s in ('1cell','8cell') else f"{s}\n(n=6/6)" if s=='2cell' else f"{s}\n(n=9/9)" for s in stages])
axes[0].set_ylabel('DEGs (padj<0.05, |log2FC|>1)')
axes[0].set_title('Fig 9A. IVF vs PA DEGs per cleavage stage (42-sample bulk RNA-seq)')
axes[0].legend()
for i, (d, u) in enumerate(zip(down, up)):
    axes[0].text(i - w/2, d + 8, str(d), ha='center', fontsize=8, color='#d62728')
    axes[0].text(i + w/2, u + 8, str(u), ha='center', fontsize=8, color='#1f77b4')

# Panel B: 4cell volcano
df = pd.read_csv(os.path.join(BASE, "m15_deseq2_4cell_IVFvsPA.csv"))
padj = df['padj'].clip(lower=1e-300)
sig = (df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 1)
cols = np.where(sig & (df['log2FoldChange'] > 0), '#1f77b4',
        np.where(sig & (df['log2FoldChange'] < 0), '#d62728', 'grey'))
axes[1].scatter(df['log2FoldChange'], -np.log10(padj), c=cols, s=5, alpha=0.4, rasterized=True)
axes[1].axhline(-np.log10(0.05), color='grey', ls='--', alpha=0.5)
axes[1].axvline(1, color='grey', ls=':', alpha=0.3); axes[1].axvline(-1, color='grey', ls=':', alpha=0.3)
axes[1].set_xlabel('log2 FC (PA - IVF)'); axes[1].set_ylabel('-log10(padj)')
axes[1].set_title('Fig 9B. Volcano: IVF vs PA @ 4-cell (n=9/9, best-powered)')

fig.tight_layout()
fig.savefig(OUT, dpi=140, bbox_inches='tight')
print("saved", OUT)
print("per-stage DEGs (down/up):", dict(zip(stages, zip(down, up))))
