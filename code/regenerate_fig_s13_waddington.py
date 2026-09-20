#!/usr/bin/env python
"""重绘 Fig S13（Waddington optimal-transport），用含 GSE164812 的 m5_ot 数据。

数据源：m5_ot_vivo_fate.csv（in_vivo E0-E8 细胞的 fate entropy + UMAP）
        m5_ot_displacement.csv（in_vivo/IVF/PA 的 OT 发育位移）
逻辑与 rebuild_supp_svgs_real.py 的 Waddington 部分一致。
"""
import os, pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'Arial',
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
})

P = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
OUT = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/submission/figures"

fate = pd.read_csv(f"{P}/m5_ot_vivo_fate.csv")
disp = pd.read_csv(f"{P}/m5_ot_displacement.csv")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
sc = ax.scatter(fate['UMAP1'], fate['UMAP2'], c=fate['fate_entropy'],
                cmap='viridis', s=6, alpha=0.6)
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2'); ax.set_title('Fate entropy')
plt.colorbar(sc, ax=ax, label='entropy', shrink=0.8)

ax = axes[1]
for cond in disp['condition'].dropna().unique():
    sub = disp[disp['condition'] == cond]
    tmp = sub.groupby('stage')['displacement'].mean()
    ax.plot(range(len(tmp)), tmp.values, 'o-', label=cond)
ax.set_xticks(range(len(tmp)))
ax.set_xticklabels(tmp.index, fontsize=7, rotation=45)
ax.set_ylabel('Displacement'); ax.set_title('OT displacement'); ax.legend(fontsize=7)

# --- panel labels A/B (bold, top-left, white bbox) ---
for ch, a in zip('AB', axes):
    a.text(0.01, 0.97, ch, transform=a.transAxes, fontsize=12, fontweight='bold',
           va='top', ha='left', zorder=20,
           bbox=dict(boxstyle='square,pad=0.18', facecolor='white',
                     edgecolor='none', alpha=0.85))

fig.savefig(f"{OUT}/Fig9_waddington_ot.png", dpi=300, bbox_inches='tight')
fig.savefig(f"{OUT}/Fig9_waddington_ot.svg", format='svg', bbox_inches='tight')
fig.savefig(f"{OUT}/Fig9_waddington_ot.pdf", format='pdf', bbox_inches='tight')
plt.close()
print("done: Fig9_waddington_ot.{png,svg,pdf}")
