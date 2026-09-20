"""
m15_fig11_replot.py — 从已有的 GO 富集结果 CSV 重绘 Figure 11。
不重扫 1.3GB gene2go，直接从 m15_fig11_go_enrich.csv 读 p 值。

诚实呈现策略（避免 p-hacking）：
- BH 校正后无 term 达 q<0.05（富集信号弱）。
- 改用以 raw p<0.05 的 top 富集 term（明确标注「未校正、提示性探索」）。
- 图注与正文一致说明：经多重校正后无显著 term，代谢/TCA 相关通路未系统性富集。
"""
import os, pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed/deseq2'
CSV  = os.path.join(BASE, 'm15_fig11_go_enrich.csv')
PNG  = os.path.join(BASE, 'm15_fig11_go_enrich.png')

df = pd.read_csv(CSV)
for c in ['p', 'q', 'deg_n', 'bg_n']:
    df[c] = pd.to_numeric(df[c], errors='coerce')
df = df.dropna(subset=['p'])
df['log10p'] = -np.log10(df['p'].clip(lower=1e-300))

domains = [('Process', 'Biological Process (BP)'),
           ('Function', 'Molecular Function (MF)'),
           ('Component', 'Cellular Component (CC)')]
colors  = {'Process': '#2E7D32', 'Function': '#1565C0', 'Component': '#E65100'}

fig, axes = plt.subplots(1, 3, figsize=(15, 6.2))
for ax, (dom, title) in zip(axes, domains):
    sub = df[df['domain'] == dom].sort_values('p').head(8).sort_values('p')
    if len(sub) == 0:
        ax.text(0.5, 0.5, 'no term (p<0.05)', ha='center', va='center')
        ax.set_title(title); continue
    y = np.arange(len(sub))
    ax.barh(y, sub['log10p'], color=colors[dom])
    ax.set_yticks(y)
    ax.set_yticklabels([f"{t[:42]}" for t in sub['term']], fontsize=7.5)
    for yi, (lp, p) in enumerate(zip(sub['log10p'], sub['p'])):
        ax.text(lp + 0.05, yi, f"{p:.1e}", va='center', fontsize=6.5, color='#333')
    ax.set_xlabel('-log10(raw p)', fontsize=8)
    ax.set_title(f"{title}\n(top 8, p<0.05, n={len(sub)})", fontsize=9)
    ax.invert_yaxis()
fig.suptitle("Figure 11 — Top GO terms enriched in stage-conserved PA-down/up DEGs\n"
             "raw p-values (NOT BH-adjusted; exploratory). After BH correction no term reached q<0.05.",
             fontsize=10, y=1.02)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(PNG, dpi=150, bbox_inches='tight')
print('saved', PNG, os.path.getsize(PNG), 'bytes')

# 同时导出 raw-p top 表供正文引用
sig = df[df['p'] < 0.05].sort_values('p')
out = sig[['GO', 'domain', 'term', 'bg_n', 'deg_n', 'p', 'q']].head(30)
out.to_csv(os.path.join(BASE, 'm15_fig11_top_rawp.tsv'), sep='\t', index=False)
print('top raw-p terms exported:', len(out))
print(out[['domain', 'term', 'p', 'q']].to_string(index=False))
