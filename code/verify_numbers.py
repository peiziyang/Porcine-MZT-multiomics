#!/usr/bin/env python
import os, pandas as pd
BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
P = lambda *x: os.path.join(BASE, *x)

print("=== Fig2D: IVF vs PA (GSE164812 FPKM Wilcoxon) ===")
d = pd.read_csv(P("m3_ivf_vs_pa_results.csv"))
sig = d[d['significant'] == True] if 'significant' in d.columns else d[d['padj'] < 0.05]
print("  total rows:", len(d))
print("  significant (True):", int((d['significant']==True).sum()))
print("  padj<0.05:", int((d['padj']<0.05).sum()))
# PA-down = log2FC_PA_vs_IVF < 0
print("  PA-down (log2FC<0 & sig):", int(((d['log2FC_PA_vs_IVF']<0)&(d['significant']==True)).sum()))
print("  PA-up   (log2FC>0 & sig):", int(((d['log2FC_PA_vs_IVF']>0)&(d['significant']==True)).sum()))
# with |log2FC|>0.5
m = (d['significant']==True)&(d['log2FC_PA_vs_IVF'].abs()>0.5)
print("  |log2FC|>0.5 & sig:", int(m.sum()), " down:", int((m&(d['log2FC_PA_vs_IVF']<0)).sum()), " up:", int((m&(d['log2FC_PA_vs_IVF']>0)).sum()))

print("\n=== Fig3: Early vs Late pseudobulk DESeq2 (GSE168106) ===")
d = pd.read_csv(P("m3_deseq2_results_early_vs_late.csv"))
deg = (d['padj']<0.05)&(d['log2FoldChange'].abs()>1)
print("  genes tested:", len(d))
print("  DEGs (padj<0.05, |log2FC|>1):", int(deg.sum()))
print("   up (late>early):", int((deg&(d['log2FoldChange']>0)).sum()))
print("   down:", int((deg&(d['log2FoldChange']<0)).sum()))

print("\n=== Fig5: pySCENIC regulons ===")
try:
    a = pd.read_csv(P("pyscenic_out","aucell_scores_full.csv"), index_col=0)
    print("  AUCell matrix shape (cells x regulons):", a.shape)
except Exception as e:
    print("  err", e)
# regulons count from regulons_ctx
for fn in ["regulons_ctx.csv","aucell_scores.csv"]:
    fp = P("pyscenic_out", fn)
    if os.path.exists(fp):
        r = pd.read_csv(fp, index_col=0) if fn.endswith('csv') else None
        print(f"  {fn} shape:", r.shape if r is not None else '?')

print("\n=== Fig8: MOFA+ variance ===")
for fn in ["m8_mofa_r2.csv","m8_mofa_loadings.csv"]:
    fp = P(fn)
    if os.path.exists(fp):
        r = pd.read_csv(fp, index_col=0)
        print(f"  {fn} shape:", r.shape)
        print("    cols:", list(r.columns)[:6])

print("\n=== Fig2/Fig7 trajectory means (m5) ===")
for fn in ["m5_top50_divergent_genes.csv","m5_trajectory_comparison.png","m5_waddington_ot.png"]:
    print("  ", "OK" if os.path.exists(P(fn)) else "MISSING", fn)
# try rss results
fp = P("m9_lineage_rss_results.csv")
if os.path.exists(fp):
    r = pd.read_csv(fp)
    print("  m9 rss shape:", r.shape, "cols:", list(r.columns)[:8])
