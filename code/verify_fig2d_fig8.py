#!/usr/bin/env python
import os, pandas as pd
BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
P = lambda *x: os.path.join(BASE, *x)

print("=== Fig2D precise recompute (GSE164812 FPKM Wilcoxon) ===")
d = pd.read_csv(P("m3_ivf_vs_pa_results.csv"))
print("columns:", list(d.columns))
padj = d['padj']
lfc = d['log2FC_PA_vs_IVF']
n_sig = int((padj < 0.05).sum())
n_sig_abs05 = int(((padj < 0.05) & (lfc.abs() > 0.5)).sum())
down = int(((padj < 0.05) & (lfc < 0) & (lfc.abs() > 0.5)).sum())
up   = int(((padj < 0.05) & (lfc > 0) & (lfc.abs() > 0.5)).sum())
print(f"  padj<0.05 total            : {n_sig}")
print(f"  padj<0.05 & |log2FC|>0.5   : {n_sig_abs05}")
print(f"     of which PA-down (lfc<0): {down}")
print(f"     of which PA-up   (lfc>0): {up}")
# also check the 'significant' column meaning
print("  'significant'==True count   :", int((d['significant']==True).sum()))
print("  significant col unique vals  :", d['significant'].unique()[:5])
# what threshold gives 186?
for thr in [0, 0.25, 0.5, 1.0]:
    c = int(((padj<0.05)&(lfc.abs()>thr)).sum())
    print(f"  padj<0.05 & |lfc|>{thr} : {c}")

print("\n=== Fig8 MOFA+ variance check ===")
r = pd.read_csv(P("m8_mofa_r2.csv"), index_col=0)
print(r.to_string())
tot = r['R2_pct'].sum()
print("  total variance explained (sum R2_pct):", round(tot,2))
if 'F1' in r.index:
    print("  F1 R2_pct:", round(r.loc['F1','R2_pct'],2))

print("\nDONE")
