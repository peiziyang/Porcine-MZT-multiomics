#!/usr/bin/env python
import os, pandas as pd
BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
P = lambda *x: os.path.join(BASE, *x)
d = pd.read_csv(P("m3_ivf_vs_pa_results.csv"))
padj, lfc = d['padj'], d['log2FC_PA_vs_IVF']
print("=== significant column value_counts ===")
print(d['significant'].value_counts())
nd = (d['significant']=='Down (PA<IVF)').sum()
nu = (d['significant']=='Up (PA>IVF)').sum()
print(f"  Down label count={nd}, Up label count={nu}, sum={nd+nu}")
# infer threshold of significant column
sub = d[d['significant']!='NS']
if len(sub):
    print("  |lfc| min among labeled:", round(sub['log2FC_PA_vs_IVF'].abs().min(),3))
    print("  padj max among labeled:", round(sub['padj'].max(),4))
# count at |lfc|>1.0
m1 = (padj<0.05)&(lfc.abs()>1.0)
print(f"\n  padj<0.05 & |lfc|>1.0 : {int(m1.sum())}  down={int((m1&(lfc<0)).sum())} up={int((m1&(lfc>0)).sum())}")
# m3c_validation outputs
import glob
for f in glob.glob(P("m3c_validation","*")):
    print("  m3c file:", os.path.basename(f))
fp = P("m3c_validation","peg_meg_overlap.csv")
if os.path.exists(fp):
    print(pd.read_csv(fp).to_string())
