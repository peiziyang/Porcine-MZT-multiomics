#!/usr/bin/env python
"""Verify the headline numbers cited in v6 against the actual output CSVs."""
import os, sys
import pandas as pd

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
P = lambda *x: os.path.join(BASE, *x)

print("=" * 70)
print("VERIFY v6 headline claims against real output files")
print("=" * 70)

# ---- Fig 2D: IVF vs PA on GSE164812 (Wilcoxon on FPKM) ----
f = P("m3_ivf_vs_pa_results.csv")
if os.path.exists(f):
    d = pd.read_csv(f)
    print(f"\n[Fig2D] m3_ivf_vs_pa_results.csv  shape={d.shape}")
    print("  columns:", list(d.columns)[:12])
    # try to find significance / log2fc columns
    for col in d.columns:
        if 'log2' in col.lower() or 'lfc' in col.lower():
            lfc = col
        if 'padj' in col.lower() or 'pval' in col.lower() or 'adj' in col.lower():
            padj = col
    # fallback: just count rows
    print(f"  rows={len(d)}")
else:
    print("\n[Fig2D] MISSING:", f)

# ---- Fig 3: Early vs Late pseudobulk DESeq2 GSE168106 ----
f = P("m3_deseq2_results_early_vs_late.csv")
if os.path.exists(f):
    d = pd.read_csv(f)
    print(f"\n[Fig3] m3_deseq2_results_early_vs_late.csv  shape={d.shape}")
    print("  columns:", list(d.columns)[:12])
    print(f"  rows={len(d)}")
else:
    print("\n[Fig3] MISSING:", f)

# ---- Fig 5: pySCENIC regulons ----
f = P("pyscenic_out")
if os.path.isdir(f):
    files = os.listdir(f)
    print(f"\n[Fig5] pyscenic_out/ files={len(files)}")
    for x in sorted(files)[:8]:
        print("   ", x)
else:
    print("\n[Fig5] MISSING dir:", f)

# ---- Fig 8: MOFA+ ----
for f in [P("m8_mofa_figure.png"), P("m8_mofa_loadings.csv"), P("m8_mofa_r2.csv")]:
    print(f"[Fig8] {'OK ' if os.path.exists(f) else 'MISSING'} {os.path.basename(f)}")

# ---- Fig 7: Waddington-OT ----
for f in [P("m5_waddington_ot.png"), P("m9_lineage_rss_results.csv")]:
    print(f"[Fig7] {'OK ' if os.path.exists(f) else 'MISSING'} {os.path.basename(f)}")

# ---- Image-link existence check for every figure cited in v6 md ----
import re
md = open("E:/Workbuddy/2026-07-27-11-58-27/manuscript/manuscript_full_with_figures_v6.md", encoding="utf-8").read()
imgs = re.findall(r'!\[[^\]]*\]\(([^)]+)\)', md)
print(f"\n[IMG LINKS] total={len(imgs)}")
miss = 0
for im in imgs:
    if not os.path.exists(im):
        miss += 1
        print("  MISSING:", im)
print(f"  missing images: {miss}")

print("\nDONE")
