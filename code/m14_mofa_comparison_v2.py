#!/usr/bin/env python
"""MOFA+ multi-omics v2 downstream analysis.

Inputs (from m13_v2):
  - mofa_multiomics_factors_v2.csv       (n_cells x n_factors)
  - mofa_multiomics_weights_RNA_v2.csv   (n_genes x n_factors)
  - mofa_multiomics_weights_METH_v2.csv  (n_genes x n_factors)
  - mofa_multiomics_r2_per_view_v2.csv   (factor, view, r2_pct)

Outputs:
  - factor_top_genes_v2.csv      top N genes per factor per view (RNA & METH)
  - factor_cross_view_v2.csv     per-factor cross-view concordance, meth-specificity
  - factor_summary_v2.csv        factor-level summary table (R2, theme tags, top genes)
  - figure_multiomics_v2.png     6-panel summary figure
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, mannwhitneyu, fisher_exact

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa"
META = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/m2_metadata.csv"
OUT  = BASE

mf   = pd.read_csv(f"{OUT}/mofa_multiomics_factors_v2.csv", index_col=0)
mw_r = pd.read_csv(f"{OUT}/mofa_multiomics_weights_RNA_v2.csv", index_col=0)
mw_m = pd.read_csv(f"{OUT}/mofa_multiomics_weights_METH_v2.csv", index_col=0)
r2   = pd.read_csv(f"{OUT}/mofa_multiomics_r2_per_view_v2.csv")

factors = [c for c in mf.columns if c.startswith("F")]
n_factors = len(factors)

# donor mapping
meta = pd.read_csv(META)
donor_map = meta.set_index("cell")["donor"].to_dict()
donors = [donor_map.get(c, "unknown") for c in mf.index]
donor_colors = {"AF1":"#E41A1C","AF2":"#377EB8","AF3":"#4DAF4A","AF4":"#984EA3","AF5":"#FF7F00"}

# ============== 1. Per-factor top genes (each view) ==============
TOP_N = 30
top_rows = []
for fi in factors:
    # RNA view: top positive & top negative
    rna_w = mw_r[fi]
    meth_w = mw_m[fi]
    rna_top_pos = rna_w.sort_values(ascending=False).head(TOP_N).index.tolist()
    rna_top_neg = rna_w.sort_values(ascending=True).head(TOP_N).index.tolist()
    meth_top_pos = meth_w.sort_values(ascending=False).head(TOP_N).index.tolist()
    meth_top_neg = meth_w.sort_values(ascending=True).head(TOP_N).index.tolist()
    top_rows.append({"factor": fi, "view": "RNA_pos",
                     "genes": ";".join(rna_top_pos),
                     "abs_mean_weight": float(np.abs(rna_w.loc[rna_top_pos]).mean())})
    top_rows.append({"factor": fi, "view": "RNA_neg",
                     "genes": ";".join(rna_top_neg),
                     "abs_mean_weight": float(np.abs(rna_w.loc[rna_top_neg]).mean())})
    top_rows.append({"factor": fi, "view": "METH_pos",
                     "genes": ";".join(meth_top_pos),
                     "abs_mean_weight": float(np.abs(meth_w.loc[meth_top_pos]).mean())})
    top_rows.append({"factor": fi, "view": "METH_neg",
                     "genes": ";".join(meth_top_neg),
                     "abs_mean_weight": float(np.abs(meth_w.loc[meth_top_neg]).mean())})
top_df = pd.DataFrame(top_rows)
top_df.to_csv(f"{OUT}/factor_top_genes_v2.csv", index=False)
print(f"Top genes: {top_df.shape}", flush=True)

# ============== 2. Cross-view concordance per factor ==============
# For each factor: Pearson r between RNA weight and METH weight across genes.
# Also: |METH weight| - |RNA weight| distribution => methylation-specific genes.
cross_rows = []
for fi in factors:
    rna_w = mw_r[fi].values
    meth_w = mw_m[fi].values
    r, pv = pearsonr(rna_w, meth_w)
    # genes where |METH weight| > |RNA weight| + 1 sd => "methylation-driven"
    rna_abs = np.abs(rna_w); meth_abs = np.abs(meth_w)
    delta = meth_abs - rna_abs
    thr = delta.mean() + delta.std()
    meth_driven = (delta > thr).sum()
    rna_driven = (-delta > thr).sum()
    cross_rows.append({
        "factor": fi,
        "rna_meth_corr_r": float(r),
        "rna_meth_corr_p": float(pv),
        "rna_mean_abs_w": float(rna_abs.mean()),
        "meth_mean_abs_w": float(meth_abs.mean()),
        "meth_specificity_delta": float(delta.mean()),
        "n_meth_driven_genes": int(meth_driven),
        "n_rna_driven_genes": int(rna_driven),
    })
cross_df = pd.DataFrame(cross_rows)
cross_df.to_csv(f"{OUT}/factor_cross_view_v2.csv", index=False)
print(f"Cross-view: {cross_df.shape}", flush=True)

# ============== 3. Factor-level summary ==============
# Pull R2 per view per factor
r2_pivot = r2.pivot(index="factor", columns="view", values="r2_pct").reset_index()
r2_pivot.columns.name = None
r2_pivot = r2_pivot.rename(columns={"RNA":"r2_RNA_pct", "METH":"r2_METH_pct"})
r2_pivot["r2_total_pct"] = r2_pivot["r2_RNA_pct"].fillna(0) + r2_pivot["r2_METH_pct"].fillna(0)

summary = r2_pivot.merge(cross_df, on="factor")
# Add top 5 gene labels (RNA positive) for each factor
top_label = []
for fi in factors:
    rna_top = mw_r[fi].sort_values(ascending=False).head(5).index.tolist()
    meth_top = mw_m[fi].sort_values(ascending=False).head(5).index.tolist()
    top_label.append({
        "factor": fi,
        "top_RNA_pos_genes": ";".join(rna_top),
        "top_METH_pos_genes": ";".join(meth_top),
    })
top_label_df = pd.DataFrame(top_label)
summary = summary.merge(top_label_df, on="factor")
summary.to_csv(f"{OUT}/factor_summary_v2.csv", index=False)
print(f"Summary: {summary.shape}", flush=True)
print(summary[["factor","r2_RNA_pct","r2_METH_pct","rna_meth_corr_r","n_meth_driven_genes","top_RNA_pos_genes","top_METH_pos_genes"]].to_string(index=False))

# ============== 4. Visualization ==============
fig, axes = plt.subplots(3, 2, figsize=(15, 17))
fig.suptitle("MOFA+ Multi-omics v2 (RNA + CpG Methylation): 32 GV Oocytes\n"
             "ARD + spike-slab weights, 7-factor model", fontsize=13, fontweight="bold")

# Panel A: Factor scores heatmap
ax = axes[0,0]
vmin, vmax = -2.5, 2.5
im = ax.imshow(mf.values.T, aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)
ax.set_yticks(range(n_factors))
ax.set_yticklabels(factors)
ax.set_xticks([])
ax.set_title("A: Factor scores (z-scaled)")
plt.colorbar(im, ax=ax, shrink=0.7)

# Panel B: R2 per view per factor (stacked bar)
ax = axes[0,1]
xpos = np.arange(n_factors)
rna_r2 = summary["r2_RNA_pct"].values
meth_r2 = summary["r2_METH_pct"].fillna(0).values
ax.bar(xpos, rna_r2, color="#377EB8", label="RNA", edgecolor="black")
ax.bar(xpos, meth_r2, bottom=rna_r2, color="#FF7F00", label="METH", edgecolor="black")
ax.set_xticks(xpos); ax.set_xticklabels(factors)
ax.set_ylabel("Variance explained (%)")
ax.set_title("B: Per-factor R² (stacked: RNA + METH)")
ax.legend(loc="upper right", fontsize=9)

# Panel C: Scores x Donor (top 2 factors)
ax = axes[1,0]
for di, donor in enumerate(["AF1","AF2","AF3","AF4","AF5"]):
    x_jit = np.random.normal(di, 0.08, size=sum([d==donor for d in donors]))
    vals = [mf.loc[c, "F1"] for c, d in zip(mf.index, donors) if d == donor]
    ax.scatter(x_jit, vals, c=donor_colors[donor], s=18, alpha=0.8, edgecolors="none")
ax.set_xticks(range(5)); ax.set_xticklabels(["AF1","AF2","AF3","AF4","AF5"])
ax.axhline(0, color="gray", ls="--", lw=0.5)
ax.set_xlabel("Donor"); ax.set_ylabel("F1 score")
ax.set_title("C: F1 scores × donor (biological variability)")

# Panel D: Cross-view correlation heatmap
ax = axes[1,1]
im2 = ax.imshow(cross_df["rna_meth_corr_r"].values.reshape(1,-1),
                aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(n_factors)); ax.set_xticklabels(factors)
ax.set_yticks([0]); ax.set_yticklabels(["RNA-METH r"])
for i, r in enumerate(cross_df["rna_meth_corr_r"].values):
    ax.text(i, 0, f"{r:+.2f}", ha="center", va="center", fontsize=9,
            color="white" if abs(r)>0.5 else "black")
ax.set_title("D: Cross-view weight correlation (Pearson r)")
plt.colorbar(im2, ax=ax, shrink=0.6)

# Panel E: Methylation-specificity per factor
ax = axes[2,0]
delta = cross_df["meth_specificity_delta"].values
colors = ["#FF7F00" if d > 0 else "#377EB8" for d in delta]
ax.barh(factors, delta, color=colors, edgecolor="black")
ax.axvline(0, color="gray", lw=0.8)
ax.set_xlabel("Δ |METH weight| − |RNA weight|  (positive = methylation-driven)")
ax.set_title("E: Methylation specificity per factor")
ax.invert_yaxis()

# Panel F: Top meth-driven genes per factor
ax = axes[2,1]
ax.axis("off")
text_lines = ["Methylation-driven top-3 genes per factor:", ""]
for fi in factors:
    rna_w = mw_r[fi]
    meth_w = mw_m[fi]
    delta = np.abs(meth_w) - np.abs(rna_w)
    top = delta.sort_values(ascending=False).head(3)
    text_lines.append(f"{fi}: {', '.join(top.index)}")
ax.text(0.05, 0.95, "\n".join(text_lines), transform=ax.transAxes,
        fontsize=8.5, family="monospace", va="top",
        bbox=dict(boxstyle="round", facecolor="lightyellow", edgecolor="gray"))

plt.tight_layout()
plt.savefig(f"{OUT}/figure_multiomics_v2.png", dpi=200, bbox_inches="tight")
print(f"Saved: {OUT}/figure_multiomics_v2.png", flush=True)

# ============== 5. Print key findings ==============
print("\n=== KEY FINDINGS (v2) ===")
for _, row in summary.iterrows():
    print(f"  {row['factor']}: R² RNA={row['r2_RNA_pct']:.1f}%  R² METH={row['r2_METH_pct']:.1f}%  "
          f"RNA-METH r={row['rna_meth_corr_r']:+.2f}  "
          f"meth-specific genes={int(row['n_meth_driven_genes'])}")
    print(f"     top RNA+: {row['top_RNA_pos_genes']}")
    print(f"     top METH+: {row['top_METH_pos_genes']}")

print("\nDONE.")