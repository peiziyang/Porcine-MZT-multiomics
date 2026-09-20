#!/usr/bin/env python
"""MZT trajectory docking: project GV multi-omics MOFA+ F1 onto E0-E10 development.

Hypothesis:
  The GV oocyte multi-omics F1 (shared RNA+METH program) reflects maternal
  preparation for MZT. If true, projecting F1 weights onto E0-E10 single-cell
  RNA should show:
    1. F1 projection score DECLINES during MZT (maternal clearance)
    2. F2 (RNA-only cytoskeletal) follows a different trajectory
    3. The decline rate correlates with known ZGA timing

Inputs:
  data/processed/m7_expr_cellsxgenes.csv         (1833 cells x 12548 genes)
  data/processed/m2_metadata.csv                  (stage, condition per cell)
  data/processed/m4_mofa/mofa_multiomics_weights_RNA_v2.csv   (MOFA+ W)
  data/processed/m4_mofa/mofa_multiomics_factors_v2.csv       (MOFA+ Z)

Outputs:
  data/processed/m4_mofa/mzt_trajectory_projection.csv
  data/processed/m4_mofa/mzt_trajectory_figure.png
"""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"

# ====== 1. Load data ======
print("Loading data ...", flush=True)
expr = pd.read_csv(f"{BASE}/m7_expr_cellsxgenes.csv", index_col=0)  # (1833, 12548)
meta = pd.read_csv(f"{BASE}/m2_metadata.csv")
w_rna = pd.read_csv(f"{BASE}/m4_mofa/mofa_multiomics_weights_RNA_v2.csv", index_col=0)
z_gv = pd.read_csv(f"{BASE}/m4_mofa/mofa_multiomics_factors_v2.csv", index_col=0)

fcols = [c for c in z_gv.columns if c.startswith("F")]
print(f"  Expression: {expr.shape}", flush=True)
print(f"  Metadata: {meta.shape}", flush=True)
print(f"  MOFA+ W_RNA: {w_rna.shape}, {len(fcols)} factors", flush=True)

# ====== 2. Harmonize genes between atlas and MOFA+ ======
mo_genes = list(w_rna.index)  # 2552 gene symbols
atlas_genes = list(expr.columns)  # 12548 gene symbols

common_genes = sorted(set(mo_genes) & set(atlas_genes))
print(f"  Common genes: {len(common_genes)} of {len(mo_genes)} MOFA+ genes", flush=True)

# Subset both
expr = expr[common_genes]  # (1833, ~2500)
w_rna = w_rna.loc[common_genes]  # (~2500, k)

# Add stage info
cell_stage = meta.set_index("cell")["stage"].to_dict()
expr["stage"] = [cell_stage.get(c, "Unknown") for c in expr.index]

# ====== 3. Filter to E0-E10 developmental cells ======
valid_stages = ["E0","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"]
e_stage_mask = expr["stage"].isin(valid_stages)
expr_dev = expr[e_stage_mask].copy()
stage_order = expr_dev.pop("stage")
print(f"  E0-E10 cells: {len(expr_dev)}", flush=True)

# ====== 4. Project: Z_proj = X @ W (in standardized space) ======
X = expr_dev.values
X_mean = X.mean(axis=0)
X_std = X.std(axis=0).clip(min=1e-10)
Xs = (X - X_mean) / X_std
W = w_rna[fcols].values

Z_proj = Xs @ W  # (n_cells, n_factors) — projected scores

proj_df = pd.DataFrame(Z_proj, index=expr_dev.index, columns=fcols)
proj_df["stage"] = list(stage_order)

# Per-stage mean
stage_mean = proj_df.groupby("stage")[fcols].mean()

# Build output
proj_df.to_csv(f"{BASE}/m4_mofa/mzt_trajectory_projection.csv")
print(f"  Saved projection: {proj_df.shape}", flush=True)
stage_mean.to_csv(f"{BASE}/m4_mofa/mzt_stage_mean_projection.csv")
print(f"  Stages: {sorted(stage_mean.index.tolist())}", flush=True)

# ====== 5. Correlation: projection score vs developmental ordinal ======
# Map stage to numeric ordinal
stage_ordinal = {}
for i, s in enumerate(["E0","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"]):
    stage_ordinal[s] = i

proj_df["ordinal"] = proj_df["stage"].map(stage_ordinal)

print("\nProjected factor scores vs developmental time:", flush=True)
for fc in fcols:
    r, p = spearmanr(proj_df[fc], proj_df["ordinal"])
    direction = "↑" if r > 0 else "↓"
    print(f"  {fc}: r={r:+.3f}, p={p:.1e} {direction}  (mean E0={stage_mean.loc['E0',fc]:+.2f}, E10={stage_mean.loc['E10',fc]:+.2f})",
          flush=True)

# ====== 6. Visualization ======
# Figure: trajectory of F1-F7 projection scores across E0-E10
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle("GV Oocyte Multi-omics MOFA+ F1-F7 Projected onto E0-E10 Atlas",
             fontsize=14, fontweight="bold")

# Panel A: F1-F4 line plots across stages
ax = axes[0, 0]
colors = ["#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#999999"]
for fi, fc in enumerate(fcols):
    vals = stage_mean[fc]
    stages = stage_mean.index.tolist()
    ax.plot(stages, vals, "o-", color=colors[fi], label=fc, linewidth=2, markersize=6)
ax.set_xlabel("Developmental Stage")
ax.set_ylabel("Mean projected factor score")
ax.set_title("A: Projected MOFA+ factor scores across E0-E10")
ax.legend(fontsize=8)
ax.axhline(0, color="gray", ls="--", lw=0.5)

# Panel B: F1 decline — highlight ZGA window
ax = axes[0, 1]
f1_vals = stage_mean["F1"]
stages = stage_mean.index.tolist()
# Color by zone: E0-E2 (maternal), E3/E4 (ZGA), E5-E10 (post-ZGA)
zone_colors = []
for s in stages:
    if int(s[1:]) <= 2: zone_colors.append("#E41A1C")
    elif int(s[1:]) <= 4: zone_colors.append("#FF7F00")
    else: zone_colors.append("#377EB8")
ax.bar(stages, f1_vals, color=zone_colors, edgecolor="black")
# Annotation
r_f1, p_f1 = spearmanr(proj_df["F1"], proj_df["ordinal"])
ax.text(0.5, 0.95, f"Spearman r={r_f1:.3f}, p={p_f1:.1e}",
        transform=ax.transAxes, ha="center", fontsize=11, fontweight="bold")
ax.set_xlabel("Developmental Stage")
ax.set_ylabel("F1 projected score")
ax.set_title("B: F1 (shared RNA+METH) — maternal program clearance")
# Legend
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#E41A1C",label="Maternal (E0-E2)"),
                   Patch(color="#FF7F00",label="ZGA (E3-E4)"),
                   Patch(color="#377EB8",label="Post-ZGA (E5-E10)")],
          fontsize=8)

# Panel C: F1 vs F2 divergence over stages
ax = axes[1, 0]
f1 = stage_mean["F1"].values
f2 = stage_mean["F2"].values
for i, s in enumerate(stages):
    ax.scatter(f1[i], f2[i], s=80+int(s[1:])*10, c="#E41A1C", alpha=0.7, edgecolors="black")
    ax.annotate(s, (f1[i], f2[i]), fontsize=8, ha="left")
ax.set_xlabel("F1 score (shared maternal)")
ax.set_ylabel("F2 score (RNA-specific)")
ax.axhline(0,color="gray",ls="--",lw=0.5)
ax.axvline(0,color="gray",ls="--",lw=0.5)
ax.set_title("C: F1 vs F2 — maternal → ZGA transition")

# Panel D: F1 projection heatmap (cells ordered by stage)
ax = axes[1, 1]
# Sort cells by stage ordinal
cell_order = proj_df.sort_values("ordinal").index
stage_labels = [s for s in cell_order for _ in range(1)]  # just cell labels
Z_sorted = proj_df.loc[cell_order, fcols].values.T
im = ax.imshow(Z_sorted, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
# Add stage boundaries
prev = ""
boundary_x = []
for i, c in enumerate(cell_order):
    s = proj_df.loc[c, "stage"]
    if s != prev:
        prev = s
        boundary_x.append(i)
for bx in boundary_x:
    ax.axvline(bx, color="black", lw=0.3, alpha=0.5)
ax.set_yticks(range(len(fcols)))
ax.set_yticklabels(fcols)
ax.set_xticks([])
ax.set_title(f"D: Per-cell projection (E0-E10, n={len(cell_order)})")
plt.colorbar(im, ax=ax, shrink=0.7)

plt.tight_layout()
plt.savefig(f"{BASE}/m4_mofa/mzt_trajectory_figure.png", dpi=200, bbox_inches="tight")
print(f"\nSaved: {BASE}/m4_mofa/mzt_trajectory_figure.png", flush=True)

# ====== 7. Summary ======
print("\n=== MZT TRAJECTORY SUMMARY ===")
print(f"E0-E10 cells projected: {len(proj_df)}")
for fc in fcols:
    r, p = spearmanr(proj_df[fc], proj_df["ordinal"])
    print(f"  {fc}: Spearman r={r:+.3f} p={p:.1e} — "
          f"{'declines with development' if r<-0.2 else 'rises' if r>0.2 else 'no clear trend'}")
print("\nDONE.")