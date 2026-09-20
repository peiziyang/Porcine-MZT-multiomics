#!/usr/bin/env python
"""
Lineage-level regulon specificity analysis v2.
Strategy: annotate GSE168106 blastocyst cells by lineage using
canonical marker expression, then compute a simpler specificity score:
  specificity = mean_AUC_in_lineage / mean_AUC_in_other_lineages (fold-change)
Including a more reliable cluster annotation step.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/"
OUT = BASE + "m9_lineage_rss_figure.png"
OUT_TABLE = BASE + "m9_lineage_rss_results.csv"

# Load data
harmony = pd.read_csv(BASE + "m2_harmony_corrected.csv")
meta = pd.read_csv(BASE + "m2_metadata.csv")
auc = pd.read_csv(BASE + "pyscenic_out/aucell_scores_full.csv")
expr = pd.read_csv(BASE + "m7_hvg_symbol.csv", index_col=0)
stage_map = pd.read_csv(BASE + "m7_cell_stage_map.csv")

# Align cells
auc_cells = set(auc["Cell"].values)
harmony_cells = set(harmony["cell"].values)
common = list(auc_cells & harmony_cells)
print(f"Common cells: {len(common)}")

# Build label maps
cond_map = dict(zip(meta["cell"], meta["condition"]))
ds_map = dict(zip(meta["cell"], meta["dataset"]))
sl_map = dict(zip(stage_map["cell_id"], stage_map["stage_label"]))

# Subset to GSE168106 blastocyst cells (E6-E14, in vivo)
blast = [c for c in common 
         if cond_map.get(c) == "in_vivo"
         and sl_map.get(c, "").startswith("E")
         and int(sl_map[c].replace("E","")) >= 6]
print(f"Blastocyst cells (E6-E14, in vivo): {len(blast)}")

if len(blast) < 30:
    # Use all GSE168106 cells
    blast = [c for c in common if ds_map.get(c) == "GSE168106"]
    print(f"Fallback: all GSE168106: {len(blast)}")

# Cluster on Harmony HPCs
hpc_cols = [f"HPC{i}" for i in range(1, 51)]
blast_idx = [list(harmony["cell"]).index(c) for c in blast]
blast_hpc = harmony.iloc[blast_idx][hpc_cols].values.astype(float)

k = min(4, len(blast) // 10)
k = max(k, 2)
kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
clusters = kmeans.fit_predict(blast_hpc)

# Check marker expression per cluster
MARKERS = {"ICM": ["POU5F1","NANOG","SOX2","DPPA3"],
           "TE": ["CDX2","GATA3","KRT8","KRT18","ELF5"],
           "PrE": ["GATA4","GATA6","SOX17","PDGFRA"],
           "EPI": ["POU5F1","NANOG"]}
avail_m = [m for ms in MARKERS.values() for m in ms if m in expr.index]
avail_m = list(set(avail_m))
print(f"Available markers: {avail_m}")

expr_sub = expr.loc[avail_m, blast].T.values  # cells × markers
ex_scaled = StandardScaler().fit_transform(expr_sub)
marker_df = pd.DataFrame(ex_scaled, index=blast, columns=avail_m)
marker_df["cluster"] = clusters
clus_marker = marker_df.groupby("cluster").mean()
print("\nMarker profiles per cluster (z-scored):")
print(clus_marker.round(2))

# Auto-annotate
lin_scores = {}
for lin, mkrs in MARKERS.items():
    mk = [m for m in mkrs if m in clus_marker.columns]
    if mk:
        lin_scores[lin] = clus_marker[mk].mean(axis=1)
lin_df = pd.DataFrame(lin_scores)
anno = lin_df.idxmax(axis=1)
print("Annotation:", anno.to_dict())

# Alternative annotation: look at individual strong markers
# If cluster has high GATA4/GATA6 → PrE; high CDX2/GATA3 → TE; high POU5F1/NANOG → ICM
clust_anno = {}
for c in range(k):
    row = clus_marker.loc[c]
    # ICM: POU5F1 high
    if row.get("POU5F1", -2) > 0.3 and row.get("GATA4", -2) < 0:
        clust_anno[c] = "ICM"
    elif row.get("GATA4", -2) > 0.3 or row.get("GATA6", -2) > 0.3:
        clust_anno[c] = "PrE"
    elif row.get("CDX2", -2) > 0.3 or row.get("GATA3", -2) > 0.3:
        clust_anno[c] = "TE"
    elif row.get("KRT8", -2) > 0.3 and row.get("POU5F1", -2) > 0.3:
        clust_anno[c] = "ICM/TE mixed"
    else:
        clust_anno[c] = "Uncertain"
print("Rule-based annotation:", clust_anno)

# Use rule-based if possible, else use auto-annotation
use_anno = clust_anno if any(v != "Uncertain" for v in clust_anno.values()) else anno
cluster_df = pd.DataFrame({"cell": blast, "cluster": clusters})
cluster_df["lineage"] = cluster_df["cluster"].map(use_anno)

# AUCell join
regulon_cols = [c for c in auc.columns if c != "Cell"]
auc_clust = auc.merge(cluster_df, left_on="Cell", right_on="cell", how="inner")
print(f"\nAUC + cluster: {auc_clust.shape[0]} cells")

# Per-lineage mean AUC
mean_reg = auc_clust.groupby("lineage")[regulon_cols].mean()
lineages = list(mean_reg.index)
print(f"Lineages: {lineages}")

# Specificity score: for each regulon, max lineage mean / mean of others
spec_scores = {}
for reg in regulon_cols:
    vals = mean_reg[reg].values
    if vals.sum() > 0 and len(vals) > 1:
        max_idx = np.argmax(vals)
        max_val = vals[max_idx]
        other_mean = np.mean([v for i, v in enumerate(vals) if i != max_idx])
        fc = max_val / max(other_mean, 1e-10)
        spec_scores[reg] = (lineages[max_idx], fc)

# Top lineage-specific regulons
top_by_lineage = {}
for lin in lineages:
    lin_regs = [(reg, fc) for reg, (l, fc) in spec_scores.items() if l == lin]
    lin_regs.sort(key=lambda x: -x[1])
    top_by_lineage[lin] = lin_regs[:6]
    print(f"\nTop regulons in {lin}:")
    for reg, fc in lin_regs[:6]:
        print(f"  {reg}: FC={fc:.2f}")

# Figure: 2-panel
fig, axes = plt.subplots(1, 2, figsize=(14, max(5, len(avail_m)*0.35)))

# Panel 1: Marker heatmap
im1 = axes[0].imshow(clus_marker.values, cmap="RdBu_r", aspect="auto",
                      vmin=-2, vmax=2)
axes[0].set_yticks(range(len(clus_marker.index)))
axes[0].set_yticklabels([f"Cluster {c} ({use_anno.get(c,'?')})" for c in clus_marker.index], fontsize=8)
axes[0].set_xticks(range(len(clus_marker.columns)))
axes[0].set_xticklabels(clus_marker.columns, fontsize=8, rotation=45, ha="right")
axes[0].set_title("A. Cluster marker profiles (z-scored expression)", fontsize=10, fontweight="bold")
fig.colorbar(im1, ax=axes[0], shrink=0.8, label="z-score")

# Panel 2: Top lineage-specific regulons
all_top = []
seen = set()
for lin in lineages:
    for reg, fc in top_by_lineage.get(lin, []):
        if reg not in seen:
            all_top.append((reg, lin))
            seen.add(reg)

if all_top:
    top_regs = [t[0] for t in all_top][:20]
    # Build matrix: mean AUC × (lineage × regulon)
    plot_data = mean_reg[top_regs].T
    im2 = axes[1].imshow(plot_data.values, cmap="YlOrRd", aspect="auto")
    axes[1].set_yticks(range(len(plot_data.index)))
    axes[1].set_yticklabels([r.replace("(+)","") for r in plot_data.index], fontsize=7)
    axes[1].set_xticks(range(len(plot_data.columns)))
    axes[1].set_xticklabels(plot_data.columns, fontsize=9)
    axes[1].set_title("B. Lineage-specific regulon mean AUC", fontsize=10, fontweight="bold")
    fig.colorbar(im2, ax=axes[1], shrink=0.8, label="Mean AUCell score")

fig.suptitle("Blastocyst lineage regulon specificity in the porcine preimplantation atlas",
             fontsize=11, fontweight="bold")
plt.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"\nSaved {OUT}")

# Save results
rss_out = pd.DataFrame([(reg, lin, fc) for reg, (lin, fc) in spec_scores.items()],
                       columns=["regulon", "lineage", "FC"])
rss_out.to_csv(OUT_TABLE, index=False)
print(f"Saved {OUT_TABLE}")
