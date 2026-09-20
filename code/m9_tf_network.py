#!/usr/bin/env python
"""
Core TF regulatory network plot for Fig.5 (new panel or supplementary).
Reads pySCENIC regulon co-expression edges and plots a network of top TFs
with their target genes.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/"
OUT = BASE + "m9_tf_network.png"

# 1. Load regulon co-expression data
rc = pd.read_csv(BASE + "pyscenic_out/regulons_coexpr.csv")
print(f"Loaded {rc.shape[0]} TFs with co-expression targets")

# 2. Select top TFs by AUCell mean activity
au = pd.read_csv(BASE + "pyscenic_out/aucell_scores_full.csv")
regulon_cols = [c for c in au.columns if c != "Cell"]
mean_auc = au[regulon_cols].mean().sort_values(ascending=False)
print("\nTop regulons by mean AUC:")
for r, v in mean_auc.head(10).items():
    print(f"  {r}: {v:.4f}")

# Pick diverse TFs for the network (covering different biological programs)
TOP_TFS = ["ESRRA", "POU5F1", "MXD3", "GATA4", "PAX6", "TFAP2C"]
available_tfs = [t for t in TOP_TFS if t in rc["TF"].values]
print(f"\nAvailable core TFs: {available_tfs}")

# 3. Build graph
G = nx.DiGraph()
for _, row in rc[rc["TF"].isin(available_tfs)].iterrows():
    tf = row["TF"]
    targets = str(row["targets"]).split(",")[:20]  # top 20 targets per TF
    G.add_node(tf, type="TF", size=800)
    for tgt in targets:
        tgt = tgt.strip()
        if tgt:
            G.add_node(tgt, type="target", size=100)
            G.add_edge(tf, tgt)

# 4. Plot
plt.figure(figsize=(13, 10))
# Use a more spread out layout (kamada_kawai often works better)
pos = nx.spring_layout(G, k=2.5, seed=42, iterations=100)

# Node colors
tf_nodes = [n for n, attr in G.nodes(data=True) if attr["type"] == "TF"]
target_nodes = [n for n, attr in G.nodes(data=True) if attr["type"] == "target"]

nx.draw_networkx_nodes(G, pos, nodelist=tf_nodes, node_size=1800,
                       node_color="#E74C3C", edgecolors="black", linewidths=1.5)
nx.draw_networkx_nodes(G, pos, nodelist=target_nodes, node_size=60,
                       node_color="#85C1E9", alpha=0.7)
nx.draw_networkx_edges(G, pos, alpha=0.3, arrows=True, arrowsize=10,
                       edge_color="gray", width=0.5)
# Labels with white background for readability
nx.draw_networkx_labels(G, pos, labels={n: n for n in tf_nodes},
                        font_size=11, font_weight="bold", font_color="white")

# Legend
import matplotlib.patches as mpatches
legend_TF = mpatches.Patch(color="#E74C3C", label="Transcription factor (regulon master)")
legend_target = mpatches.Patch(color="#85C1E9", label="Target gene")
plt.legend(handles=[legend_TF, legend_target], loc="lower left", fontsize=10, framealpha=0.9)

plt.title("Core transcription factor regulatory networks in the porcine preimplantation embryo\n(6 TFs × 20 targets each; co-expression edges from pySCENIC GRNBoost2)",
          fontsize=12, fontweight="bold")
plt.axis("off")
plt.tight_layout()
plt.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"Saved {OUT}")
print(f"Graph: {len(tf_nodes)} TFs, {len(target_nodes)} targets, {G.number_of_edges()} edges")
