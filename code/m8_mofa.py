#!/usr/bin/env python
"""MOFA+ (mofapy2) multi-group factor analysis of the porcine preimplantation atlas.

Design (respecting mofapy2's data model):
  - 1 RNA view, 3 groups = the 3 studies with retained gene-level HVG expression:
      GSE168106 (in-vivo reference), GSE112380 (pre-gastrulation), GSE234116 (oocyte).
  - IVF/PA (GSE164812, n=21/20) are EXCLUDED: their per-cell expression was not
    retained in the HVG matrix; their divergence is quantified separately in
    Fig.2 (M5), Fig.5 (regulons) and Fig.7 (Waddington-OT displacement).
  - Features = 3000 HVGs (shared gene order across groups).
  - Factors = 10 (ARD-pruned). Seed fixed for reproducibility.
Outputs:
  m8_mofa_factors.csv, m8_mofa_loadings.csv, m8_mofa_r2.csv,
  m8_mofa_topgenes.csv, m8_mofa_factor_by_stage.csv, m8_mofa_factor_by_dataset.csv
"""
import os, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/"
SEED = 42
GENES_USED = 3000
FACTORS = 10

# ---- 1. load HVG expression (genes x cells) and cell->dataset map ----
hvg = pd.read_csv(BASE + "m7_hvg_symbol.csv", index_col=0)   # 3000 genes x 1833 cells
genes = list(hvg.index)
cells = list(hvg.columns)
X = hvg.values.T                                               # cells x genes (1833 x 3000)
X = np.asarray(X, dtype=np.float64)

stage_map = pd.read_csv(BASE + "m7_cell_stage_map.csv")        # cell_id, dataset, stage_label, stage_order
cell2ds = dict(zip(stage_map["cell_id"], stage_map["dataset"]))
cell2stage = dict(zip(stage_map["cell_id"], stage_map["stage_label"]))
cell2order = dict(zip(stage_map["cell_id"], stage_map["stage_order"]))

GROUPS = ["GSE168106", "GSE112380", "GSE234116"]
GROUP_LABEL = {"GSE168106": "in_vivo (GSE168106)",
               "GSE112380": "GSE112380 (pre-gastrulation)",
               "GSE234116": "GSE234116 (oocyte)"}

# Build per-group matrices with identical gene order
mats, names_per_group = [], []
for g in GROUPS:
    gcells = [c for c in cells if cell2ds.get(c) == g]
    mats.append(X[[cells.index(c) for c in gcells], :])
    names_per_group.append(gcells)
    print(f"  group {g}: n={len(gcells)} cells")

# mofapy2 expects data[view][group] = (samples, features)
data = [mats]                       # 1 view, 3 groups
views_names = ["RNA"]
groups_names = [GROUP_LABEL[g] for g in GROUPS]
features_names = [genes]            # one per view
samples_names = [list(np.concatenate(names_per_group))] if False else names_per_group
# samples_names must be a nested list length G, each = sample names for that group
samples_names = names_per_group

# ---- 2. run MOFA+ ----
from mofapy2.run.entry_point import entry_point
ent = entry_point()
ent.set_data_options(scale_groups=True, scale_views=False, center_groups=True)
ent.set_data_matrix(data, likelihoods=["gaussian"] * len(views_names),
                    views_names=views_names, groups_names=groups_names,
                    samples_names=samples_names, features_names=features_names)
ent.set_model_options(factors=FACTORS, ard_weights=True, ard_factors=False, spikeslab_weights=True)
ent.set_train_options(iter=1000, convergence_mode="fast", seed=SEED, quiet=True)
ent.build()
ent.run()

# ---- 3. extract factors & loadings ----
Z = ent.model.nodes["Z"].getExpectations()["E"]                 # (N, K)
Wexp = ent.model.nodes["W"].getExpectations()
W = Wexp[0]["E"]                                                 # (D, K) shared across groups
N, K = Z.shape
print(f"  factors Z: {Z.shape}, loadings W: {W.shape}")

# global cell order = concatenation of groups in the order provided
global_cells = []
for gnames in names_per_group:
    global_cells.extend(gnames)
assert len(global_cells) == N, f"cell order mismatch {len(global_cells)} vs {N}"

# ---- 4. variance explained per factor (proportion of total variance) ----
# total variance of data (per-gene, summed across genes), using centered Y
Yc = X - X.mean(axis=0, keepdims=True)
total_var = float(np.sum(Yc.var(axis=0)))
r2 = []
for k in range(K):
    comp = np.outer(Z[:, k], W[:, k])        # N x D reconstruction component
    comp_var = float(np.sum(comp.var(axis=0)))
    r2.append(comp_var / total_var)
r2 = np.array(r2)
print("  R2 per factor:", np.round(r2, 4))

# ---- 5. top genes per factor ----
top_genes = {}
for k in range(K):
    order = np.argsort(-np.abs(W[:, k]))
    tg = [(genes[i], float(W[i, k])) for i in order[:25]]
    top_genes[k] = tg

# ---- 6. factor score by dataset & by stage (effect sizes) ----
ds_of_cell = [cell2ds.get(c) for c in global_cells]
stage_of_cell = [cell2stage.get(c) for c in global_cells]
order_of_cell = [cell2order.get(c) for c in global_cells]

df_factors = pd.DataFrame(Z, columns=[f"F{k+1}" for k in range(K)])
df_factors.insert(0, "cell", global_cells)
df_factors["dataset"] = ds_of_cell
df_factors["stage_label"] = stage_of_cell
df_factors["stage_order"] = order_of_cell

# factor means per dataset (group offset) -> effect size
fac_by_ds = df_factors.groupby("dataset")[[f"F{k+1}" for k in range(K)]].mean()
# factor means per stage (developmental axis)
fac_by_stage = df_factors.dropna(subset=["stage_order"]).groupby("stage_label")[
    [f"F{k+1}" for k in range(K)]].mean()
fac_by_stage_order = df_factors.dropna(subset=["stage_order"]).groupby("stage_order")[
    [f"F{k+1}" for k in range(K)]].mean().sort_index()

# ---- 7. save ----
df_factors.to_csv(BASE + "m8_mofa_factors.csv", index=False)
pd.DataFrame(W, index=genes, columns=[f"F{k+1}" for k in range(K)]).to_csv(
    BASE + "m8_mofa_loadings.csv")
pd.DataFrame({"factor": [f"F{k+1}" for k in range(K)], "R2": r2,
             "R2_pct": r2 * 100}).to_csv(BASE + "m8_mofa_r2.csv", index=False)
tg_rows = []
for k in range(K):
    for rank, (g, w) in enumerate(top_genes[k], 1):
        tg_rows.append({"factor": f"F{k+1}", "rank": rank, "gene": g, "weight": w})
pd.DataFrame(tg_rows).to_csv(BASE + "m8_mofa_topgenes.csv", index=False)
fac_by_ds.to_csv(BASE + "m8_mofa_factor_by_dataset.csv")
fac_by_stage.to_csv(BASE + "m8_mofa_factor_by_stage.csv")

# ---- 8. diagnostics to stdout ----
print("\n=== Factor R2 (variance explained, %) ===")
for k in range(K):
    print(f"  F{k+1}: {r2[k]*100:.2f}%")
print("\n=== Top 10 genes on highest-R2 factor ===")
topk = int(np.argmax(r2))
for g, w in top_genes[topk][:10]:
    print(f"  {g:12s} {w:+.3f}")
print("\n=== Factor mean by dataset (group offset) ===")
print(fac_by_ds.round(3).to_string())
print("\n=== Factor mean by stage (top 3 factors) ===")
print(fac_by_stage_order[[f"F{k+1}" for k in range(min(3, K))]].round(3).to_string())
print("\nDONE. Seed=%d  N=%d  genes=%d  factors=%d" % (SEED, N, len(genes), K))
