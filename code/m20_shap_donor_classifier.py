#!/usr/bin/env python
"""XGBoost donor classifier + SHAP interpretation on 32 GV oocytes.

Rationale: We have 5 donors (AF1-AF5), 32 cells total. With small sample
size, leave-one-out cross-validation gives honest estimate. SHAP values
identify which genes drive donor classification.

KEY COMPARISON: Are MOFA+ F1 top genes also SHAP top genes?
If yes -> independent validation of MOFA+ factor interpretation.
If no -> donors carry different signal than MOFA+ factors capture.

Inputs:
  data/processed/m4_mofa/rna_view_for_mofa.csv   (32 cells x 2552 genes)
  data/processed/m2_metadata.csv                 (donor labels)
  data/processed/m4_mofa/factor_top_genes_v2.csv (MOFA+ top genes per factor)

Outputs:
  data/processed/m4_mofa/shap_analysis_results.csv
  data/processed/m4_mofa/shap_vs_mofa_comparison.csv
  data/processed/m4_mofa/figure_shap.png
"""
import os, sys
import pandas as pd, numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import xgboost as xgb
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import accuracy_score
import shap

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
MFA  = f"{BASE}/m4_mofa"
OUT  = MFA

# ============== 1. Load data ==============
print("Loading data ...", flush=True)
rna = pd.read_csv(f"{MFA}/rna_view_for_mofa.csv", index_col=0)
meta = pd.read_csv(f"{BASE}/m2_metadata.csv")
donor_map = meta.set_index("cell")["donor"].to_dict()

cells = rna.index.tolist()
X = rna.values.astype(np.float32)
y = np.array([donor_map[c] for c in cells])
print(f"  X: {X.shape}, y labels: {np.unique(y)}", flush=True)

donor_count = pd.Series(y).value_counts().to_dict()
print(f"  Donor distribution: {donor_count}", flush=True)

# ============== 2. Train XGBoost with LOO-CV ==============
print("\nTraining XGBoost with Leave-One-Out CV ...", flush=True)
loo = LeaveOneOut()
preds = np.zeros_like(y, dtype=object)
all_models = []

for fold, (train_idx, test_idx) in enumerate(loo.split(X)):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train = y[train_idx]
    
    # Label encoding
    classes = np.unique(y_train)
    y_enc = np.searchsorted(classes, y_train)
    
    model = xgb.XGBClassifier(
        n_estimators=100,
        max_depth=3,
        learning_rate=0.1,
        objective="multi:softprob",
        num_class=len(classes),
        random_state=42,
        verbosity=0,
    )
    model.fit(X_train, y_enc)
    all_models.append(model)
    
    # Predict test cell
    test_pred_enc = model.predict(X_test)[0]
    preds[test_idx[0]] = classes[test_pred_enc]

acc = accuracy_score(y, preds)
print(f"  LOO-CV accuracy: {acc*100:.1f}% ({sum(p==t for p,t in zip(preds,y))}/{len(y)})", flush=True)

# Confusion
print("  Predictions vs truth:")
for cell, true, pred in zip(cells, y, preds):
    if true != pred:
        print(f"    {cell}: true={true}, pred={pred}")

# ============== 3. SHAP analysis on full model ==============
print("\nTraining full model for SHAP ...", flush=True)
classes = np.unique(y)
y_enc = np.searchsorted(classes, y)

full_model = xgb.XGBClassifier(
    n_estimators=100,
    max_depth=3,
    learning_rate=0.1,
    objective="multi:softprob",
    num_class=len(classes),
    random_state=42,
    verbosity=0,
)
full_model.fit(X, y_enc)

# SHAP values
explainer = shap.TreeExplainer(full_model)
shap_values = explainer.shap_values(X)  # shape (n_cells, n_genes, n_classes) or (n_classes, n_cells, n_genes)
print(f"  SHAP shape: {np.array(shap_values).shape}", flush=True)

# For multi-class: take mean |SHAP| across classes
if len(shap_values.shape) == 3:  # (n_classes, n_cells, n_genes)
    mean_abs_shap = np.mean(np.abs(shap_values), axis=(0, 2))  # (n_genes,)
elif len(shap_values.shape) == 2:  # (n_cells, n_genes) - binary
    mean_abs_shap = np.mean(np.abs(shap_values), axis=0)

shap_df = pd.DataFrame({
    "gene": rna.columns,
    "mean_abs_shap": mean_abs_shap,
})
shap_df = shap_df.sort_values("mean_abs_shap", ascending=False)
shap_df.to_csv(f"{OUT}/shap_analysis_results.csv", index=False)
print(f"  SHAP top 10 genes: {shap_df.head(10)['gene'].tolist()}", flush=True)

# ============== 4. SHAP vs MOFA+ comparison ==============
print("\nComparing SHAP top genes with MOFA+ top genes ...", flush=True)
top_df = pd.read_csv(f"{MFA}/factor_top_genes_v2.csv")
# Get all unique top genes from MOFA+
mofa_genes = set()
for _, row in top_df.iterrows():
    mofa_genes.update(str(row["genes"]).split(";"))
print(f"  MOFA+ unique top genes: {len(mofa_genes)}", flush=True)

# SHAP top 50
shap_top50 = set(shap_df.head(50)["gene"])
print(f"  SHAP top 50: {len(shap_top50)}", flush=True)

# Overlap
overlap = shap_top50 & mofa_genes
print(f"  Overlap: {len(overlap)} genes", flush=True)
print(f"  Overlapping genes: {sorted(overlap)[:20]}", flush=True)

# Per-factor overlap
factor_overlap = []
for fac in top_df["factor"].unique():
    fac_sub = top_df[top_df["factor"]==fac]
    fac_genes = set()
    for _, r in fac_sub.iterrows():
        fac_genes.update(str(r["genes"]).split(";"))
    ov = shap_top50 & fac_genes
    factor_overlap.append({
        "factor": fac,
        "n_mofa_top_genes": len(fac_genes),
        "n_overlap_shap_top50": len(ov),
        "overlap_genes": ";".join(sorted(ov)),
    })
pd.DataFrame(factor_overlap).to_csv(f"{OUT}/shap_vs_mofa_comparison.csv", index=False)
print("\nPer-factor overlap with SHAP top 50:")
for r in factor_overlap:
    print(f"  {r['factor']}: {r['n_overlap_shap_top50']}/{r['n_mofa_top_genes']} -- {r['overlap_genes'][:80]}")

# ============== 5. Visualization ==============
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle("XGBoost Donor Classifier + SHAP on Pig GV Oocytes (n=32)",
             fontsize=14, fontweight="bold")

# A: SHAP summary barplot (top 30)
ax = axes[0, 0]
top30 = shap_df.head(30).sort_values("mean_abs_shap")
ax.barh(range(len(top30)), top30["mean_abs_shap"], color="#377EB8", edgecolor="black")
ax.set_yticks(range(len(top30)))
ax.set_yticklabels(top30["gene"], fontsize=8)
ax.set_xlabel("Mean |SHAP value| across cells")
ax.set_title("A: Top 30 SHAP feature importance")
ax.invert_yaxis()

# B: Confusion matrix (predicted vs actual donor)
ax = axes[0, 1]
classes = sorted(set(y))
n_cls = len(classes)
cm = np.zeros((n_cls, n_cls), dtype=int)
for t, p in zip(y, preds):
    ti = classes.index(t); pi = classes.index(p)
    cm[ti, pi] += 1
im = ax.imshow(cm, cmap="Blues")
ax.set_xticks(range(n_cls)); ax.set_xticklabels(classes)
ax.set_yticks(range(n_cls)); ax.set_yticklabels(classes)
ax.set_xlabel("Predicted"); ax.set_ylabel("True")
ax.set_title(f"B: LOO-CV Confusion (acc={acc*100:.1f}%)")
for i in range(n_cls):
    for j in range(n_cls):
        if cm[i,j] > 0:
            ax.text(j, i, str(cm[i,j]), ha="center", va="center", color="white" if cm[i,j]>1 else "black")
plt.colorbar(im, ax=ax, shrink=0.6)

# C: Per-factor overlap barplot
ax = axes[1, 0]
fac = pd.DataFrame(factor_overlap)
xpos = np.arange(len(fac))
ax.bar(xpos - 0.2, fac["n_mofa_top_genes"], 0.4, color="#377EB8", label="MOFA+ top genes")
ax.bar(xpos + 0.2, fac["n_overlap_shap_top50"], 0.4, color="#E41A1C", label="Overlap with SHAP top 50")
ax.set_xticks(xpos); ax.set_xticklabels(fac["factor"])
ax.set_ylabel("Gene count")
ax.set_title("C: MOFA+ vs SHAP top gene overlap per factor")
ax.legend()

# D: SHAP vs MOFA+ rank comparison (top 30 of each)
ax = axes[1, 1]
shap_top30_genes = set(shap_df.head(30)["gene"])
mofa_factor_top = {}
for _, row in top_df.iterrows():
    if row["factor"] == "F1" and row["view"] == "RNA_pos":
        for g in str(row["genes"]).split(";"):
            mofa_factor_top[g] = "F1_RNA+"
    elif row["factor"] == "F1" and row["view"] == "METH_pos":
        for g in str(row["genes"]).split(";"):
            mofa_factor_top[g] = "F1_METH+"

# Venn-like scatter
shap_ranks = {g: i for i, g in enumerate(shap_df["gene"])}
mofa_ranks = {g: i for i, g in enumerate(shap_df["gene"]) if g in mofa_factor_top}
common = list(set(shap_ranks.keys()) & set(mofa_factor_top.keys()))
if common:
    x = [shap_ranks[g] for g in common]
    y_ = [mofa_ranks[g] for g in common]
    ax.scatter(x, y_, color="#E41A1C", s=30, alpha=0.7, edgecolors="black")
    ax.set_xlabel("SHAP rank")
    ax.set_ylabel("MOFA+ F1 rank")
    ax.set_title("D: Genes in both SHAP top-30 and F1 top")
else:
    ax.text(0.5, 0.5, "No overlap", ha="center", va="center", transform=ax.transAxes)

plt.tight_layout()
plt.savefig(f"{OUT}/figure_shap.png", dpi=200, bbox_inches="tight")
print(f"\nSaved: {OUT}/figure_shap.png", flush=True)

# ============== 6. Summary ==============
print("\n=== SUMMARY ===")
print(f"LOO-CV accuracy: {acc*100:.1f}%")
print(f"Top 10 SHAP genes: {', '.join(shap_df.head(10)['gene'])}")
print(f"\nFactor-wise overlap (SHAP top 50 vs MOFA+):")
for r in factor_overlap:
    if r["n_overlap_shap_top50"] > 0:
        print(f"  {r['factor']}: {r['n_overlap_shap_top50']} overlapping genes")
print(f"\nTotal SHAP top 50 genes also in MOFA+ top: {len(shap_top50 & mofa_genes)}")

print("\nDONE.")