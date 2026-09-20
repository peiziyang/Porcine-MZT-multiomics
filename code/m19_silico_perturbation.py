#!/usr/bin/env python
"""In silico TF perturbation in pig GV oocytes.

Uses pySCENIC co-expression regulons (regulons_coexpr.csv) for TF-target edges,
estimates TF activity per cell via mean(standardized targets), and simulates
knockdown/overexpression via proportional target-gene adjustment.

Then projects perturbed cells into MOFA+ latent space and compares
with known PA embryo deviation direction.

Inputs:
  data/processed/pyscenic_out/regulons_coexpr.csv   (246 TFs x targets)
  data/processed/m4_mofa/rna_view_for_mofa.csv       (32 cells x 2552 genes)
  data/processed/m4_mofa/mofa_multiomics_factors_v2.csv
  data/processed/m4_mofa/mofa_multiomics_weights_RNA_v2.csv
"""
import os
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE   = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
SCENIC = f"{BASE}/pyscenic_out"
MFA    = f"{BASE}/m4_mofa"
OUT    = f"{BASE}/perturbation"
os.makedirs(OUT, exist_ok=True)

# ========== 1. Load regulons ==========
print("Loading regulons ...", flush=True)
reg_df = pd.read_csv(f"{SCENIC}/regulons_coexpr.csv")
regulons = {}
for _, r in reg_df.iterrows():
    tf = r["TF"]
    targets = [t.strip() for t in str(r["targets"]).split(",") if t.strip()]
    regulons[tf] = targets
print(f"  {len(regulons)} TFs loaded", flush=True)

# ========== 2. Load expression + MOFA ==========
print("Loading expression ...", flush=True)
rna = pd.read_csv(f"{MFA}/rna_view_for_mofa.csv", index_col=0)
cells = rna.index.tolist()
gene_list = rna.columns.tolist()
gene_idx = {g: i for i, g in enumerate(gene_list)}
n, g = len(cells), len(gene_list)
X = rna.values  # log-normalized

Z_df = pd.read_csv(f"{MFA}/mofa_multiomics_factors_v2.csv", index_col=0)
W_df = pd.read_csv(f"{MFA}/mofa_multiomics_weights_RNA_v2.csv", index_col=0)
fcols = [c for c in Z_df.columns if c.startswith("F")]
W = W_df[fcols].values  # (g, k)
Z = Z_df[fcols].values   # (n, k)

# ========== 3. Identify key TFs ==========
# TFs from MOFA+ F1 top genes that are actually TFs with regulons
key_tfs = ["DNMT1", "POU5F1", "ESRRA", "TFAP2C", "ATF3", "GATA4", "NFE2L3", "XBP1"]
key_tfs = [t for t in key_tfs if t in regulons]
print(f"Key TFs with regulons: {key_tfs}", flush=True)
for tf in key_tfs:
    in_data = sum(1 for t in regulons[tf] if t in gene_idx)
    print(f"  {tf}: {len(regulons[tf])} regulon targets, {in_data} in expression data", flush=True)

# ========== 4. Compute per-cell TF activity ==========
# Activity = mean(z-scored expression of all target genes)
print("Computing TF activities ...", flush=True)
tf_act = {}
for tf in key_tfs:
    targets_in = [t for t in regulons[tf] if t in gene_idx]
    if len(targets_in) < 3:
        print(f"  {tf}: insufficient targets ({len(targets_in)}), skip", flush=True)
        continue
    idxs = [gene_idx[t] for t in targets_in]
    Xt = X[:, idxs]
    Xt_z = (Xt - Xt.mean(axis=0)) / (Xt.std(axis=0) + 1e-10)
    tf_act[tf] = Xt_z.mean(axis=1)
    print(f"  {tf}: activity range [{tf_act[tf].min():.3f}, {tf_act[tf].max():.3f}]", flush=True)

# ========== 5. In silico perturbation ==========
def perturb(X_orig, tf, act_vec, effect_factor, regulon_targets, gene_index):
    """
    TF perturbation: for cell i, adjust target gene j by:
      Δ = correlation(TF, target) * effect_factor * act_vec[i] * σ_j
    where correlation is computed from the 32 cells, and σ_j is gene std.
    This produces a proportional shift in target expression.
    """
    targets_in = [t for t in regulon_targets if t in gene_index]
    if not targets_in:
        return X_orig.copy()
    
    # Compute TF-target correlation from data
    tf_idx = gene_index.get(tf)
    if tf_idx is None:
        return X_orig.copy()
    
    tf_expr = X_orig[:, tf_idx]
    corrs = []
    for t in targets_in:
        t_idx = gene_index[t]
        corr = np.corrcoef(tf_expr, X_orig[:, t_idx])[0, 1]
        corrs.append(max(abs(corr), 0.1) * np.sign(corr) if np.isfinite(corr) else 0.1)
    
    X_mod = X_orig.copy()
    for i in range(n):
        scale = effect_factor * act_vec[i]
        for j, t in enumerate(targets_in):
            t_idx = gene_index[t]
            delta = corrs[j] * scale * X_orig[i, t_idx]
            X_mod[i, t_idx] = np.clip(X_orig[i, t_idx] + delta, 0, 20)
    return X_mod

def project_mofa(X_mod, X_train_mean, X_train_std):
    """Project (n, g) into MOFA+ latent space (n, k)."""
    Xs = (X_mod - X_train_mean) / (X_train_std + 1e-10)
    # Z ≈ Xs @ W @ (W.T @ W)^-1
    WTW_inv = np.linalg.pinv(W.T @ W)
    return Xs @ W @ WTW_inv

X_mean = X.mean(axis=0)
X_std = X.std(axis=0)

print("Running perturbations ...", flush=True)
all_rows = []

EFFECTS = [("100%_KD", -1.0), ("50%_KD", -0.5), ("50%_OE", +0.5), ("100%_OE", +1.0)]

for tf in key_tfs:
    if tf not in tf_act:
        continue
    act = tf_act[tf]
    for label, factor in EFFECTS:
        X_mod = perturb(X, tf, act, factor, regulons[tf], gene_idx)
        Z_mod = project_mofa(X_mod, X_mean, X_std)
        shift = Z_mod - Z
        for ci, cell in enumerate(cells):
            row = {"cell": cell, "tf": tf, "perturbation": label,
                   "shift_magnitude": float(np.linalg.norm(shift[ci]))}
            for fi, fc in enumerate(fcols):
                row[f"shift_{fc}"] = shift[ci, fi]
            all_rows.append(row)
        avg = shift.mean(axis=0)
        print(f"  {tf:8s} {label:8s}: avg shift F1={avg[0]:+.3f} F2={avg[1]:+.3f} "
              f"|shift|={np.linalg.norm(avg):.3f}", flush=True)

shift_df = pd.DataFrame(all_rows)
shift_df.to_csv(f"{OUT}/perturbation_vectors.csv", index=False)
print(f"\nSaved {len(shift_df)} perturbation vectors", flush=True)

# ========== 6. PA deviation comparison ==========
# Load conserved PA-downregulated genes from earlier bulk validation
print("\nLoading PA deviation reference ...", flush=True)
try:
    deg_f = f"{BASE}/deseq2/m15_conserved_degs.csv"
    if os.path.exists(deg_f):
        degs = pd.read_csv(deg_f)
        if "gene" in degs.columns:
            pa_genes = degs[degs["padj"] < 0.05]["gene"]
            if pd.api.types.is_numeric_dtype(degs["log2FoldChange"]):
                pa_lfc = degs[degs["padj"] < 0.05].set_index("gene")["log2FoldChange"]
                # PA direction in gene space: genes that are DOWN in PA have negative LFC
                pa_dir = np.zeros(g)
                n_in = 0
                for gene in pa_genes:
                    if gene in gene_idx:
                        pa_dir[gene_idx[gene]] = -1.0 if pa_lfc.get(gene, 0) < 0 else 1.0
                        n_in += 1
                pa_dir /= (np.linalg.norm(pa_dir) + 1e-10)
                print(f"  PA direction: {n_in} genes mapped", flush=True)
            else:
                pa_dir = None
        else:
            pa_dir = None
    else:
        pa_dir = None
except Exception as e:
    print(f"  PA DEG load failed: {e}", flush=True)
    pa_dir = None

# Fallback
if pa_dir is None:
    pa_genes_fallback = ["KLHL15", "MT-ATP6", "MT-ND5", "IDH2", "PKM", "NOP9", "GNL3", "EXOSC9"]
    pa_dir = np.zeros(g)
    for gene in pa_genes_fallback:
        if gene in gene_idx:
            pa_dir[gene_idx[gene]] = -1.0
    pa_dir /= (np.linalg.norm(pa_dir) + 1e-10)
    print(f"  Using fallback PA genes", flush=True)

pa_z = W.T @ pa_dir  # (k,)

pd.DataFrame({"factor": fcols, "PA_direction": pa_z}).to_csv(
    f"{OUT}/pa_deviation_direction.csv", index=False)

# ========== 7. Cosine similarity: perturbation vs PA ==========
cos_vals = []
for tf in key_tfs:
    if tf not in tf_act: continue
    for label, _ in EFFECTS:
        sub = shift_df[(shift_df["tf"]==tf) & (shift_df["perturbation"]==label)]
        pert_vec = sub[[f"shift_{fc}" for fc in fcols]].mean().values
        cos = np.dot(pert_vec, pa_z) / (np.linalg.norm(pert_vec)*np.linalg.norm(pa_z)+1e-10)
        cos_vals.append({"tf": tf, "perturbation": label, "cos_sim": float(cos),
                         "shift_mag": sub["shift_magnitude"].mean()})

cos_df = pd.DataFrame(cos_vals)
cos_df.to_csv(f"{OUT}/perturbation_vs_pa.csv", index=False)
cos_df_sorted = cos_df.sort_values("cos_sim", ascending=False)
print("\n=== Perturbation vs PA alignment (cosine similarity) ===")
for _, r in cos_df_sorted.head(10).iterrows():
    print(f"  {r['tf']:8s} {r['perturbation']:8s}: cos_sim={r['cos_sim']:+.4f}")

# ========== 8. Visualization ==========
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle("In Silico TF Perturbation in Pig GV Oocytes (MOFA+ Latent Space)",
             fontsize=14, fontweight="bold")

# A: DNMT1 perturbation → shift in F1-F2
ax = axes[0,0]
dnmt1 = shift_df[shift_df["tf"]=="DNMT1"]
cc = {"100%_KD":"#E41A1C","50%_KD":"#FF7F00","50%_OE":"#377EB8","100%_OE":"#4DAF4A"}
mm = {"100%_KD":"v","50%_KD":"s","50%_OE":"^","100%_OE":"o"}
for label in EFFECTS:
    l = label[0]
    sub = dnmt1[dnmt1["perturbation"]==l]
    ax.scatter(sub["shift_F1"], sub["shift_F2"], c=cc[l], marker=mm[l],
               s=40, alpha=0.85, label=l, edgecolors="black", linewidth=0.3)
ax.axhline(0,color="gray",ls="--",lw=0.5); ax.axvline(0,color="gray",ls="--",lw=0.5)
ax.set_xlabel("Δ F1 (RNA+METH shared axis)"); ax.set_ylabel("Δ F2 (RNA-specific)")
ax.set_title("A: DNMT1 perturbation → MOFA+ F1-F2 shift")
ax.legend(fontsize=8, loc="lower right")

# B: Perturbation vs PA alignment
ax = axes[0,1]
top10 = cos_df_sorted.head(10)
colors = ["#E41A1C" if "KD" in r["perturbation"] else "#377EB8"
          for _, r in top10.iterrows()]
labels = [f"{r['tf']}_{r['perturbation']}" for _, r in top10.iterrows()]
ax.barh(range(len(top10)), top10["cos_sim"], color=colors, edgecolor="black", height=0.6)
ax.set_yticks(range(len(top10)))
ax.set_yticklabels(labels, fontsize=7.5)
ax.axvline(0, color="gray", lw=0.8)
ax.set_xlabel("Cosine similarity with PA deviation vector")
ax.set_title("B: Top perturbation vs PA deviation alignment")
ax.invert_yaxis()

# C: Factor-level shift heatmap (all TFs, 100% KD)
ax = axes[1,0]
kd_tfs = [t for t in key_tfs if t in tf_act]
kd_data = np.zeros((len(kd_tfs), len(fcols)))
for ti, tf in enumerate(kd_tfs):
    sub = shift_df[(shift_df["tf"]==tf) & (shift_df["perturbation"]=="100%_KD")]
    for fi, fc in enumerate(fcols):
        kd_data[ti, fi] = sub[f"shift_{fc}"].mean()
im = ax.imshow(kd_data.T, aspect="auto", cmap="RdBu_r", vmin=-0.8, vmax=0.8)
ax.set_xticks(range(len(kd_tfs)))
ax.set_xticklabels(kd_tfs, rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(len(fcols))); ax.set_yticklabels(fcols)
ax.set_title("C: 100% KD effect per factor per TF")
for i in range(len(kd_tfs)):
    for j in range(len(fcols)):
        ax.text(i, j, f"{kd_data[i,j]:.2f}", ha="center", va="center", fontsize=6)
plt.colorbar(im, ax=ax, shrink=0.7)

# D: Summary magnitude per TF
ax = axes[1,1]
for ti, tf in enumerate(kd_tfs):
    sub_k = shift_df[(shift_df["tf"]==tf)&(shift_df["perturbation"]=="100%_KD")]
    sub_o = shift_df[(shift_df["tf"]==tf)&(shift_df["perturbation"]=="100%_OE")]
    mag_k = sub_k["shift_magnitude"].mean()
    mag_o = sub_o["shift_magnitude"].mean()
    ax.bar(ti-0.15, mag_k, 0.3, color="#E41A1C", alpha=0.8, label="KD" if ti==0 else "")
    ax.bar(ti+0.15, mag_o, 0.3, color="#377EB8", alpha=0.8, label="OE" if ti==0 else "")
ax.set_xticks(range(len(kd_tfs)))
ax.set_xticklabels(kd_tfs, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Mean |shift| in MOFA+ space")
ax.set_title("D: Overall perturbation magnitude")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(f"{OUT}/figure_perturbation.png", dpi=200, bbox_inches="tight")
print(f"\nSaved: {OUT}/figure_perturbation.png", flush=True)

# ========== 9. Summary ==========
print("\n=== SUMMARY ===")
best = cos_df_sorted.iloc[0]
print(f"Best alignment: {best['tf']} {best['perturbation']} — cos_sim = {best['cos_sim']:.4f}")
print(f"This means: {best['tf']} perturbation shifts cells in a direction that is")
print(f"  {'aligned with' if best['cos_sim']>0 else 'opposite to'} the PA deviation direction.")
if best['cos_sim'] > 0.5:
    print("  → STRONG SUPPORT for in silico prediction matching biological observation!")
elif best['cos_sim'] > 0.3:
    print("  → MODERATE SUPPORT — directionally consistent.")
else:
    print("  → WEAK or NO alignment — perturbation may affect different pathways than PA.")
print("\nDONE.")