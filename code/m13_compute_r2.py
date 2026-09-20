#!/usr/bin/env python
"""Compute MOFA+ per-view variance explained with proper view scaling.

MOFA+ uses scale_views=True internally, so the per-view scaling is:
  X_internal = (X - mu_v) / sigma_v   per view
  W matrices are in the scaled units
  Reconstruction in original units: X_hat = (Z @ W^T) * sigma_v + mu_v

This script computes per-factor R² as:
  R²_k = Var(recon_k) / Var(Xc)   where Xc = X - mu_v
  And total R² per view: 1 - SS_res / SS_tot.
"""
import pandas as pd, numpy as np

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa"

rna = pd.read_csv(f"{BASE}/rna_view_for_mofa.csv", index_col=0)
meth = pd.read_csv(f"{BASE}/meth_view_for_mofa.csv", index_col=0)
factors = pd.read_csv(f"{BASE}/mofa_multiomics_factors_v2.csv", index_col=0)
w_rna = pd.read_csv(f"{BASE}/mofa_multiomics_weights_RNA_v2.csv", index_col=0)
w_meth = pd.read_csv(f"{BASE}/mofa_multiomics_weights_METH_v2.csv", index_col=0)

fcols = [c for c in factors.columns if c.startswith("F")]
Z = factors[fcols].values
W_rna = w_rna[fcols].values
W_meth = w_meth[fcols].values

def per_view_r2(X, Z, W):
    """X: (n, d), Z: (n, k), W: (d, k) — W in scaled units."""
    n, d = X.shape
    mu = X.mean(axis=0)
    sigma = X.std(axis=0).clip(min=1e-12)
    # Internal scaled data and reconstruction
    Xs = (X - mu) / sigma
    Xs_hat = Z @ W.T
    # Per-factor contribution
    ss_tot = (Xs ** 2).sum()
    ss_res_total = ((Xs - Xs_hat) ** 2).sum()
    total_r2 = 1 - ss_res_total / ss_tot
    per_factor = []
    for k in range(W.shape[1]):
        recon_k = np.outer(Z[:, k], W[:, k])
        # Variance contribution
        var_k = (recon_k ** 2).sum()
        per_factor.append({
            "k": k,
            "factor": fcols[k],
            "var_contribution_pct": float(var_k / ss_tot * 100),
        })
    return total_r2, per_factor

r2_rna, pf_rna = per_view_r2(rna.values, Z, W_rna)
r2_meth, pf_meth = per_view_r2(meth.values, Z, W_meth)

print(f"Total R² (all 7 factors):")
print(f"  RNA:  {r2_rna*100:.2f}%")
print(f"  METH: {r2_meth*100:.2f}%")
print()

# Per-factor per-view R² share
rows = []
for p in pf_rna:
    rows.append({"factor": p["factor"], "view": "RNA", "r2_pct": p["var_contribution_pct"]})
for p in pf_meth:
    rows.append({"factor": p["factor"], "view": "METH", "r2_pct": p["var_contribution_pct"]})
# Add totals
rows.append({"factor": "TOTAL", "view": "RNA",  "r2_pct": r2_rna * 100})
rows.append({"factor": "TOTAL", "view": "METH", "r2_pct": r2_meth * 100})

r2_df = pd.DataFrame(rows)
r2_df.to_csv(f"{BASE}/mofa_multiomics_r2_per_view_v2.csv", index=False)
print(r2_df.pivot(index="factor", columns="view", values="r2_pct").round(2).to_string())
print(f"\nSaved: {BASE}/mofa_multiomics_r2_per_view_v2.csv", flush=True)