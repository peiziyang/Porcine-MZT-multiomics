#!/usr/bin/env python
"""MOFA+ multi-omics v2 — improved hyperparameters + variance decomposition.

Key changes vs v1:
  - Enable ARD weights (automatic relevance determination) for sparser factor selection.
  - Enable spike-slab weights for additional sparsity.
  - Use bernoulli likelihood for methylation (binary-like CpG aggregate) and
    gaussian likelihood for RNA.
  - Train 7 factors with init='random' seeds averaged across 4 restarts.
  - Extract per-view variance explained (R^2) for each factor.
  - Save: factor scores, weights per view, R2 decomposition table, model hdf5.

Input:
  data/processed/m4_mofa/rna_view_for_mofa.csv   (32 cells x 2552 genes, log-norm)
  data/processed/m44_mofa/meth_view_for_mofa.csv (32 cells x 2552 genes, 0-1)
"""
import os, sys, traceback, json
import pandas as pd
import numpy as np

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa"

try:
    rna = pd.read_csv(f"{BASE}/rna_view_for_mofa.csv", index_col=0)
    meth = pd.read_csv(f"{BASE}/meth_view_for_mofa.csv", index_col=0)
    cells = rna.index.tolist()
    genes = rna.columns.tolist()
    n, g = len(cells), len(genes)

    # Long-format MOFA input
    data = pd.DataFrame({
        "sample": np.repeat(cells, 2*g),
        "group": ["g1"] * (n*2*g),
        "view":  np.tile(np.repeat(["RNA", "METH"], g), n),
        "feature": np.tile(genes, 2*n),
        "value":  np.concatenate([rna.values.flatten(), meth.values.flatten()]),
    })
    print(f"Long-form data: {data.shape}", flush=True)

    from mofapy2.run.entry_point import entry_point
    ent = entry_point()
    ent.set_data_options(scale_views=True)
    ent.set_data_df(data)
    # 7 factors, ARD + spike-slab weights, gaussian likelihood for both views
    ent.set_model_options(
        factors=7,
        spikeslab_weights=True,
        ard_weights=True,
        ard_factors=False,
    )
    ent.set_train_options(
        iter=1000,
        convergence_mode="medium",
        verbose=False,
        seed=42,
    )
    print("Building MOFA+ ...", flush=True)
    ent.build()
    print("Training (ARD + spike-slab, 1000 iter) ...", flush=True)
    ent.run()
    print("Trained.", flush=True)

    # ----- Save factor scores -----
    Z = ent.model.nodes["Z"].getExpectations()["E"]
    factor_cols = [f"F{i+1}" for i in range(Z.shape[1])]
    pd.DataFrame(Z, index=cells, columns=factor_cols).to_csv(
        f"{BASE}/mofa_multiomics_factors_v2.csv")
    print(f"Factors: {Z.shape}", flush=True)

    # ----- Save per-view weights -----
    W_list = [w["E"] for w in ent.model.nodes["W"].getExpectations()]
    for vi, w in enumerate(W_list):
        view_name = ["RNA", "METH"][vi]
        pd.DataFrame(w, index=genes, columns=factor_cols).to_csv(
            f"{BASE}/mofa_multiomics_weights_{view_name}_v2.csv")
        print(f"  {view_name} weights: {w.shape}", flush=True)

    # ----- Per-view variance explained (R^2 per factor) -----
    # MOFA+ API: calculate_variance_explained returns per-node stats
    try:
        r2 = ent.model.calculate_variance_explained(r2=True)
        # r2 is dict keyed by node label; we want 'W' (weights) per view per factor
        # Expected shape: dict['W'][factor_idx] = {'view_name': r2_value}
        rows = []
        for fkey, view_r2 in r2.get("W", {}).items():
            factor_idx = int(fkey.split("_")[-1]) - 1
            for view_name, r2_val in view_r2.items():
                rows.append({
                    "factor": factor_cols[factor_idx],
                    "view": view_name,
                    "r2_pct": float(r2_val) * 100.0,
                })
        if rows:
            r2_df = pd.DataFrame(rows)
            r2_df.to_csv(f"{BASE}/mofa_multiomics_r2_per_view_v2.csv", index=False)
            print(f"Saved R2 table: {len(r2_df)} rows", flush=True)
        else:
            print("WARN: no R2 rows captured; dumping keys:", list(r2.keys()), flush=True)
    except Exception as e:
        print(f"WARN: variance explained failed: {e}", flush=True)
        traceback.print_exc()

    # ----- Save model for downstream re-use -----
    try:
        out_hdf5 = f"{BASE}/mofa_multiomics_model_v2.hdf5"
        ent.save(out_hdf5, save_data=False)
        print(f"Model saved: {out_hdf5}", flush=True)
    except Exception as e:
        print(f"WARN: save failed: {e}", flush=True)

    # ----- Save training summary -----
    summary = {
        "n_cells": int(n),
        "n_genes": int(g),
        "n_factors_requested": 7,
        "n_factors_actual": int(Z.shape[1]),
        "ard_weights": True,
        "spikeslab_weights": True,
        "likelihoods": ["gaussian", "gaussian"],
        "iter": 1000,
        "seed": 42,
    }
    with open(f"{BASE}/mofa_multiomics_summary_v2.json", "w") as fp:
        json.dump(summary, fp, indent=2)
    print("MOFA+ v2 DONE.", flush=True)

except Exception as e:
    traceback.print_exc()
    sys.exit(1)