#!/usr/bin/env python
"""MOFA+ multi-omics: RNA + CpG methylation on GV oocytes."""
import pandas as pd, numpy as np
import sys, traceback

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa"

try:
    rna = pd.read_csv(f"{BASE}/rna_view_for_mofa.csv", index_col=0)
    meth = pd.read_csv(f"{BASE}/meth_view_for_mofa.csv", index_col=0)
    cells = rna.index.tolist()
    genes = rna.columns.tolist()
    n, g = len(cells), len(genes)
    
    # Build long-format data efficiently
    data = pd.DataFrame({
        "sample": np.repeat(cells, 2*g),
        "group": ["group1"] * (n * 2 * g),
        "view": np.tile(np.repeat(["RNA", "METH"], g), n),
        "feature": np.tile(genes, 2*n),
        "value": np.concatenate([rna.values.flatten(), meth.values.flatten()])
    })
    print(f"Data: {data.shape}", flush=True)
    
    # Train MOFA+
    from mofapy2.run.entry_point import entry_point
    ent = entry_point()
    ent.set_data_options(scale_views=True)
    ent.set_data_df(data)
    ent.set_model_options(factors=5)
    ent.set_train_options(iter=500, convergence_mode="fast", verbose=True, seed=42)
    print("Building...", flush=True)
    ent.build()
    print("Training...", flush=True)
    ent.run()
    
    # Extract factors: Z has shape (n_cells, n_factors)
    Z = ent.model.nodes["Z"].getExpectations()["E"]
    print(f"Factors Z: {Z.shape}", flush=True)
    pd.DataFrame(Z, index=cells, columns=[f"F{i+1}" for i in range(Z.shape[1])]).to_csv(
        f"{BASE}/mofa_multiomics_factors.csv")
    
    # Extract weights: W is list of (n_features, n_factors) per view
    W = [w["E"] for w in ent.model.nodes["W"].getExpectations()]
    for vi, w in enumerate(W):
        view_name = ["RNA", "METH"][vi]
        pd.DataFrame(w, index=genes, columns=[f"F{i+1}" for i in range(w.shape[1])]).to_csv(
            f"{BASE}/mofa_multiomics_weights_{view_name}.csv")
        print(f"  {view_name} weights: {w.shape}", flush=True)
    
    print("MOFA+ DONE!", flush=True)
    
except Exception as e:
    traceback.print_exc()
    sys.exit(1)
