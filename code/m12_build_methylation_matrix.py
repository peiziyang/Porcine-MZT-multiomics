#!/usr/bin/env python
"""
Build methylation matrix from CGmap gene_methyl intermediate files
and run MOFA+ multi-view (RNA + Methylation).
Run this after cgmap_wsl_v2.sh completes.
"""
import os, sys, glob, gzip
import pandas as pd
import numpy as np

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
OUT_DIR = os.path.join(BASE, "data/processed/m4_mofa")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Read all gene_methyl.txt files (rebuilt from cpg.bed + gene BED if needed)
# For now, we'll process the cpg.bed files and run bedtools ourselves
# or read pre-computed results

# 2. Build methylation matrix by reading gene_methyl.txt files
gene_methyl_txts = sorted(glob.glob(os.path.join(OUT_DIR, "_tmp_*_gene_methyl.txt")))
print(f"Found {len(gene_methyl_txts)} gene_methyl files")

all_data = {}
for fpath in gene_methyl_txts:
    fname = os.path.basename(fpath)
    cell_id = "_".join(fname.split("_")[1:]).replace("_gene_methyl.txt", "")
    cell_data = {}
    with open(fpath) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5 and parts[4] not in (".", "0", "") and float(parts[4]) > 0:
                cell_data[parts[3]] = float(parts[4])
    all_data[cell_id] = cell_data

if not all_data:
    print("No gene_methyl files found. Run wsl_cgmap_v2.sh first.")
    sys.exit(1)

# 3. Build methylation matrix
cells = list(all_data.keys())
all_genes = sorted(set(g for c in all_data.values() for g in c))
print(f"  {len(cells)} cells, {len(all_genes)} genes")

meth_df = pd.DataFrame(np.nan, index=cells, columns=all_genes)
for cell_id, gene_vals in all_data.items():
    for gid, methyl in gene_vals.items():
        meth_df.loc[cell_id, gid] = methyl

meth_df = meth_df.fillna(meth_df.median())
print(f"  Methylation matrix: {meth_df.shape}")
meth_df.to_csv(os.path.join(OUT_DIR, "methylation_matrix.csv"))

# 4. MOFA+ multi-view
print(f"\nMOFA+ multi-view...")
rna = pd.read_csv(os.path.join(BASE, "data/processed/m7_hvg_symbol.csv"), index_col=0).T
common = sorted(set(meth_df.index) & set(rna.index))
print(f"  Common cells: {len(common)}")

if len(common) >= 5:
    cg = sorted(set(meth_df.columns) & set(rna.columns))
    top = (rna[cg].var().sort_values(ascending=False).head(3000).index
           if len(cg) > 3000 else cg)
    
    ma = pd.DataFrame(0.5, index=rna.index, columns=top)
    for g in meth_df.columns:
        if g in ma.columns:
            ma[g] = meth_df[g]
    
    pd.concat([rna.loc[common, top], ma.loc[common, top]], axis=1).to_csv(
        os.path.join(OUT_DIR, "multiomics_views.csv"))
    
    try:
        from mofapy2.run.entry_point import entry_point
        ent = entry_point()
        ent.set_data_options(scale_views=True)
        vd = pd.concat([rna.loc[common, top].add_prefix("RNA_"),
                        ma.loc[common, top].add_prefix("METH_")], axis=1)
        ent.set_data_df(vd)
        ent.set_model_options(factors=5, likelihoods="gaussian")
        ent.set_train_options(iter=500, convergence_mode="fast", verbose=True, seed=42)
        ent.build()
        ent.run()
        pd.DataFrame(ent.get_factors(), index=common).to_csv(
            os.path.join(OUT_DIR, "mofa_multiomics_factors.csv"))
        print(f"  MOFA+ done!")
    except Exception as e:
        print(f"  MOFA+ failed: {type(e).__name__}")
        # Save the views anyway
        print(f"  Views saved for manual analysis")

print("\nDONE!")
