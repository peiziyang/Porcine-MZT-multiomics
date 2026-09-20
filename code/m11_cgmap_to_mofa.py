#!/usr/bin/env python
"""
CGmap → gene-level methylation matrix → MOFA+ multi-view.
Reads individual .CGmap.gz files, computes per-gene average CpG methylation,
builds gene × cell methylation matrix matching RNA data, runs MOFA+.
"""
import os, sys, gzip, glob, time, bisect
import numpy as np
import pandas as pd

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
CGMAP_DIR = os.path.join(BASE, "data/raw/cgmap")
GTF_FILE = os.path.join(BASE, "data/raw/pig_genes.gtf.gz")
OUT_DIR = os.path.join(BASE, "data/processed/m4_mofa")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Parse GTF
print("Parsing GTF for gene annotations...")
gene_annot = {}  # gene_id → (chrom, start, end)
with gzip.open(GTF_FILE, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "gene": continue
        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        attr = parts[8]
        gid = None
        for a in attr.split(";"):
            a = a.strip()
            if a.startswith("gene_id "):
                gid = a.split('"')[1] if '"' in a else a.split()[1]
                break
        if gid:
            gene_annot[gid] = (chrom, start, end)

print(f"  {len(gene_annot)} genes loaded")

# Build chromosome-specific gene lists for interval search
# {chrom: [(start, end, gene_id), ...]}
from collections import defaultdict
chrom_genes = defaultdict(list)
for gid, (chrom, start, end) in gene_annot.items():
    chrom_genes[chrom].append((start, end, gid))

# Sort each chrom by start position
for chrom in chrom_genes:
    chrom_genes[chrom].sort(key=lambda x: x[0])

print(f"  {len(chrom_genes)} chromosomes")

# 2. Process CGmap files
cgmap_files = sorted(glob.glob(os.path.join(CGMAP_DIR, "*.CGmap.gz")))
print(f"\nProcessing {len(cgmap_files)} CGmap files...")

# Gene-level methylation accumulator: {gene_id: {cell_id: [methyl_vals]}}
gene_methyl = defaultdict(lambda: defaultdict(list))

for fpath in cgmap_files:
    fname = os.path.basename(fpath)
    cell_id = "_".join(fname.split("_")[1:]).replace(".CGmap.gz", "")
    
    fstart = time.time()
    n_cpg = 0
    n_mapped = 0
    
    with gzip.open(fpath, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8: continue
            context = parts[3]
            if context != "CG": continue  # only CpG sites
            
            chrom = parts[0]
            pos = int(parts[2])
            methyl = float(parts[5])
            n_cpg += 1
            
            # Find containing gene via binary search
            if chrom in chrom_genes:
                genes = chrom_genes[chrom]
                idx = bisect.bisect_right([g[0] for g in genes], pos) - 1
                if idx >= 0:
                    g_start, g_end, gid = genes[idx]
                    if g_start <= pos <= g_end:
                        gene_methyl[gid][cell_id].append(methyl)
                        n_mapped += 1
            
            if n_cpg % 500000 == 0:
                pass  # progress
    
    print(f"  {cell_id}: {n_cpg:,} CpG sites, {n_mapped:,} mapped to genes ({time.time()-fstart:.0f}s)")

# 3. Build methylation matrix
print(f"\nBuilding methylation matrix...")
methylation_matrix = {}
for gid, cell_dict in gene_methyl.items():
    for cell_id, vals in cell_dict.items():
        if cell_id not in methylation_matrix:
            methylation_matrix[cell_id] = {}
        methylation_matrix[cell_id][gid] = np.mean(vals)

# Convert to DataFrame
meth_df = pd.DataFrame(methylation_matrix).T  # cells × genes
meth_df = meth_df.fillna(meth_df.median())  # fill missing with median
print(f"  Methylation matrix: {meth_df.shape}")

# Save
meth_df.to_csv(os.path.join(OUT_DIR, "methylation_matrix.csv"))
print(f"  Saved to {OUT_DIR}/methylation_matrix.csv")

# 4. Build MOFA+ input
print(f"\nBuilding MOFA+ multi-view input...")
# Load RNA expression (the 3000 HVGs)
rna = pd.read_csv(os.path.join(BASE, "data/processed/m7_hvg_symbol.csv"), index_col=0)
rna = rna.T  # cells × genes
print(f"  RNA matrix: {rna.shape}")

# Find common cells
common_cells = sorted(set(meth_df.index) & set(rna.index))
print(f"  Common cells: {len(common_cells)}")

if len(common_cells) >= 5:
    meth_sub = meth_df.loc[common_cells]
    rna_sub = rna.loc[common_cells]
    
    # Shared genes (intersect methylation and RNA gene sets)
    common_genes = sorted(set(meth_sub.columns) & set(rna_sub.columns))
    print(f"  Common genes: {len(common_genes)}")
    
    # Select top 3000 most variable genes for MOFA
    if len(common_genes) > 3000:
        gene_vars = rna_sub[common_genes].var().sort_values(ascending=False)
        top_genes = list(gene_vars.head(3000).index)
    else:
        top_genes = common_genes
    
    # Prepare MOFA+ data
    # mofapy2 input: samples (cells) × features (genes), with view annotation
    from mofapy2.run.entry_point import entry_point
    ent = entry_point()
    
    rna_view = rna_sub[top_genes].values
    meth_view = meth_sub[[g for g in top_genes if g in meth_sub.columns]]
    # Fill missing methylation genes with 0.5 (unmethylated baseline)
    meth_view_aligned = pd.DataFrame(0.5, index=rna_sub.index, columns=top_genes)
    for g in meth_view.columns:
        if g in meth_view_aligned.columns:
            meth_view_aligned[g] = meth_view[g]
    
    print(f"  RNA view: {rna_view.shape}, Meth view: {meth_view_aligned.shape}")
    
    # Build data
    data = np.zeros((len(common_cells), len(top_genes), 2))
    data[:, :, 0] = rna_view
    data[:, :, 1] = meth_view_aligned.values
    
    ent.set_data_options(scale_views=True)
    ent.set_data_df(
        pd.DataFrame(data.reshape(len(common_cells), -1),
                     index=common_cells,
                     columns=[f"RNA_{g}" for g in top_genes] + [f"METH_{g}" for g in top_genes])
    )
    # Add view metadata
    # (MOFA+ expects specific format)
    
    ent.set_model_options(factors=5, likelihoods=["gaussian", "gaussian"])
    ent.set_train_options(iter=1000, convergence_mode="fast", verbose=True)
    
    ent.build()
    ent.run()
    
    # Save factors
    factors = ent.get_factors()
    pd.DataFrame(factors, index=common_cells).to_csv(
        os.path.join(OUT_DIR, "mofa_factors.csv"))
    
    # Save weights
    weights = ent.get_weights()
    for v, w in enumerate(weights):
        pd.DataFrame(w, index=top_genes).to_csv(
            os.path.join(OUT_DIR, f"mofa_weights_view{v}.csv"))
    
    print(f"\nMOFA+ done. {len(common_cells)} cells × {5} factors")
else:
    print(f"Too few common cells ({len(common_cells)}) for MOFA+")

print("\nDone!")
