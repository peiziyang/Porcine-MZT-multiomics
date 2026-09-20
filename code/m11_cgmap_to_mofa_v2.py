#!/usr/bin/env python
"""
CGmap → gene-level methylation matrix → MOFA+ multi-view.
Optimized v2: sweep-line algorithm per chromosome, batch processing.
"""
import os, sys, gzip, glob, time
import numpy as np
import pandas as pd
from collections import defaultdict

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
CGMAP_DIR = os.path.join(BASE, "data/raw/cgmap")
GTF_FILE = os.path.join(BASE, "data/raw/pig_genes.gtf.gz")
OUT_DIR = os.path.join(BASE, "data/processed/m4_mofa")
os.makedirs(OUT_DIR, exist_ok=True)

# --- Step 1: Build gene interval index ---
print("STEP 1: Building gene interval index...", flush=True)
gene_annot = {}
with gzip.open(GTF_FILE, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "gene": continue
        chrom = parts[0]
        start, end = int(parts[3]), int(parts[4])
        attr = parts[8]
        gid = None
        for a in attr.split(";"):
            a = a.strip()
            if a.startswith("gene_id "):
                gid = a.split('"')[1] if '"' in a else a.split()[1]
                break
        if gid:
            gene_annot[gid] = (chrom, start, end)

print(f"  {len(gene_annot)} genes loaded", flush=True)

# Build per-chromosome sorted list: (gene_start, gene_end, gene_id, chrom)
# For sweep-line: we need genes sorted by start, and CpGs sorted by position
chrom_genes = defaultdict(list)
for gid, (chrom, start, end) in gene_annot.items():
    chrom_genes[chrom].append((start, end, gid))

for chrom in chrom_genes:
    chrom_genes[chrom].sort(key=lambda x: x[0])

print(f"  {len(chrom_genes)} chromosomes", flush=True)

# --- Step 2: Process CGmap files ---
cgmap_files = sorted(glob.glob(os.path.join(CGMAP_DIR, "*.CGmap.gz")))
print(f"\nSTEP 2: Processing {len(cgmap_files)} CGmap files...", flush=True)

# {gene_id: {cell_id: [methylation values]}}
gene_methyl = defaultdict(lambda: defaultdict(list))

for fpath in cgmap_files:
    fname = os.path.basename(fpath)
    cell_id = "_".join(fname.split("_")[1:]).replace(".CGmap.gz", "")
    fstart = time.time()
    
    # Read all CpG sites into per-chromosome buckets
    chrom_cpg = defaultdict(list)  # chrom -> [(pos, methyl), ...]
    with gzip.open(fpath, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6 or parts[3] != "CG":
                continue
            chrom_cpg[parts[0]].append((int(parts[2]), float(parts[5])))
    
    # For each chromosome, sweep through CpGs and assign to genes
    n_assigned = 0
    for chrom, cpgs in chrom_cpg.items():
        if chrom not in chrom_genes:
            continue
        
        cpgs.sort(key=lambda x: x[0])
        genes = chrom_genes[chrom]
        
        git = 0
        for pos, methyl in cpgs:
            # Advance gene pointer to the first gene starting at or before this CpG
            while git < len(genes) and genes[git][1] < pos:
                git += 1
            if git >= len(genes):
                break
            # Check all genes at the current position that might contain the CpG
            matched = False
            for g in range(git, len(genes)):
                g_start, g_end, gid = genes[g]
                if g_start > pos:
                    break
                if g_start <= pos <= g_end:
                    gene_methyl[gid][cell_id].append(methyl)
                    n_assigned += 1
                    matched = True
                    # Don't break - a CpG could be in overlapping genes
            # Reset git for next chromosome batch
    
    elapsed = time.time() - fstart
    total_cpg = sum(len(v) for v in chrom_cpg.values())
    print(f"  {cell_id}: {total_cpg:,} CpG sites, {n_assigned:,} assigned ({elapsed:.0f}s)", flush=True)

# --- Step 3: Build methylation matrix ---
print(f"\nSTEP 3: Building methylation matrix...", flush=True)
methylation_matrix = {}
for gid, cell_dict in gene_methyl.items():
    for cid, vals in cell_dict.items():
        if cid not in methylation_matrix:
            methylation_matrix[cid] = {}
        methylation_matrix[cid][gid] = np.mean(vals)

meth_df = pd.DataFrame(methylation_matrix).T
meth_df = meth_df.fillna(meth_df.median())
print(f"  Methylation matrix: {meth_df.shape}", flush=True)
meth_df.to_csv(os.path.join(OUT_DIR, "methylation_matrix.csv"))
print(f"  Saved", flush=True)

# --- Step 4: MOFA+ multi-view ---
print(f"\nSTEP 4: MOFA+ multi-view...", flush=True)
rna = pd.read_csv(os.path.join(BASE, "data/processed/m7_hvg_symbol.csv"), index_col=0).T
print(f"  RNA matrix: {rna.shape}", flush=True)

common_cells = sorted(set(meth_df.index) & set(rna.index))
print(f"  Common cells: {len(common_cells)}", flush=True)

if len(common_cells) >= 5:
    common_genes = sorted(set(meth_df.columns) & set(rna.columns))
    print(f"  Common genes: {len(common_genes)}", flush=True)
    
    if len(common_genes) > 3000:
        gene_vars = rna[common_genes].var().sort_values(ascending=False)
        top_genes = list(gene_vars.head(3000).index)
    else:
        top_genes = common_genes
    
    # Align methylation matrix
    meth_aligned = pd.DataFrame(0.5, index=rna.index, columns=top_genes)
    for g in meth_df.columns:
        if g in meth_aligned.columns:
            meth_aligned[g] = meth_df[g]
    
    print(f"  RNA: {rna.loc[common_cells, top_genes].shape}")
    print(f"  Meth: {meth_aligned.loc[common_cells, top_genes].shape}")
    
    # Run MOFA+
    try:
        from mofapy2.run.entry_point import entry_point
        ent = entry_point()
        ent.set_data_options(scale_views=True, scale_groups=False)
        
        # Build data for MOFA+
        rna_data = rna.loc[common_cells, top_genes].values
        meth_data = meth_aligned.loc[common_cells, top_genes].values
        
        # Create a combined data frame for MOFA
        import itertools
        sample_names = common_cells
        feature_names = [f"RNA_{g}" for g in top_genes] + [f"METH_{g}" for g in top_genes]
        all_data = np.hstack([rna_data, meth_data])
        data_df = pd.DataFrame(all_data, index=sample_names, columns=feature_names)
        
        # Add group and view metadata
        view_factor = pd.DataFrame({
            "view": ["RNA"] * len(top_genes) + ["METH"] * len(top_genes),
            "feature": feature_names,
            "group": "group1"
        })
        
        ent.set_data_df(data_df, groups_indices=[0], views_indices=[0, 1])
        ent.set_model_options(factors=5, likelihoods=["gaussian", "gaussian"])
        ent.set_train_options(iter=500, convergence_mode="fast", verbose=True, seed=42)
        ent.build()
        ent.run()
        
        # Save results
        factors = ent.get_factors()
        pd.DataFrame(factors, index=common_cells).to_csv(
            os.path.join(OUT_DIR, "mofa_multiomics_factors.csv"))
        print(f"  MOFA+ done: {len(common_cells)} cells × 5 factors", flush=True)
    except Exception as e:
        print(f"  MOFA+ failed: {type(e).__name__}: {e}", flush=True)
        # Save methylation matrix anyway for later use
        meth_aligned.loc[common_cells, top_genes].to_csv(
            os.path.join(OUT_DIR, "methylation_matrix_aligned.csv"))
        print(f"  Methylation matrix saved for manual MOFA+", flush=True)

print("\nDONE!", flush=True)
