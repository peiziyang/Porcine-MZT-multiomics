#!/usr/bin/env python
"""
CGmap → gene-level methylation matrix (memory-efficient v3).
Streaming per-chromosome: never hold all CpGs in memory at once.
"""
import os, sys, gzip, glob, time, gc
import numpy as np
import pandas as pd
from collections import defaultdict

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
CGMAP_DIR = os.path.join(BASE, "data/raw/cgmap")
GTF_FILE = os.path.join(BASE, "data/raw/pig_genes.gtf.gz")
OUT_DIR = os.path.join(BASE, "data/processed/m4_mofa")
os.makedirs(OUT_DIR, exist_ok=True)

# --- Step 1: Gene index ---
print("STEP 1: Gene interval index...", flush=True)
chrom_genes = defaultdict(list)
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
            chrom_genes[chrom].append((start, end, gid))

for chrom in chrom_genes:
    chrom_genes[chrom].sort(key=lambda x: x[0])
total_genes = sum(len(v) for v in chrom_genes.values())
print(f"  {total_genes} genes, {len(chrom_genes)} chromosomes", flush=True)

# --- Step 2: Streaming CGmap processing ---
cgmap_files = sorted(glob.glob(os.path.join(CGMAP_DIR, "*.CGmap.gz")))
print(f"\nSTEP 2: {len(cgmap_files)} files (streaming)...", flush=True)

# gene → {cell → [methyl_values]}
gene_methyl = defaultdict(lambda: defaultdict(list))

for fpath in cgmap_files:
    fname = os.path.basename(fpath)
    cell_id = "_".join(fname.split("_")[1:]).replace(".CGmap.gz", "")
    fstart = time.time()
    
    # Read entire file, streaming: per line, filter CG, find gene, accumulate
    # Use chromosome partition to avoid excess per-line lookups
    current_chrom = None
    chrom_cpgs = []  # [(pos, methyl), ...] reset per chromosome
    
    with gzip.open(fpath, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6 or parts[3] != "CG":
                continue
            chrom = parts[0]
            pos = int(parts[2])
            methyl = float(parts[5])
            
            # Chromosome changed → process previous batch
            if chrom != current_chrom and chrom_cpgs:
                if current_chrom in chrom_genes:
                    genes = chrom_genes[current_chrom]
                    # Sweep: sort CpGs by position
                    chrom_cpgs.sort(key=lambda x: x[0])
                    gi = 0
                    for p, m in chrom_cpgs:
                        while gi < len(genes) and genes[gi][1] < p:
                            gi += 1
                        if gi >= len(genes):
                            break
                        for g in range(gi, len(genes)):
                            gs, ge, gid = genes[g]
                            if gs > p:
                                break
                            if gs <= p <= ge:
                                gene_methyl[gid][cell_id].append(m)
                chrom_cpgs = []
                gc.collect()
            
            current_chrom = chrom
            chrom_cpgs.append((pos, methyl))
    
    # Process last chromosome
    if chrom_cpgs and current_chrom in chrom_genes:
        genes = chrom_genes[current_chrom]
        chrom_cpgs.sort(key=lambda x: x[0])
        gi = 0
        for p, m in chrom_cpgs:
            while gi < len(genes) and genes[gi][1] < p:
                gi += 1
            if gi >= len(genes):
                break
            for g in range(gi, len(genes)):
                gs, ge, gid = genes[g]
                if gs > p:
                    break
                if gs <= p <= ge:
                    gene_methyl[gid][cell_id].append(m)
    
    n_assigned = sum(1 for g, c in gene_methyl.items() if c.get(cell_id))
    elapsed = time.time() - fstart
    print(f"  {cell_id}: done ({elapsed:.0f}s)", flush=True)
    chrom_cpgs = None
    gc.collect()

# --- Step 3: Methylation matrix ---
print(f"\nSTEP 3: Building methylation matrix...", flush=True)
meth_data = {}
for gid, cells in gene_methyl.items():
    for cid, vals in cells.items():
        if cid not in meth_data:
            meth_data[cid] = {}
        meth_data[cid][gid] = np.mean(vals)

meth_df = pd.DataFrame(meth_data).T
meth_df = meth_df.fillna(meth_df.median())
print(f"  {meth_df.shape}", flush=True)
meth_df.to_csv(os.path.join(OUT_DIR, "methylation_matrix.csv"))
print(f"  Saved methylation_matrix.csv", flush=True)

# --- Step 4: MOFA+ ---
print(f"\nSTEP 4: MOFA+ multi-view...", flush=True)
rna = pd.read_csv(os.path.join(BASE, "data/processed/m7_hvg_symbol.csv"), index_col=0).T
common = sorted(set(meth_df.index) & set(rna.index))
print(f"  Common cells: {len(common)}", flush=True)

if len(common) >= 5:
    common_genes = sorted(set(meth_df.columns) & set(rna.columns))
    top = (rna[common_genes].var().sort_values(ascending=False).head(3000).index
           if len(common_genes) > 3000 else common_genes)
    
    meth_aligned = pd.DataFrame(0.5, index=rna.index, columns=top)
    for g in meth_df.columns:
        if g in meth_aligned.columns:
            meth_aligned[g] = meth_df[g]
    
    # Save aligned matrices for reference
    pd.concat([rna.loc[common, top], meth_aligned.loc[common, top]], axis=1).to_csv(
        os.path.join(OUT_DIR, "multiomics_views.csv"))
    print(f"  Views saved ({len(common)} cells × {len(top)} genes × 2 views)", flush=True)
    
    try:
        from mofapy2.run.entry_point import entry_point
        ent = entry_point()
        ent.set_data_options(scale_views=True)
        # Build MOFA data as separate views
        rna_data = rna.loc[common, top]
        meth_data_sub = meth_aligned.loc[common, top]
        # Stack views
        view_df = pd.concat([
            rna_data.add_prefix("RNA_"),
            meth_data_sub.add_prefix("METH_")
        ], axis=1)
        ent.set_data_df(view_df)
        ent.set_model_options(factors=5, likelihoods="gaussian")
        ent.set_train_options(iter=500, convergence_mode="fast", verbose=True, seed=42)
        ent.build()
        ent.run()
        
        pd.DataFrame(ent.get_factors(), index=common).to_csv(
            os.path.join(OUT_DIR, "mofa_multiomics_factors.csv"))
        print(f"  MOFA+ done!", flush=True)
    except Exception as e:
        print(f"  MOFA+ skipped: {type(e).__name__}", flush=True)

print("\nDONE!", flush=True)
