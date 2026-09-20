#!/usr/bin/env python
"""
CGmap → gene-level methylation matrix → MOFA+ multi-view.
Run this after user downloads individual CGmap.gz files.
"""
import os, sys, gzip
import pandas as pd
import numpy as np

CGMAP_DIR = "E:/Workbuddy/2026-07-27-11-58-27/data/raw/cgmap_indiv"
BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
GENE_ANNOT = "E:/Workbuddy/2026-07-27-11-58-27/data/raw/pig_genes.bed"  # need to create

# 1. List all downloaded CGmap.gz files
cgmap_files = [f for f in os.listdir(CGMAP_DIR) if f.endswith(".CGmap.gz")]
print(f"Found {len(cgmap_files)} CGmap files")

# 2. For each file, compute gene-level methylation
# CGmap format (from BSseeker2):
#   chr, strand, pos, context(CG/CHG/CHH), dinuc, methyl_level, C_count, CT_count
# Gene-level methylation = mean(methyl_level for CpG sites in gene body +/- 2kb)

# 3. Build gene x cell methylation matrix
# Match cells using the AF prefix naming (AF1_12, AF2_51, etc.)
cell_groups = {}
for f in cgmap_files:
    # e.g., GSM7508586_AF1_12.CGmap.gz -> AF1_12
    cell_id = "_".join(f.split("_")[1:]).replace(".CGmap.gz", "")
    cell_groups[cell_id] = f

print(f"Cells found: {list(cell_groups.keys())[:5]}...")

# 4. Load expression matrix to get the same gene set
expr = pd.read_csv(BASE + "/m7_expr_log1p_genesxcells.csv", nrows=5)
genes = list(expr.iloc[:, 0])[:100]

# 5. MOFA+ setup (placeholder - runs after methylation matrix is built)
print("\nAfter CGmap processing complete, run:")
print("  python <mofa_multiomics.py>  # MOFA+ RNA + methylation views")
