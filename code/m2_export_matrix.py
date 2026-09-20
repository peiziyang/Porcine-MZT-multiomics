"""
M7 prep: Export raw-counts expression matrix (Ensembl IDs) for pySCENIC.
Uses only the 3 datasets with true raw counts:
  GSE168106 (Han 2022, E0-E14, ~1458 cells)
  GSE112380 (Ramos-Ibeas 2019, D5-D11, ~220 cells)
  GSE234116 (Yuan 2023, GV oocytes, 62 cells)
Mixing FPKM datasets is avoided for GRN strictness.
"""
import pandas as pd
import numpy as np
import os

DATA_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw'
OUT_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
os.makedirs(OUT_DIR, exist_ok=True)

parts = []

# GSE168106
df = pd.read_csv(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_genes_counts.txt.gz', sep='\t', index_col=0)
si = pd.read_excel(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_SampleInfo.xlsx')
si_valid = si[si['Vaild'] == 'Yes'].copy()
valid = set(si_valid['ID'])
df = df[[c for c in df.columns if c in valid]]
print(f"GSE168106: {df.shape[1]} cells, {df.shape[0]} genes")
parts.append(('GSE168106', df))

# GSE112380
df = pd.read_csv(f'{DATA_DIR}/GSE112380/GSE112380_Pig_RawCounts.txt.gz', sep='\t', index_col=0)
print(f"GSE112380: {df.shape[1]} cells, {df.shape[0]} genes")
parts.append(('GSE112380', df))

# GSE234116
df = pd.read_csv(f'{DATA_DIR}/GSE234116/GSE234116_raw_counts.txt.gz', sep='\t', index_col=0)
print(f"GSE234116: {df.shape[1]} cells, {df.shape[0]} genes")
parts.append(('GSE234116', df))

# Keep only Ensembl genes present in all
ensembl_sets = []
for name, df in parts:
    gs = {g for g in df.index if str(g).startswith('ENSSSCG')}
    ensembl_sets.append(gs)
common = set.intersection(*ensembl_sets)
print(f"Common Ensembl genes: {len(common)}")

# Build unified matrix (genes x cells), log1p
expr_list = []
for name, df in parts:
    sub = df.loc[[g for g in df.index if g in common]]
    sub = sub.apply(pd.to_numeric, errors='coerce').fillna(0)
    expr_list.append(sub)
expr = pd.concat(expr_list, axis=1)
print(f"Unified raw matrix: {expr.shape[0]} genes x {expr.shape[1]} cells")

# Save raw counts (for reference)
expr.to_csv(f'{OUT_DIR}/m7_raw_counts.csv.gz', compression='gzip')
print("Saved m7_raw_counts.csv.gz")

# Save log1p (pySCENIC grn input convention: genes as index, cells as columns)
expr_log = np.log1p(expr)
expr_log.to_csv(f'{OUT_DIR}/m7_expr_log1p_genesxcells.csv')
print("Saved m7_expr_log1p_genesxcells.csv (genes x cells)")

# Also save metadata
meta_rows = []
for name, df in parts:
    for c in df.columns:
        meta_rows.append({'cell': c, 'dataset': name})
meta = pd.DataFrame(meta_rows).set_index('cell')
meta.to_csv(f'{OUT_DIR}/m7_cells_metadata.csv')
print("Saved m7_cells_metadata.csv")
print("DONE")
