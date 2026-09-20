#!/usr/bin/env python3
"""
Map pig Ensembl IDs -> gene symbols (single-step, no homologene).
For the co-expression GRN (GRNBoost2 + AUCell) we work in pig gene space,
so human ortholog mapping is NOT required here.
Also selects top highly-variable genes to speed up GRNBoost2.
Output: m7_expr_symbol.csv (symbols x cells), m7_hvg_symbol.csv (HVG subset)
"""
import pandas as pd
import numpy as np
import mygene

BASE = '/mnt/e/Workbuddy/2026-07-27-11-58-27'
PROC = f'{BASE}/data/processed'
expr = pd.read_csv(f'{PROC}/m7_expr_log1p_genesxcells.csv', index_col=0)
print(f"Pig matrix: {expr.shape[0]} genes x {expr.shape[1]} cells")

gene_ids = expr.index.tolist()
mg = mygene.MyGeneInfo()
print("Querying symbols...")
res = mg.querymany(gene_ids, scopes='ensembl.gene', species='pig',
                   fields='symbol,name', as_dataframe=True)
res = res[~res.index.duplicated(keep='first')]
sym_map = res['symbol'].dropna().to_dict()
print(f"  {len(sym_map)} genes mapped to symbols")

# Map and aggregate duplicate symbols
mapped = expr.loc[[g for g in expr.index if g in sym_map]]
mapped.index = [sym_map[g] for g in mapped.index]
mapped.index.name = 'symbol'
mapped = mapped.groupby(level=0).mean()
mapped.to_csv(f'{PROC}/m7_expr_symbol.csv')
print(f"Symbol matrix: {mapped.shape[0]} genes x {mapped.shape[1]} cells")
print("Saved m7_expr_symbol.csv")

# Select top HVGs for faster GRNBoost2
gene_var = mapped.var(axis=1).sort_values(ascending=False)
topn = min(3000, len(gene_var))
hvg = mapped.loc[gene_var.head(topn).index]
hvg.to_csv(f'{PROC}/m7_hvg_symbol.csv')
print(f"HVG subset: {hvg.shape[0]} genes -> m7_hvg_symbol.csv")
print("DONE")
