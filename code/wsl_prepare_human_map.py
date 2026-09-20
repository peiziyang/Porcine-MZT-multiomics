#!/usr/bin/env python3
"""
Prepare human-symbol expression matrix for pySCENIC.
Maps pig Ensembl IDs -> human gene symbols (via homologene) so that
the cisTarget ranking DB (human) can be used.
Input : /mnt/e/.../m7_expr_log1p_genesxcells.csv (pig Ensembl genes x cells)
Output: /mnt/e/.../m7_expr_human_symbol.csv (human symbols x cells)
"""
import pandas as pd
import numpy as np
import mygene

BASE = '/mnt/e/Workbuddy/2026-07-27-11-58-27'
expr = pd.read_csv(f'{BASE}/data/processed/m7_expr_log1p_genesxcells.csv', index_col=0)
print(f"Pig matrix: {expr.shape[0]} genes x {expr.shape[1]} cells")

gene_ids = expr.index.tolist()
mg = mygene.MyGeneInfo()

# Step 1: pig Ensembl -> homologene ID
print("Querying homologene for pig genes...")
hres = mg.querymany(gene_ids, scopes='ensembl.gene', species='pig',
                    fields='homologene', as_dataframe=True, df_index=True)
hres = hres[~hres.index.duplicated(keep='first')]
hom = hres['homologene'].dropna()

# Extract homologene IDs that are dicts/lists
hom_ids = []
pig_to_hom = {}
for pig_id, val in hom.items():
    if isinstance(val, dict) and 'id' in val:
        hid = val['id']
    elif isinstance(val, list) and len(val) > 0 and isinstance(val[0], dict) and 'id' in val[0]:
        hid = val[0]['id']
    else:
        continue
    pig_to_hom[pig_id] = hid
    hom_ids.append(hid)
print(f"  {len(pig_to_hom)} pig genes have homologene IDs")

# Step 2: homologene ID -> human symbol
print("Querying human symbols for homologene IDs...")
ures = mg.querymany(hom_ids, scopes='homologene', species='human',
                    fields='symbol', as_dataframe=True, df_index=True)
ures = ures[~ures.index.duplicated(keep='first')]
hom_to_symbol = ures['symbol'].dropna().to_dict()

# Build pig -> human symbol map
pig_to_symbol = {}
for pig_id, hid in pig_to_hom.items():
    sym = hom_to_symbol.get(hid)
    if isinstance(sym, str) and sym:
        pig_to_symbol[pig_id] = sym

print(f"  {len(pig_to_symbol)} pig genes mapped to human symbols")

# Map expression matrix (keep genes x cells: rows=genes, cols=cells for pySCENIC)
mapped = expr.loc[[g for g in expr.index if g in pig_to_symbol]]
mapped.index = [pig_to_symbol[g] for g in mapped.index]
mapped.index.name = 'human_symbol'
# Aggregate duplicate symbols (mean across matched pig genes)
mapped = mapped.groupby(level=0).mean()
mapped.to_csv(f'{BASE}/data/processed/m7_expr_human_symbol.csv')
print(f"Human-symbol matrix: {mapped.shape[0]} genes x {mapped.shape[1]} cells (genes as rows)")
print(f"Saved: m7_expr_human_symbol.csv")
