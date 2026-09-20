#!/usr/bin/env python
"""Fast GO enrichment for MOFA+ factor top genes (optimized v2).

Difference from m16_factor_go_enrichment.py:
  - Reads gene2go.gz line-by-line (streaming), filtering to only pig entries.
  - Pre-filters gene_info to only pig genes with Ensembl mapping.
  - Much faster than reading 1.3GB with pandas.

Inputs:
  factor_top_genes_v2.csv
  rna_view_for_mofa.csv (background genes)

Outputs:
  factor_go_enrichment_v2.csv
  factor_go_summary_v2.csv
"""
import os, gzip, sys
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

BASE   = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa"
ANNO   = "E:/迅雷下载"

# ============== 1. Load background genes + symbol->entrez mapping ==============
print("Loading background + gene_info ...", flush=True)
rna_view = pd.read_csv(f"{BASE}/rna_view_for_mofa.csv", index_col=0)
bg_genes = set(rna_view.columns.str.upper())

# Stream gene_info: only pig (tax_id=9606 is human, 9823 is pig)
sym2entrez = {}
with gzip.open(f"{ANNO}/Sus_scrofa.gene_info.gz", "rt", encoding="utf-8") as f:
    header = f.readline()
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        gene_id = parts[1]     # GeneID
        symbol = parts[2]       # Symbol
        dbxrefs = parts[5] if len(parts) > 5 else ""
        sym2entrez[symbol.upper()] = gene_id
print(f"  {len(sym2entrez)} pig symbols -> Entrez", flush=True)

bg_entrez = set(sym2entrez[g] for g in bg_genes if g in sym2entrez)
print(f"  {len(bg_entrez)} of {len(bg_genes)} bg genes have Entrez", flush=True)

# ============== 2. Stream gene2go: only pig entries (tax_id=9823) ==============
print("Streaming gene2go.gz (pig only) ...", flush=True)
entrez2go = defaultdict(set)
go_names = {}
go_cats = {}
go_genes = defaultdict(set)  # GO_ID -> set of entrez genes in our bg
line_count = 0
with gzip.open(f"{ANNO}/gene2go.gz", "rt", encoding="utf-8") as f:
    for line in f:
        if line_count == 0:
            line_count += 1
            continue  # skip header
        if not line.startswith("9823\t"):  # Sus scrofa tax_id
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        entrez = parts[1]; go_id = parts[2]
        go_name = parts[5] if len(parts) > 5 else ""   # GO_term name
        go_cat = parts[7] if len(parts) > 7 else ""     # Category (BP/MF/CC)
        entrez2go[entrez].add(go_id)
        go_names[go_id] = go_name
        go_cats[go_id] = go_cat
        line_count += 1
        if line_count % 100000 == 0:
            print(f"  ... {line_count} pig GOs", flush=True)
print(f"  {line_count} pig GO rows, {len(entrez2go)} genes, {len(go_names)} terms", flush=True)

# Build reverse: GO -> genes
for eid, gos in entrez2go.items():
    for g in gos:
        go_genes[g].add(eid)
# Filter to bg
for g in list(go_genes.keys()):
    go_genes[g] = go_genes[g] & bg_entrez
    if not go_genes[g]:
        del go_genes[g]
print(f"  {len(go_genes)} GO terms with bg genes", flush=True)

# ============== 3. Load factor top genes ==============
top_df = pd.read_csv(f"{BASE}/factor_top_genes_v2.csv")

# ============== 4. Enrichment ==============
print("\nRunning enrichment ...", flush=True)
all_rows = []

for _, row in top_df.iterrows():
    factor = row["factor"]; view = row["view"]
    sym_list = [s.strip() for s in str(row["genes"]).split(";") if s]
    sig_entrez = set()
    for s in sym_list:
        if s.upper() in sym2entrez:
            sig_entrez.add(sym2entrez[s.upper()])
    sig_in_bg = sig_entrez & bg_entrez
    if len(sig_in_bg) < 3:
        print(f"  {factor} {view}: only {len(sig_in_bg)} mapped, skip", flush=True)
        continue

    for go_id, go_entrez in go_genes.items():
        a = len(sig_in_bg & go_entrez)
        if a == 0:
            continue
        b = len(go_entrez) - a
        c = len(sig_in_bg) - a
        d = len(bg_entrez) - a - b - c
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        all_rows.append({
            "factor": factor, "view": view,
            "GO_ID": go_id, "GO_term": go_names.get(go_id, ""),
            "GO_category": go_cats.get(go_id, ""),
            "k_sig": a, "k_bg": len(go_entrez),
            "n_sig": len(sig_in_bg), "n_bg": len(bg_entrez),
            "p_raw": p, "n_genes_sig_total": len(sig_in_bg),
        })
    print(f"  {factor} {view}: {len(sig_in_bg)} sig genes, tested {len(go_genes)} GO terms",
          flush=True)

if not all_rows:
    print("FATAL: no enrichment rows produced")
    sys.exit(1)

enr = pd.DataFrame(all_rows)
# BH correction
_, q, _, _ = multipletests(enr["p_raw"], method="fdr_bh")
enr["q_bh"] = q
enr = enr.sort_values("p_raw")
enr.to_csv(f"{BASE}/factor_go_enrichment_v2.csv", index=False)
print(f"\nSaved: {len(enr)} rows to factor_go_enrichment_v2.csv", flush=True)

# Summary: per factor per view
summary = []
for (fac, view), grp in enr.groupby(["factor", "view"]):
    n_q05 = (grp["q_bh"] < 0.05).sum()
    best = grp.iloc[0]
    summary.append({
        "factor": fac, "view": view,
        "n_q05": n_q05, "n_raw_p005": (grp["p_raw"] < 0.005).sum(),
        "n_raw_p001": (grp["p_raw"] < 0.001).sum(),
        "best_term": f"{best['GO_ID']} {best['GO_term']}",
        "best_p": best["p_raw"], "best_q": best["q_bh"],
    })
summary_df = pd.DataFrame(summary)
summary_df.to_csv(f"{BASE}/factor_go_summary_v2.csv", index=False)

print("\n=== SUMMARY ===")
for _, r in summary_df.iterrows():
    print(f"  {r['factor']} {r['view']}: q<0.05={r['n_q05']} "
          f"raw<0.005={r['n_raw_p005']} raw<0.001={r['n_raw_p001']} "
          f"best={r['best_term'][:60]}")
print("\nDONE.")