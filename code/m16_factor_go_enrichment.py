#!/usr/bin/env python
"""Factor-level GO enrichment for MOFA+ multi-omics v2.

For each MOFA+ factor:
  - Top 30 positive + 30 negative genes per view (RNA + METH)
  - Map symbol -> Entrez via Sus_scrofa.gene_info
  - Run Fisher's exact test against background (the 2552 MOFA+ genes)
  - BH correction across all GO terms
  - Output: factor -> top enriched GO terms

Inputs:
  data/processed/m4_mofa/factor_top_genes_v2.csv     (from m14_v2)
  E:\\迅雷下载\\Sus_scrofa.gene_info.gz                (Entrez/Symbol/Ensembl mapping)
  E:\\迅雷下载\\gene2go.gz                              (Entrez -> GO terms)

Outputs:
  data/processed/m4_mofa/factor_go_enrichment.csv
  data/processed/m4_mofa/factor_go_dotplot.png
"""
import os, gzip, sys
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE   = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa"
ANNO   = "E:/迅雷下载"

# ============== 1. Load annotations ==============
print("Loading Sus scrofa gene_info ...", flush=True)
info = pd.read_csv(f"{ANNO}/Sus_scrofa.gene_info.gz", sep="\t",
                   dtype=str, na_values=["-", ""])
info = info[["GeneID", "Symbol", "dbXrefs"]].dropna(subset=["GeneID"])
# Extract Ensembl from dbXrefs (Ensembl:ENSSSCG...)
def get_ensembl(x):
    if pd.isna(x): return np.nan
    for tok in str(x).split("|"):
        if tok.startswith("Ensembl:"):
            return tok.replace("Ensembl:", "")
    return np.nan
info["Ensembl"] = info["dbXrefs"].map(get_ensembl)
info = info.dropna(subset=["Symbol"])
sym2entrez = dict(zip(info["Symbol"].str.upper(), info["GeneID"]))
print(f"  {len(sym2entrez)} symbols with Entrez ID", flush=True)

print("Loading gene2go ...", flush=True)
g2go = pd.read_csv(f"{ANNO}/gene2go.gz", sep="\t",
                   dtype=str, na_values=["-", ""])
g2go = g2go[g2go["#tax_id"] == "9823"]   # Sus scrofa tax_id
g2go = g2go[["GeneID", "GO_ID", "GO_term", "GO_category"]].drop_duplicates()
print(f"  {len(g2go)} pig GO associations", flush=True)

entrez2go = g2go.groupby("GeneID")["GO_ID"].apply(set).to_dict()
all_go_terms = set(g2go["GO_ID"])
print(f"  {len(all_go_terms)} distinct pig GO terms", flush=True)

# GO category lookup
go_cat_map = g2go.set_index("GO_ID")["GO_category"].to_dict()

# ============== 2. Load factor top genes + build background ==============
top_df = pd.read_csv(f"{BASE}/factor_top_genes_v2.csv")
# Background = all 2552 genes that went into MOFA+
rna_view = pd.read_csv(f"{BASE}/rna_view_for_mofa.csv", index_col=0)
bg_genes = set(rna_view.columns.str.upper())
bg_entrez = set(sym2entrez[g] for g in bg_genes if g in sym2entrez)
print(f"Background: {len(bg_genes)} symbols, {len(bg_entrez)} Entrez", flush=True)

# ============== 3. For each factor view, run GO enrichment ==============
def run_enrichment(symbol_list, bg_entrez, entrez2go):
    """Fisher exact per GO term."""
    # Map symbols to Entrez
    sig_entrez = set()
    for s in symbol_list:
        s_u = str(s).upper()
        if s_u in sym2entrez:
            sig_entrez.add(sym2entrez[s_u])
    # Filter to background
    sig_in_bg = sig_entrez & bg_entrez
    if len(sig_in_bg) < 3:
        return pd.DataFrame(), sig_in_bg
    # Per GO term: contingency
    rows = []
    for go_term, all_genes_with_go in [(g, entrez2go.get(g, set())) for g in all_go_terms]:
        genes_with_go_in_bg = all_genes_with_go & bg_entrez
        a = len(sig_in_bg & genes_with_go_in_bg)
        b = len(genes_with_go_in_bg) - a
        c = len(sig_in_bg) - a
        d = len(bg_entrez) - a - b - c
        if a == 0:
            continue
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "GO_ID": go_term,
            "GO_category": go_cat_map.get(go_term, ""),
            "k_sig": a,
            "k_bg_with": len(genes_with_go_in_bg),
            "n_sig": len(sig_in_bg),
            "n_bg": len(bg_entrez),
            "p_raw": p,
        })
    if not rows:
        return pd.DataFrame(), sig_in_bg
    df = pd.DataFrame(rows)
    # BH
    _, q, _, _ = multipletests(df["p_raw"], method="fdr_bh")
    df["q_bh"] = q
    df = df.sort_values("p_raw")
    return df, sig_in_bg

# ============== 4. Run for every factor × view ==============
all_enrich = []
summary_rows = []
print("\nRunning GO enrichment ...", flush=True)
for _, row in top_df.iterrows():
    factor = row["factor"]; view = row["view"]
    sym_list = [s for s in row["genes"].split(";") if s]
    enr, sig_in_bg = run_enrichment(sym_list, bg_entrez, entrez2go)
    if enr.empty:
        print(f"  {factor} {view}: no enrichment", flush=True)
        continue
    enr["factor"] = factor
    enr["view"]   = view
    # Annotate term names
    name_map = g2go.drop_duplicates("GO_ID").set_index("GO_ID")["GO_term"].to_dict()
    enr["GO_term_name"] = enr["GO_ID"].map(name_map).fillna("")
    all_enrich.append(enr)
    n_sig = (enr["q_bh"] < 0.05).sum()
    summary_rows.append({
        "factor": factor, "view": view, "n_sig_genes": len(sig_in_bg),
        "n_sig_q05": n_sig, "best_p": enr["p_raw"].min(),
        "best_term": enr.iloc[0]["GO_ID"] + " " + enr.iloc[0]["GO_term_name"],
    })
    print(f"  {factor} {view}: {n_sig} BH<0.05 GO terms ({len(sig_in_bg)} sig genes mapped)", flush=True)

if not all_enrich:
    print("ERROR: no enrichment produced.", flush=True)
    sys.exit(1)

enrich_df = pd.concat(all_enrich, ignore_index=True)
enrich_df.to_csv(f"{BASE}/factor_go_enrichment_full.csv", index=False)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f"{BASE}/factor_go_summary.csv", index=False)
print(f"\nSaved: {BASE}/factor_go_enrichment_full.csv ({len(enrich_df)} rows)", flush=True)
print(f"Saved: {BASE}/factor_go_summary.csv", flush=True)

# ============== 5. Visualization: top GO terms per factor-view (raw-p < 0.01) ==============
# Pick top 5 terms per factor-view combo
display_df = enrich_df[enrich_df["p_raw"] < 0.01].copy()
display_df["neglog10p"] = -np.log10(display_df["p_raw"].clip(lower=1e-12))
display_df["label"] = display_df["GO_term_name"].str[:40] + " (" + display_df["GO_ID"] + ")"
top_per_group = (display_df.sort_values("p_raw")
                  .groupby(["factor", "view"])
                  .head(5)
                  .reset_index(drop=True))
print(f"Top terms for plot: {len(top_per_group)}", flush=True)

if len(top_per_group) > 0:
    fig, ax = plt.subplots(figsize=(13, max(6, 0.32 * len(top_per_group))))
    # y axis: factor-view-term
    ypos = np.arange(len(top_per_group))
    colors = {"RNA_pos":"#1f77b4","RNA_neg":"#aec7e8","METH_pos":"#ff7f0e","METH_neg":"#ffbb78"}
    bar_colors = top_per_group["view"].map(colors).fillna("gray")
    ax.barh(ypos, top_per_group["neglog10p"], color=bar_colors, edgecolor="black", height=0.7)
    ax.set_yticks(ypos)
    ylabels = top_per_group["factor"] + " | " + top_per_group["view"] + " | " + top_per_group["label"]
    ax.set_yticklabels(ylabels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("-log10(raw p-value)")
    ax.set_title("Factor-level GO enrichment: top 5 terms per factor-view (raw p < 0.01)\n"
                 "Background: 2552 MOFA+ genes; Fisher exact; BH q in CSV")
    # Legend
    from matplotlib.patches import Patch
    legend_handles = [Patch(color=c, label=v) for v, c in colors.items()]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=8)
    ax.axvline(-np.log10(0.05), color="red", ls="--", lw=0.8, alpha=0.5, label="p=0.05")
    plt.tight_layout()
    plt.savefig(f"{BASE}/factor_go_dotplot.png", dpi=200, bbox_inches="tight")
    print(f"Saved: {BASE}/factor_go_dotplot.png", flush=True)
else:
    print("No terms passed raw p<0.01; skipping dotplot.", flush=True)

# ============== 6. Print top finding per factor ==============
print("\n=== TOP ENRICHMENT PER FACTOR ===")
for (fac, view), grp in enrich_df.groupby(["factor", "view"]):
    best = grp.iloc[0]
    print(f"  {fac} {view}: {best['GO_term_name']} ({best['GO_ID']}, cat={best['GO_category']}) "
          f"k={int(best['k_sig'])}/{int(best['n_sig'])}  p={best['p_raw']:.2e}  q={best['q_bh']:.2e}")

print("\nDONE.")