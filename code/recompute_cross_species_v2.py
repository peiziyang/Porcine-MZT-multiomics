#!/usr/bin/env python
"""
Recompute cross-species comparison with the CORRECT F1 gene set.

Fix: previous cross-species analysis used F1_combined (union of F1 RNA_pos U
F1 METH_pos = 44 genes from factor_top_genes_v2.csv), which DIFFERS from the
"F1-top44" gene set (top 44 by |F1 RNA weight|, Table S4) used in the PA
gene-set reanalysis. Also, the "all detectable symbol-matched gene pairs"
background was actually just those 44 genes, not genome-wide.

This script:
  1. Defines F1 genes = top 44 by |F1 RNA weight| (matches Table S4).
  2. Maps ALL pig genes (2,552) to human (exact) / mouse (case-insensitive).
  3. Genome-wide background = all pig genes with a symbol match.
  4. F1 overlap (X/44 in both species), F1 cross-species Spearman correlation.
  5. Genome-wide cross-species correlation.
  6. Expression-decile-matched permutation (10,000 iters, without replacement).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
MF = os.path.join(BASE, "data/processed/m4_mofa")
ANNO = "E:/迅雷下载"
OUT = MF

rng = np.random.default_rng(42)
N_PERM = 10000

# ── 1. F1-top44 genes (top 44 by |F1 RNA weight|) ──
w_rna = pd.read_csv(os.path.join(MF, "mofa_multiomics_weights_RNA_v2.csv"), index_col=0)
f1_top44 = list(w_rna["F1"].abs().sort_values(ascending=False).head(44).index)
print(f"F1-top44 genes: {len(f1_top44)}")

# ── 2. Pig GV oocyte mean expression (across 32 cells) ──
pig = pd.read_csv(os.path.join(MF, "rna_view_for_mofa.csv"), index_col=0)
pig_mean = pig.mean(axis=0)  # Series indexed by gene symbol

# ── 3. Human / mouse expression matrices ──
hu = pd.read_csv(os.path.join(ANNO, "GSE44183_human_expression_mat.txt.gz"),
                 sep="\t", index_col=0)
mo = pd.read_csv(os.path.join(ANNO, "GSE44183_mouse_expression_mat.txt.gz"),
                 sep="\t", index_col=0)
hu_mean = hu.mean(axis=1)  # across 29 human samples
mo_mean = mo.mean(axis=1)  # across 18 mouse samples
hu_genes = set(hu.index)
mo_upper = {g.upper(): g for g in mo.index}  # case-insensitive map

# ── 4. Symbol matching ──
def human_match(g):
    return g in hu_genes

def mouse_match(g):
    return g.upper() in mo_upper

def mouse_symbol(g):
    return mo_upper.get(g.upper())

# F1-top44 mapping
f1_in_hu = [g for g in f1_top44 if human_match(g)]
f1_in_mo = [g for g in f1_top44 if mouse_match(g)]
f1_both = [g for g in f1_top44 if human_match(g) and mouse_match(g)]
print(f"\nF1-top44: human match n={len(f1_in_hu)}, mouse match n={len(f1_in_mo)}, both n={len(f1_both)}")

# Genome-wide background (all pig genes with symbol match)
all_pig_genes = list(pig.columns)
gw_hu = [g for g in all_pig_genes if human_match(g)]
gw_mo = [g for g in all_pig_genes if mouse_match(g)]
print(f"Genome-wide: human match n={len(gw_hu)}, mouse match n={len(gw_mo)}")

# ── 5. Cross-species Spearman correlations ──
def spearman_hu(genes):
    pv = [pig_mean[g] for g in genes]
    hv = [hu_mean[g] for g in genes]
    return spearmanr(pv, hv)

def spearman_mo(genes):
    pv = [pig_mean[g] for g in genes]
    mv = [mo_mean[mouse_symbol(g)] for g in genes]
    return spearmanr(pv, mv)

r_hu_f1, p_hu_f1 = spearman_hu(f1_in_hu)
r_mo_f1, p_mo_f1 = spearman_mo(f1_in_mo)
r_hu_gw, p_hu_gw = spearman_hu(gw_hu)
r_mo_gw, p_mo_gw = spearman_mo(gw_mo)

print(f"\nF1 pig-human: rho={r_hu_f1:.4f}, p={p_hu_f1:.4e}, n={len(f1_in_hu)}")
print(f"F1 pig-mouse: rho={r_mo_f1:.4f}, p={p_mo_f1:.4e}, n={len(f1_in_mo)}")
print(f"Genome-wide pig-human: rho={r_hu_gw:.4f}, p={p_hu_gw:.4e}, n={len(gw_hu)}")
print(f"Genome-wide pig-mouse: rho={r_mo_gw:.4f}, p={p_mo_gw:.4e}, n={len(gw_mo)}")

# ── 6. Expression-decile-matched permutation (without replacement) ──
def _bucketize(f1_genes, gw_genes):
    """Return (f1_buckets, bucket_index_map) where bucket_index_map[b] is a
    numpy array of positional indices into gw_genes for decile b."""
    gw_expr = pig_mean[gw_genes].values
    f1_expr = pig_mean[f1_genes].values
    edges = np.percentile(gw_expr, np.arange(0, 101, 10))
    f1_b = np.clip(np.digitize(f1_expr, edges) - 1, 0, 9)
    gw_b = np.clip(np.digitize(gw_expr, edges) - 1, 0, 9)
    from collections import defaultdict
    bmap = defaultdict(list)
    for i, b in enumerate(gw_b):
        bmap[int(b)].append(i)
    bmap = {b: np.array(v, dtype=np.int64) for b, v in bmap.items()}
    return f1_b, bmap

def _perm(f1_genes, gw_genes, get_val, obs_rho, n_perm=N_PERM):
    f1_b, bmap = _bucketize(f1_genes, gw_genes)
    used = np.zeros(len(gw_genes), dtype=bool)
    rho = np.zeros(n_perm)
    valid = 0
    for i in range(n_perm):
        used[:] = False
        picked = np.empty(len(f1_genes), dtype=np.int64)
        ok = True
        for k, fb in enumerate(f1_b):
            j = -1
            for span in range(10):
                lo, hi = max(0, fb - span), min(9, fb + span)
                cands = []
                for b in range(lo, hi + 1):
                    arr = bmap.get(b)
                    if arr is not None:
                        avail = arr[~used[arr]]
                        if len(avail):
                            cands.append(avail)
                if cands:
                    pool = np.concatenate(cands)
                    j = pool[rng.integers(len(pool))]
                    break
            if j < 0:
                ok = False
                break
            picked[k] = j
            used[j] = True
        if ok and len(picked) >= 5:
            vals = np.array([get_val(gw_genes[j]) for j in picked])
            pv = np.array([pig_mean[gw_genes[j]] for j in picked])
            if np.all(np.isfinite(vals)) and np.all(np.isfinite(pv)):
                rho[valid], _ = spearmanr(pv, vals)
                valid += 1
    rho = rho[:valid]
    p = float(np.mean(np.abs(rho) >= abs(obs_rho))) if len(rho) else np.nan
    return p

print(f"\nRunning expression-matched permutation ({N_PERM} iters)...", flush=True)
p_hu_perm = _perm(f1_in_hu, gw_hu, lambda g: hu_mean[g], r_hu_f1)
print(f"Pig-human: empirical p={p_hu_perm:.4f}", flush=True)
p_mo_perm = _perm(f1_in_mo, gw_mo, lambda g: mo_mean[mouse_symbol(g)], r_mo_f1)
print(f"Pig-mouse:  empirical p={p_mo_perm:.4f}", flush=True)

# ── 7. Save new ortholog table (Table S3) ──
records = []
for g in sorted(f1_top44):
    rec = {"pig_gene": g, "factor": "F1"}
    if g in pig_mean.index:
        rec["pig_GV_mean"] = float(pig_mean[g])
    if human_match(g):
        rec["human_mean"] = float(hu_mean[g])
        rec["human_var"] = float(hu.loc[g].var())
    if mouse_match(g):
        rec["mouse_mean"] = float(mo_mean[mouse_symbol(g)])
        rec["mouse_var"] = float(mo.loc[mouse_symbol(g)].var())
    rec["in_human"] = human_match(g)
    rec["in_mouse"] = mouse_match(g)
    records.append(rec)
ortho_df = pd.DataFrame(records)
ortho_df.to_csv(os.path.join(OUT, "cross_species_orthologs_v2.csv"), index=False)
print(f"\nSaved cross_species_orthologs_v2.csv: {len(ortho_df)} rows")

# ── 8. Save permutation results ──
summary = pd.DataFrame({
    "test": ["F1_observed", "F1_observed", "genome_wide", "genome_wide"],
    "species": ["human", "mouse", "human", "mouse"],
    "rho": [r_hu_f1, r_mo_f1, r_hu_gw, r_mo_gw],
    "p_value": [p_hu_f1, p_mo_f1, p_hu_gw, p_mo_gw],
    "expr_matched_perm_p": [p_hu_perm, p_mo_perm, np.nan, np.nan],
    "n_genes": [len(f1_in_hu), len(f1_in_mo), len(gw_hu), len(gw_mo)],
})
summary.to_csv(os.path.join(OUT, "cross_species_permutation_results_v2.csv"), index=False)
print("Saved cross_species_permutation_results_v2.csv")
print("\n=== FINAL NUMBERS ===")
print(f"F1 overlap both species: {len(f1_both)}/{len(f1_top44)}")
print(f"F1 human n={len(f1_in_hu)}, mouse n={len(f1_in_mo)}")
print(f"F1 pig-human rho={r_hu_f1:.3f}, nominal p={p_hu_f1:.2e}, perm p={p_hu_perm:.3f}")
print(f"F1 pig-mouse rho={r_mo_f1:.3f}, nominal p={p_mo_f1:.2e}, perm p={p_mo_perm:.3f}")
print(f"GW pig-human rho={r_hu_gw:.3f}, p={p_hu_gw:.2e}, n={len(gw_hu)}")
print(f"GW pig-mouse rho={r_mo_gw:.3f}, p={p_mo_gw:.2e}, n={len(gw_mo)}")
print("DONE")
