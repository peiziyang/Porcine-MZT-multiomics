#!/usr/bin/env python
"""Verify Table 2 empirical P-values by re-running the gene-set permutation test."""
import os
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control as bh_fdr

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
DESEQ = os.path.join(BASE, "data/processed/deseq2/m15_deseq2_overall_IVFvsPA_stageAdj.csv")
MF = os.path.join(BASE, "data/processed/m4_mofa")
N_PERM = 10000
rng = np.random.default_rng(42)

# ── stage-adjusted DESeq2 log2FC by symbol ──
de = pd.read_csv(DESEQ)
de = de[de['gene_symbol'].notna() & ~de['gene_symbol'].astype(str).str.startswith('ENSSSCG')]
lfc_by_sym = {}
for _, r in de.iterrows():
    lfc_by_sym[str(r['gene_symbol'])] = float(r['log2FoldChange'])

# ── pig GV mean expression (decile matching) ──
pig = pd.read_csv(os.path.join(MF, "rna_view_for_mofa.csv"), index_col=0)
pig_mean = pig.mean(axis=0)

# background = mappable DESeq2 genes that also have pig GV expression
bg_genes = [g for g in lfc_by_sym if g in pig_mean.index]
print(f"background n = {len(bg_genes)}")

bg_expr = pig_mean[bg_genes].values
bg_lfc = np.array([lfc_by_sym[g] for g in bg_genes])
edges = np.percentile(bg_expr, np.arange(0, 101, 10))
bg_bucket = np.clip(np.digitize(bg_expr, edges) - 1, 0, 9)
from collections import defaultdict
bmap = defaultdict(list)
for i, b in enumerate(bg_bucket):
    bmap[int(b)].append(i)
bmap = {b: np.array(v) for b, v in bmap.items()}

# ── gene sets ──
s4 = pd.read_csv(os.path.join(BASE, "manuscript/submission/tables/Table_S4_gene_set_membership.csv"))
sets = {}
for f in ['F1', 'F3', 'F4', 'F6']:
    sets[f] = list(s4[s4[f + '_member'] == True]['gene'])

def gene_set_test(genes, obs_mean):
    genes_in = [g for g in genes if g in lfc_by_sym]
    # expression deciles of the observed genes (based on pig GV expression where available)
    obs_buckets = []
    obs_expr = []
    for g in genes_in:
        if g in pig_mean.index:
            obs_expr.append(pig_mean[g])
        else:
            obs_expr.append(bg_expr.mean())  # fallback to mean
    obs_expr = np.array(obs_expr)
    obs_b = np.clip(np.digitize(obs_expr, edges) - 1, 0, 9)
    n = len(genes_in)
    rho = np.zeros(N_PERM)
    used = np.zeros(len(bg_genes), dtype=bool)
    valid = 0
    for i in range(N_PERM):
        used[:] = False
        picked = np.empty(n, dtype=np.int64)
        ok = True
        for k, fb in enumerate(obs_b):
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
        if ok:
            rho[valid] = np.mean(bg_lfc[picked])
            valid += 1
    rho = rho[:valid]
    p = (1 + np.sum(np.abs(rho) >= abs(obs_mean))) / (1 + valid)
    return p, len(genes_in), obs_mean

print(f"\n{'set':4s} {'mapped':>6s} {'obs_mean':>10s} {'empirical_P':>12s}")
rows = []
for f in ['F1', 'F3', 'F4', 'F6']:
    genes_in = [g for g in sets[f] if g in lfc_by_sym]
    obs_mean = np.mean([lfc_by_sym[g] for g in genes_in])
    p, mapped, om = gene_set_test(sets[f], obs_mean)
    rows.append({'gene_set': f, 'mapped_n': mapped, 'defined_n': 44,
                 'obs_mean_log2FC': om, 'empirical_p': p})
    print(f"{f:4s} {mapped:6d} {om:10.3f} {p:12.4f}")

# BH correction
ps = [r['empirical_p'] for r in rows]
bh = bh_fdr(ps)
for r, b in zip(rows, bh):
    r['bh_adj_p'] = b
    print(f"{r['gene_set']:4s}  BH-adjusted P = {b:.4f}")

# ── 落盘结果文件（数据可复现铁律：permutation P 值必须保存结果文件）──
out_df = pd.DataFrame(rows)
out_df['background_n'] = len(bg_genes)
out_df['n_perm'] = N_PERM
out_df['seed'] = 42
out_dir = os.path.join(MF, 'geneset_test')
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, 'gene_set_permutation_results.csv')
out_df.to_csv(out_path, index=False)
print(f"\n已保存结果文件: {out_path}")
print("（此文件是 Table 2 empirical P / BH 值的权威来源）")
