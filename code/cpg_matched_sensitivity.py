#!/usr/bin/env python
"""
CpG coverage matched sensitivity analysis.

Question: are F1 genes' higher informative promoter CpG counts explained by
gene length, expression level, or promoter CpG density?

Approach:
  1. Compute gene length (merged exon span) from the GTF.
  2. Get GV-oocyte mean expression from rna_view_for_mofa.csv.
  3. Bin background genes by expression decile x gene-length quartile.
  4. For each F1 gene with CpG data, sample matched background genes from the
     same expression x length bin; compare informative CpG counts.
  5. Also compare promoter CpG DENSITY (informative CpG / gene length) directly.

Saves to data/processed/m4_mofa/donor_analysis/../cpg_matched_analysis/
"""
import os, re, gzip, warnings
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
warnings.filterwarnings('ignore')

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
MF = os.path.join(BASE, "data", "processed", "m4_mofa")
OUT = os.path.join(MF, "cpg_matched_analysis")
os.makedirs(OUT, exist_ok=True)
GTF = os.path.join(BASE, "data", "raw", "pig_genes.gtf.gz")

# 1. Gene length (merged exon span) from GTF
print("Computing gene lengths from GTF (merged exon span)...")
gene_exons = {}
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 9 or parts[2] != 'exon':
            continue
        gid_m = re.search(r'gene_id "([^"]+)"', parts[8])
        if not gid_m:
            continue
        gid = gid_m.group(1)
        start, end = int(parts[3]), int(parts[4])
        gene_exons.setdefault(gid, []).append((start, end))

def merged_span(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total, cur_s, cur_e = 0, intervals[0][0], intervals[0][1]
    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    total += cur_e - cur_s + 1
    return total

gene_len = {g: merged_span(exons) for g, exons in gene_exons.items()}
print(f"  {len(gene_len)} genes with length")

# 2. Per-gene CpG stats
pg = pd.read_csv(os.path.join(MF, "per_gene_cpg_stats.csv"))
pg['gene_len'] = pg['gene_id'].map(gene_len)
print(f"  per_gene_cpg_stats: {len(pg)} rows, {pg['gene_len'].notna().sum()} with length")

# 3. GV expression
rna = pd.read_csv(os.path.join(MF, "rna_view_for_mofa.csv"), index_col=0)
rna_mean = rna.mean(axis=0)
pg['gv_expr'] = pg['gene_name'].map(rna_mean)
print(f"  {pg['gv_expr'].notna().sum()} genes with GV expression")

# 4. F1 genes from Table S1
s1 = pd.read_csv(os.path.join(BASE, "manuscript", "submission", "tables", "Table_S1_cpg_stats.csv"))
f1_genes = set(s1['gene_name'])
pg['is_f1'] = pg['gene_name'].isin(f1_genes)

# Restrict to genes with length + expression for matching
df = pg.dropna(subset=['gene_len', 'gv_expr']).copy()
df = df[df['gene_len'] > 0]
print(f"\n  Genes available for matching: {len(df)} (F1: {df['is_f1'].sum()})")

# Observed difference
f1 = df[df['is_f1']]
bg = df[~df['is_f1']]
obs_med_f1 = f1['mean_promoter_ncpg'].median()
obs_med_bg = bg['mean_promoter_ncpg'].median()
print(f"  Observed median informative CpG: F1={obs_med_f1:.1f}, background={obs_med_bg:.1f}")
u, p_obs = mannwhitneyu(f1['mean_promoter_ncpg'], bg['mean_promoter_ncpg'], alternative='two-sided')
print(f"  Unmatched MWU P = {p_obs:.2e}")

# 5. Expression x length matched comparison
df['expr_bin'] = pd.qcut(df['gv_expr'], 10, labels=False, duplicates='drop')
df['len_bin'] = pd.qcut(df['gene_len'], 10, labels=False, duplicates='drop')
f1 = df[df['is_f1']]
bg = df[~df['is_f1']]

rng = np.random.default_rng(42)
n_perm = 1000
# For each permutation, for each F1 gene, sample 1 matched background gene from same (expr_bin, len_bin)
null_medians = np.zeros(n_perm)
for i in range(n_perm):
    matched = []
    for _, row in f1.iterrows():
        pool = bg[(bg['expr_bin'] == row['expr_bin']) & (bg['len_bin'] == row['len_bin'])]
        if len(pool) > 0:
            matched.append(pool.sample(1, random_state=rng).iloc[0])
    if len(matched) >= 20:
        null_medians[i] = np.median([m['mean_promoter_ncpg'] for m in matched])
obs_match_med = np.median(null_medians)
p_expr_len = (np.sum(null_medians >= obs_med_f1) + 1) / (n_perm + 1)
print(f"\n  Expression x length matched background median: {obs_match_med:.1f}")
print(f"  F1 observed median vs matched background: {obs_med_f1:.1f} vs {obs_match_med:.1f}")
print(f"  Permutation P (F1 >= matched null): {p_expr_len:.4f}")

# 6. Promoter CpG DENSITY (per kb gene)
df['cpg_density_kb'] = df['mean_promoter_ncpg'] / (df['gene_len'] / 1000)
f1_d = df[df['is_f1']]['cpg_density_kb']
bg_d = df[~df['is_f1']]['cpg_density_kb']
u2, p_dens = mannwhitneyu(f1_d, bg_d, alternative='two-sided')
print(f"\n  CpG density per kb gene: F1 median={f1_d.median():.2f}, background median={bg_d.median():.2f}")
print(f"  MWU P (density) = {p_dens:.2e}")

# 7. Correlation of informative CpG with expression and length (confounder check)
r_expr, p_expr = spearmanr(df['gv_expr'], df['mean_promoter_ncpg'])
r_len, p_len = spearmanr(df['gene_len'], df['mean_promoter_ncpg'])
print(f"\n  Spearman: informative CpG vs expression r={r_expr:.3f} (P={p_expr:.1e})")
print(f"  Spearman: informative CpG vs gene length r={r_len:.3f} (P={p_len:.1e})")

# Save
pd.DataFrame({
    'metric': ['F1_median_infoCpG', 'bg_median_infoCpG', 'matched_bg_median_infoCpG',
               'expr_len_matched_permP', 'F1_density_per_kb', 'bg_density_per_kb',
               'density_MWU_P', 'corr_expr_rho', 'corr_expr_P', 'corr_len_rho', 'corr_len_P'],
    'value': [obs_med_f1, obs_med_bg, obs_match_med, p_expr_len,
              f1_d.median(), bg_d.median(), p_dens,
              r_expr, p_expr, r_len, p_len]
}).to_csv(os.path.join(OUT, "cpg_matched_sensitivity.csv"), index=False)

# matched null distribution
np.save(os.path.join(OUT, "matched_null_medians.npy"), null_medians)
print("\nSaved to", OUT)
print("DONE")
