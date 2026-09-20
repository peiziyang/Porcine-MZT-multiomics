#!/usr/bin/env python
"""
Reviewer fix #5: Cross-species background permutation.
Test whether the F1 cross-species correlation exceeds random expectation
using genome-wide background and expression-matched random gene sets.
"""
import pandas as pd, numpy as np, os
from scipy.stats import spearmanr
import warnings; warnings.filterwarnings('ignore')

BASE = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
MF = os.path.join(BASE, 'm4_mofa')
OUT = os.path.join(MF)

# ── Load F1 genes ──
f1_w = pd.read_csv(os.path.join(MF, 'mofa_multiomics_weights_RNA_v2.csv'), index_col=0)
f1_genes_44 = list(f1_w['F1'].abs().sort_values(ascending=False).head(44).index)

# ── Load cross-species data ──
ortho = pd.read_csv(os.path.join(MF, 'cross_species_orthologs.csv'))
print(f"Total ortholog records: {len(ortho)}")

# ── Build genome-wide expression data ──
# Use pig GV mean and human/mouse mean where available
hu_genes_all = ortho[ortho['in_human']].dropna(subset=['pig_GV_mean', 'human_mean'])
mo_genes_all = ortho[ortho['in_mouse']].dropna(subset=['pig_GV_mean', 'mouse_mean'])

# F1 genes with orthologs
hu_f1 = [g for g in f1_genes_44 if g in set(hu_genes_all['pig_gene'])]
mo_f1 = [g for g in f1_genes_44 if g in set(mo_genes_all['pig_gene'])]

hu_f1_df = hu_genes_all[hu_genes_all['pig_gene'].isin(hu_f1)]
mo_f1_df = mo_genes_all[mo_genes_all['pig_gene'].isin(mo_f1)]

# Observed correlations
r_hu_obs, p_hu_obs = spearmanr(hu_f1_df['pig_GV_mean'], hu_f1_df['human_mean'])
r_mo_obs, p_mo_obs = spearmanr(mo_f1_df['pig_GV_mean'], mo_f1_df['mouse_mean'])

print(f"\nObserved F1 correlations:")
print(f"  Pig-Human: rho={r_hu_obs:.3f}, p={p_hu_obs:.2e}, n={len(hu_f1_df)}")
print(f"  Pig-Mouse:  rho={r_mo_obs:.3f}, p={p_mo_obs:.2e}, n={len(mo_f1_df)}")

# ── Permutation test 1: Random gene sets of same size ──
print(f"\nPermutation 1: Random gene sets (n={len(hu_f1)}, 10,000 iterations)...")
n_perm = 10000
hu_rand_rho = np.zeros(n_perm)
mo_rand_rho = np.zeros(n_perm)

hu_genes_pool = list(set(hu_genes_all['pig_gene']))
mo_genes_pool = list(set(mo_genes_all['pig_gene']))

for i in range(n_perm):
    rand_genes_hu = np.random.choice(hu_genes_pool, len(hu_f1), replace=False)
    rand_genes_mo = np.random.choice(mo_genes_pool, len(mo_f1), replace=False)

    sub_hu = hu_genes_all[hu_genes_all['pig_gene'].isin(rand_genes_hu)]
    sub_mo = mo_genes_all[mo_genes_all['pig_gene'].isin(rand_genes_mo)]

    if len(sub_hu) >= 5:
        hu_rand_rho[i], _ = spearmanr(sub_hu['pig_GV_mean'], sub_hu['human_mean'])
    if len(sub_mo) >= 5:
        mo_rand_rho[i], _ = spearmanr(sub_mo['pig_GV_mean'], sub_mo['mouse_mean'])

p_hu_perm = np.mean(hu_rand_rho >= r_hu_obs)
p_mo_perm = np.mean(mo_rand_rho >= r_mo_obs)

print(f"  Pig-Human: empirical p={p_hu_perm:.4f}, "
      f"null mean={np.mean(hu_rand_rho):.3f}, null std={np.std(hu_rand_rho):.3f}")
print(f"  Pig-Mouse:  empirical p={p_mo_perm:.4f}, "
      f"null mean={np.mean(mo_rand_rho):.3f}, null std={np.std(mo_rand_rho):.3f}")

# ── Permutation test 2: Expression-matched random gene sets ──
print(f"\nPermutation 2: Expression-matched random genes...")
f1_hu_expr = hu_f1_df['pig_GV_mean'].values
f1_mo_expr = mo_f1_df['pig_GV_mean'].values

all_hu_expr = hu_genes_all['pig_GV_mean'].values
all_mo_expr = mo_genes_all['pig_GV_mean'].values

hu_rank_rho = np.zeros(n_perm)
mo_rank_rho = np.zeros(n_perm)

for i in range(n_perm):
    # Sample genes with similar expression distribution (decile-matched)
    hu_deciles = np.percentile(all_hu_expr, np.arange(0, 101, 10))
    matched_hu = []
    for f1val in f1_hu_expr:
        d = np.digitize(f1val, hu_deciles) - 1
        d = max(0, min(9, d))
        candidates = hu_genes_all[
            (hu_genes_all['pig_GV_mean'] >= hu_deciles[d]) &
            (hu_genes_all['pig_GV_mean'] < hu_deciles[min(d+1, 10)])
            ]
        if len(candidates) > 0:
            matched_hu.append(candidates.sample(1)['pig_gene'].values[0])
    if len(matched_hu) >= 5:
        sub = hu_genes_all[hu_genes_all['pig_gene'].isin(matched_hu)]
        hu_rank_rho[i], _ = spearmanr(sub['pig_GV_mean'], sub['human_mean'])

    mo_deciles = np.percentile(all_mo_expr, np.arange(0, 101, 10))
    matched_mo = []
    for f1val in f1_mo_expr:
        d = np.digitize(f1val, mo_deciles) - 1
        d = max(0, min(9, d))
        candidates = mo_genes_all[
            (mo_genes_all['pig_GV_mean'] >= mo_deciles[d]) &
            (mo_genes_all['pig_GV_mean'] < mo_deciles[min(d+1, 10)])
            ]
        if len(candidates) > 0:
            matched_mo.append(candidates.sample(1)['pig_gene'].values[0])
    if len(matched_mo) >= 5:
        sub = mo_genes_all[mo_genes_all['pig_gene'].isin(matched_mo)]
        mo_rank_rho[i], _ = spearmanr(sub['pig_GV_mean'], sub['mouse_mean'])

p_hu_match = np.mean(hu_rank_rho[hu_rank_rho != 0] >= r_hu_obs)
p_mo_match = np.mean(mo_rank_rho[mo_rank_rho != 0] >= r_mo_obs)

print(f"  Pig-Human: empirical p={p_hu_match:.4f} (expression-matched)")
print(f"  Pig-Mouse:  empirical p={p_mo_match:.4f} (expression-matched)")

# ── Genome-wide correlation ──
print(f"\nGenome-wide Pig-Human correlation (all {len(hu_genes_all)} orthologs):")
r_hu_gw, p_hu_gw = spearmanr(hu_genes_all['pig_GV_mean'], hu_genes_all['human_mean'])
r_mo_gw, p_mo_gw = spearmanr(mo_genes_all['pig_GV_mean'], mo_genes_all['mouse_mean'])
print(f"  Pig-Human: rho={r_hu_gw:.3f}, p={p_hu_gw:.2e}")
print(f"  Pig-Mouse:  rho={r_mo_gw:.3f}, p={p_mo_gw:.2e}")

# ── Save results ──
summary = pd.DataFrame({
    'test': ['F1_observed', 'F1_observed', 'genome_wide', 'genome_wide',
             'random_perm', 'random_perm', 'expr_matched', 'expr_matched'],
    'species': ['human', 'mouse', 'human', 'mouse', 'human', 'mouse', 'human', 'mouse'],
    'rho': [r_hu_obs, r_mo_obs, r_hu_gw, r_mo_gw,
            np.mean(hu_rand_rho), np.mean(mo_rand_rho),
            np.mean(hu_rank_rho[hu_rank_rho!=0]), np.mean(mo_rank_rho[mo_rank_rho!=0])],
    'p_value': [p_hu_obs, p_mo_obs, p_hu_gw, p_mo_gw,
                p_hu_perm, p_mo_perm, p_hu_match, p_mo_match],
    'n_genes': [len(hu_f1_df), len(mo_f1_df), len(hu_genes_all), len(mo_genes_all),
                len(hu_f1), len(mo_f1), len(hu_f1), len(mo_f1)],
})
summary.to_csv(os.path.join(MF, 'cross_species_permutation_results.csv'), index=False)
print(f"\nResults saved to: {os.path.join(MF, 'cross_species_permutation_results.csv')}")
print("DONE")
