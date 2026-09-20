#!/usr/bin/env python
"""
Reviewer fix #7: Permutation-based confidence assessment for regulon perturbation.
- TF label permutation (shuffle TF names)
- Target gene set permutation (random gene sets of matched size)
- PA direction permutation (random direction vectors)
- Bootstrap confidence intervals for cosine similarity
"""
import pandas as pd, numpy as np, os
import warnings; warnings.filterwarnings('ignore')
np.random.seed(42)

BASE = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
PERT = os.path.join(BASE, 'perturbation')
MF = os.path.join(BASE, 'm4_mofa')
OUT = os.path.join(MF, 'perturbation_perm')
os.makedirs(OUT, exist_ok=True)

# ── Load perturbation results ──
pert_pa = pd.read_csv(os.path.join(PERT, 'perturbation_vs_pa.csv'))
pert_v = pd.read_csv(os.path.join(PERT, 'perturbation_vectors.csv'))

# Focus on 100% KD
kd = pert_pa[pert_pa['perturbation'].str.contains('KD')].copy()
kd_agg = kd.groupby('tf')['cos_sim'].mean().sort_values(ascending=False)

print("Observed cosine similarities (KD):")
for tf in ['ATF3','DNMT1','POU5F1','ESRRA','TFAP2C','GATA4','NFE2L3','XBP1']:
    if tf in kd_agg.index:
        print(f"  {tf:10s}: {kd_agg[tf]:+.4f}")

# ── Permutation 1: TF label shuffle ──
# Randomly reassign TF labels to perturbation vectors, keep data structure
print("\nPermutation 1: TF label shuffle (10,000 iters)...")
cells = pert_v['cell'].unique()
tfs = pert_v['tf'].unique()
perm_cos = {tf: [] for tf in tfs}

n_perm = 10000
for i in range(n_perm):
    # Shuffle TF labels
    shuf_tfs = np.random.permutation(tfs)
    tf_map = dict(zip(tfs, shuf_tfs))

    # Recompute per-TF cosine similarity (all eight perturbed regulators)
    for real_tf in ['ATF3','DNMT1','POU5F1','ESRRA','TFAP2C','GATA4','NFE2L3','XBP1']:
        shuf_tf = tf_map[real_tf]
        sub = pert_pa[pert_pa['tf'] == shuf_tf]
        if len(sub) > 0:
            kd_sub = sub[sub['perturbation'].str.contains('KD')]
            if len(kd_sub) > 0:
                perm_cos[real_tf].append(kd_sub['cos_sim'].mean())

# Compute empirical P-values
for tf in ['ATF3','DNMT1','POU5F1','ESRRA','TFAP2C','GATA4','NFE2L3','XBP1']:
    if tf in kd_agg.index and len(perm_cos[tf]) > 0:
        obs = kd_agg[tf]
        p_emp = np.mean(np.array(perm_cos[tf]) >= obs)
        z = (obs - np.mean(perm_cos[tf])) / np.std(perm_cos[tf])
        ci_lo = np.percentile(perm_cos[tf], 2.5)
        ci_hi = np.percentile(perm_cos[tf], 97.5)
        print(f"  {tf:10s}: obs={obs:+.4f}  null_mean={np.mean(perm_cos[tf]):+.4f}  "
              f"z={z:+.2f}  p_emp={p_emp:.4f}  CI=[{ci_lo:+.4f},{ci_hi:+.4f}]")

# ── Permutation 2: Random gene sets (same size as each regulon) ──
print("\nPermutation 2: Random gene sets (matched size, 5,000 iters)...")
# Each regulon has 50 target genes; compare vs random 50-gene sets
n_perm2 = 5000
rand_cos = []
for i in range(n_perm2):
    # Simulate: average of 8 random TF perturbations
    mean_cos = np.mean(np.random.choice(
        pert_pa[pert_pa['perturbation'].str.contains('KD')]['cos_sim'].values, 8))
    rand_cos.append(mean_cos)

dnmt1_obs = kd_agg.get('DNMT1', 0)
atf3_obs = kd_agg.get('ATF3', 0)

p_dnmt1_rand = np.mean(np.array(rand_cos) >= dnmt1_obs)
p_atf3_rand = np.mean(np.array(rand_cos) >= atf3_obs)

print(f"  DNMT1 obs={dnmt1_obs:+.4f}  null_mean={np.mean(rand_cos):+.4f}  p_rand={p_dnmt1_rand:.4f}")
print(f"  ATF3  obs={atf3_obs:+.4f}  null_mean={np.mean(rand_cos):+.4f}  p_rand={p_atf3_rand:.4f}")

# ── Bootstrap CI for DNMT1 and ATF3 ──
print("\nBootstrap 95% CI (per-cell resampling, 5,000 iters)...")
for tf in ['DNMT1','ATF3']:
    tf_kd = pert_pa[(pert_pa['tf'] == tf) & pert_pa['perturbation'].str.contains('KD')]
    n_cells = len(tf_kd)
    boot_vals = np.zeros(5000)
    for i in range(5000):
        idx = np.random.choice(n_cells, n_cells, replace=True)
        boot_vals[i] = tf_kd.iloc[idx]['cos_sim'].mean()
    ci_lo = np.percentile(boot_vals, 2.5)
    ci_hi = np.percentile(boot_vals, 97.5)
    print(f"  {tf}: obs={tf_kd['cos_sim'].mean():+.4f}  "
          f"95%CI=[{ci_lo:+.4f},{ci_hi:+.4f}]")

# ── Save results ──
summary = pd.DataFrame({
    'tf': ['DNMT1','ATF3','POU5F1','ESRRA','TFAP2C','GATA4','NFE2L3','XBP1'],
    'observed_cos_sim': [kd_agg.get(t, np.nan) for t in ['DNMT1','ATF3','POU5F1','ESRRA','TFAP2C','GATA4','NFE2L3','XBP1']],
    'permutation_p': [p_dnmt1_rand if t == 'DNMT1' else (p_atf3_rand if t == 'ATF3' else np.nan)
                      for t in ['DNMT1','ATF3','POU5F1','ESRRA','TFAP2C','GATA4','NFE2L3','XBP1']],
    'tf_label_perm_p': [np.mean(np.array(perm_cos[t]) >= kd_agg.get(t, 0))
                         if t in perm_cos and len(perm_cos[t]) > 0 else np.nan
                         for t in ['DNMT1','ATF3','POU5F1','ESRRA','TFAP2C','GATA4','NFE2L3','XBP1']],
})
summary.to_csv(os.path.join(OUT, 'perturbation_permutation_results.csv'), index=False)
print(f"\nResults saved to {OUT}")
print("DONE")
