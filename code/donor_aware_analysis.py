#!/usr/bin/env python
"""
Donor-aware sensitivity analysis for the GV-oocyte MOFA+ model.

1. Per-donor oocyte counts (AF1-AF5)
2. Donor-wise factor score distributions + Kruskal-Wallis test per factor
3. Donor-label permutation: shuffle donor labels, recompute KW H for F1-F7,
   permutation P = fraction of null H >= observed H
4. Leave-one-donor-out core gene stability (archived p1_lo_donor_stability.csv)

Saves to data/processed/m4_mofa/donor_analysis/.
"""
import os, warnings
import numpy as np
import pandas as pd
from scipy.stats import kruskal
warnings.filterwarnings('ignore')

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
MF = os.path.join(BASE, "data", "processed", "m4_mofa")
OUT = os.path.join(MF, "donor_analysis")
os.makedirs(OUT, exist_ok=True)

factors = pd.read_csv(os.path.join(MF, "mofa_multiomics_factors_v2.csv"), index_col=0)
cells = list(factors.index)
donor = pd.Series([c.split('_')[0] for c in cells], index=cells)
print("Per-donor oocyte counts:")
print(donor.value_counts().sort_index().to_string())

# 1. donor summary
summary = pd.DataFrame({'donor': donor.value_counts().index,
                        'n_oocytes': donor.value_counts().values}).sort_values('donor')
summary.to_csv(os.path.join(OUT, "donor_oocyte_counts.csv"), index=False)

# 2. donor-wise factor scores
donor_score = factors.copy()
donor_score['donor'] = donor
donor_mean = donor_score.groupby('donor').mean()
donor_mean.to_csv(os.path.join(OUT, "factor_scores_by_donor.csv"))

# 3. Observed KW + donor-label permutation
print("\nDonor-wise factor score Kruskal-Wallis test + donor-label permutation:")
rng = np.random.default_rng(42)
n_perm = 1000
rows = []
for fc in factors.columns:
    obs_groups = [factors.loc[cells, fc][donor == d].values for d in donor.unique()]
    h_obs, p_nom = kruskal(*obs_groups)
    null_h = np.zeros(n_perm)
    for i in range(n_perm):
        perm_donor = donor.values.copy()
        rng.shuffle(perm_donor)
        g = [factors.loc[cells, fc].values[perm_donor == d] for d in donor.unique()]
        null_h[i], _ = kruskal(*g)
    p_perm = (np.sum(null_h >= h_obs) + 1) / (n_perm + 1)
    rows.append({'factor': fc, 'KW_H': h_obs, 'KW_nominal_p': p_nom,
                 'KW_permutation_p': p_perm})
    print(f"  {fc}: H={h_obs:.2f}, nominal p={p_nom:.4f}, perm p={p_perm:.4f}")
kw_df = pd.DataFrame(rows)
kw_df.to_csv(os.path.join(OUT, "factor_donor_kruskal.csv"), index=False)

# 4. Leave-one-donor-out stability (archived re-trained model results)
print("\nLeave-one-donor-out core gene stability (p1_lo_donor_stability.csv):")
loo = pd.read_csv(os.path.join(MF, "p1_lo_donor_stability.csv"))
print(loo[['excluded_donor', 'core_f1_present', 'core_f1_total']].to_string(index=False))
print(f"  Mean recovery: {loo['core_f1_present'].mean():.1f}/{loo['core_f1_total'].iloc[0]}")

print("\nDONE. Outputs in", OUT)
