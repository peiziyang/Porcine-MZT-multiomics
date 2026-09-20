#!/usr/bin/env python
"""LO-dataset won't work (832 cells all from GSE168106).
Instead: bootstrap CI on stage correlations + leave-one-stage-out."""
import pandas as pd, numpy as np, os
from scipy.stats import spearmanr
warnings = __import__('warnings'); warnings.filterwarnings('ignore')

BASE = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
MF = os.path.join(BASE, 'm4_mofa')
OUT = os.path.join(MF, 'lo_dataset')
os.makedirs(OUT, exist_ok=True)

tproj = pd.read_csv(os.path.join(MF, 'mzt_trajectory_projection.csv'), index_col=0)
meta = pd.read_csv(os.path.join(BASE, 'm2_metadata.csv'))
meta_idx = meta.set_index('cell')
tproj['dataset'] = [meta_idx.loc[c,'dataset'] if c in meta_idx.index else '?' for c in tproj.index]

# All data from GSE168106 - note this as limitation
print("NOTE: All 832 E0-E10 projection cells from GSE168106 (single dataset).")
print("Cannot perform cross-dataset validation.")
print("Instead: bootstrap CI + leave-one-stage-out.\n")

FCOLS = ['F1','F2','F3','F4','F5','F6','F7']
stage_order = ['E0','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10']

# ── Bootstrap CI on stage correlations ──
print("BOOTSTRAP 95% CI (resample cells within each stage, 5000 iter):")
n_boot = 5000
boot_results = {fc: np.zeros(n_boot) for fc in ['F1','F3','F4','F6']}

for i in range(n_boot):
    # Resample cells within each stage
    boot_cells = []
    for s in stage_order:
        sub = tproj[tproj['stage'] == s]
        if len(sub) > 0:
            idx = np.random.choice(len(sub), len(sub), replace=True)
            boot_cells.append(sub.iloc[idx])
    if not boot_cells:
        continue
    boot_df = pd.concat(boot_cells)
    
    for fc in ['F1','F3','F4','F6']:
        means = [boot_df[boot_df['stage'] == s][fc].mean() 
                for s in stage_order if s in boot_df['stage'].values]
        ords = np.arange(len(means))
        if len(means) >= 5:
            r, _ = spearmanr(ords, means)
            boot_results[fc][i] = r

for fc in ['F1','F3','F4','F6']:
    vals = boot_results[fc]
    ci_lo = np.percentile(vals, 2.5)
    ci_hi = np.percentile(vals, 97.5)
    print(f"  {fc}: rho={np.median(vals):+.3f}  95%CI=[{ci_lo:+.3f},{ci_hi:+.3f}]")

# ── Leave-one-stage-out ──
print("\nLEAVE-ONE-STAGE-OUT (drop each stage, recompute correlation):")
lo_results = []
full_rho = {}
for fc in ['F1','F3','F4','F6']:
    means = []
    ords = []
    for i, s in enumerate(stage_order):
        sub = tproj[tproj['stage'] == s]
        if len(sub) > 0:
            means.append(sub[fc].mean())
            ords.append(i)
    full_rho[fc], _ = spearmanr(ords, means)

for drop_s in stage_order:
    remaining = [s for s in stage_order if s != drop_s]
    if len(remaining) < 5:
        continue
    for fc in ['F1','F3','F4','F6']:
        means = []
        ords = []
        for i, s in enumerate(stage_order):
            if s == drop_s:
                continue
            sub = tproj[tproj['stage'] == s]
            if len(sub) > 0:
                means.append(sub[fc].mean())
                ords.append(i)
        if len(means) >= 5:
            r, p = spearmanr(ords, means)
            delta = r - full_rho[fc]
            lo_results.append({'dropped_stage': drop_s, 'factor': fc,
                              'full_rho': full_rho[fc], 'lo_rho': r,
                              'delta': delta, 'n_stages': len(means)})

lo_df = pd.DataFrame(lo_results)
for fc in ['F1','F3','F4','F6']:
    sub = lo_df[lo_df['factor'] == fc]
    max_delta = sub['delta'].abs().max()
    max_delta_stage = sub.loc[sub['delta'].abs().idxmax(), 'dropped_stage']
    print(f"  {fc}: full={full_rho[fc]:+.3f}  max|delta|={max_delta:.3f} (drop {max_delta_stage})")

lo_df.to_csv(os.path.join(OUT, 'lo_stage_sensitivity.csv'), index=False)
print(f"\nSaved to {OUT}")
print("DONE")
