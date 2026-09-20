#!/usr/bin/env python
"""
Reviewer fix #6: Leave-one-dataset-out projection sensitivity.
Tests whether MOFA+ factor projections are driven by a single dataset
by excluding each dataset in turn and recomputing stage correlations.
"""
import pandas as pd, numpy as np, os
from scipy.stats import spearmanr
import warnings; warnings.filterwarnings('ignore')

BASE = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
MF = os.path.join(BASE, 'm4_mofa')
OUT = os.path.join(MF, 'lo_dataset')

os.makedirs(OUT, exist_ok=True)

# ── Load projection data ──
tproj = pd.read_csv(os.path.join(MF, 'mzt_trajectory_projection.csv'), index_col=0)
meta = pd.read_csv(os.path.join(BASE, 'm2_metadata.csv'))

# Merge dataset info
meta_idx = meta.set_index('cell')
tproj['dataset'] = [meta_idx.loc[c, 'dataset'] if c in meta_idx.index else 'unknown'
                    for c in tproj.index]
tproj['stage'] = tproj['stage'].astype(str)

datasets = sorted([d for d in tproj['dataset'].unique() if d != 'unknown'])
print(f"Datasets in projection: {datasets}")
n_cells_per_ds = tproj['dataset'].value_counts()
print(f"Cells per dataset:\n{n_cells_per_ds}")

FCOLS = ['F1','F2','F3','F4','F5','F6','F7']
stage_order = ['E0','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10']

# ── Full model (all datasets) ──
print("\nFULL MODEL (all datasets):")
full_means = {}
for fc in FCOLS:
    vals = []
    for s in stage_order:
        sub = tproj[tproj['stage'] == s]
        if len(sub) > 0:
            vals.append(sub[fc].mean())
        else:
            vals.append(np.nan)
    full_means[fc] = np.array(vals)

full_results = {}
for fc in FCOLS:
    valid = ~np.isnan(full_means[fc])
    if valid.sum() >= 5:
        r, p = spearmanr(np.arange(len(stage_order))[valid], full_means[fc][valid])
        full_results[fc] = (r, p)
        print(f"  {fc}: rho={r:+.3f}, p={p:.4f}")

# ── LO-dataset ──
print("\nLEAVE-ONE-DATASET-OUT:")
lo_results = []
for drop_ds in datasets:
    sub = tproj[tproj['dataset'] != drop_ds]
    print(f"\n  Excluding {drop_ds} ({len(sub)} cells remain):")

    for fc in ['F1','F3','F4','F6']:
        vals = []
        for s in stage_order:
            sub_s = sub[sub['stage'] == s]
            if len(sub_s) > 0:
                vals.append(sub_s[fc].mean())
            else:
                vals.append(np.nan)
        arr = np.array(vals)
        valid = ~np.isnan(arr)
        if valid.sum() >= 5:
            r, p = spearmanr(np.arange(len(stage_order))[valid], arr[valid])
            full_r, full_p = full_results[fc]
            delta_r = r - full_r
            print(f"    {fc}: rho={r:+.3f} (full={full_r:+.3f}, delta={delta_r:+.3f})  p={p:.4f}")
            lo_results.append({
                'excluded_dataset': drop_ds,
                'factor': fc,
                'n_cells_remain': len(sub),
                'full_rho': full_r,
                'lo_rho': r,
                'delta_rho': delta_r,
                'lo_p': p,
                'n_stages': valid.sum(),
            })

lo_df = pd.DataFrame(lo_results)
lo_df.to_csv(os.path.join(OUT, 'lo_dataset_sensitivity.csv'), index=False)

# Summary: max absolute deviation per factor
print("\nMAX |delta rho| per factor across all dataset exclusions:")
for fc in ['F1','F3','F4','F6']:
    sub_lo = lo_df[lo_df['factor'] == fc]
    max_delta = sub_lo['delta_rho'].abs().max()
    print(f"  {fc}: max |delta| = {max_delta:.3f}")

print(f"\nResults saved to {OUT}")
print("DONE")
