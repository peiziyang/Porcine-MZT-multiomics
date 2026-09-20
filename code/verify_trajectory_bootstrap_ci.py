#!/usr/bin/env python
"""验证正文轨迹 bootstrap CI 并落盘可复现结果（seed=42）。

逻辑同 reviewer_fix_lo_stage.py（每个 stage 内 resample 细胞 → 重算 per-stage pseudobulk
mean → Spearman rho vs stage 顺序），但向量化 + 固定 seed + 落盘。
"""
import pandas as pd, numpy as np, os
from scipy.stats import spearmanr
import warnings; warnings.filterwarnings('ignore')

BASE = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
MF = os.path.join(BASE, 'm4_mofa')
OUT = os.path.join(MF, 'lo_dataset')
os.makedirs(OUT, exist_ok=True)

tproj = pd.read_csv(os.path.join(MF, 'mzt_trajectory_projection.csv'), index_col=0)
stage_order = ['E0','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10']
FCOLS = ['F1','F3','F4','F6']
# 预分组（向量化，避免循环内重复布尔索引）
groups = {s: tproj[tproj['stage'] == s][FCOLS].values
          for s in stage_order if len(tproj[tproj['stage'] == s]) > 0}

n_boot = 5000
seed = 42
rng = np.random.default_rng(seed)
boot = {fc: np.zeros(n_boot) for fc in FCOLS}

for i in range(n_boot):
    sm = np.vstack([groups[s][rng.integers(0, len(groups[s]), len(groups[s]))].mean(axis=0)
                    for s in stage_order])
    ords = np.arange(sm.shape[0])
    for k, fc in enumerate(FCOLS):
        boot[fc][i], _ = spearmanr(ords, sm[:, k])

rows = []
for fc in FCOLS:
    v = boot[fc]
    lo, hi = np.percentile(v, 2.5), np.percentile(v, 97.5)
    rows.append({'factor': fc, 'median_rho': np.median(v),
                 'ci_lo': lo, 'ci_hi': hi, 'n_boot': n_boot, 'seed': seed})
    print(f"  {fc}: median_rho={np.median(v):+.3f}  95%CI=[{lo:+.3f},{hi:+.3f}]")

pd.DataFrame(rows).to_csv(os.path.join(OUT, 'trajectory_bootstrap_ci.csv'), index=False)
print('已落盘:', os.path.join(OUT, 'trajectory_bootstrap_ci.csv'))
print('正文声称: F1 CI=[-0.58,-0.15](全负), F3 CI=[-0.93,-0.86]')
