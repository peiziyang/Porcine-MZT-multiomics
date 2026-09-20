"""
M5 upgrade: Waddington-OT lineage reconstruction on the Harmony-corrected atlas.

Design:
- Use GSE168106 in_vivo stages E0-E8 as the developmental backbone.
- Standardise the 50-D Harmony space to unit variance so every PC contributes equally.
- Compute entropic-OT (Sinkhorn) couplings between consecutive stages.
- Cluster terminal E8 cells into K=3 fate clusters; propagate coupling mass forward
  to obtain, for each cell, a probability distribution over the 3 terminal fates.
- Fate entropy = plasticity; max fate probability = commitment.
- Map IVF/PA 1C/2C/4C/8C to E2/E3/E4/E5 and compute their OT developmental
  displacement as RMS distance to the corresponding in_vivo stage distribution
  in the standardised Harmony space.

Outputs:
  m5_waddington_ot.png            4-panel figure
  m5_ot_vivo_fate.csv             in_vivo E0-E8 cells + fate metrics
  m5_ot_displacement.csv          per-cell displacement for in_vivo/IVF/PA
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ot
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import os, warnings
warnings.filterwarnings('ignore')

# 注意：Waddington-OT 需要 IVF/PA 细胞（GSE164812）作为映射对象，
# 因此使用"含 GSE164812"的 Harmony embedding（1,955 细胞版），
# 而非图谱重跑后移除 GSE164812 的 1,910 细胞版。
OUT = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
Z = pd.read_csv(f'{OUT}/m2_harmony_corrected_with_gse164812.csv', index_col=0)
meta = pd.read_csv(f'{OUT}/m2_metadata_with_gse164812.csv', index_col=0)
common = Z.index.intersection(meta.index)
Z = Z.loc[common].values
meta = meta.loc[common].copy()

early_stages = ['E0','E1','E2','E3','E4','E5','E6','E7','E8']
stage_order = {s:i for i,s in enumerate(early_stages)}

# ---------- 1. Standardise HPC space on the in_vivo backbone ----------
vivo_mask = (meta['condition']=='in_vivo') & meta['stage'].isin(early_stages)
idx_vivo = meta[vivo_mask].index
Z_vivo_raw = Z[vivo_mask]
stages_vivo = meta.loc[idx_vivo, 'stage'].values

scaler = StandardScaler()
Z_vivo = scaler.fit_transform(Z_vivo_raw)
Z_all = scaler.transform(Z)  # standardised coordinates for all cells

# replace in-vivo rows in Z_all with the already-standardised version (identical)
print(f"in_vivo E0-E8 cells: {Z_vivo.shape[0]}")
stage_cells = {s: Z_vivo[stages_vivo==s] for s in early_stages}
for s in early_stages:
    print(f"  {s}: {stage_cells[s].shape[0]} cells")

# ---------- 2. Terminal fate clusters (E8) ----------
K = 3
kmeans_e8 = KMeans(n_clusters=K, random_state=42, n_init=10).fit(stage_cells['E8'])
e8_cluster_labels = kmeans_e8.labels_
print(f"E8 fate clusters: {np.bincount(e8_cluster_labels)}")

# ---------- 3. Entropic OT couplings between consecutive stages ----------
couplings = {}
for i in range(len(early_stages)-1):
    s0, s1 = early_stages[i], early_stages[i+1]
    X0, X1 = stage_cells[s0], stage_cells[s1]
    n0, n1 = X0.shape[0], X1.shape[0]
    if n0 < 2 or n1 < 2:
        print(f"  skip {s0}->{s1}: too few cells")
        continue
    C = ot.dist(X0, X1, metric='sqeuclidean')
    a = np.ones(n0) / n0
    b = np.ones(n1) / n1
    epsilon = 0.05 * float(np.median(C))
    P = ot.sinkhorn(a, b, C, reg=epsilon, numItermax=1000, stopThr=1e-5)
    T = P / P.sum(axis=1, keepdims=True)   # row-stochastic transition matrix
    couplings[(s0,s1)] = T
    print(f"  {s0}->{s1}: coupling {P.shape}, eps={epsilon:.3f}, transport cost={np.sum(P*C):.3f}")

# ---------- 4. Forward fate probabilities to E8 fate clusters ----------
# For each stage, compute cells x E8_cells forward probabilities, then aggregate to clusters.
fate_to_clusters = {}
# E8: identity over its own cells
fate_to_clusters['E8'] = np.eye(stage_cells['E8'].shape[0])
for s in reversed(early_stages[:-1]):
    s_next = early_stages[early_stages.index(s)+1]
    if (s,s_next) not in couplings:
        continue
    # cells_s x cells_next  ->  cells_s x cells_E8
    fate_to_clusters[s] = couplings[(s,s_next)] @ fate_to_clusters[s_next]

entropy = np.full(Z_vivo.shape[0], np.nan)
maxprob = np.full(Z_vivo.shape[0], np.nan)
for s in early_stages:
    if s not in fate_to_clusters:
        continue
    P_cells_to_e8 = fate_to_clusters[s]  # cells_at_s x E8_cells
    # aggregate to clusters: cells_at_s x K
    P_clusters = np.zeros((P_cells_to_e8.shape[0], K))
    for k in range(K):
        P_clusters[:, k] = P_cells_to_e8[:, e8_cluster_labels==k].sum(axis=1)
    ep = -np.sum(P_clusters * np.log(P_clusters + 1e-12), axis=1)
    mp = P_clusters.max(axis=1)
    mask = stages_vivo == s
    entropy[mask] = ep
    maxprob[mask] = mp

meta_ot = meta.loc[idx_vivo].copy()
meta_ot['fate_entropy'] = entropy
meta_ot['max_fate_prob'] = maxprob

# ---------- 5. IVF/PA developmental displacement ----------
stage_map = {'1C':'E2','2C':'E3','4C':'E4','8C':'E5'}
displacements = []

# in_vivo self-reference (same-stage cells vs same-stage distribution)
for s in ['E2','E3','E4','E5']:
    X_ref = stage_cells[s]
    for i, z in enumerate(X_ref):
        d2 = np.sum((X_ref - z)**2, axis=1)
        displacements.append({
            'condition': 'in_vivo', 'stage': s, 'ref_stage': s,
            'cell': f'{s}_{i}', 'displacement': float(np.sqrt(d2.mean()))
        })

# IVF / PA
for cond in ['IVF','PA']:
    mask = meta['condition'] == cond
    idx_cond = meta[mask].index
    Z_cond = Z_all[mask]
    stages_cond = meta.loc[idx_cond, 'stage'].values
    for cell, z, s in zip(idx_cond, Z_cond, stages_cond):
        if s not in stage_map:
            continue
        ref = stage_map[s]
        X_ref = stage_cells[ref]
        d2 = np.sum((X_ref - z)**2, axis=1)
        displacements.append({
            'condition': cond, 'stage': s, 'ref_stage': ref,
            'cell': cell, 'displacement': float(np.sqrt(d2.mean()))
        })

disp_df = pd.DataFrame(displacements)

# ---------- 6. Plots ----------
umap_df = pd.read_csv(f'{OUT}/m2_harmony_umap_with_gse164812.csv', index_col=0)
coords = umap_df.loc[idx_vivo, ['UMAP1_harmony','UMAP2_harmony']].values

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# A: UMAP colored by fate entropy
ax = axes[0,0]
sc = ax.scatter(coords[:,0], coords[:,1], c=meta_ot['fate_entropy'],
                cmap='viridis_r', s=5, alpha=0.8, rasterized=True)
ax.set_title('A. Waddington-OT fate entropy (in vivo E0–E8)')
ax.set_xticks([]); ax.set_yticks([])
plt.colorbar(sc, ax=ax, shrink=0.6, label='Fate entropy (nats)')

# B: entropy & max prob vs stage
ax = axes[0,1]
summ = meta_ot.groupby('stage').agg({'fate_entropy':['mean','sem'],'max_fate_prob':['mean','sem']}).reindex(early_stages)
x = np.arange(len(early_stages))
ax.errorbar(x, summ[('fate_entropy','mean')], yerr=summ[('fate_entropy','sem')],
            marker='o', color='#2ca02c', linewidth=2, markersize=8, label='Fate entropy')
ax.set_ylabel('Mean fate entropy (nats)', color='#2ca02c')
ax.tick_params(axis='y', labelcolor='#2ca02c')
ax2 = ax.twinx()
ax2.errorbar(x, summ[('max_fate_prob','mean')], yerr=summ[('max_fate_prob','sem')],
             marker='s', color='#d62728', linewidth=2, markersize=7, label='Max fate probability')
ax2.set_ylabel('Mean max fate probability', color='#d62728')
ax2.tick_params(axis='y', labelcolor='#d62728')
ax.set_xticks(x); ax.set_xticklabels(early_stages)
ax.set_xlabel('Developmental stage')
ax.set_title('B. Fate plasticity decreases; commitment increases')

# C: displacement boxplot
ax = axes[1,0]
positions, data_list, labels = [], [], []
pos = 0
colors = []
for s in ['E2','E3','E4','E5']:
    for c, color in [('in_vivo','#2ca02c'), ('IVF','#d62728'), ('PA','#1f77b4')]:
        vals = disp_df[(disp_df['condition']==c)&(disp_df['ref_stage']==s)]['displacement'].values
        data_list.append(vals)
        positions.append(pos)
        labels.append(f'{c}\n{s}')
        colors.append(color)
        pos += 1
    pos += 0.5

bp = ax.boxplot(data_list, positions=positions, widths=0.5, patch_artist=True,
                showfliers=False, medianprops={'color':'black','linewidth':1.5})
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color); patch.set_alpha(0.7)
for vals, p in zip(data_list, positions):
    jitter = np.random.normal(p, 0.04, size=len(vals))
    ax.scatter(jitter, vals, c='black', s=8, alpha=0.3, zorder=3)
ax.set_xticks(positions)
ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel('OT developmental displacement (std. HPC units)')
ax.set_title('C. PA deviates most from the in vivo reference')

# D: UMAP colored by max fate probability
ax = axes[1,1]
sc2 = ax.scatter(coords[:,0], coords[:,1], c=meta_ot['max_fate_prob'],
                 cmap='plasma', s=5, alpha=0.8, rasterized=True)
ax.set_title('D. Max probability toward a terminal E8 fate')
ax.set_xticks([]); ax.set_yticks([])
plt.colorbar(sc2, ax=ax, shrink=0.6, label='Max fate probability')

plt.tight_layout()
fig.savefig(f'{OUT}/m5_waddington_ot.png', dpi=150, bbox_inches='tight')
plt.close()
print("saved m5_waddington_ot.png")

# ---------- 7. Save outputs ----------
meta_ot.to_csv(f'{OUT}/m5_ot_vivo_fate.csv')
disp_df.to_csv(f'{OUT}/m5_ot_displacement.csv')
print("saved m5_ot_vivo_fate.csv, m5_ot_displacement.csv")

print("\nMean displacement (std. HPC units) by condition & reference stage:")
print(disp_df.groupby(['condition','ref_stage'])['displacement'].agg(['mean','std','count']).round(3))

print("\nOverall mean displacement across matched stages E2-E5:")
print(disp_df[disp_df['ref_stage'].isin(['E2','E3','E4','E5'])].groupby('condition')['displacement'].agg(['mean','std','count']).round(3))

print("\nFate entropy by stage:")
print(meta_ot.groupby('stage')['fate_entropy'].agg(['mean','std','count']).round(3))
print("DONE")
