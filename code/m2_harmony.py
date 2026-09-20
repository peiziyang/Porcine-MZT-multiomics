"""
M2 + Harmony: 用 Harmony 替换原 PCA 空间 dataset 均值中心化做批次校正，
重做 Fig.1 UMAP，并定量对比"批次消除效果 vs 发育阶段主轴保留度"。

输入（均为已存产物，无需重跑原始整合）：
  data/processed/m2_pca.csv       (1955 cells x 50 PC, 未校正)
  data/processed/m2_metadata.csv  (含 dataset 批次标签 / stage_group 发育标签)

输出：
  m2_umap_harmony_by_dataset.png / by_stage_group.png / by_source.png
  m2_umap_by_dataset_pca.png     (原 PCA 中心化版，作对照)
  m2_harmony_comparison.png      (定量对比柱状图)
  m2_harmony_metrics.csv         (数字依据)
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from sklearn.linear_model import LinearRegression
import os

OUT = '/mnt/e/Workbuddy/2026-07-27-11-58-27/data/processed'
pca = pd.read_csv(f'{OUT}/m2_pca.csv', index_col=0)
meta = pd.read_csv(f'{OUT}/m2_metadata.csv', index_col=0)
common = pca.index.intersection(meta.index)
X = pca.loc[common].values.astype(float)
meta = meta.loc[common].copy()
print(f"cells={X.shape[0]}, pcs={X.shape[1]}")

# ---- Harmony on the 50-PC embedding, batch = dataset ----
# harmonypy 2.0 run_harmony expects data_mat as PCs x cells; returns Z_corr as cells x PCs
from harmonypy import run_harmony
ho = run_harmony(X.T, meta, vars_use='dataset', verbose=False)
Z = np.asarray(ho.Z_corr, dtype=float)   # cells x PCs
print(f"harmony corrected embedding: {Z.shape}")

# Save corrected embedding for downstream OT/GRN etc.
pca_corrected = pd.DataFrame(Z, index=meta.index, columns=[f'HPC{i+1}' for i in range(Z.shape[1])])
pca_corrected.to_csv(f'{OUT}/m2_harmony_corrected.csv')
print("saved m2_harmony_corrected.csv")

# ---- UMAP (raw PCA vs Harmony) ----
import umap
red = umap.UMAP(n_components=2, random_state=42, n_neighbors=30, min_dist=0.3)
Xu = red.fit_transform(X)
Zu = red.fit_transform(Z)

# Save UMAP coords
umap_df = pd.DataFrame(np.column_stack([Xu, Zu]), index=meta.index,
                       columns=['UMAP1_raw','UMAP2_raw','UMAP1_harmony','UMAP2_harmony'])
umap_df.to_csv(f'{OUT}/m2_harmony_umap.csv')
print("saved m2_harmony_umap.csv")

def make_plot(coords, color_col, title, fname, cmap='tab20'):
    fig, ax = plt.subplots(figsize=(10, 8))
    cats = meta[color_col].astype(str)
    uc = sorted(cats.unique())
    cmap_obj = plt.get_cmap(cmap, len(uc))
    for i, c in enumerate(uc):
        m = (cats == c).values
        ax.scatter(coords[m, 0], coords[m, 1], s=3, alpha=0.7,
                   color=cmap_obj(i), label=c, rasterized=True)
    ax.legend(markerscale=5, fontsize=7, loc='upper left',
              bbox_to_anchor=(1.02, 1), frameon=False)
    ax.set_title(title, fontsize=14); ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout(); fig.savefig(f'{OUT}/{fname}', dpi=150, bbox_inches='tight'); plt.close()
    print("saved", fname)

make_plot(Xu, 'dataset', 'PCA-centering (current method)\ncolored by dataset', 'm2_umap_by_dataset_pca.png')
make_plot(Zu, 'dataset', 'Harmony batch correction\ncolored by dataset', 'm2_umap_harmony_by_dataset.png')
make_plot(Zu, 'stage_group', 'Harmony batch correction\ncolored by developmental stage', 'm2_umap_harmony_by_stage_group.png')
make_plot(Zu, 'source', 'Harmony batch correction\ncolored by reproductive source', 'm2_umap_harmony_by_source.png')

# ---- quantitative comparison ----
def var_partition(emb, cat):
    dummies = pd.get_dummies(meta[cat].astype(str))
    if dummies.shape[1] < 2:
        return 0.0
    r2 = [LinearRegression().fit(dummies.values, emb[:, d]).score(dummies.values, emb[:, d])
          for d in range(emb.shape[1])]
    return float(np.mean(r2))

def sil(emb, cat):
    return float(silhouette_score(emb, meta[cat].astype(str).values, metric='euclidean'))

batch_raw, batch_har = sil(X, 'dataset'), sil(Z, 'dataset')
stage_raw, stage_har = sil(X, 'stage_group'), sil(Z, 'stage_group')
vb_raw, vb_har = var_partition(X, 'dataset'), var_partition(Z, 'dataset')
vs_raw, vs_har = var_partition(X, 'stage_group'), var_partition(Z, 'stage_group')

print(f"Silhouette dataset: raw={batch_raw:.3f} -> harmony={batch_har:.3f}")
print(f"Silhouette stage:   raw={stage_raw:.3f} -> harmony={stage_har:.3f}")
print(f"VarPartition batch R2: raw={vb_raw:.4f} -> harmony={vb_har:.4f}")
print(f"VarPartition stage R2: raw={vs_raw:.4f} -> harmony={vs_har:.4f}")

pd.DataFrame({
    'metric': ['sil_dataset', 'sil_stage', 'varp_batch', 'varp_stage'],
    'raw': [batch_raw, stage_raw, vb_raw, vs_raw],
    'harmony': [batch_har, stage_har, vb_har, vs_har],
}).to_csv(f'{OUT}/m2_harmony_metrics.csv', index=False)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].bar(['raw\n(PCA-center)', 'Harmony'], [batch_raw, batch_har], color=['#bbb', '#2ca02c'])
axes[0].set_title('Dataset (batch) silhouette\nlower = better mixing'); axes[0].set_ylim(0, 1)
for i, v in enumerate([batch_raw, batch_har]):
    axes[0].text(i, v + 0.02, f'{v:.3f}', ha='center', fontsize=10)
axes[1].bar(['raw\n(PCA-center)', 'Harmony'], [stage_raw, stage_har], color=['#bbb', '#1f77b4'])
axes[1].set_title('Stage-group silhouette\nhigher = better separation'); axes[1].set_ylim(0, 1)
for i, v in enumerate([stage_raw, stage_har]):
    axes[1].text(i, v + 0.02, f'{v:.3f}', ha='center', fontsize=10)
plt.tight_layout(); fig.savefig(f'{OUT}/m2_harmony_comparison.png', dpi=150); plt.close()
print("saved m2_harmony_comparison.png")
print("DONE")
