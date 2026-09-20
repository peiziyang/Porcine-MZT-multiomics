#!/usr/bin/env python
"""完整重跑猪着床前胚胎图谱（移除 GSE164812，得到 1,914 细胞）。

数据流（与 m2_integration.py + m2_harmony.py 一致，但移除 GSE164812）：
  5 数据集加载 → 统一 Ensembl ID → 交集 common genes → QC → StandardScaler
  → PCA(50) → Harmony(batch=dataset) → UMAP(30, min_dist 0.3)

输出（覆盖 data/processed/ 下对应文件）：
  m2_pca.csv               (1914 x 50, 未校正 PCA)
  m2_harmony_corrected.csv (1914 x 50, Harmony 校正后)
  m2_harmony_umap.csv      (1914 x 4, UMAP1/2_raw + UMAP1/2_harmony)
  m2_metadata.csv          (1914 x N, 含 stage_group)
"""
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
import harmonypy
import os
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw'
OUT_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'

# ============ PART 1: 加载 5 个数据集（移除 GSE164812；GSE242553 RefSeq 跳过） ============
print("=" * 60)
print("PART 1: Loading 5 datasets (GSE164812 removed)")
print("=" * 60)

datasets = {}

# --- GSE168106: Han 2022 Cell Res, E0-E14 ---
print("\n[1/5] GSE168106 ...")
df = pd.read_csv(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_genes_counts.txt.gz',
                 sep='\t', index_col=0)
si = pd.read_excel(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_SampleInfo.xlsx')
si_valid = si[si['Vaild'] == 'Yes'].copy()
valid_barcodes = set(si_valid['ID'].values)
df = df[[c for c in df.columns if c in valid_barcodes]]
cell_to_stage = dict(zip(si_valid['ID'], si_valid['Stage_II']))
cell_to_lane = dict(zip(si_valid['ID'], si_valid['Lane']))
metadata = pd.DataFrame({
    'stage': [cell_to_stage.get(c, 'Unknown') for c in df.columns],
    'donor': [cell_to_lane.get(c, 'Unknown') for c in df.columns],
    'dataset': 'GSE168106', 'source': 'in_vivo', 'condition': 'in_vivo',
    'study': 'Han2022',
}, index=df.columns)
datasets['GSE168106'] = {'expr': df, 'meta': metadata, 'format': 'counts'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE112380: Ramos-Ibeas 2019 ---
print("\n[2/5] GSE112380 ...")
df = pd.read_csv(f'{DATA_DIR}/GSE112380/GSE112380_Pig_RawCounts.txt.gz', sep='\t', index_col=0)
metadata = pd.DataFrame({
    'stage': ['preimplantation'] * df.shape[1],
    'donor': ['pooled'] * df.shape[1],
    'dataset': 'GSE112380', 'source': 'in_vivo', 'condition': 'in_vivo',
    'study': 'RamosIbeas2019',
}, index=df.columns)
datasets['GSE112380'] = {'expr': df, 'meta': metadata, 'format': 'counts'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE160334: Du 2021, GV/MII oocytes ---
print("\n[3/5] GSE160334 ...")
df = pd.read_csv(f'{DATA_DIR}/GSE160334/GSE160334_gene_FPKM_matrix.txt.gz', sep='\t', index_col=0)
meta_list = []
for col in df.columns:
    if col.startswith('GV'):
        stage, rep = 'GV', col[2:]
    elif col.startswith('MII'):
        stage, rep = 'MII', col[3:]
    else:
        stage, rep = 'Unknown', '1'
    meta_list.append({'stage': stage, 'donor': f'donor_{rep}', 'condition': 'in_vivo',
                      'source': 'in_vivo', 'dataset': 'GSE160334', 'study': 'Du2021'})
metadata = pd.DataFrame(meta_list, index=df.columns)
datasets['GSE160334'] = {'expr': df, 'meta': metadata, 'format': 'fpkm'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE139512: Kong 2019, oocyte-blastocyst ---
print("\n[4/5] GSE139512 ...")
df = pd.read_csv(f'{DATA_DIR}/GSE139512/GSE139512_fpkm.txt.gz', sep='\t', index_col=0)
meta_list = []
for col in df.columns:
    cl = col.lower()
    if 'oo' in cl or 'oocyte' in cl: stage = 'Oocyte'
    elif '1-cell' in cl or '1_cell' in cl: stage = '1C'
    elif '2-cell' in cl or '2_cell' in cl: stage = '2C'
    elif '4-cell' in cl or '4_cell' in cl: stage = '4C'
    elif '8-cell' in cl or '8_cell' in cl: stage = '8C'
    elif 'mor' in cl: stage = 'Morula'
    elif 'blasto' in cl or 'bl' in cl: stage = 'Blastocyst'
    else: stage = 'Unknown'
    meta_list.append({'stage': stage, 'donor': f'donor_{col.split("_")[0] if "_" in col else "pool"}',
                      'condition': 'in_vivo', 'source': 'in_vivo',
                      'dataset': 'GSE139512', 'study': 'Kong2019'})
metadata = pd.DataFrame(meta_list, index=df.columns)
datasets['GSE139512'] = {'expr': df, 'meta': metadata, 'format': 'fpkm'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE234116: Yuan 2023, GV oocytes (multi-omics RNA) ---
print("\n[5/5] GSE234116 ...")
df = pd.read_csv(f'{DATA_DIR}/GSE234116/GSE234116_raw_counts.txt.gz', sep='\t', index_col=0)
meta_list = []
for col in df.columns:
    animal = col.split('_')[0] if '_' in col else 'unknown'
    meta_list.append({'stage': 'GV', 'donor': animal, 'condition': 'in_vivo',
                      'source': 'in_vivo', 'dataset': 'GSE234116', 'study': 'Yuan2023'})
metadata = pd.DataFrame(meta_list, index=df.columns)
datasets['GSE234116'] = {'expr': df, 'meta': metadata, 'format': 'counts'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# ============ PART 2: common genes (5 数据集交集) ============
print("\n" + "=" * 60)
print("PART 2: Common Ensembl genes")
print("=" * 60)
ensembl_sets = []
for name in ['GSE168106', 'GSE112380', 'GSE160334', 'GSE139512', 'GSE234116']:
    genes = {g for g in datasets[name]['expr'].index if str(g).startswith('ENSSSCG')}
    ensembl_sets.append(genes)
    print(f"  {name}: {len(genes)} Ensembl genes")
common_genes = sorted(set.intersection(*ensembl_sets))
print(f"\n  Common Ensembl genes across 5 datasets: {len(common_genes)}")

# ============ PART 3: 统一矩阵 ============
print("\n" + "=" * 60)
print("PART 3: Unified matrix")
print("=" * 60)
all_expr_parts, all_meta_parts = [], []
for name, d in datasets.items():
    expr, meta = d['expr'], d['meta']
    avail = [g for g in common_genes if g in expr.index]
    sub = expr.loc[avail]
    sub = sub.apply(pd.to_numeric, errors='coerce').fillna(0)
    sub_log = np.log1p(sub) if d['format'] == 'counts' else np.log1p(sub.clip(lower=0))
    all_expr_parts.append(sub_log.T)
    all_meta_parts.append(meta)
expr_matrix = pd.concat(all_expr_parts, axis=0)
meta_matrix = pd.concat(all_meta_parts, axis=0)
print(f"  Unified matrix: {expr_matrix.shape[0]} cells x {expr_matrix.shape[1]} genes")
for ds in meta_matrix['dataset'].unique():
    print(f"    {ds}: {(meta_matrix['dataset'] == ds).sum()} cells")

# ============ PART 4: QC ============
cell_totals = expr_matrix.sum(axis=1)
keep = cell_totals > cell_totals.quantile(0.01)
expr_f = expr_matrix.loc[keep]
meta_f = meta_matrix.loc[keep]
print(f"\n  After QC: {expr_f.shape[0]} cells (removed {sum(~keep)})")

# ============ PART 5: 标准化 + PCA ============
X = expr_f.values.astype(np.float32)
X_scaled = StandardScaler().fit_transform(X)
pca = PCA(n_components=50, random_state=42)
X_pca = pca.fit_transform(X_scaled)
print(f"\n  PCA: {X_pca.shape[1]} PCs, top5 var={pca.explained_variance_ratio_[:5].round(4)}")

# 保存未校正 PCA
pca_df = pd.DataFrame(X_pca, index=meta_f.index, columns=[f'PC{i+1}' for i in range(50)])
pca_df.to_csv(f'{OUT_DIR}/m2_pca.csv')
print("  saved m2_pca.csv")

# ============ PART 6: Harmony (batch=dataset) ============
print("\n  Running Harmony (batch=dataset) ...")
ho = harmonypy.run_harmony(X_pca.T, meta_f, vars_use='dataset', verbose=False)
Z = np.asarray(ho.Z_corr, dtype=float)  # cells x PCs
print(f"  Harmony corrected: {Z.shape}")
harm_df = pd.DataFrame(Z, index=meta_f.index, columns=[f'HPC{i+1}' for i in range(Z.shape[1])])
harm_df.to_csv(f'{OUT_DIR}/m2_harmony_corrected.csv')
print("  saved m2_harmony_corrected.csv")

# ============ PART 7: UMAP (raw centering + Harmony) ============
red = umap.UMAP(n_components=2, random_state=42, n_neighbors=30, min_dist=0.3)
# raw 对照：对未校正 PCA 做 dataset 中心化（与 m2_integration 一致）
X_c = X_pca.copy()
for ds in meta_f['dataset'].unique():
    mask = (meta_f['dataset'] == ds).values
    if mask.sum() > 1:
        X_c[mask] = X_c[mask] - X_c[mask].mean(axis=0)
Xu = red.fit_transform(X_c)
Zu = red.fit_transform(Z)
umap_df = pd.DataFrame(np.column_stack([Xu, Zu]), index=meta_f.index,
                       columns=['UMAP1_raw', 'UMAP2_raw', 'UMAP1_harmony', 'UMAP2_harmony'])
umap_df.to_csv(f'{OUT_DIR}/m2_harmony_umap.csv')
print("  saved m2_harmony_umap.csv")

# ============ PART 8: stage_group + metadata ============
def consolidate_stage(row):
    s = str(row['stage'])
    if s in ['Oocyte', 'GV']: return 'Oocyte'
    if s in ['1C', 'MII']: return '1C/MII'
    if s == '2C': return '2C'
    if s == '4C': return '4C'
    if s == '8C': return '8C'
    if s in ['E0', 'E1']: return 'Zygote/E0-E1'
    if s == 'E2': return 'E2 (2C)'
    if s == 'E3': return 'E3 (4-8C)'
    if s in ['E4', 'E5']: return 'Morula (E4-5)'
    if s in ['E6', 'E7', 'E8']: return 'Early blastocyst (E6-8)'
    if s in ['E9', 'E10']: return 'Hatched blastocyst (E9-10)'
    if s in ['E11', 'E12', 'E13', 'E14']: return 'Post-impl (E11-14)'
    if s in ['p10', 'p54', 'p64']: return 'pgEpiSC (stem cell)'
    if 'Mor' in s or 'mor' in s: return 'Morula'
    if 'blasto' in s.lower() or 'Blast' in s: return 'Blastocyst'
    if s == 'preimplantation': return 'Preimplant (D5-11)'
    return f'Other ({s[:8]})'

meta_f['stage_group'] = meta_f.apply(consolidate_stage, axis=1)
meta_f['UMAP1'] = Xu[:, 0]
meta_f['UMAP2'] = Xu[:, 1]
meta_f.index.name = 'cell'
meta_f.to_csv(f'{OUT_DIR}/m2_metadata.csv')
print("  saved m2_metadata.csv")

print("\n=== 数据集分布 ===")
print(meta_f['dataset'].value_counts().to_string())
print("\n=== stage_group 分布 ===")
print(meta_f['stage_group'].value_counts().to_string())
print(f"\nDONE: 总细胞数 = {len(meta_f)}")
