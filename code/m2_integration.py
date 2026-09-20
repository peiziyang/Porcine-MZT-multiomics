"""
M2: 猪着床前胚胎单细胞参考图谱 — 数据整合
==============================================
读取 7 个公开 scRNA-seq 数据集，统一 Ensembl 基因 ID，
质控、归一化、批次校正、降维可视化。

所有数据均来自已发表公开数据集（非自产）。
"""
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw'
OUT_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'

# ============================================================
# PART 1: Load all datasets
# ============================================================
print("=" * 60)
print("PART 1: Loading datasets")
print("=" * 60)

datasets = {}

# --- GSE168106: Han 2022 Cell Res, E0-E14 ---
print("\n[1/7] GSE168106 (Han 2022) ...")
df = pd.read_csv(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_genes_counts.txt.gz',
                 sep='\t', index_col=0)
si = pd.read_excel(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_SampleInfo.xlsx')
# Filter valid cells
si_valid = si[si['Vaild'] == 'Yes'].copy()
valid_barcodes = set(si_valid['ID'].values)
df_cols = [c for c in df.columns if c in valid_barcodes]
df = df[df_cols]
# Map metadata
cell_to_stage = dict(zip(si_valid['ID'], si_valid['Stage_II']))
cell_to_lane = dict(zip(si_valid['ID'], si_valid['Lane']))
metadata_168106 = pd.DataFrame({
    'cell': df.columns,
    'stage': [cell_to_stage.get(c, 'Unknown') for c in df.columns],
    'donor': [cell_to_lane.get(c, 'Unknown') for c in df.columns],
    'dataset': 'GSE168106',
    'source': 'in_vivo',
    'condition': 'in_vivo',
    'study': 'Han2022'
}).set_index('cell')
datasets['GSE168106'] = {'expr': df, 'meta': metadata_168106, 'format': 'counts'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE112380: Ramos-Ibeas 2019 Nat Commun, D5-D11 ---
print("\n[2/7] GSE112380 (Ramos-Ibeas 2019) ...")
df = pd.read_csv(f'{DATA_DIR}/GSE112380/GSE112380_Pig_RawCounts.txt.gz',
                 sep='\t', index_col=0)
# Parse stage from the paper: cells from D5-D11 embryos
# Column names are sequencing IDs, need to infer stage from known metadata
# For now, assign "preimplantation" and refine later
metadata_112380 = pd.DataFrame({
    'cell': df.columns,
    'stage': ['preimplantation'] * df.shape[1],
    'donor': ['pooled'] * df.shape[1],  # 28 embryos pooled
    'dataset': 'GSE112380',
    'source': 'in_vivo',
    'condition': 'in_vivo',
    'study': 'RamosIbeas2019'
}).set_index('cell')
# Note: Ramos-Ibeas 2019 has sex info (15M, 13F), could add later
datasets['GSE112380'] = {'expr': df, 'meta': metadata_112380, 'format': 'counts'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE164812: Du 2021 Sci Rep, IVF/PA 1-8 cell ---
print("\n[3/7] GSE164812 (Du 2021) ...")
df = pd.read_csv(f'{DATA_DIR}/GSE164812/GSE164812_gene_FPKM_matrix.txt.gz',
                 sep='\t', index_col=0)
# Parse metadata from column names: "IVF 1-cell-1", "PA 2-cell-1-1"
meta_list = []
for col in df.columns:
    parts = col.replace('  ', ' ').split(' ')
    if len(parts) >= 2:
        cond = parts[0]  # IVF or PA
        stage_raw = parts[1]  # 1-cell, 2-cell, 4-cell, 8-cell
        rep = parts[2] if len(parts) > 2 else '1'
    else:
        cond, stage_raw, rep = 'Unknown', col, '1'
    meta_list.append({
        'condition': cond,
        'stage': stage_raw.replace('-cell', 'C'),
        'stage_clean': stage_raw,
        'source': cond,
        'dataset': 'GSE164812',
        'study': 'Du2021'
    })
metadata_164812 = pd.DataFrame(meta_list, index=df.columns)
# Assign donor based on condition+replicate grouping (proxy)
metadata_164812['donor'] = metadata_164812['condition'] + '_' + metadata_164812.index
# Use clean stage label (e.g. '1C-1' -> '1C')
metadata_164812['stage'] = metadata_164812['stage'].str.replace(r'-\d.*', '', regex=True)
datasets['GSE164812'] = {'expr': df, 'meta': metadata_164812, 'format': 'fpkm'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE160334: Du 2021, GV/MII oocytes ---
print("\n[4/7] GSE160334 (Du 2021 GV/MII) ...")
df = pd.read_csv(f'{DATA_DIR}/GSE160334/GSE160334_gene_FPKM_matrix.txt.gz',
                 sep='\t', index_col=0)
meta_list = []
for col in df.columns:
    if col.startswith('GV'):
        stage, rep = 'GV', col[2:]
    elif col.startswith('MII'):
        stage, rep = 'MII', col[3:]
    else:
        stage, rep = 'Unknown', '1'
    meta_list.append({
        'stage': stage,
        'donor': f'donor_{rep}',
        'condition': 'in_vivo',
        'source': 'in_vivo',
        'dataset': 'GSE160334',
        'study': 'Du2021'
    })
metadata_160334 = pd.DataFrame(meta_list, index=df.columns)
datasets['GSE160334'] = {'expr': df, 'meta': metadata_160334, 'format': 'fpkm'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE139512: Kong 2019 FASEB J, oocyte-blastocyst ---
print("\n[5/7] GSE139512 (Kong 2019) ...")
df = pd.read_csv(f'{DATA_DIR}/GSE139512/GSE139512_fpkm.txt.gz',
                 sep='\t', index_col=0)
meta_list = []
for col in df.columns:
    col_lower = col.lower()
    if 'oo' in col_lower or 'oocyte' in col_lower:
        stage = 'Oocyte'
    elif '1-cell' in col_lower or '1_cell' in col_lower:
        stage = '1C'
    elif '2-cell' in col_lower or '2_cell' in col_lower:
        stage = '2C'
    elif '4-cell' in col_lower or '4_cell' in col_lower:
        stage = '4C'
    elif '8-cell' in col_lower or '8_cell' in col_lower:
        stage = '8C'
    elif 'mor' in col_lower:
        stage = 'Morula'
    elif 'blasto' in col_lower or 'bl' in col_lower:
        stage = 'Blastocyst'
    else:
        stage = 'Unknown'
    meta_list.append({
        'stage': stage,
        'donor': f'donor_{col.split("_")[0] if "_" in col else "pool"}',
        'condition': 'in_vivo',
        'source': 'in_vivo',
        'dataset': 'GSE139512',
        'study': 'Kong2019'
    })
metadata_139512 = pd.DataFrame(meta_list, index=df.columns)
# Clean up Unknown labels with a more aggressive parse
metadata_139512['stage'] = metadata_139512['stage'].astype(str)
datasets['GSE139512'] = {'expr': df, 'meta': metadata_139512, 'format': 'fpkm'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE234116: Yuan 2023 Cell Mol Life Sci, GV oocytes (multi-omics RNA) ---
print("\n[6/7] GSE234116 (Yuan 2023 multi-omics RNA) ...")
df = pd.read_csv(f'{DATA_DIR}/GSE234116/GSE234116_raw_counts.txt.gz',
                 sep='\t', index_col=0)
meta_list = []
for col in df.columns:
    # AF{animal}_{oocyte_id}
    parts = col.split('_')
    animal = parts[0] if len(parts) > 0 else 'unknown'
    meta_list.append({
        'stage': 'GV',
        'donor': animal,
        'condition': 'in_vivo',
        'source': 'in_vivo',
        'dataset': 'GSE234116',
        'study': 'Yuan2023'
    })
metadata_234116 = pd.DataFrame(meta_list, index=df.columns)
datasets['GSE234116'] = {'expr': df, 'meta': metadata_234116, 'format': 'counts'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes")

# --- GSE242553: pig blastoid ---
print("\n[7/7] GSE242553 (pig blastoid) ...")
df = pd.read_csv(f'{DATA_DIR}/GSE242553/GSE242553_Transcript_fpkm.txt.gz',
                 sep='\t', index_col=0)
# GSE242553 uses RefSeq IDs (XM_/NM_), not Ensembl. Map or skip for now.
# For M2, include with a note that genes are RefSeq and won't overlap well.
metadata_242553 = pd.DataFrame({
    'stage': ['Blastoid'] * df.shape[1],
    'donor': ['blastoid_1'] * df.shape[1],
    'condition': ['in_vitro'] * df.shape[1],
    'source': ['blastoid_ESC'] * df.shape[1],
    'dataset': ['GSE242553'] * df.shape[1],
    'study': ['Blastoid']
}, index=df.columns)
datasets['GSE242553'] = {'expr': df, 'meta': metadata_242553, 'format': 'fpkm'}
print(f"  -> {df.shape[1]} cells, {df.shape[0]} genes (RefSeq IDs - limited overlap)")

# ============================================================
# PART 2: Harmonize gene IDs (Ensembl)
# ============================================================
print("\n" + "=" * 60)
print("PART 2: Harmonizing gene IDs")
print("=" * 60)

# All 6 main datasets use ENSSSCG... IDs. Find common genes.
ensembl_sets = []
for name in ['GSE168106', 'GSE112380', 'GSE164812', 'GSE160334',
             'GSE139512', 'GSE234116']:
    genes = set(datasets[name]['expr'].index)
    genes = {g for g in genes if g.startswith('ENSSSCG')}
    ensembl_sets.append(genes)
    print(f"  {name}: {len(genes)} Ensembl genes")

common_genes = set.intersection(*ensembl_sets)
print(f"\n  Common Ensembl genes across 6 datasets: {len(common_genes)}")

# For GSE242553, keep separate (RefSeq) - won't overlap well
# but we'll include it for completeness in the final plots

# ============================================================
# PART 3: Build unified expression matrix
# ============================================================
print("\n" + "=" * 60)
print("PART 3: Building unified matrix")
print("=" * 60)

gene_list = sorted(common_genes)
all_expr_parts = []
all_meta_parts = []

for name, d in datasets.items():
    if name == 'GSE242553':
        continue  # Skip for now (RefSeq IDs)
    
    expr = d['expr']
    meta = d['meta']
    
    # Subset to common genes
    avail = [g for g in gene_list if g in expr.index]
    sub = expr.loc[avail]
    
    # Convert to numeric (handle any string values like 'NA')
    sub = sub.apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # Log-transform
    if d['format'] == 'counts':
        sub_log = np.log1p(sub)
    else:  # fpkm
        sub_log = np.log1p(sub.clip(lower=0))
    
    sub_log.index.name = 'gene'
    sub_log.columns.name = 'cell'
    all_expr_parts.append(sub_log.T)  # cells x genes
    all_meta_parts.append(meta)

# Concatenate
expr_matrix = pd.concat(all_expr_parts, axis=0)
meta_matrix = pd.concat(all_meta_parts, axis=0)

n_cells = expr_matrix.shape[0]
n_genes = expr_matrix.shape[1]
print(f"  Unified matrix: {n_cells} cells x {n_genes} genes")
print(f"  Dataset breakdown:")
for ds in meta_matrix['dataset'].unique():
    n = (meta_matrix['dataset'] == ds).sum()
    print(f"    {ds}: {n} cells")

# ============================================================
# PART 4: QC filtering
# ============================================================
print("\n" + "=" * 60)
print("PART 4: QC filtering")
print("=" * 60)

# Simple QC: filter cells with very low total expression
cell_totals = expr_matrix.sum(axis=1)
q01 = cell_totals.quantile(0.01)
keep_cells = cell_totals > q01
expr_filtered = expr_matrix.loc[keep_cells]
meta_filtered = meta_matrix.loc[keep_cells]
print(f"  After QC: {expr_filtered.shape[0]} cells (removed {sum(~keep_cells)} low-quality)")

# ============================================================
# PART 5: Normalize & scale
# ============================================================
print("\n" + "=" * 60)
print("PART 5: Normalization")
print("=" * 60)

X = expr_filtered.values.astype(np.float32)

# Center each gene to mean=0 (per-gene standardization for PCA)
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
print(f"  Data shape after scaling: {X_scaled.shape}")

# ============================================================
# PART 6: PCA
# ============================================================
print("\n" + "=" * 60)
print("PART 6: PCA")
print("=" * 60)

pca = PCA(n_components=50, random_state=42)
X_pca = pca.fit_transform(X_scaled)
print(f"  PCA: {X_pca.shape[1]} components")
print(f"  Explained variance (top 5): {pca.explained_variance_ratio_[:5].round(4)}")

# ============================================================
# PART 7: Simple batch correction (center each dataset in PCA)
# ============================================================
print("\n" + "=" * 60)
print("PART 7: Dataset-level centering (batch correction)")
print("=" * 60)

X_corrected = X_pca.copy()
for ds in meta_filtered['dataset'].unique():
    mask = (meta_filtered['dataset'] == ds).values
    if mask.sum() > 1:
        center = X_pca[mask].mean(axis=0)
        X_corrected[mask] = X_pca[mask] - center
        print(f"  {ds}: centered {mask.sum()} cells")

print(f"  Post-correction std per component: {X_corrected.std(axis=0)[:5].round(4)}")

# ============================================================
# PART 8: UMAP
# ============================================================
print("\n" + "=" * 60)
print("PART 8: UMAP")
print("=" * 60)

reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=30, min_dist=0.3)
X_umap = reducer.fit_transform(X_corrected)
print(f"  UMAP shape: {X_umap.shape}")

# ============================================================
# PART 9: Visualization
# ============================================================
print("\n" + "=" * 60)
print("PART 9: Generating plots")
print("=" * 60)

import os
os.makedirs(f'{OUT_DIR}', exist_ok=True)

def make_plot(meta, coords, color_col, title, filename, cmap='tab20'):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    categories = meta[color_col].astype(str)
    unique_cats = sorted(categories.unique())
    
    colors = plt.cm.get_cmap(cmap, len(unique_cats))
    for i, cat in enumerate(unique_cats):
        mask = categories == cat
        ax.scatter(coords[mask, 0], coords[mask, 1],
                  s=3, alpha=0.7, label=cat,
                  color=colors(i), rasterized=True)
    
    ax.legend(markerscale=5, fontsize=7, loc='upper left',
             bbox_to_anchor=(1.02, 1), frameon=False)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout()
    fig.savefig(f'{OUT_DIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")

# Plot by dataset source
make_plot(meta_filtered, X_umap, 'dataset',
          'Pig Preimplantation Embryo scRNA-seq Atlas\n(Colored by Dataset)',
          'm2_umap_by_dataset.png')

# Plot by developmental stage (cleaner consolidated labels)
def consolidate_stage(row):
    s = str(row['stage'])
    if s in ['Oocyte', 'GV']: return 'Oocyte'
    if s in ['1C', 'MII']: return '1C/MII'
    if s == '2C': return '2C'
    if s == '4C': return '4C'
    if s == '8C': return '8C'
    if s in ['E0','E1']: return 'Zygote/E0-E1'
    if s == 'E2': return 'E2 (2C)'
    if s == 'E3': return 'E3 (4-8C)'
    if s in ['E4','E5']: return 'Morula (E4-5)'
    if s in ['E6','E7','E8']: return 'Early blastocyst (E6-8)'
    if s in ['E9','E10']: return 'Hatched blastocyst (E9-10)'
    if s in ['E11','E12','E13','E14']: return 'Post-impl (E11-14)'
    if s in ['p10','p54','p64']: return 'pgEpiSC (stem cell)'
    if 'Mor' in s or 'mor' in s: return 'Morula'
    if 'blasto' in s.lower() or 'Blast' in s: return 'Blastocyst'
    if s == 'preimplantation': return 'Preimplant (D5-11)'
    return f'Other ({s[:8]})'

meta_filtered['stage_group'] = meta_filtered.apply(consolidate_stage, axis=1)
make_plot(meta_filtered, X_umap, 'stage_group',
          'Pig Preimplantation Embryo scRNA-seq Atlas\n(Colored by Consolidated Stage)',
          'm2_umap_by_stage_group.png')

# Plot by source (IVF/PA/in_vivo/blastoid)
make_plot(meta_filtered, X_umap, 'source',
          'Pig Preimplantation Embryo scRNA-seq Atlas\n(Colored by Source)',
          'm2_umap_by_source.png')

# ============================================================
# PART 10: Summary statistics
# ============================================================
print("\n" + "=" * 60)
print("PART 10: Summary")
print("=" * 60)

summary = meta_filtered.groupby(['dataset', 'stage']).size().reset_index(name='n_cells')
print("\nCell counts per dataset and stage:")
print(summary.to_string())

# Save metadata and coordinates for downstream use
meta_filtered['UMAP1'] = X_umap[:, 0]
meta_filtered['UMAP2'] = X_umap[:, 1]
meta_filtered.to_csv(f'{OUT_DIR}/m2_metadata.csv')

# Save PCA coordinates
pca_df = pd.DataFrame(X_pca, index=meta_filtered.index,
                      columns=[f'PC{i+1}' for i in range(X_pca.shape[1])])
pca_df.to_csv(f'{OUT_DIR}/m2_pca.csv')

print(f"\nDone! All outputs saved to {OUT_DIR}/")
print(f"  - m2_metadata.csv: cell metadata + UMAP coordinates")
print(f"  - m2_pca.csv: PCA coordinates")
print(f"  - m2_umap_by_dataset.png")
print(f"  - m2_umap_by_stage.png")
print(f"  - m2_umap_by_source.png")
