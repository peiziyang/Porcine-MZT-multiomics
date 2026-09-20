"""
M5: 跨条件伪时序轨迹比较 — IVF vs PA vs in_vivo
=================================================
用 GSE168106 (体内 E0-E8) + GSE164812 (IVF/PA 1C-8C)：
PCA pseudotime 对齐 → per-gene loess 平滑 → 条件特异偏离基因。

LAVDC 护栏:
- #1 条件标签来自 GEO metadata
- #2 轨迹比较不报细胞级 p，报效应量
- #3 方法来自文献回顾
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.interpolate import UnivariateSpline
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw'
OUT_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'

print("="*60)
print("M5: Cross-condition Pseudotime Trajectory Comparison")
print("="*60)

# ============================================================
# 1. Load GSE168106 (in_vivo, E0-E8)
# ============================================================
print("\n[1/5] Loading GSE168106 (in vivo E0-E8)...")
counts_vivo = pd.read_csv(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_genes_counts.txt.gz',
                          sep='\t', index_col=0)
si = pd.read_excel(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_SampleInfo.xlsx')
si_valid = si[si['Vaild'] == 'Yes'].copy()
cell_to_day = dict(zip(si_valid['ID'], si_valid['Stage_II']))
cell_to_lane = dict(zip(si_valid['ID'], si_valid['Lane']))

# Select E0-E8 stages only (early embryo, matching IVF/PA stage range)
early_days = ['E0','E1','E2','E3','E4','E5','E6','E7','E8']
vivo_cells = [c for c in counts_vivo.columns
              if c in cell_to_day and pd.notna(cell_to_day[c])
              and cell_to_day[c] in early_days]
counts_vivo = counts_vivo[vivo_cells]

# Create metadata
vivo_meta = pd.DataFrame({
    'cell': vivo_cells,
    'day': [cell_to_day[c] for c in vivo_cells],
    'condition': 'in_vivo',
    'dataset': 'GSE168106'
}).set_index('cell')

print(f"  {len(vivo_cells)} cells, {counts_vivo.shape[0]} genes")

# ============================================================
# 2. Load GSE164812 (IVF/PA, 1C-8C)
# ============================================================
print("\n[2/5] Loading GSE164812 (IVF/PA 1C-8C)...")
fpkm_paivf = pd.read_csv(f'{DATA_DIR}/GSE164812/GSE164812_gene_FPKM_matrix.txt.gz',
                         sep='\t', index_col=0)
fpkm_paivf = fpkm_paivf.apply(pd.to_numeric, errors='coerce').fillna(0)

# Parse condition+stage
def parse_gse164812(col):
    parts = col.replace('  ',' ').split(' ')
    cond = parts[0]
    stage_raw = parts[1] if len(parts)>1 else ''
    if '1-cell' in stage_raw: stage = '1C'
    elif '2-cell' in stage_raw: stage = '2C'
    elif '4-cell' in stage_raw: stage = '4C'
    elif '8-cell' in stage_raw: stage = '8C'
    else: stage = 'unknown'
    return cond, stage

parsed = [parse_gse164812(c) for c in fpkm_paivf.columns]
paivf_meta = pd.DataFrame({
    'cell': fpkm_paivf.columns,
    'condition': [p[0] for p in parsed],
    'stage': [p[1] for p in parsed],
    'dataset': 'GSE164812'
}).set_index('cell')

print(f"  {fpkm_paivf.shape[1]} cells (IVF={sum(paivf_meta['condition']=='IVF')}, PA={sum(paivf_meta['condition']=='PA')})")

# ============================================================
# 3. Common genes, normalize, combine
# ============================================================
print("\n[3/5] Finding common genes and combining...")
common = sorted(set(counts_vivo.index) & set(fpkm_paivf.index))
common = [g for g in common if g.startswith('ENSSSCG')]
print(f"  Common Ensembl genes: {len(common)}")

# Normalize each dataset separately (per-gene z-score)
vivo_mat = np.log1p(counts_vivo.loc[common])
paivf_mat = np.log2(fpkm_paivf.loc[common] + 1)

# Z-score per gene within each dataset
vivo_z = (vivo_mat - vivo_mat.mean(axis=1).values.reshape(-1,1)) / (vivo_mat.std(axis=1).values.reshape(-1,1) + 1e-8)
paivf_z = (paivf_mat - paivf_mat.mean(axis=1).values.reshape(-1,1)) / (paivf_mat.std(axis=1).values.reshape(-1,1) + 1e-8)

# Combine: rows = genes, cols = cells
combined = pd.concat([pd.DataFrame(vivo_z, index=common, columns=vivo_cells),
                       pd.DataFrame(paivf_z, index=common, columns=fpkm_paivf.columns)], axis=1)
combined_meta = pd.concat([vivo_meta, paivf_meta])
print(f"  Combined matrix: {combined.shape[1]} cells x {combined.shape[0]} genes")

# ============================================================
# 4. PCA pseudotime
# ============================================================
print("\n[4/5] PCA pseudotime...")
X = combined.T.values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=20, random_state=42)
X_pca = pca.fit_transform(X_scaled)
print(f"  PC1 var: {pca.explained_variance_ratio_[0]:.3f}")

# PC1 as pseudotime
pseudotime = X_pca[:, 0]
# Normalize to [0,1]
pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())

combined_meta['pseudotime'] = pseudotime
combined_meta['PC1'] = X_pca[:, 0]
combined_meta['PC2'] = X_pca[:, 1]

# Check: does pseudotime correlate with developmental stage?
stage_order = {'E0':0,'E1':1,'E2':2,'E3':3,'E4':4,'E5':5,'E6':6,'E7':7,'E8':8,
               '1C':2,'2C':3,'4C':4,'8C':5}
vivo_meta_sub = combined_meta[combined_meta['dataset']=='GSE168106']
estimated_order = vivo_meta_sub['day'].map(stage_order)
r_spearman, p_spearman = spearmanr(estimated_order, vivo_meta_sub['pseudotime'])
print(f"  Pseudotime vs stage Spearman r = {r_spearman:.3f}, p = {p_spearman:.2e}")

# Also for IVF/PA
paivf_meta_sub = combined_meta[combined_meta['dataset']=='GSE164812']
estimated_order2 = paivf_meta_sub['stage'].map(stage_order)
r2, p2 = spearmanr(estimated_order2, paivf_meta_sub['pseudotime'])
print(f"  IVF/PA: pseudotime vs stage Spearman r = {r2:.3f}, p = {p2:.2e}")

# ============================================================
# 5. Per-gene trajectory fitting + condition divergence
# ============================================================
print("\n[5/5] Per-gene trajectory comparison...")

# For each condition, fit a loess-like spline along pseudotime
conditions = ['in_vivo', 'IVF', 'PA']
n_pt = 100  # evaluation points
pt_grid = np.linspace(0, 1, n_pt)

# Get per-condition expression matrices
cond_mats = {}
cond_pts = {}
for cond in conditions:
    cells = combined_meta[combined_meta['condition'] == cond].index
    if len(cells) >= 5:
        cond_mats[cond] = combined[cells]
        cond_pts[cond] = combined_meta.loc[cells, 'pseudotime'].values

n_genes = len(common)

# Sample genes for analysis (run on all but report top)
# For speed, use every 10th gene for full trajectory analysis
# and compute divergence for all genes using simplified metric
print("  Computing gene-level divergence scores...")

# Metrics: per-gene correlation with pseudotime, per-condition
divergence_scores = {}
for cond in conditions:
    if cond not in cond_mats: continue
    mat = cond_mats[cond]
    pt = cond_pts[cond]
    # For each gene, compute spline fit to pseudotime
    corrs = []
    rms_curves = []
    for i, gene in enumerate(mat.index):
        vals = mat.iloc[i].values
        # Simple linear correlation with pseudotime
        if len(set(vals)) > 1 and len(set(pt)) > 1:
            r, _ = spearmanr(vals, pt)
            corrs.append(abs(r))
        else:
            corrs.append(0)
        # Spline for trajectory curve
        try:
            spl = UnivariateSpline(pt, vals, s=len(pt)*0.5, k=3)
            curve = spl(pt_grid)
            rms_curves.append(curve)
        except:
            rms_curves.append(np.zeros(n_pt))
    divergence_scores[cond] = {
        'corr_pseudotime': np.array(corrs),
        'curves': np.array(rms_curves)
    }

# Compute pairwise trajectory divergence between conditions
# Euclidean distance between smoothed curves
gene_divergence = {}
for pair in [('IVF', 'PA'), ('IVF', 'in_vivo'), ('PA', 'in_vivo')]:
    c1, c2 = pair
    if c1 in divergence_scores and c2 in divergence_scores:
        curves1 = divergence_scores[c1]['curves']
        curves2 = divergence_scores[c2]['curves']
        # Per-gene curve distance (L2 norm)
        dists = np.sqrt(np.sum((curves1 - curves2)**2, axis=1))
        gene_divergence[pair] = dists

# Find top condition-divergent genes
top_genes = {}
for pair, dists in gene_divergence.items():
    top_idx = np.argsort(dists)[-30:]  # top 30
    top_genes[pair] = [(common[i], dists[i]) for i in reversed(top_idx)]

print(f"\n  Top divergent genes (IVF vs PA):")
for g, d in top_genes[('IVF', 'PA')][:10]:
    print(f"    {g}: divergence={d:.4f}")

print(f"\n  Top divergent genes (IVF vs in_vivo):")
for g, d in top_genes[('IVF', 'in_vivo')][:10]:
    print(f"    {g}: divergence={d:.4f}")

print(f"\n  Top divergent genes (PA vs in_vivo):")
for g, d in top_genes[('PA', 'in_vivo')][:10]:
    print(f"    {g}: divergence={d:.4f}")

# ============================================================
# 6. Visualization: Fig.2 panels
# ============================================================
print("\n[6/6] Generating Fig.2 panels...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# 6a: PC1 vs PC2, colored by condition
ax = axes[0, 0]
colors_cond = {'in_vivo': '#2ca02c', 'IVF': '#d62728', 'PA': '#1f77b4'}
for cond, c in colors_cond.items():
    mask = combined_meta['condition'] == cond
    if mask.sum() > 0:
        ax.scatter(combined_meta.loc[mask, 'PC1'].values,
                  combined_meta.loc[mask, 'PC2'].values,
                  c=c, s=25, alpha=0.7, label=f'{cond} ({mask.sum()})',
                  edgecolors='white', linewidth=0.3, rasterized=True)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%}) — Pseudotime')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
ax.set_title('A. Common Embedding: IVF, PA, and in vivo\n(Pseudotime = PC1)')
ax.legend(fontsize=8)

# 6b: Pseudotime density per condition
ax = axes[0, 1]
for cond, c in colors_cond.items():
    pts = combined_meta[combined_meta['condition'] == cond]['pseudotime']
    if len(pts) > 0:
        ax.hist(pts, bins=20, alpha=0.5, color=c, label=cond)
ax.set_xlabel('Pseudotime')
ax.set_ylabel('Cell Count')
ax.set_title('B. Pseudotime Distribution by Condition')
ax.legend(fontsize=8)

# 6c: Spearman r (correlation with pseudotime) for top genes
ax = axes[0, 2]
genes_to_plot = [g for g, _ in top_genes[('IVF', 'PA')][:5]]
for gene in genes_to_plot:
    gene_idx = list(common).index(gene)
    for cond, c in colors_cond.items():
        if cond in cond_mats:
            cells_c = [c for c in cond_mats[cond].columns if gene in cond_mats[cond].index]
            if len(cells_c) > 0:
                vals = cond_mats[cond].loc[gene].values
                pts = cond_pts[cond]
                ax.scatter(pts, vals, c=c, alpha=0.4, s=15)
ax.set_xlabel('Pseudotime')
ax.set_ylabel('Expression (z-score)')
ax.set_title('C. Top Divergent Genes (IVF vs PA)')
ax.legend(fontsize=7)

# 6d: Example trajectory (one gene, fitted splines)
ax = axes[1, 0]
example_gene = top_genes[('IVF', 'PA')][0][0]
gene_idx = list(common).index(example_gene)
for cond, c in colors_cond.items():
    if cond in cond_mats:
        pts = cond_pts[cond]
        vals = cond_mats[cond].loc[example_gene].values
        ax.scatter(pts, vals, c=c, alpha=0.5, s=12, label=f'{cond}')
        try:
            spl = UnivariateSpline(pts, vals, s=len(pts)*0.5, k=3)
            ax.plot(pt_grid, spl(pt_grid), color=c, linewidth=2, alpha=0.8)
        except:
            pass
ax.set_xlabel('Pseudotime')
ax.set_ylabel('Expression (z-score)')
ax.set_title(f'D. Trajectory: {example_gene[:20]}')
ax.legend(fontsize=8)

# 6e: Divergence heatmap (top genes x condition pairs)
ax = axes[1, 1]
div_df = pd.DataFrame({
    'IVF vs PA': gene_divergence[('IVF', 'PA')],
    'IVF vs in_vivo': gene_divergence[('IVF', 'in_vivo')],
    'PA vs in_vivo': gene_divergence[('PA', 'in_vivo')]
}, index=common)
top_overall = div_df.sum(axis=1).nlargest(30).index
div_top = div_df.loc[top_overall]
# Normalize per row for heatmap
div_norm = div_top.apply(lambda x: (x-x.min())/(x.max()-x.min() if x.max()>x.min() else 1), axis=1)
im = ax.imshow(div_norm.values, aspect='auto', cmap='RdBu_r')
ax.set_xticks(range(3))
ax.set_xticklabels(div_norm.columns, fontsize=8, rotation=30)
ax.set_yticks(range(30))
ax.set_yticklabels([g[:15]+'...' for g in div_norm.index], fontsize=5)
ax.set_title('E. Gene-level Trajectory Divergence\n(Top 30 Overall)')
plt.colorbar(im, ax=ax, shrink=0.7, label='Normalized divergence')

# 6f: Summary statistics
ax = axes[1, 2]
ax.axis('off')
summary = [
    "Cross-condition Trajectory Comparison",
    "======================================",
    f"Genes analyzed: {n_genes}",
    f"Pseudotime: PCA-PC1 (var={pca.explained_variance_ratio_[0]:.1%})",
    f"Spearman r (pseudotime vs stage):",
    f"  in_vivo:  r={r_spearman:.3f}",
    f"  IVF/PA:   r={r2:.3f}",
    "",
    "Top Divergence (IVF vs PA):",
    f"  1. {top_genes[('IVF','PA')][0][0][:20]}",
    f"  2. {top_genes[('IVF','PA')][1][0][:20]}",
    f"  3. {top_genes[('IVF','PA')][2][0][:20]}",
    "",
    "Top Divergence (PA vs in_vivo):",
    f"  1. {top_genes[('PA','in_vivo')][0][0][:20]}",
    "",
    "Method: PCA embedding + per-gene",
    "spline fitting + L2 curve distance",
    "[LAVDC #2: effect sizes reported,",
    " cell-level p-values avoided]"
]
for i, line in enumerate(summary):
    ax.text(0.02, 0.98-i*0.052, line, fontsize=8, family='monospace',
           transform=ax.transAxes, verticalalignment='top')

plt.suptitle('Fig.2: Cross-condition Developmental Trajectory Comparison\n(IVF vs PA vs in_vivo, GSE164812 + GSE168106)',
            fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(f'{OUT_DIR}/m5_trajectory_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: m5_trajectory_comparison.png")

# ============================================================
# 7. Save results
# ============================================================
combined_meta.to_csv(f'{OUT_DIR}/m5_trajectory_metadata.csv')

pd.DataFrame({'gene': common,
              'IVF_vs_PA_divergence': gene_divergence[('IVF', 'PA')],
              'IVF_vs_invivo_divergence': gene_divergence[('IVF', 'in_vivo')],
              'PA_vs_invivo_divergence': gene_divergence[('PA', 'in_vivo')]
}).to_csv(f'{OUT_DIR}/m5_gene_divergence.csv', index=False)

print(f"\nDone! Outputs: {OUT_DIR}/")
print(f"  m5_trajectory_comparison.png  — Fig.2 complete")
print(f"  m5_trajectory_metadata.csv")
print(f"  m5_gene_divergence.csv")
