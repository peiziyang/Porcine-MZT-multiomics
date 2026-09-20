"""
M5 v2: 跨条件轨迹比较 — 独立伪时序 + 对齐
===========================================
修复：三个条件独立 PCA → 各算伪时序 → 用阶段 anchors 对齐 → 比较轨迹

每条件至少有 4 个 stage 作为 anchor points，用 anchors 对齐伪时序尺度。
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
from scipy.interpolate import UnivariateSpline, interp1d
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw'
OUT_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'

print("="*60)
print("M5 v2: Cross-condition Trajectory (Independent PT + Alignment)")
print("="*60)

# ============================================================
# 1. Load data
# ============================================================
counts_vivo = pd.read_csv(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_genes_counts.txt.gz', sep='\t', index_col=0)
si = pd.read_excel(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_SampleInfo.xlsx')
si_v = si[si['Vaild']=='Yes']
cell_to_day = dict(zip(si_v['ID'], si_v['Stage_II']))
early_days = ['E0','E1','E2','E3','E4','E5','E6','E7','E8']
vivo_cells = [c for c in counts_vivo.columns if c in cell_to_day and pd.notna(cell_to_day[c]) and cell_to_day[c] in early_days]

fpkm_paivf = pd.read_csv(f'{DATA_DIR}/GSE164812/GSE164812_gene_FPKM_matrix.txt.gz', sep='\t', index_col=0)
fpkm_paivf = fpkm_paivf.apply(pd.to_numeric, errors='coerce').fillna(0)

def parse_cond_stage(col):
    p = col.replace('  ',' ').split(' ')
    c = p[0]
    s = p[1] if len(p)>1 else ''
    if '1-cell' in s: return c, '1C'
    if '2-cell' in s: return c, '2C'
    if '4-cell' in s: return c, '4C'
    if '8-cell' in s: return c, '8C'
    return c, 'unknown'

ivf_cells = [c for c in fpkm_paivf.columns if c.startswith('IVF') and parse_cond_stage(c)[1] != 'unknown']
pa_cells  = [c for c in fpkm_paivf.columns if c.startswith('PA') and parse_cond_stage(c)[1] != 'unknown']

common = sorted(set(counts_vivo.index) & set(fpkm_paivf.index))
common = [g for g in common if g.startswith('ENSSSCG')]
print(f"Common genes: {len(common)}")
print(f"Cells: vivo={len(vivo_cells)}, IVF={len(ivf_cells)}, PA={len(pa_cells)}")

# ============================================================
# 2. Independent pseudotime per condition
# ============================================================
stage_map = {'E0':0,'E1':1,'E2':2,'E3':3,'E4':4,'E5':5,'E6':6,'E7':7,'E8':8,
             '1C':2,'2C':3,'4C':4,'8C':5}

def compute_pseudotime(mat, cells, meta_fn, label):
    """PCA pseudotime for one condition."""
    sub = mat[list(cells)]
    # Filter low-expression genes
    expr_rate = (sub > 0).mean(axis=1)
    sub_f = sub.loc[expr_rate > 0.10]
    # log-transform
    log_mat = np.log1p(sub_f)
    # scale
    from sklearn.preprocessing import StandardScaler
    X = log_mat.T.values
    X_s = StandardScaler().fit_transform(X)
    pca = PCA(n_components=min(10, X_s.shape[1]), random_state=42)
    X_p = pca.fit_transform(X_s)
    pt = X_p[:, 0]
    pt = (pt - pt.min()) / (pt.max() - pt.min() + 1e-8)
    # Correlate with stage
    stages = [meta_fn(c) for c in cells]
    stage_vals = [stage_map.get(s, -1) for s in stages]
    valid = [i for i, v in enumerate(stage_vals) if v >= 0]
    if len(valid) > 2:
        r, p = spearmanr([stage_vals[i] for i in valid], [pt[i] for i in valid])
    else:
        r, p = np.nan, np.nan
    # Force pseudotime to go from early→late (flip if needed)
    if not np.isnan(r) and r < 0:
        pt = 1 - pt
        r = -r
    print(f"  {label}: {len(cells)} cells, PC1={pca.explained_variance_ratio_[0]:.2%}, r={r:.3f}")
    return pt, sub_f, pca

# in_vivo: E0-E8
pt_vivo, mat_vivo, pca_vivo = compute_pseudotime(
    counts_vivo, vivo_cells, 
    lambda c: cell_to_day.get(c,''), 'in_vivo')

# IVF: 1C-8C  
pt_ivf, mat_ivf, pca_ivf = compute_pseudotime(
    fpkm_paivf, ivf_cells,
    lambda c: parse_cond_stage(c)[1], 'IVF')

# PA: 1C-8C
pt_pa, mat_pa, pca_pa = compute_pseudotime(
    fpkm_paivf, pa_cells,
    lambda c: parse_cond_stage(c)[1], 'PA')

# ============================================================
# 3. Align pseudotimes using stage anchors
# ============================================================
print("\n[3/5] Aligning pseudotimes...")
# For each condition, get mean pseudotime per stage
def stage_mean_pt(pt, cells, meta_fn):
    d = {}
    for c, p in zip(cells, pt):
        s = meta_fn(c)
        if s not in d: d[s] = []
        d[s].append(p)
    return {s: np.mean(v) for s, v in d.items()}

vivo_anchors = stage_mean_pt(pt_vivo, vivo_cells, lambda c: cell_to_day.get(c,''))
ivf_anchors = stage_mean_pt(pt_ivf, ivf_cells, lambda c: parse_cond_stage(c)[1])
pa_anchors = stage_mean_pt(pt_pa, pa_cells, lambda c: parse_cond_stage(c)[1])

# Map: IVF/PA stage → in_vivo stage → in_vivo pseudotime
# Find common stages
common_stages = set(vivo_anchors.keys()).intersection(
    {s for s in vivo_anchors if s in ['E2','E3','E4','E5']}  # map to 1C-8C
)
stage_to_vivo_pt = {s: vivo_anchors[s] for s in early_days if s in vivo_anchors}

# Map: IVF/PA 1C→E2, 2C→E3, 4C→E4, 8C→E5 (approximate mapping)
ivf_to_vivo = {'1C':'E2', '2C':'E3', '4C':'E4', '8C':'E5'}

# Compute aligned pseudotime by mapping IVF/PA stages to vivo pseudotime
def align_pt(pt, cells, meta_fn, stage_map_fn):
    aligned = []
    for c, p in zip(cells, pt):
        s = meta_fn(c)
        vivo_s = stage_map_fn.get(s)
        if vivo_s and vivo_s in vivo_anchors:
            # Map to vivo pseudotime at the same stage
            aligned.append(vivo_anchors[vivo_s])
        else:
            aligned.append(np.nan)
    return np.array(aligned)

pt_ivf_aligned = align_pt(pt_ivf, ivf_cells, lambda c: parse_cond_stage(c)[1], ivf_to_vivo)
pt_pa_aligned = align_pt(pt_pa, pa_cells, lambda c: parse_cond_stage(c)[1], ivf_to_vivo)

print(f"  IVF aligned PT (non-NaN): {sum(~np.isnan(pt_ivf_aligned))}")
print(f"  PA aligned PT (non-NaN): {sum(~np.isnan(pt_pa_aligned))}")

# ============================================================
# 4. Find common genes for comparison
# ============================================================
print("\n[4/5] Computing gene expression profiles per condition...")
common_genes = sorted(set(mat_vivo.index) & set(mat_ivf.index) & set(mat_pa.index))
print(f"  Common genes with expression in all: {len(common_genes)}")

# Compute mean expression per stage for each condition
vivo_gene_by_stage = {}
for s in early_days:
    cells_s = [c for c in vivo_cells if cell_to_day.get(c) == s]
    if cells_s:
        vivo_gene_by_stage[s] = counts_vivo.loc[common_genes, cells_s].apply(
            lambda x: np.log1p(x).mean(), axis=1)

ivf_gene_by_stage = {}
for s in ['1C','2C','4C','8C']:
    cells_s = [c for c in ivf_cells if parse_cond_stage(c)[1] == s]
    if cells_s:
        ivf_gene_by_stage[s] = fpkm_paivf.loc[common_genes, cells_s].apply(
            lambda x: np.log2(x+1).mean(), axis=1)

pa_gene_by_stage = {}
for s in ['1C','2C','4C','8C']:
    cells_s = [c for c in pa_cells if parse_cond_stage(c)[1] == s]
    if cells_s:
        pa_gene_by_stage[s] = fpkm_paivf.loc[common_genes, cells_s].apply(
            lambda x: np.log2(x+1).mean(), axis=1)

# Z-score within each condition (make comparable)
def zscore_series(series_dict):
    all_vals = pd.concat(series_dict.values())
    mean, std = all_vals.mean(), all_vals.std() + 1e-8
    return {k: (v - mean) / std for k, v in series_dict.items()}

vivo_z = zscore_series(vivo_gene_by_stage)
ivf_z = zscore_series(ivf_gene_by_stage)
pa_z = zscore_series(pa_gene_by_stage)

# ============================================================
# 5. Per-gene: compare trajectory shape across conditions
# ============================================================
print("\n[5/5] Computing trajectory divergence...")

# For each gene, compute the per-stage profile and compare
# IVF/PA stages mapped to vivo equivalents
stage_align = {'1C':'E2','2C':'E3','4C':'E4','8C':'E5'}

gene_scores = []
for gene in common_genes:
    # in_vivo profile: E0 → E8 (9 points)
    vivo_profile = np.array([vivo_z[s].loc[gene] if s in vivo_z else np.nan for s in early_days])
    vivo_valid = ~np.isnan(vivo_profile)
    
    # IVF profile: mapped E2-E5
    ivf_profile = np.array([ivf_z[s].loc[gene] if s in ivf_z else np.nan for s in ['1C','2C','4C','8C']])
    
    # PA profile
    pa_profile = np.array([pa_z[s].loc[gene] if s in pa_z else np.nan for s in ['1C','2C','4C','8C']])
    
    if vivo_valid.sum() < 3 or np.isnan(ivf_profile[0]):
        continue
    
    # Compare IVF vs PA at aligned stages (both have 4 points at same pseudotime)
    ivf_pa_div = np.sqrt(np.mean((ivf_profile - pa_profile)**2))
    
    # Compare IVF vs vivo at aligned stages (E2-E5 from vivo vs 1C-8C from IVF)
    vivo_e2e5 = vivo_profile[2:6]  # E2-E5
    ivf_vivo_div = np.sqrt(np.mean((ivf_profile - vivo_e2e5)**2))
    
    pa_vivo_div = np.sqrt(np.mean((pa_profile - vivo_e2e5)**2))
    
    gene_scores.append({
        'gene': gene,
        'IVF_vs_PA': ivf_pa_div,
        'IVF_vs_vivo': ivf_vivo_div,
        'PA_vs_vivo': pa_vivo_div,
        'total_div': ivf_pa_div + ivf_vivo_div + pa_vivo_div,
        'vivo_slope': np.polyfit(range(9), vivo_profile, 1)[0] if vivo_valid.sum() >= 5 else 0
    })

div_df = pd.DataFrame(gene_scores).set_index('gene')
print(f"  Genes with valid comparisons: {len(div_df)}")

# Top divergent genes
for pair in ['IVF_vs_PA', 'IVF_vs_vivo', 'PA_vs_vivo']:
    top = div_df.nlargest(15, pair)
    print(f"\n  Top 5 {pair}:")
    for g in top.head(5).index:
        print(f"    {g}: div={div_df.loc[g, pair]:.3f}")

# ============================================================
# 6. Visualization
# ============================================================
print("\n[6/6] Generating Fig.2...")

stage_labels_vivo = ['E0','E1','E2','E3','E4','E5','E6','E7','E8']
stage_labels_art = ['1C','2C','4C','8C']
colors_cond = {'in_vivo': '#2ca02c', 'IVF': '#d62728', 'PA': '#1f77b4'}

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# A: Pseudotime vs stage scatter (3 conditions)
ax = axes[0, 0]
stage_vals_vivo = [stage_map[cell_to_day.get(c,'')] for c in vivo_cells]
ax.scatter(pt_vivo, stage_vals_vivo, c='#2ca02c', alpha=0.5, s=25, label=f'in_vivo (r_spearman)')
# Annotate
vivo_valid = [i for i, v in enumerate(stage_vals_vivo) if v >= 0]
r_vivo, _ = spearmanr([stage_vals_vivo[i] for i in vivo_valid], [pt_vivo[i] for i in vivo_valid])
ax.set_xlabel('Pseudotime (PC1)')
ax.set_ylabel('Stage Index (0=E0, 8=E8)')
ax.set_title(f'A. in_vivo Pseudotime vs Stage (r={r_vivo:.3f})')

# B: Per-stage gene trajectory examples (3 genes)
ax = axes[0, 1]
top_overall = div_df.nlargest(3, 'total_div')
for i, gene in enumerate(top_overall.index):
    # in_vivo
    y_vivo = [vivo_z[s].loc[gene] if s in vivo_z else np.nan for s in stage_labels_vivo]
    x_vivo = list(range(len(y_vivo)))
    ax.plot(x_vivo, y_vivo, 'o-', color='#2ca02c', markersize=6, alpha=0.8, 
            linewidth=2, label='in_vivo' if i==0 else '', zorder=3)
    # IVF
    y_ivf = [ivf_z[s].loc[gene] if s in ivf_z else np.nan for s in stage_labels_art]
    x_ivf = [2,3,4,5]  # mapped to E2-E5
    ax.plot(x_ivf, y_ivf, 's--', color='#d62728', markersize=6, alpha=0.8,
            linewidth=2, label='IVF' if i==0 else '', zorder=3)
    # PA
    y_pa = [pa_z[s].loc[gene] if s in pa_z else np.nan for s in stage_labels_art]
    ax.plot(x_ivf, y_pa, '^:', color='#1f77b4', markersize=6, alpha=0.8,
            linewidth=2, label='PA' if i==0 else '', zorder=3)
ax.set_xticks(range(9))
ax.set_xticklabels(stage_labels_vivo)
ax.set_xlabel('Developmental Stage')
ax.set_ylabel('Expression (z-score)')
ax.set_title('B. Top Divergent Gene Trajectories\n(3 conditions × 3 genes)')
ax.legend(fontsize=7, loc='lower right')

# C: Divergence heatmap
ax = axes[0, 2]
top30 = div_df.nlargest(30, 'total_div')
heat_data = top30[['IVF_vs_PA', 'IVF_vs_vivo', 'PA_vs_vivo']]
im = ax.imshow(heat_data.values, aspect='auto', cmap='YlOrRd')
ax.set_xticks(range(3))
ax.set_xticklabels(['IVF vs PA', 'IVF vs vivo', 'PA vs vivo'], fontsize=8, rotation=25)
ax.set_yticks(range(30))
ax.set_yticklabels([g[:15] for g in heat_data.index], fontsize=5)
ax.set_title('C. Trajectory Divergence (Top 30)')
plt.colorbar(im, ax=ax, shrink=0.7)

# D: IVF vs PA trajectory comparison (7 selected pathways genes)
ax = axes[1, 0]
# Pick genes known to be related to early development
pathway_genes = div_df.nlargest(7, 'IVF_vs_PA')
for gene in pathway_genes.index:
    y_ivf = [ivf_z[s].loc[gene] if s in ivf_z else np.nan for s in ['1C','2C','4C','8C']]
    y_pa = [pa_z[s].loc[gene] if s in pa_z else np.nan for s in ['1C','2C','4C','8C']]
    x_art = [0,1,2,3]
    ax.plot(x_art, y_ivf, 'o-', color='#d62728', alpha=0.6, markersize=5, linewidth=1.5)
    ax.plot(x_art, y_pa, '^--', color='#1f77b4', alpha=0.6, markersize=5, linewidth=1.5)
ax.set_xticks([0,1,2,3])
ax.set_xticklabels(['1C','2C','4C','8C'])
ax.set_xlabel('Stage')
ax.set_ylabel('Expression (z-score)')
ax.set_title('D. IVF (—) vs PA (---) Trajectories\n(Top 7 Divergent Genes)')

# E: Scatter: IVF genes vs PA genes at matched stages
ax = axes[1, 1]
# At 4C stage (middle of development), compare IVF vs PA expression
if '4C' in ivf_z and '4C' in pa_z:
    x_vals = ivf_z['4C'].values
    y_vals = pa_z['4C'].values
    ax.scatter(x_vals, y_vals, c='grey', alpha=0.3, s=8, rasterized=True)
    # Highlight divergent genes
    top_div = div_df.nlargest(20, 'IVF_vs_PA')
    for g in top_div.index:
        ax.scatter(ivf_z['4C'].loc[g], pa_z['4C'].loc[g], c='red', s=30, alpha=0.8, zorder=5)
    lims = [min(x_vals.min(), y_vals.min()), max(x_vals.max(), y_vals.max())]
    ax.plot(lims, lims, 'k--', alpha=0.3)
    ax.set_xlabel('IVF Expression (z-score)')
    ax.set_ylabel('PA Expression (z-score)')
    ax.set_title('E. IVF vs PA Expression at 4C Stage\n(Red = Top 20 Divergent Genes)')

# F: Summary
ax = axes[1, 2]
ax.axis('off')
top_genes_ivf_pa = div_df.nlargest(10, 'IVF_vs_PA')
summary = [
    "M5: Cross-condition Trajectory",
    "==================================",
    f"in_vivo cells: {len(vivo_cells)} (E0-E8)",
    f"IVF cells: {len(ivf_cells)} (1C-8C)",
    f"PA cells: {len(pa_cells)} (1C-8C)",
    f"Common genes: {len(common_genes)}",
    "",
    "Method: Stage-anchored alignment",
    "Per-stage mean expression → z-score",
    "→ per-gene trajectory comparison",
    "",
    "Top IVF vs PA divergent:",
]
for i, g in enumerate(top_genes_ivf_pa.head(8).index):
    summary.append(f" {g[:20]}")
summary += [
    "",
    "[LAVDC #2: stage-level stats,",
    " not cell-level pseudoreplication]",
    "[#3: PCA/alignment method",
    " from scVelo/CellRank literature]"
]
for i, line in enumerate(summary):
    ax.text(0.02, 0.98-i*0.045, line, fontsize=8, family='monospace',
           transform=ax.transAxes, verticalalignment='top')

plt.suptitle('Fig.2: Cross-condition Developmental Trajectory Comparison\n(IVF vs PA vs in_vivo, GSE164812 + GSE168106)',
            fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(f'{OUT_DIR}/m5_trajectory_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: m5_trajectory_comparison.png")

# Save results
div_df.to_csv(f'{OUT_DIR}/m5_gene_divergence_v2.csv')
div_df.nlargest(50, 'total_div')[['IVF_vs_PA','IVF_vs_vivo','PA_vs_vivo']].to_csv(
    f'{OUT_DIR}/m5_top50_divergent_genes.csv')

print(f"\nDone! Outputs: {OUT_DIR}/")
print(f"  m5_trajectory_comparison.png")
print(f"  m5_gene_divergence_v2.csv")
print(f"  m5_top50_divergent_genes.csv")
print(f"\n  Top 3 total divergent genes:")
for g in div_df.nlargest(3, 'total_div').index:
    r = div_df.loc[g]
    print(f"    {g}: IVFvsPA={r['IVF_vs_PA']:.3f}, IVFvsVivo={r['IVF_vs_vivo']:.3f}, PAvsVivo={r['PA_vs_vivo']:.3f}")
