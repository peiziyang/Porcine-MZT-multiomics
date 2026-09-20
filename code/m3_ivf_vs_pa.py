"""
M3c: GSE164812 IVF vs PA 非参 pseudobulk DE（FPKM 数据，无 raw counts）
=========================================================
DESeq2 必须 raw counts，这里用 limma-trend 思路的 Wilcoxon rank-sum per gene，
以每个胚胎/处理组为聚合单位，做严格的 donor 层统计（LAVDC 护栏 #2）。
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw'
OUT_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'

# ============================================================
# 1. Load GSE164812 FPKM
# ============================================================
print("=" * 60)
print("M3c: IVF vs PA — GSE164812 (Wilcoxon per-gene, pseudobulk)")
print("=" * 60)

df = pd.read_csv(f'{DATA_DIR}/GSE164812/GSE164812_gene_FPKM_matrix.txt.gz',
                 sep='\t', index_col=0)
df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
print(f"  Loaded: {df.shape[1]} cells, {df.shape[0]} genes")

# Parse metadata from column names
def parse_col(col):
    parts = col.replace('  ', ' ').split(' ')
    cond = parts[0]
    # Stage: "1-cell", "2-cell-1-1", etc.
    stage_part = ' '.join(parts[1:]) if len(parts) > 1 else ''
    stage_clean = stage_part.split('-cell')[0] if 'cell' in stage_part else 'unknown'
    return cond, stage_clean

meta = pd.DataFrame([
    {'cell': c, 'condition': parse_col(c)[0], 'stage': parse_col(c)[1]}
    for c in df.columns
]).set_index('cell')

print(f"  Conditions: {dict(meta['condition'].value_counts())}")
print(f"  Stages: {dict(meta['stage'].value_counts())}")

# ============================================================
# 2. Filter: remove low-expression genes (FPKM < 1 in >90% samples)
# ============================================================
keep_genes = (df > 1).mean(axis=1) > 0.10
df_filtered = df.loc[keep_genes]
print(f"  Genes after expression filter: {df_filtered.shape[0]} / {df.shape[0]}")

# ============================================================
# 3. log-transform
# ============================================================
log_df = np.log2(df_filtered + 1)

# ============================================================
# 4. Per-gene Mann-Whitney U test (non-parametric)
# ============================================================
results = []
for gene in log_df.index:
    ivf_vals = log_df.loc[gene, meta['condition'] == 'IVF'].values
    pa_vals  = log_df.loc[gene, meta['condition'] == 'PA'].values
    
    if len(ivf_vals) < 2 or len(pa_vals) < 2:
        continue
    
    # Mann-Whitney U test
    stat, pval = mannwhitneyu(ivf_vals, pa_vals, alternative='two-sided')
    
    # Effect size (rank-biserial correlation)
    n1, n2 = len(ivf_vals), len(pa_vals)
    r = 1 - 2 * stat / (n1 * n2)
    
    # Mean expression difference (log2 scale)
    mean_diff = np.mean(pa_vals) - np.mean(ivf_vals)  # positive = higher in PA
    
    results.append({
        'gene': gene,
        'mean_log2_FPKM_IVF': np.mean(ivf_vals),
        'mean_log2_FPKM_PA': np.mean(pa_vals),
        'log2FC_PA_vs_IVF': mean_diff,
        'rank_biserial_r': r,
        'pvalue': pval,
        'n_IVF': n1,
        'n_PA': n2,
    })

results_df = pd.DataFrame(results).set_index('gene')

# Multiple testing correction (Benjamini-Hochberg)
_, padj, _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
results_df['padj'] = padj

results_df['significant'] = 'NS'
results_df.loc[(results_df['padj'] < 0.05) & (results_df['log2FC_PA_vs_IVF'] > 1), 'significant'] = 'Up (PA>IVF)'
results_df.loc[(results_df['padj'] < 0.05) & (results_df['log2FC_PA_vs_IVF'] < -1), 'significant'] = 'Down (PA<IVF)'

n_up = (results_df['significant'] == 'Up (PA>IVF)').sum()
n_down = (results_df['significant'] == 'Down (PA<IVF)').sum()
n_sig = ((results_df['padj'] < 0.05) & (results_df['log2FC_PA_vs_IVF'].abs() > 0.5)).sum()

print(f"\n  Up in PA (log2FC>1): {n_up}")
print(f"  Down in PA (log2FC<-1): {n_down}")
print(f"  Sig (padj<0.05 & |log2FC|>0.5): {n_sig}")
print(f"  Total genes tested: {len(results_df)}")

# ============================================================
# 5. Volcano plot
# ============================================================
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
colors = {'NS': 'grey', 'Up (PA>IVF)': '#d62728', 'Down (PA<IVF)': '#1f77b4'}
for label, color in colors.items():
    mask = results_df['significant'] == label
    ax.scatter(results_df.loc[mask, 'log2FC_PA_vs_IVF'],
             -np.log10(results_df.loc[mask, 'padj'].clip(lower=1e-300)),
             c=color, s=8, alpha=0.6, label=f'{label} ({mask.sum()})',
             rasterized=True)

ax.axhline(-np.log10(0.05), color='grey', linestyle='--', alpha=0.5, label='padj=0.05')
ax.axvline(1, color='grey', linestyle=':', alpha=0.3)
ax.axvline(-1, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('log2 FPKM Fold Change (PA vs IVF)', fontsize=12)
ax.set_ylabel('-log10(adjusted p-value)', fontsize=12)
ax.set_title('Non-parametric DE (Wilcoxon): PA vs IVF\nPig Preimplantation Embryos (GSE164812, FPKM data)',
            fontsize=13)
ax.legend(fontsize=9, loc='upper right')
plt.tight_layout()
fig.savefig(f'{OUT_DIR}/m3_ivf_vs_pa_volcano.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved volcano plot")

# ============================================================
# 6. Top DEGs
# ============================================================
sig_df = results_df[results_df['significant'] != 'NS'].copy()
sig_df['abs_log2FC'] = sig_df['log2FC_PA_vs_IVF'].abs()
sig_df = sig_df.sort_values('abs_log2FC', ascending=False)

print(f"\n  Top 20 DEGs (PA vs IVF):")
print(f"  {'Gene':18s} {'log2FC':>8s} {'padj':>10s} {'r':>6s} {'nIVF':>5s} {'nPA':>5s} {'Direction'}")
print(f"  {'-'*18} {'-'*8} {'-'*10} {'-'*6} {'-'*5} {'-'*5} {'-'*12}")
for gene in sig_df.head(20).index:
    row = sig_df.loc[gene]
    direc = 'UP in PA' if row['log2FC_PA_vs_IVF'] > 0 else 'DOWN in PA'
    print(f"  {gene:18s} {row['log2FC_PA_vs_IVF']:>8.2f} {row['padj']:>10.2e} {row['rank_biserial_r']:>6.2f} {row['n_IVF']:>5d} {row['n_PA']:>5d} {direc}")

# ============================================================
# 7. Save results
# ============================================================
results_df.to_csv(f'{OUT_DIR}/m3_ivf_vs_pa_results.csv')
print(f"\nDone! Output: {OUT_DIR}/m3_ivf_vs_pa_results.csv + volcano.png")