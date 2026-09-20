"""
M3: Pseudobulk 差异表达分析 — Donor/样本层统计
==================================================
用 GSE168106 原始 counts，按 E-day 做 pseudobulk 聚合，
DESeq2 比较 early (E0-E3) vs late (E11-E14) 胚胎发育阶段。

严格遵循 LAVDC 护栏 #2：以生物学重复为单位做 pseudobulk，
报效应量+adjusted p-value，功效不足标注。
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw'
OUT_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'

# ============================================================
# 1. Load data
# ============================================================
print("=" * 60)
print("M3: Pseudobulk DE — GSE168106 E0-E3 vs E11-E14")
print("=" * 60)

counts = pd.read_csv(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_genes_counts.txt.gz',
                     sep='\t', index_col=0)
si = pd.read_excel(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_SampleInfo.xlsx')
si_valid = si[si['Vaild'] == 'Yes'].copy()

# Map each cell barcode to its developmental day
cell_to_day = dict(zip(si_valid['ID'], si_valid['Stage_II']))
cell_to_lane = dict(zip(si_valid['ID'], si_valid['Lane']))

# Only keep cells with valid day annotation
valid_cells = [c for c in counts.columns if c in cell_to_day and pd.notna(cell_to_day[c])]
counts = counts[valid_cells]
print(f"  Loaded: {len(valid_cells)} cells with day annotation")

# ============================================================
# 2. Define groups
# ============================================================
early_days = ['E0', 'E1', 'E2', 'E3']
late_days  = ['E11', 'E12', 'E13', 'E14']

# Group cells by E-day → pseudobulk
def make_pseudobulk(counts, cell_to_day, days, label):
    """Aggregate cells by E-day into pseudobulk samples."""
    pb_samples = {}
    for day in days:
        day_cells = [c for c in counts.columns if cell_to_day.get(c) == day]
        if len(day_cells) > 0:
            pb = counts[day_cells].sum(axis=1)
            pb.name = f'{label}_{day}'
            pb_samples[f'{label}_{day}'] = pb
    return pd.DataFrame(pb_samples)

pb_early = make_pseudobulk(counts, cell_to_day, early_days, 'Early')
pb_late  = make_pseudobulk(counts, cell_to_day, late_days, 'Late')

# Combine
pb_all = pd.concat([pb_early, pb_late], axis=1)
print(f"  Early pseudobulk samples: {list(pb_early.columns)}")
print(f"  Late pseudobulk samples:  {list(pb_late.columns)}")

# ============================================================
# 3. Metadata for DESeq2
# ============================================================
metadata = pd.DataFrame({
    'sample': pb_all.columns,
    'condition': ['Early'] * pb_early.shape[1] + ['Late'] * pb_late.shape[1]
}).set_index('sample')

# ============================================================
# 4. Filter low-expressed genes
# ============================================================
# Keep genes with at least 10 counts total across all pseudobulk samples
gene_sums = pb_all.sum(axis=1)
pb_filtered = pb_all.loc[gene_sums >= 10]
print(f"  Genes after filtering (min 10 total counts): {pb_filtered.shape[0]} / {pb_all.shape[0]}")

# ============================================================
# 5. Run DESeq2
# ============================================================
print("\n  Running DESeq2...")
dds = DeseqDataSet(
    counts=pb_filtered.T.astype(int),
    metadata=metadata,
    design="~condition",
)
dds.deseq2()

# Get results for Early vs Late
stat_res = DeseqStats(dds, contrast=["condition", "Early", "Late"])
stat_res.summary()
results = stat_res.results_df

# ============================================================
# 6. Process results
# ============================================================
results['log2FC'] = results['log2FoldChange']
results['-log10(padj)'] = -np.log10(results['padj'].clip(lower=1e-300))
results['significant'] = 'NS'
results.loc[(results['padj'] < 0.05) & (results['log2FC'] > 1), 'significant'] = 'Up (Late>Early)'
results.loc[(results['padj'] < 0.05) & (results['log2FC'] < -1), 'significant'] = 'Down (Late<Early)'

n_up = (results['significant'] == 'Up (Late>Early)').sum()
n_down = (results['significant'] == 'Down (Late<Early)').sum()
n_total = len(results)
print(f"\n  Results: {n_up} up-regulated, {n_down} down-regulated (padj<0.05, |log2FC|>1)")
print(f"  Total genes tested: {n_total}")
print(f"  Sig at padj<0.05 (any FC): {(results['padj'] < 0.05).sum()}")

# ============================================================
# 7. Volcano plot
# ============================================================
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
colors = {'NS': 'grey', 'Up (Late>Early)': '#d62728', 'Down (Late<Early)': '#1f77b4'}
for label, color in colors.items():
    mask = results['significant'] == label
    ax.scatter(results.loc[mask, 'log2FC'],
              results.loc[mask, '-log10(padj)'],
              c=color, s=5, alpha=0.6, label=f'{label} ({mask.sum()})',
              rasterized=True)

# Label top genes
if n_up + n_down > 0:
    top_genes = results[results['significant'] != 'NS'].nlargest(10, 'log2FC')
    top_genes = pd.concat([top_genes, results[results['significant'] != 'NS'].nsmallest(10, 'log2FC')])
    for gene in top_genes.index[:20]:
        row = results.loc[gene]
        ax.annotate(gene.replace('ENSSSCG', 'E'),
                   (row['log2FC'], row['-log10(padj)']),
                   fontsize=6, alpha=0.8)

ax.axhline(-np.log10(0.05), color='grey', linestyle='--', alpha=0.5, label='padj=0.05')
ax.axvline(1, color='grey', linestyle=':', alpha=0.3)
ax.axvline(-1, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('log2 Fold Change (Late vs Early)', fontsize=12)
ax.set_ylabel('-log10(adjusted p-value)', fontsize=12)
ax.set_title('Pseudobulk DESeq2: Late (E11-E14) vs Early (E0-E3)\nPig Preimplantation Embryo Development',
            fontsize=14)
ax.legend(fontsize=8, loc='upper right')
plt.tight_layout()
fig.savefig(f'{OUT_DIR}/m3_volcano_early_vs_late.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved volcano plot")

# ============================================================
# 8. MA plot
# ============================================================
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
base_mean = results['baseMean']
log2fc = results['log2FC']
colors_arr = ['grey'] * len(results)
for i, sig in enumerate(results['significant']):
    if sig.startswith('Up'): colors_arr[i] = '#d62728'
    elif sig.startswith('Down'): colors_arr[i] = '#1f77b4'
ax.scatter(np.log10(base_mean + 1), log2fc, c=colors_arr, s=3, alpha=0.5, rasterized=True)
ax.axhline(0, color='grey', linestyle='-', alpha=0.5)
ax.axhline(1, color='grey', linestyle=':', alpha=0.3)
ax.axhline(-1, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('log10(Mean Expression + 1)', fontsize=12)
ax.set_ylabel('log2 Fold Change', fontsize=12)
ax.set_title('MA Plot: Late vs Early', fontsize=14)
plt.tight_layout()
fig.savefig(f'{OUT_DIR}/m3_maplot_early_vs_late.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved MA plot")

# ============================================================
# 9. Identify and annotate top genes (map Ensembl to symbols using pig annotation)
# ============================================================
print(f"\n  Top DEGs:")
top_df = results[results['significant'] != 'NS'].copy()
top_df = top_df.sort_values('padj')
top_df['abs_log2FC'] = top_df['log2FC'].abs()
top_df = top_df.sort_values('abs_log2FC', ascending=False)

print(f"\n  {'Gene ID':20s} {'log2FC':>8s} {'padj':>10s}  {'Direction'}")
print(f"  {'-'*20} {'-'*8} {'-'*10}  {'-'*10}")
for gene in top_df.head(20).index:
    row = top_df.loc[gene]
    direc = 'UP (Late)' if row['log2FC'] > 0 else 'DOWN (Late)'
    print(f"  {gene:20s} {row['log2FC']:>8.2f} {row['padj']:>10.2e}  {direc}")

# ============================================================
# 10. Save results
# ============================================================
results_clean = results[['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'significant']].copy()
results_clean.index.name = 'Gene'
results_clean.to_csv(f'{OUT_DIR}/m3_deseq2_results_early_vs_late.csv')

# Per-donor pseudobulk counts
pb_all.index.name = 'Gene'
pb_all.to_csv(f'{OUT_DIR}/m3_pseudobulk_counts.csv')

# Summary stats
summary = {
    'test': 'Early (E0-E3) vs Late (E11-E14)',
    'method': 'Pseudobulk DESeq2 per E-day',
    'design': '~ condition',
    'n_early_samples': pb_early.shape[1],
    'n_late_samples': pb_late.shape[1],
    'n_genes_tested': n_total,
    'n_sig_up': n_up,
    'n_sig_down': n_down,
    'n_sig_total': n_up + n_down,
    'notes': 'E-day used as replicate unit. For donor-level stats, individual embryo IDs needed.',
    'warning': 'Small n per group (n=4 each). Effect sizes reported; power limited.'
}
for k, v in summary.items():
    print(f"  {k}: {v}")

print(f"\nDone! Outputs: {OUT_DIR}/")
print(f"  - m3_volcano_early_vs_late.png")
print(f"  - m3_maplot_early_vs_late.png")
print(f"  - m3_deseq2_results_early_vs_late.csv")
print(f"  - m3_pseudobulk_counts.csv")
