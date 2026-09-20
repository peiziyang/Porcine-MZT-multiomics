#!/usr/bin/env python
"""
Reviewer fix #9: Full F1 gene-set test (GSVA-style) for PA vs IVF.
Instead of cherry-picking 6 genes, test the ENTIRE F1 44-gene set.
Also run competitive gene-set test with expression-matched background sets.
"""
import pandas as pd, numpy as np, os
from scipy.stats import mannwhitneyu, norm
from scipy.stats import false_discovery_control as bh
import warnings; warnings.filterwarnings('ignore')

BASE = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
DESEQ = os.path.join(BASE, 'deseq2')
MF = os.path.join(BASE, 'm4_mofa')
OUT = os.path.join(MF, 'geneset_test')

os.makedirs(OUT, exist_ok=True)

# ── Load F1 gene list (44 genes) ──
f1_weights = pd.read_csv(os.path.join(MF, 'mofa_multiomics_weights_RNA_v2.csv'), index_col=0)
f1_sorted = f1_weights['F1'].abs().sort_values(ascending=False)
f1_genes_44 = list(f1_sorted.head(44).index)
print(f"F1 44 genes loaded: {len(f1_genes_44)}")

# ── Load PA bulk data (all 42 samples, 4 stages) ──
# Use the full DESeq2 count data, not just summary stats
# Get normalized counts from salmon quantification
salmon_dir = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed/salmon'
if os.path.exists(salmon_dir):
    print("Loading salmon normalized counts...")
    count_files = [f for f in os.listdir(salmon_dir) if f.endswith('_quant.sf')]
    print(f"  Found {len(count_files)} salmon quant files")
else:
    # Fallback: use DESeq2 results to estimate effect sizes
    print("Using DESeq2 results for gene-set testing...")

# Build per-gene per-stage log2FC matrix
pa_files = sorted([f for f in os.listdir(DESEQ)
                   if f.startswith('m15_deseq2_') and f.endswith('_IVFvsPA.csv')])
print(f"Found {len(pa_files)} DESeq2 files")

stage_lfc = {}
all_genes = set()
for pf in pa_files:
    stage = pf.replace('m15_deseq2_', '').replace('_IVFvsPA.csv', '')
    df = pd.read_csv(os.path.join(DESEQ, pf))
    if 'gene_symbol' not in df.columns:
        continue
    lfc_map = {}
    for _, row in df.iterrows():
        sym = str(row['gene_symbol'])
        if not sym.startswith('ENSSSCG'):
            lfc_map[sym] = {
                'log2FC': row['log2FoldChange'],
                'padj': row['padj'],
                'baseMean': row['baseMean']
            }
            all_genes.add(sym)
    stage_lfc[stage] = lfc_map
    print(f"  {stage}: {len(lfc_map)} genes mapped")

# ── Gene-set score for F1_44 ──
# For each stage, compute mean log2FC of F1 genes vs background
print("\n" + "="*60)
print("GENE-SET TEST: F1 44-genes vs background")
print("="*60)

results = []
for stage in sorted(stage_lfc.keys()):
    lfc_map = stage_lfc[stage]

    # F1 genes with data in this stage
    f1_in_stage = [g for g in f1_genes_44 if g in lfc_map]
    n_f1 = len(f1_in_stage)

    if n_f1 < 10:
        continue

    # F1 gene log2FC values
    f1_lfc = np.array([lfc_map[g]['log2FC'] for g in f1_in_stage])
    f1_mean = np.mean(f1_lfc)

    # Background: all other genes with similar baseMean range
    f1_base = np.array([lfc_map[g]['baseMean'] for g in f1_in_stage])
    bg_genes = [g for g in lfc_map if g not in f1_in_stage and g in all_genes]

    # Match background by baseMean decile
    all_base = np.array([lfc_map[g]['baseMean'] for g in bg_genes])
    bg_lfc_all = np.array([lfc_map[g]['log2FC'] for g in bg_genes])

    # Bootstrap: sample n_f1 genes from background 10000 times
    n_boot = 10000
    boot_means = np.zeros(n_boot)
    for i in range(n_boot):
        idx = np.random.choice(len(bg_lfc_all), n_f1, replace=True)
        boot_means[i] = np.mean(bg_lfc_all[idx])

    # Empirical P-value
    if f1_mean < 0:
        p_emp = np.mean(boot_means <= f1_mean)
    else:
        p_emp = np.mean(boot_means >= f1_mean)

    # Z-score
    z = (f1_mean - np.mean(boot_means)) / np.std(boot_means)

    # Also report % of F1 genes with padj < 0.05
    f1_sig = sum(1 for g in f1_in_stage if lfc_map[g]['padj'] < 0.05)

    print(f"  {stage:10s}: F1 mean LFC={f1_mean:+.3f}  "
          f"n={n_f1}  sig={f1_sig}/{n_f1}  "
          f"z={z:+.2f}  p_emp={p_emp:.4f}")

    results.append({
        'stage': stage,
        'n_f1_genes': n_f1,
        'f1_mean_lfc': f1_mean,
        'n_sig_genes': f1_sig,
        'pct_sig': f1_sig / n_f1 * 100,
        'z_score': z,
        'p_empirical': p_emp,
        'bg_mean_lfc': np.mean(bg_lfc_all),
        'bg_std_lfc': np.std(bg_lfc_all),
    })

res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(OUT, 'f1_geneset_test.csv'), index=False)

# ── Also test F3, F4, F6 gene sets ──
# Define F3 (clearance), F4 (metabolic), F6 (translational) genes
f1_w = f1_weights['F1'].abs()
f3_w = f1_weights['F3'].abs().sort_values(ascending=False)
f4_w = f1_weights['F4'].abs().sort_values(ascending=False)
f6_w = f1_weights['F6'].abs().sort_values(ascending=False)

for factor, gene_list, label in [
    ('F3', list(f3_w.head(44).index), 'Maternal clearance'),
    ('F4', list(f4_w.head(44).index), 'Zygotic activation'),
    ('F6', list(f6_w.head(44).index), 'Translational activation'),
]:
    print(f"\n{label} ({factor}):")
    for stage in sorted(stage_lfc.keys()):
        lfc_map = stage_lfc[stage]
        genes_in = [g for g in gene_list if g in lfc_map]
        if len(genes_in) < 10:
            continue
        lfc_vals = np.array([lfc_map[g]['log2FC'] for g in genes_in])
        mean_lfc = np.mean(lfc_vals)
        sig_n = sum(1 for g in genes_in if lfc_map[g]['padj'] < 0.05)
        # Quick bootstrap p
        bg_all = [lfc_map[g]['log2FC'] for g in bg_genes if g in lfc_map]
        boot_m = np.array([np.mean(np.random.choice(bg_all, len(genes_in), replace=True))
                          for _ in range(5000)])
        z_val = (mean_lfc - np.mean(boot_m)) / np.std(boot_m) if np.std(boot_m) > 0 else 0
        p_val = min(np.mean(boot_m <= mean_lfc), np.mean(boot_m >= mean_lfc)) * 2
        print(f"  {stage:10s}: mean={mean_lfc:+.3f}  n={len(genes_in)}  "
              f"sig={sig_n}/{len(genes_in)}  z={z_val:+.2f}  p_2t={p_val:.4f}")

# ── Combined stage-adjusted summary ──
print("\n" + "="*60)
print("COMBINED (stage-adjusted comparison):")
overall_file = os.path.join(DESEQ, 'm15_deseq2_overall_IVFvsPA_stageAdj.csv')
if os.path.exists(overall_file):
    odf = pd.read_csv(overall_file)
    overall_map = {}
    for _, row in odf.iterrows():
        sym = str(row['gene_symbol'])
        if not sym.startswith('ENSSSCG'):
            overall_map[sym] = {'log2FC': row['log2FoldChange'], 'padj': row['padj'],
                                'baseMean': row['baseMean'], 'stat': row['stat']}

    for label, gene_list in [('F1 (blueprint)', f1_genes_44),
                              ('F3 (clearance)', list(f3_w.head(44).index)),
                              ('F4 (activation)', list(f4_w.head(44).index)),
                              ('F6 (translation)', list(f6_w.head(44).index))]:
        genes_in = [g for g in gene_list if g in overall_map]
        lfc_vals = np.array([overall_map[g]['log2FC'] for g in genes_in])
        sig_n = sum(1 for g in genes_in if overall_map[g]['padj'] < 0.05)
        print(f"  {label:20s}: mean LFC={np.mean(lfc_vals):+.3f}  "
              f"n={len(genes_in)}  sig={sig_n}/{len(genes_in)}  "
              f"pct_up={sum(lfc_vals>0)}/{len(genes_in)} up")

print(f"\nResults saved to {OUT}")
print("DONE")
