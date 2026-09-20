#!/usr/bin/env python
"""Fix S3, S5, S6, S12 with correct column names from actual data."""
import matplotlib; matplotlib.use('Agg')
matplotlib.rcParams.update({'font.family':'Arial','svg.fonttype':'none','savefig.dpi':300})
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd, numpy as np, os

ROOT = r'E:/Workbuddy/2026-07-27-11-58-27'
OUT = os.path.join(ROOT, 'manuscript', 'submission', 'figures')

def save(fig, name):
    fig.savefig(os.path.join(OUT, f'{name}.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT, f'{name}.svg'), format='svg', bbox_inches='tight')
    print(f'  {name}.svg ({os.path.getsize(os.path.join(OUT, f"{name}.svg"))//1024} KB)')

# ============================================================
# S3: DESeq2 Summary (fix NaN alignment)
# ============================================================
print('S3: DESeq2 Summary')
df = pd.read_csv(f'{ROOT}/data/processed/deseq2/m15_deseq2_overall_IVFvsPA_stageAdj.csv')
lfc = df['log2FoldChange'].values
padj = df['padj'].values
valid = ~np.isnan(padj) & ~np.isnan(lfc)
lfc_c = lfc[valid]
padj_c = padj[valid]

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.suptitle('S3: DESeq2 Stage-Adjusted PA-vs-IVF Summary (SRP301735, n=42, 8,726 genes with padj)', fontsize=11, fontweight='bold')

ax = axes[0]
ax.hist(lfc_c, bins=60, color='#377EB8', edgecolor='black', lw=0.3)
ax.axvline(0, color='gray', ls='--', lw=1)
ax.axvline(np.mean(lfc_c), color='#E41A1C', lw=2, ls='-', label=f'mean={np.mean(lfc_c):.2f}')
ax.set_xlabel('log2 FC (PA vs IVF)'); ax.set_ylabel('Gene count')
ax.set_title('A: LFC Distribution (genes with padj)'); ax.legend(fontsize=8)

ax = axes[1]
sig = padj_c < 0.05
ax.scatter(lfc_c[~sig], -np.log10(np.clip(padj_c[~sig], 1e-300, None)), c='gray', s=5, alpha=0.3)
ax.scatter(lfc_c[sig], -np.log10(np.clip(padj_c[sig], 1e-300, None)), c='red', s=8, alpha=0.5)
ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
ax.set_xlabel('log2 FC'); ax.set_ylabel('-log10(padj)')
ax.set_title(f'B: Volcano ({sum(sig)} significant)')

ax = axes[2]
bm = df['baseMean'].dropna().values
ax.hist(np.log10(bm + 1), bins=50, color='#4DAF4A', edgecolor='black', lw=0.3)
ax.set_xlabel('log10(baseMean)'); ax.set_ylabel('Gene count')
ax.set_title('C: BaseMean Distribution')
save(fig, 'supp_fig_deseq2_summary')
plt.close('all')

# ============================================================
# S5: GO Enrichment (correct columns)
# ============================================================
print('S5: GO Enrichment')
go = pd.read_csv(f'{ROOT}/data/processed/m4_mofa/factor_go_enrichment_v2.csv')
# Use p_raw column, GO_term column
top_terms = go.groupby('factor').apply(lambda x: x.nsmallest(5, 'p_raw')).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(12, 8))
fig.suptitle('S5: GO Enrichment (Top 5 Terms per Factor)', fontsize=11, fontweight='bold')
terms = []
factors_list = []
pvals = []
for _, row in top_terms.iterrows():
    terms.append(str(row['GO_term'])[:60])
    factors_list.append(str(row['factor']))
    pvals.append(-np.log10(row['p_raw'] + 1e-300))

colors = {'F1':'#377EB8','F3':'#E41A1C','F4':'#4DAF4A','F6':'#984EA3',
          'F2':'#FF7F00','F5':'#A65628','F7':'#F781BF'}
bar_colors = [colors.get(f, '#999') for f in factors_list]

ax.barh(range(len(terms)), pvals, color=bar_colors)
ax.set_yticks(range(len(terms)))
ax.set_yticklabels([f'{f} | {t}' for f, t in zip(factors_list, terms)], fontsize=8)
ax.set_xlabel('-log10(p-value)'); ax.invert_yaxis()
handles = [mpatches.Patch(color=c, label=f) for f, c in colors.items() if f in set(factors_list)]
ax.legend(handles=handles, fontsize=8, title='Factor', ncol=2)
save(fig, 'supp_fig_go_enrich')
plt.close('all')

# ============================================================
# S6: IVF vs PA Volcano (correct columns)
# ============================================================
print('S6: IVF vs PA Volcano')
df_ivf = pd.read_csv(f'{ROOT}/data/processed/m3_ivf_vs_pa_results.csv')
lfc = df_ivf['log2FC_PA_vs_IVF'].values
pv = df_ivf['padj'].values
valid = ~np.isnan(pv) & ~np.isnan(lfc)
lfc_c = lfc[valid]
pv_c = pv[valid]
sig = (pv_c < 0.05) & (np.abs(lfc_c) > 0.5)

fig, ax = plt.subplots(figsize=(8, 6))
fig.suptitle('S6: IVF vs PA Differential Expression (GSE164812)', fontsize=11, fontweight='bold')
ax.scatter(lfc_c[~sig], -np.log10(np.clip(pv_c[~sig], 1e-300, None)), c='gray', s=6, alpha=0.3)
ax.scatter(lfc_c[sig], -np.log10(np.clip(pv_c[sig], 1e-300, None)), c='red', s=8, alpha=0.5)
ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
ax.set_xlabel('log2 FC (PA vs IVF)'); ax.set_ylabel('-log10(padj)')
ax.set_title(f'Volcano: {sum(sig)} significant genes (|FC|>0.5, padj<0.05)', fontsize=10)
save(fig, 'supp_fig_ivf_vs_pa_volcano')
plt.close('all')

# ============================================================
# S12: SHAP (correct columns)
# ============================================================
print('S12: SHAP Analysis')
shap = pd.read_csv(f'{ROOT}/data/processed/m4_mofa/shap_analysis_results.csv')
top = shap.nlargest(15, 'mean_abs_shap')

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('S12: SHAP Feature Importance (Donor Classification)', fontsize=11, fontweight='bold')

ax = axes[0]
clrs = plt.cm.viridis(np.linspace(0, 1, len(top)))
ax.barh(range(len(top)), top['mean_abs_shap'].values, color=clrs)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top['gene'].values, fontsize=9)
ax.set_xlabel('Mean |SHAP|'); ax.set_title('A: Top 15 Features by SHAP Importance'); ax.invert_yaxis()

ax = axes[1]
ax.scatter(top['mean_abs_shap'].values, range(len(top)), c=clrs, s=80, edgecolors='black', lw=0.5)
for i, (g, v) in enumerate(zip(top['gene'], top['mean_abs_shap'])):
    ax.text(v + 0.01, i, f'  {g}', fontsize=8, va='center')
ax.set_yticks([]); ax.set_xlabel('Mean |SHAP|')
ax.set_title('B: Feature Importance Lollipop')
save(fig, 'supp_fig_shap')
plt.close('all')

print('\n=== 4 FIXED SVGs DONE ===')
