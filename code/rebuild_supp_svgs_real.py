#!/usr/bin/env python
"""重建 14 个补充图 SVG（真实数据，fonttype='none'，可编辑）。
S1 (cpg_informative) 已有正确 SVG，跳过。"""
import os, pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

plt.rcParams.update({
    'font.family': 'Arial',
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
})

ROOT = r"E:/Workbuddy/2026-07-27-11-58-27"
P = os.path.join(ROOT, "data/processed")
OUT = os.path.join(ROOT, "manuscript/submission/figures")
os.makedirs(OUT, exist_ok=True)

def save(fig, name):
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT, f"{name}.svg"), format='svg', bbox_inches='tight')
    print(f"  ✓ {name}.svg")

# ============ S2: atlas UMAP ============
print("S2 atlas UMAP")
meta = pd.read_csv(f"{P}/m2_metadata.csv")
umap = pd.read_csv(f"{P}/m2_harmony_umap.csv")
meta = meta.merge(umap, on='cell', how='left')
meta['UMAP1'] = meta['UMAP1_harmony']
meta['UMAP2'] = meta['UMAP2_harmony']
fig, ax = plt.subplots(figsize=(9, 7))
stage_order = sorted(meta['stage_group'].dropna().unique(), key=str)
cmap = plt.cm.tab10
for i, st in enumerate(stage_order):
    sub = meta[meta['stage_group'] == st]
    ax.scatter(sub['UMAP1'], sub['UMAP2'], s=8, alpha=0.6, c=[cmap(i % 10)], label=st)
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.set_title('Integrated atlas UMAP (1,955 cells, Harmony-corrected)')
ax.legend(fontsize=6, ncol=3, title='Stage', loc='best', markerscale=3)
save(fig, 'supp_fig_atlas_umap'); plt.close()

# ============ S3: GO enrichment ============
print("S3 GO enrichment")
go = pd.read_csv(f"{P}/m4_mofa/factor_go_enrichment_v2.csv")
go['log10p'] = -np.log10(go['p_raw'].clip(lower=1e-300))
top = go.groupby('factor').apply(lambda x: x.nsmallest(5, 'p_raw')).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(11, 8))
colors = {'F1':'#377EB8','F2':'#FF7F00','F3':'#E41A1C','F4':'#4DAF4A','F5':'#984EA3','F6':'#A65628','F7':'#F781BF'}
for f in sorted(top['factor'].unique()):
    sub = top[top['factor']==f]
    ax.barh(sub.index[::-1], sub['log10p'][::-1], color=colors.get(f,'#999'), label=f)
terms = [f"{r['factor']} | {str(r['GO_term'])[:45]}" for _, r in top.iterrows()]
ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms, fontsize=7.5)
ax.set_xlabel('-log10(raw p)'); ax.invert_yaxis()
ax.set_title('GO enrichment (top 5 terms per factor, raw p)')
ax.legend(fontsize=8, ncol=4)
save(fig, 'supp_fig_go_enrich'); plt.close()

# ============ S4: conserved DEG heatmap ============
print("S4 conserved DEG heatmap")
stages = ['1cell','2cell','4cell','8cell']
stage_dfs = {}
for st in stages:
    d = pd.read_csv(f"{P}/deseq2/m15_deseq2_{st}_IVFvsPA.csv")
    d['gene_symbol'] = d['gene_symbol'].astype(str)
    stage_dfs[st] = d
deg_by_gene = {}
for st in stages:
    d = stage_dfs[st]
    sig = d[(d['padj']<0.05) & (d['log2FoldChange'].abs()>1)]
    for _, r in sig.iterrows():
        g = r['gene_symbol']
        if g.startswith('ENSSSCG'): continue
        deg_by_gene.setdefault(g, {})[st] = r['log2FoldChange']
cons = {g: v for g, v in deg_by_gene.items() if len(v) >= 3}
top_genes = sorted(cons, key=lambda g: -abs(np.mean(list(cons[g].values()))))[:30]
mat = np.array([[cons[g].get(st, 0) for st in stages] for g in top_genes])
fig, ax = plt.subplots(figsize=(6, max(6, len(top_genes)*0.28)))
im = ax.imshow(mat, aspect='auto', cmap='RdBu_r', vmin=-6, vmax=6)
ax.set_yticks(range(len(top_genes))); ax.set_yticklabels(top_genes, fontsize=8)
ax.set_xticks(range(4)); ax.set_xticklabels(['1C','2C','4C','8C'])
ax.set_title('Conserved DEGs (|log2FC|>1, padj<0.05, ≥3 stages)')
plt.colorbar(im, ax=ax, label='log2FC', shrink=0.8)
save(fig, 'supp_fig_conserved_heatmap'); plt.close()

# ============ S5: DESeq2 summary ============
print("S5 DESeq2 summary")
deseq = pd.read_csv(f"{P}/deseq2/m15_deseq2_overall_IVFvsPA_stageAdj.csv")
lfc = deseq['log2FoldChange'].values; padj = deseq['padj'].values
valid = ~np.isnan(padj) & ~np.isnan(lfc)
lfc_c, padj_c = lfc[valid], padj[valid]
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
ax = axes[0]; ax.hist(lfc_c, bins=60, color='#377EB8', edgecolor='black', lw=0.3)
ax.axvline(np.mean(lfc_c), color='#E41A1C', lw=2, label=f'mean={np.mean(lfc_c):.2f}')
ax.set_xlabel('log2FC (PA vs IVF)'); ax.set_ylabel('Genes'); ax.set_title('LFC distribution'); ax.legend(fontsize=8)
ax = axes[1]; sig = padj_c < 0.05
ax.scatter(lfc_c[~sig], -np.log10(padj_c[~sig].clip(1e-300)), c='gray', s=4, alpha=0.3)
ax.scatter(lfc_c[sig], -np.log10(padj_c[sig].clip(1e-300)), c='red', s=8, alpha=0.5)
ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
ax.set_xlabel('log2FC'); ax.set_ylabel('-log10(padj)'); ax.set_title(f'Volcano ({sum(sig)} sig)')
ax = axes[2]; bm = deseq['baseMean'].dropna().values
ax.hist(np.log10(bm+1), bins=50, color='#4DAF4A', edgecolor='black', lw=0.3)
ax.set_xlabel('log10(baseMean)'); ax.set_ylabel('Genes'); ax.set_title('baseMean')
save(fig, 'supp_fig_deseq2_summary'); plt.close()

# ============ S6: factor correlation matrix (15×15) ============
print("S6 factor correlation matrix")
mf = pd.read_csv(f"{P}/m4_mofa/mofa_multiomics_factors.csv", index_col=0)
rf = pd.read_csv(f"{P}/m8_mofa_factors.csv", index_col=0)
rf_factors = rf[[f"F{i}" for i in range(1,11)]]
common = sorted(set(mf.index) & set(rf.index))
rf_c = rf_factors.loc[common]
mo = mf.loc[common]
corr = pd.concat([rf_c.add_prefix("RNA_"), mo.add_prefix("MO_")], axis=1).corr()
fig, ax = plt.subplots(figsize=(9, 8))
im = ax.imshow(corr.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
ax.set_xticks(range(len(corr.columns))); ax.set_yticks(range(len(corr.columns)))
ax.set_xticklabels(corr.columns, fontsize=6, rotation=90)
ax.set_yticklabels(corr.columns, fontsize=6)
ax.set_title('Factor correlation: RNA-only (10) vs Multi-omics (5)')
plt.colorbar(im, ax=ax, shrink=0.7)
save(fig, 'supp_fig_factor_correlation'); plt.close()

# ============ S7: IVF vs PA volcano ============
print("S7 IVF vs PA volcano")
v = pd.read_csv(f"{P}/m3_ivf_vs_pa_results.csv")
lfc = v['log2FC_PA_vs_IVF'].values; pv = v['padj'].values
valid = ~np.isnan(pv) & ~np.isnan(lfc)
lfc_c, pv_c = lfc[valid], pv[valid]
sig = (pv_c < 0.05) & (np.abs(lfc_c) > 0.5)
fig, ax = plt.subplots(figsize=(7.5, 6))
ax.scatter(lfc_c[~sig], -np.log10(pv_c[~sig].clip(1e-300)), c='gray', s=5, alpha=0.3)
ax.scatter(lfc_c[sig], -np.log10(pv_c[sig].clip(1e-300)), c='red', s=8, alpha=0.5)
ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
ax.set_xlabel('log2FC (PA vs IVF)'); ax.set_ylabel('-log10(padj)')
ax.set_title(f'IVF vs PA (GSE164812): {sum(sig)} significant')
save(fig, 'supp_fig_ivf_vs_pa_volcano'); plt.close()

# ============ S8: metabolism ============
print("S8 metabolism")
vivo = pd.read_csv(f"{P}/m6_metabolism_vivo.csv")
paivf = pd.read_csv(f"{P}/m6_metabolism_paivf.csv")
ratio = pd.read_csv(f"{P}/m6_metabolism_ratio.csv")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
ax = axes[0]
for st in sorted(vivo['stage'].dropna().unique(), key=str):
    sub = vivo[vivo['stage']==st]
    ax.scatter(sub['glycolysis'], sub['oxphos'], s=8, alpha=0.5, label=st)
ax.set_xlabel('Glycolysis score'); ax.set_ylabel('OxPhos score'); ax.set_title('Metabolic scores (in vivo)')
ax.legend(fontsize=6, ncol=2)
ax = axes[1]
r = ratio.sort_values('stage')
ax.plot(r['stage'], r['gly_mean'], 'o-', label='Glycolysis')
ax.plot(r['stage'], r['ox_mean'], 's-', label='OxPhos')
ax.set_xlabel('Stage'); ax.set_ylabel('Mean score'); ax.set_title('Scores by stage'); ax.legend(fontsize=7)
ax = axes[2]
ax.bar(r['stage'], r['ratio'], color='#4DAF4A')
ax.set_xlabel('Stage'); ax.set_ylabel('Gly/Ox ratio'); ax.set_title('Glycolysis:OxPhos ratio')
save(fig, 'supp_fig_metabolism'); plt.close()

# ============ S9: RNA-only MOFA+ ============
print("S9 RNA-only MOFA+")
rf9 = pd.read_csv(f"{P}/m8_mofa_factors.csv", index_col=0)
r2 = pd.read_csv(f"{P}/m8_mofa_r2.csv")
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
ax = axes[0,0]
r2_sorted = r2.sort_values('factor')
ax.bar(r2_sorted['factor'], r2_sorted['R2_pct'], color='#377EB8')
ax.set_xlabel('Factor'); ax.set_ylabel('R² (%)'); ax.set_title('Variance explained')
ax = axes[0,1]
stage_col = 'stage_label' if 'stage_label' in rf9.columns else 'dataset'
if stage_col in rf9.columns:
    for f in ['F1','F2','F3']:
        tmp = rf9.groupby(stage_col)[f].mean()
        ax.plot(range(len(tmp)), tmp.values, 'o-', label=f)
    ax.set_xticks(range(len(tmp))); ax.set_xticklabels(tmp.index, fontsize=6, rotation=45)
ax.set_xlabel('Stage'); ax.set_ylabel('Mean factor score'); ax.set_title('Factor by stage'); ax.legend(fontsize=7)
ax = axes[1,0]
ax.scatter(rf9['F1'], rf9['F2'], c='#4DAF4A', s=8, alpha=0.5)
ax.set_xlabel('F1'); ax.set_ylabel('F2'); ax.set_title('F1 vs F2')
ax = axes[1,1]
ax.hist(rf9['F1'].dropna(), bins=40, color='#E41A1C', edgecolor='black', lw=0.3)
ax.set_xlabel('F1 score'); ax.set_ylabel('Cells'); ax.set_title('F1 distribution')
save(fig, 'supp_fig_mofa_rna_only'); plt.close()

# ============ S10: overview ============
print("S10 overview")
otypes = pd.read_csv(f"{P}/m4_oocyte_types.csv")
odegs = pd.read_csv(f"{P}/m4_oocyte_type_degs.csv")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
ax = axes[0]
for t in sorted(otypes['oocyte_type'].dropna().unique()):
    sub = otypes[otypes['oocyte_type']==t]
    ax.hist(sub['total_transcripts'], bins=30, alpha=0.5, label=t)
ax.set_xlabel('Total transcripts'); ax.set_ylabel('Cells'); ax.set_title('Transcript count by type'); ax.legend(fontsize=7)
ax = axes[1]
lfc_o = odegs['log2FC_II_vs_I'].dropna().values
pv_o = odegs['padj'].dropna().values
sig_o = (pv_o < 0.05) & (np.abs(lfc_o) > 0.5)
ax.scatter(lfc_o[~sig_o], -np.log10(pv_o[~sig_o].clip(1e-300)), c='gray', s=4, alpha=0.3)
ax.scatter(lfc_o[sig_o], -np.log10(pv_o[sig_o].clip(1e-300)), c='red', s=7, alpha=0.5)
ax.set_xlabel('log2FC (II vs I)'); ax.set_ylabel('-log10(padj)'); ax.set_title('Oocyte type DEGs')
ax = axes[2]
top_deg = odegs[odegs['padj']<0.05].nsmallest(15, 'log2FC_II_vs_I')
ax.barh(range(len(top_deg)), top_deg['log2FC_II_vs_I'].values, color='#984EA3')
ax.set_yticks(range(len(top_deg))); ax.set_yticklabels(top_deg['gene'], fontsize=7)
ax.set_xlabel('log2FC'); ax.invert_yaxis(); ax.set_title('Top DEGs')
save(fig, 'supp_fig_overview'); plt.close()

# ============ S11: regulon by stage ============
print("S11 regulon by stage")
aucell = pd.read_csv(f"{P}/pyscenic_out/aucell_scores_full.csv", index_col=0)
sm = pd.read_csv(f"{P}/m7_cell_stage_map.csv")
top_reg = aucell.var().abs().sort_values(ascending=False).head(15).index
stage_col = 'stage_label' if 'stage_label' in sm.columns else 'stage_order'
sm = sm.set_index('cell_id')
aucell = aucell.loc[aucell.index.isin(sm.index)]
sm = sm.loc[aucell.index]
fig, ax = plt.subplots(figsize=(10, 5.5))
for reg in top_reg:
    tmp = aucell.groupby(sm[stage_col])[reg].mean()
    ax.plot(range(len(tmp)), tmp.values, 'o-', markersize=4, lw=1.2, label=reg.split('(')[0])
ax.set_xticks(range(len(tmp))); ax.set_xticklabels(tmp.index, fontsize=7, rotation=45)
ax.set_ylabel('Mean AUC'); ax.set_title('Regulon activity by stage (top 15)')
ax.legend(fontsize=6, ncol=3)
save(fig, 'supp_fig_regulon_by_stage'); plt.close()

# ============ S12: regulon heatmap ============
print("S12 regulon heatmap")
top10 = aucell.var().abs().sort_values(ascending=False).head(10).index
order = sm.sort_values(stage_col).index
mat = aucell.loc[order, top10].T.values
fig, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(mat, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
ax.set_yticks(range(len(top10))); ax.set_yticklabels([r.split('(')[0] for r in top10], fontsize=8)
ax.set_xlabel('Cells (ordered by stage)')
plt.colorbar(im, ax=ax, label='AUC', shrink=0.8)
ax.set_title('Regulon activity heatmap (top 10)')
save(fig, 'supp_fig_regulon_heatmap'); plt.close()

# ============ S13: SHAP ============
print("S13 SHAP")
shap = pd.read_csv(f"{P}/m4_mofa/shap_analysis_results.csv")
top_shap = shap.nsmallest(20, 'mean_abs_shap') if 'mean_abs_shap' in shap.columns else shap.nlargest(20, 'mean_abs_shap')
top_shap = shap.sort_values('mean_abs_shap', ascending=False).head(20)
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(range(len(top_shap)), top_shap['mean_abs_shap'].values[::-1], color='#FF7F00')
ax.set_yticks(range(len(top_shap))); ax.set_yticklabels(top_shap['gene'].values[::-1], fontsize=8)
ax.set_xlabel('Mean |SHAP|'); ax.set_title('SHAP feature importance (donor classification)')
save(fig, 'supp_fig_shap'); plt.close()

# ============ S14: trajectory divergence ============
print("S14 trajectory divergence")
gdiv = pd.read_csv(f"{P}/m5_gene_divergence_v2.csv")
top50 = pd.read_csv(f"{P}/m5_top50_divergent_genes.csv")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
ax.hist(gdiv['total_div'].clip(0, 5), bins=60, color='#377EB8', edgecolor='black', lw=0.3)
ax.set_xlabel('Total divergence'); ax.set_ylabel('Genes'); ax.set_title('Gene-wise divergence')
ax = axes[1]
t = top50.head(15)
x = np.arange(len(t))
ax.barh(x, t['IVF_vs_PA'].values[::-1], color='#E41A1C', alpha=0.7, label='IVF vs PA')
ax.set_yticks(x); ax.set_yticklabels(t['gene'].values[::-1], fontsize=7)
ax.set_xlabel('Divergence'); ax.set_title('Top divergent genes'); ax.legend(fontsize=7)
save(fig, 'supp_fig_trajectory_divergence'); plt.close()

# ============ S15: Waddington OT ============
print("S15 Waddington OT")
fate = pd.read_csv(f"{P}/m5_ot_vivo_fate.csv")
disp = pd.read_csv(f"{P}/m5_ot_displacement.csv")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
ax.scatter(fate['UMAP1'], fate['UMAP2'], c=fate['fate_entropy'], cmap='viridis', s=6, alpha=0.6)
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2'); ax.set_title('Fate entropy')
plt.colorbar(ax.collections[0], ax=ax, label='entropy', shrink=0.8)
ax = axes[1]
for cond in disp['condition'].dropna().unique():
    sub = disp[disp['condition']==cond]
    tmp = sub.groupby('stage')['displacement'].mean()
    ax.plot(range(len(tmp)), tmp.values, 'o-', label=cond)
ax.set_xticks(range(len(tmp))); ax.set_xticklabels(tmp.index, fontsize=7, rotation=45)
ax.set_ylabel('Displacement'); ax.set_title('OT displacement'); ax.legend(fontsize=7)
save(fig, 'supp_fig_waddington_ot'); plt.close()

print("\n=== 全部 14 个补充图 SVG 重建完成 ===")
