#!/usr/bin/env python
"""批量修复 14 个补充图中的重叠/排序/标签问题。

修复清单:
- S3 conserved_heatmap: 增宽 fig 防止 colorbar 标签裁切
- S7 metabolism: stage 按解析数字排序（不是字母序）+ rotation
- S8 mofa_rna_only: stage 按解析数字排序 + F1 distribution 标题加边距
- S9 overview: x 轴禁用 1e7 科学计数 + gene 名 ENSSSCG→symbol
- S10 regulon_by_stage: stage 按解析数字排序
- S13 trajectory_divergence: gene 名 ENSSSCG→symbol
"""
import os, re, gzip
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    plt.close(fig)


# ENSSSCG → gene symbol 映射
ID2SYM = {}
for path in [os.path.join(ROOT, "data/reference/Sus_scrofa.Sscrofa11.1.113.gtf.gz"),
             os.path.join(ROOT, "data/raw/pig_genes.gtf.gz")]:
    if not os.path.exists(path): continue
    try:
        with gzip.open(path, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                p = line.strip().split('\t')
                if len(p) < 9 or p[2] != 'gene': continue
                attrs = dict(re.findall(r'(\w+)\s*"([^"]+)"', p[8]))
                gid = attrs.get('gene_id')
                sym = attrs.get('gene_name', '')
                if gid and sym and sym != gid:
                    ID2SYM[gid] = sym
        break
    except Exception:
        continue
# 合并 deseq2 映射
d = pd.read_csv(os.path.join(P, "deseq2/m15_deseq2_overall_IVFvsPA_stageAdj.csv"))
for _, r in d.iterrows():
    gid = str(r['gene'])
    sym = str(r.get('gene_symbol', ''))
    if sym and not sym.startswith('ENSSSCG') and not sym.lower() == 'nan':
        ID2SYM[gid] = sym
print(f"ID2SYM 映射数: {len(ID2SYM)}")


def to_symbol(g):
    """把 ENSSSCG ID 转成 gene symbol，找不到返回原 ID."""
    if not g or not str(g).startswith('ENSSSCG'):
        return g
    return ID2SYM.get(g, g)


def parse_stage_key(s):
    """把 'E3' / 'GV' / 'in_vitro' / 'pgEpiSC' 转成可比较的 key.
    返回 (group_rank, num). E0..E14 → (0, 0..14), GV → (1, 0), in_vitro → (2, 0), pgEpiSC → (3, 0)."""
    m = re.match(r'E(\d+)', str(s))
    if m: return (0, int(m.group(1)))
    if str(s).startswith('GV'): return (1, 0)
    if 'vitro' in str(s).lower(): return (2, 0)
    if 'Epi' in str(s): return (3, 0)
    return (4, 0)


# ============ S3 conserved heatmap（增宽） ============
print("S3 conserved heatmap")
stages = ['1cell', '2cell', '4cell', '8cell']
stage_dfs = {}
for st in stages:
    df = pd.read_csv(f"{P}/deseq2/m15_deseq2_{st}_IVFvsPA.csv")
    df['gene_symbol'] = df['gene_symbol'].astype(str)
    stage_dfs[st] = df
deg_by_gene = {}
for st in stages:
    df = stage_dfs[st]
    sig = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 1)]
    for _, r in sig.iterrows():
        g = str(r['gene_symbol'])
        if g.startswith('ENSSSCG'): continue
        deg_by_gene.setdefault(g, {})[st] = r['log2FoldChange']
cons = {g: v for g, v in deg_by_gene.items() if len(v) >= 3}
top_genes = sorted(cons, key=lambda g: -abs(np.mean(list(cons[g].values()))))[:30]
mat = np.array([[cons[g].get(st, 0) for st in stages] for g in top_genes])
fig, ax = plt.subplots(figsize=(7, max(6, len(top_genes) * 0.28)))
im = ax.imshow(mat, aspect='auto', cmap='RdBu_r', vmin=-6, vmax=6)
ax.set_yticks(range(len(top_genes)))
ax.set_yticklabels(top_genes, fontsize=8)
ax.set_xticks(range(4))
ax.set_xticklabels(['1C', '2C', '4C', '8C'])
ax.set_title('Recurrent DEGs (|log2FC|>1, padj<0.05, \u22653 stages)')
plt.colorbar(im, ax=ax, label='log2FC', shrink=0.7, pad=0.02)
save(fig, 'supp_fig_conserved_heatmap')

# ============ S7 metabolism（解析 stage 排序） ============
print("S7 metabolism")
vivo = pd.read_csv(f"{P}/m6_metabolism_vivo.csv")
paivf = pd.read_csv(f"{P}/m6_metabolism_paivf.csv")
ratio = pd.read_csv(f"{P}/m6_metabolism_ratio.csv")
# 解析排序
ratio_sorted = ratio.copy()
ratio_sorted['_sort_key'] = ratio_sorted['stage'].map(parse_stage_key)
ratio_sorted = ratio_sorted.sort_values('_sort_key')
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
ax = axes[0]
# scatter 用所有 stage（字母序无所谓，scatter 没有 tick label）
stage_order_vivo = sorted(vivo['stage'].dropna().unique(), key=parse_stage_key)
cmap_m = plt.cm.tab20
for i, st in enumerate(stage_order_vivo):
    sub = vivo[vivo['stage'] == st]
    ax.scatter(sub['glycolysis'], sub['oxphos'], s=8, alpha=0.5,
               c=[cmap_m(i % 20)], label=st)
ax.set_xlabel('Glycolysis score')
ax.set_ylabel('OxPhos score')
ax.set_title('Metabolic scores (in vivo)')
ax.legend(fontsize=5.5, ncol=2, loc='best')
ax = axes[1]
r = ratio_sorted
x_pos = np.arange(len(r))
ax.plot(x_pos, r['gly_mean'], 'o-', label='Glycolysis')
ax.plot(x_pos, r['ox_mean'], 's-', label='OxPhos')
ax.set_xticks(x_pos)
ax.set_xticklabels(r['stage'].tolist(), fontsize=7, rotation=45, ha='right')
ax.set_xlabel('Stage')
ax.set_ylabel('Mean score')
ax.set_title('Scores by stage')
ax.legend(fontsize=7)
ax = axes[2]
ax.bar(x_pos, r['ratio'], color='#4DAF4A')
ax.set_xticks(x_pos)
ax.set_xticklabels(r['stage'].tolist(), fontsize=7, rotation=45, ha='right')
ax.set_xlabel('Stage')
ax.set_ylabel('Gly/Ox ratio')
ax.set_title('Glycolysis:OxPhos ratio')
fig.tight_layout()
save(fig, 'supp_fig_metabolism')

# ============ S8 mofa_rna_only（解析 stage 排序 + F1 dist 边距） ============
print("S8 mofa_rna_only")
rf9 = pd.read_csv(f"{P}/m8_mofa_factors.csv", index_col=0)
r2 = pd.read_csv(f"{P}/m8_mofa_r2.csv")
fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
ax = axes[0, 0]
r2_sorted = r2.sort_values('factor')
ax.bar(r2_sorted['factor'], r2_sorted['R2_pct'], color='#377EB8')
ax.set_xlabel('Factor')
ax.set_ylabel('R² (%)')
ax.set_title('Variance explained')
ax = axes[0, 1]
stage_col = 'stage_label' if 'stage_label' in rf9.columns else 'dataset'
if stage_col in rf9.columns:
    # 解析排序 stage
    stage_order_unique = sorted(rf9[stage_col].dropna().unique(), key=parse_stage_key)
    for f in ['F1', 'F2', 'F3']:
        tmp = rf9.groupby(stage_col)[f].mean().reindex(stage_order_unique)
        ax.plot(np.arange(len(tmp)), tmp.values, 'o-', label=f)
    ax.set_xticks(np.arange(len(stage_order_unique)))
    ax.set_xticklabels(stage_order_unique, fontsize=6, rotation=45, ha='right')
ax.set_xlabel('Stage')
ax.set_ylabel('Mean factor score')
ax.set_title('Factor by stage')
ax.legend(fontsize=7)
ax = axes[1, 0]
ax.scatter(rf9['F1'], rf9['F2'], c='#4DAF4A', s=8, alpha=0.5)
ax.set_xlabel('F1')
ax.set_ylabel('F2')
ax.set_title('F1 vs F2')
ax = axes[1, 1]
ax.hist(rf9['F1'].dropna(), bins=40, color='#E41A1C', edgecolor='black', lw=0.3)
ax.set_xlabel('F1 score')
ax.set_ylabel('Cells')
# 标题上加 padding 防止和 x 轴重叠
ax.set_title('F1 distribution', pad=15)
fig.tight_layout()
save(fig, 'supp_fig_mofa_rna_only')

# ============ S9 overview（disable 1e7 + gene symbol） ============
print("S9 overview")
otypes = pd.read_csv(f"{P}/m4_oocyte_types.csv")
odegs = pd.read_csv(f"{P}/m4_oocyte_type_degs.csv")
# 把 odegs gene 列映射到 symbol
odegs['gene_sym'] = odegs['gene'].astype(str).map(to_symbol)
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
ax = axes[0]
for t in sorted(otypes['oocyte_type'].dropna().unique()):
    sub = otypes[otypes['oocyte_type'] == t]
    ax.hist(sub['total_transcripts'], bins=30, alpha=0.5, label=t)
ax.set_xlabel('Total transcripts')
ax.set_ylabel('Cells')
ax.set_title('Transcript count by type')
ax.ticklabel_format(style='plain', axis='x')  # 禁用 1e7
# set_xticks 用 plain 数值
ax.set_xticks(np.linspace(otypes['total_transcripts'].min(),
                           otypes['total_transcripts'].max(), 5))
ax.tick_params(axis='x', labelsize=8, rotation=15)
ax.legend(fontsize=7)
ax = axes[1]
lfc_o = odegs['log2FC_II_vs_I'].dropna().values
pv_o = odegs['padj'].dropna().values
sig_o = (pv_o < 0.05) & (np.abs(lfc_o) > 0.5)
ax.scatter(lfc_o[~sig_o], -np.log10(pv_o[~sig_o].clip(1e-300)), c='gray', s=4, alpha=0.3)
ax.scatter(lfc_o[sig_o], -np.log10(pv_o[sig_o].clip(1e-300)), c='red', s=7, alpha=0.5)
ax.set_xlabel('log2FC (II vs I)')
ax.set_ylabel('-log10(padj)')
ax.set_title('Oocyte type DEGs')
ax = axes[2]
top_deg = odegs[odegs['padj'] < 0.05].nsmallest(15, 'log2FC_II_vs_I')
top_deg['gene_sym'] = top_deg['gene'].astype(str).map(to_symbol)
ax.barh(range(len(top_deg)), top_deg['log2FC_II_vs_I'].values, color='#984EA3')
ax.set_yticks(range(len(top_deg)))
ax.set_yticklabels(top_deg['gene_sym'].tolist(), fontsize=7)
ax.set_xlabel('log2FC')
ax.invert_yaxis()
ax.set_title('Top DEGs')
fig.tight_layout()
save(fig, 'supp_fig_overview')

# ============ S10 regulon by stage（解析排序） ============
print("S10 regulon by stage")
aucell = pd.read_csv(f"{P}/pyscenic_out/aucell_scores_full.csv", index_col=0)
sm = pd.read_csv(f"{P}/m7_cell_stage_map.csv")
top_reg = aucell.var().abs().sort_values(ascending=False).head(15).index
stage_col = 'stage_label' if 'stage_label' in sm.columns else 'stage_order'
sm = sm.set_index('cell_id')
aucell = aucell.loc[aucell.index.isin(sm.index)]
sm = sm.loc[aucell.index]
stage_order_unique = sorted(sm[stage_col].dropna().unique(), key=parse_stage_key)
fig, ax = plt.subplots(figsize=(10.5, 5.5))
for reg in top_reg:
    tmp = aucell.groupby(sm[stage_col])[reg].mean().reindex(stage_order_unique)
    ax.plot(np.arange(len(tmp)), tmp.values, 'o-', markersize=4, lw=1.2,
            label=reg.split('(')[0])
ax.set_xticks(np.arange(len(stage_order_unique)))
ax.set_xticklabels(stage_order_unique, fontsize=7, rotation=45, ha='right')
ax.set_ylabel('Mean AUC')
ax.set_title('Regulon activity by stage (top 15)')
ax.legend(fontsize=6, ncol=3, loc='best')
fig.tight_layout()
save(fig, 'supp_fig_regulon_by_stage')

# ============ S13 trajectory divergence（gene symbol） ============
print("S13 trajectory divergence")
gdiv = pd.read_csv(f"{P}/m5_gene_divergence_v2.csv")
top50 = pd.read_csv(f"{P}/m5_top50_divergent_genes.csv")
top50['gene_sym'] = top50['gene'].astype(str).map(to_symbol)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
ax.hist(gdiv['total_div'].clip(0, 5), bins=60, color='#377EB8', edgecolor='black', lw=0.3)
ax.set_xlabel('Total divergence')
ax.set_ylabel('Genes')
ax.set_title('Gene-wise divergence')
ax = axes[1]
t = top50.head(15)
x = np.arange(len(t))
ax.barh(x, t['IVF_vs_PA'].values[::-1], color='#E41A1C', alpha=0.7, label='IVF vs PA')
ax.set_yticks(x)
ax.set_yticklabels(t['gene_sym'].values[::-1], fontsize=7)
ax.set_xlabel('Divergence')
ax.set_title('Top divergent genes')
ax.legend(fontsize=7)
fig.tight_layout()
save(fig, 'supp_fig_trajectory_divergence')

# ============ S1 已正确（之前美化过），S2/S4/S5/S6/S11/S12/S14 已正确，不动 ============
print("\n=== 6 处修复完成（S3/S7/S8/S9/S10/S13）===")
print("其余 8 个补充图已正确，跳过")