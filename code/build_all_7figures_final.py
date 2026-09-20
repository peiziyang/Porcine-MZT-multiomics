#!/usr/bin/env python
"""
Build all 7 main figures for BOR submission.
Design rules:
- NO overlapping sub-panels (use GridSpec, never inset_axes)
- Arial font, SVG fonttype='none' (AI editable)
- Correct data column names verified from survey
- Panel labels A/B/C... in top-left
- 300 dpi PNG backup for each SVG
"""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial',
    'svg.fonttype': 'none',
    'savefig.dpi': 300,
    'font.size': 7,
    'axes.titlesize': 8,
    'axes.labelsize': 7,
})
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import FancyBboxPatch, Circle
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
from adjustText import adjust_text
import os, warnings, gzip, re
warnings.filterwarnings('ignore')

# ── Paths ──
BASE  = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
MF    = os.path.join(BASE, 'm4_mofa')
PERT  = os.path.join(BASE, 'perturbation')
DESEQ = os.path.join(BASE, 'deseq2')
SCEN  = os.path.join(BASE, 'pyscenic_out')
OUT   = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures'
os.makedirs(OUT, exist_ok=True)

# ── Load all data ──
meta    = pd.read_csv(os.path.join(BASE, 'm2_metadata.csv'))
umap_df = pd.read_csv(os.path.join(BASE, 'm2_harmony_umap.csv'))
factor  = pd.read_csv(os.path.join(MF, 'mofa_multiomics_factors_v2.csv'), index_col=0)
r2_df   = pd.read_csv(os.path.join(MF, 'mofa_multiomics_r2_per_view_v2.csv'))
cross_v = pd.read_csv(os.path.join(MF, 'factor_cross_view_v2.csv'))
wrna    = pd.read_csv(os.path.join(MF, 'mofa_multiomics_weights_RNA_v2.csv'), index_col=0)
wmeth   = pd.read_csv(os.path.join(MF, 'mofa_multiomics_weights_METH_v2.csv'), index_col=0)
ortho   = pd.read_csv(os.path.join(MF, 'cross_species_orthologs.csv'))
stg_m   = pd.read_csv(os.path.join(MF, 'mzt_stage_mean_projection.csv'), index_col=0)
tproj   = pd.read_csv(os.path.join(MF, 'mzt_trajectory_projection.csv'), index_col=0)
pert_v  = pd.read_csv(os.path.join(PERT, 'perturbation_vectors.csv'))
pert_pa = pd.read_csv(os.path.join(PERT, 'perturbation_vs_pa.csv'))
div_df  = pd.read_csv(os.path.join(BASE, 'm5_gene_divergence_v2.csv'))
met     = pd.read_csv(os.path.join(BASE, 'm6_metabolism_ratio.csv'))
# DESeq2 overall
de_overall = pd.read_csv(os.path.join(DESEQ, 'm15_deseq2_overall_IVFvsPA_stageAdj.csv'))
# Volcano data
m3_de = pd.read_csv(os.path.join(BASE, 'm3_ivf_vs_pa_results.csv'))

# ── Helper ──
def panel_label(ax, label, x=0.02, y=0.97, fontsize=6):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize,
            fontweight='bold', va='top', ha='left')

FCOLS = ['F1','F2','F3','F4','F5','F6','F7']

# ═══════════════════════════════════════════════════════════════════
# FIGURE 1: Atlas UMAP + MOFA+ Multi-omics (3x2 = 6 panels)
# ═══════════════════════════════════════════════════════════════════
print("=== Fig1: Atlas UMAP + MOFA+ ===")
fig1 = plt.figure(figsize=(5.8, 9.0))
gs1 = fig1.add_gridspec(4, 2, hspace=0.35, wspace=0.3,
                        top=0.97, bottom=0.03, left=0.06, right=0.98)

# ── A: Atlas UMAP ──
ax = fig1.add_subplot(gs1[0, 0])
panel_label(ax, 'A')
# Merge metadata with UMAP
meta_umap = meta.copy()
meta_umap['UMAP1'] = umap_df['UMAP1_harmony'].values
meta_umap['UMAP2'] = umap_df['UMAP2_harmony'].values
stage_order = sorted(meta_umap['stage_group'].dropna().unique(), key=str)
cmap = plt.cm.tab10
for i, st in enumerate(stage_order):
    sub = meta_umap[meta_umap['stage_group'] == st]
    if len(sub):
        ax.scatter(sub['UMAP1'], sub['UMAP2'], s=8, alpha=0.6,
                   c=[cmap(i % 10)], label=st, rasterized=True)
ax.set_title('Porcine preimplantation atlas (1,955 cells)', fontsize=8, fontweight='bold')
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.legend(loc='lower right', fontsize=6, ncol=2, markerscale=2, title='Stage', framealpha=0.85)
ax.tick_params(labelsize=7)

# ── B: Factor scores heatmap ──
ax = fig1.add_subplot(gs1[0, 1])
panel_label(ax, 'B')
fmat = factor[FCOLS].values.T
im = ax.imshow(fmat, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
ax.set_yticks(range(7)); ax.set_yticklabels(FCOLS, fontsize=6)
ax.set_xticks([]); ax.set_xlabel('32 GV oocytes', fontsize=7)
ax.set_title('Factor scores (z-scaled)', fontsize=8, fontweight='bold')
cbar = fig1.colorbar(im, ax=ax, shrink=0.8, aspect=20)
cbar.ax.tick_params(labelsize=7)

# ── C: Per-view variance explained ──
ax = fig1.add_subplot(gs1[1, 0])
panel_label(ax, 'C')
x = np.arange(7)
w = 0.35
r2p = r2_df.pivot(index='factor', columns='view', values='r2_pct').fillna(0)
rna_vals = [r2p.loc[f, 'RNA'] if f in r2p.index else 0 for f in FCOLS]
meth_vals = [r2p.loc[f, 'METH'] if f in r2p.index else 0 for f in FCOLS]
ax.bar(x - w/2, meth_vals, w, color='#FF7F00', edgecolor='black', lw=0.5, label='METH')
ax.bar(x + w/2, rna_vals, w, color='#377EB8', edgecolor='black', lw=0.5, label='RNA')
for i in range(7):
    if meth_vals[i] > 0.5:
        ax.text(i - w/2, meth_vals[i] + 1.5, f'{meth_vals[i]:.1f}%', ha='center', fontsize=6)
    if rna_vals[i] > 0.5:
        ax.text(i + w/2, rna_vals[i] + 1.5, f'{rna_vals[i]:.1f}%', ha='center', fontsize=6)
ax.set_xticks(x); ax.set_xticklabels(FCOLS, fontsize=6)
ax.set_ylabel('Variance contribution (%)', fontsize=7)
ax.set_title('Per-view variance contribution', fontsize=8, fontweight='bold')
ax.legend(fontsize=6, loc='upper right')
ax.set_ylim(0, max(max(rna_vals), max(meth_vals)) * 1.3)

# ── D: Cross-view weight correlation ──
ax = fig1.add_subplot(gs1[1, 1])
panel_label(ax, 'D')
r_vals = cross_v['rna_meth_corr_r'].values
colors_d = ['#E41A1C' if r > 0.5 else '#377EB8' if r < -0.3 else '#999999' for r in r_vals]
ax.bar(FCOLS, r_vals, color=colors_d, edgecolor='black', lw=0.5, width=0.6)
ax.axhline(0, color='gray', lw=0.8)
for i, r in enumerate(r_vals):
    ax.text(i, r + 0.04 if r >= 0 else r - 0.06,
            f'r={r:+.2f}', ha='center', fontsize=6, fontweight='bold')
ax.set_ylabel('RNA-METH weight Pearson r', fontsize=7)
ax.set_title('Cross-view weight correlation', fontsize=8, fontweight='bold')
ax.set_ylim(-0.15, 1.05)

# ── E: Factor scores by donor ──
ax = fig1.add_subplot(gs1[2, 0])
panel_label(ax, 'E')
donor_map = meta.set_index('cell')['donor'].to_dict()
donor_colors = {'AF1': '#E41A1C', 'AF2': '#377EB8', 'AF3': '#4DAF4A',
                'AF4': '#984EA3', 'AF5': '#FF7F00'}
offsets = np.linspace(-0.2, 0.2, 5)
for fi in range(7):
    for di, donor in enumerate(['AF1','AF2','AF3','AF4','AF5']):
        vals = []
        for c in factor.index:
            d = donor_map.get(c, 'unknown')
            if d == donor:
                vals.append(factor.loc[c, FCOLS[fi]])
        if vals:
            jitter = np.random.RandomState(42+fi*5+di).normal(0, 0.02, len(vals))
            ax.scatter(np.full(len(vals), fi+1+offsets[di]) + jitter, vals,
                      c=donor_colors[donor], s=10, alpha=0.6, edgecolors='none',
                      label=donor if fi == 0 else '')
ax.set_xticks(range(1, 8)); ax.set_xticklabels(FCOLS, fontsize=6)
ax.axhline(0, color='gray', ls='--', lw=0.5)
ax.set_xlabel('Factor', fontsize=7); ax.set_ylabel('Score', fontsize=7)
ax.set_title('Factor scores by donor', fontsize=8, fontweight='bold')
ax.legend(fontsize=6, loc='upper right', ncol=5, markerscale=0.8)
ax.tick_params(labelsize=7)

# ── F: F1 RNA vs METH weight scatter ──
ax = fig1.add_subplot(gs1[2, 1])
panel_label(ax, 'F')
r1, p1 = pearsonr(wrna['F1'].values, wmeth['F1'].values)
ax.scatter(wrna['F1'].values, wmeth['F1'].values, c='gray', alpha=0.2, s=3)
# Highlight top shared genes (top 6 only to avoid label crowding)
shared_idx = np.argsort(np.abs(wrna['F1'].values) + np.abs(wmeth['F1'].values)[::-1].argsort())[-6:] if False else np.argsort(np.abs(wrna['F1'].values) + np.abs(wmeth['F1'].values))[-6:]
texts = []
for idx in shared_idx:
    texts.append(ax.annotate(wrna.index[idx], (wrna['F1'].iloc[idx], wmeth['F1'].iloc[idx]),
               fontsize=7, fontweight='bold', color='#E41A1C'))
adjust_text(texts, ax=ax, force_text=0.5, force_points=0.3, force_static=0.2,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.4))
ax.set_xlabel('F1 RNA weight', fontsize=7); ax.set_ylabel('F1 METH weight', fontsize=7)
ax.set_title(f'F1: RNA vs METH weights (r={r1:.3f}, p={p1:.1e})', fontsize=8, fontweight='bold')
ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
ax.tick_params(labelsize=7)

# ── G: GO enrichment (top 3 terms per factor) ──
ax = fig1.add_subplot(gs1[3, :])
panel_label(ax, 'G')
go = pd.read_csv('E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa/factor_go_enrichment_v2.csv')
go['log10p'] = -np.log10(go['p_raw'].clip(lower=1e-300))
top = go.groupby('factor').apply(lambda x: x.nsmallest(3, 'p_raw')).reset_index(drop=True)
colors = {'F1':'#377EB8','F2':'#FF7F00','F3':'#E41A1C','F4':'#4DAF4A','F5':'#984EA3','F6':'#A65628','F7':'#F781BF'}
for f in sorted(top['factor'].unique()):
    sub = top[top['factor']==f]
    ax.barh(sub.index[::-1], sub['log10p'][::-1], color=colors.get(f,'#999'), label=f)
terms = [f"{r['factor']} | {str(r['GO_term'])[:48]}" for _, r in top.iterrows()]
ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms, fontsize=6)
ax.set_xlabel('-log10(raw p)', fontsize=7); ax.invert_yaxis()
ax.set_title('GO enrichment (top 3 terms per factor)', fontsize=8, fontweight='bold')
ax.legend(fontsize=6, ncol=7, loc='lower center', bbox_to_anchor=(0.5, -0.18),
           frameon=True, framealpha=0.9, edgecolor='gray')
ax.tick_params(labelsize=6)

fig1.savefig(os.path.join(OUT, 'fig1_mofa.svg'), format='svg', bbox_inches='tight')
fig1.savefig(os.path.join(OUT, 'fig1_mofa.png'), dpi=300, bbox_inches='tight')
fig1.savefig(os.path.join(OUT, 'fig1_mofa.pdf'), format='pdf', bbox_inches='tight')
    
plt.close(fig1)
print("  -> fig1_mofa.svg OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 2: Cross-species MZT conservation (2x2 = 4 panels)
# ═══════════════════════════════════════════════════════════════════
print("=== Fig2: Cross-species ===")
fig2, axes2 = plt.subplots(2, 2, figsize=(5.8, 4.8))
axes2 = axes2.flatten()

# ── A: F1 gene species overlap ──
ax = axes2[0]; panel_label(ax, 'A')
pig_only = ortho['pig_gene'].nunique()
in_human = ortho[ortho['in_human']].shape[0]
in_mouse = ortho[ortho['in_mouse']].shape[0]
in_both  = ortho[ortho['in_human'] & ortho['in_mouse']].shape[0]
human_only = in_human - in_both
mouse_only = in_mouse - in_both
pig_excl = pig_only - in_human - in_mouse + in_both

cats = ['Pig+Human\n+Mouse', 'Pig+Human\nonly', 'Pig+Mouse\nonly', 'Pig only']
vals = [in_both, human_only, mouse_only, pig_excl]
cc = ['#377EB8', '#4DAF4A', '#FF7F00', '#999999']
ax.bar(cats, vals, color=cc, edgecolor='black', lw=0.5, width=0.6)
for i, v in enumerate(vals):
    ax.text(i, v + 0.5, str(v), ha='center', fontsize=6, fontweight='bold')
ax.set_ylabel('Number of F1 genes', fontsize=7)
ax.set_title(f'F1 maternal blueprint genes (n={pig_only})', fontsize=8, fontweight='bold')

# ── B: Pig vs Human scatter ──
ax = axes2[1]; panel_label(ax, 'B')
hu_genes = ortho[ortho['in_human']].dropna(subset=['pig_GV_mean','human_mean'])
r_hu, p_hu = spearmanr(hu_genes['pig_GV_mean'], hu_genes['human_mean'])
ax.scatter(hu_genes['pig_GV_mean'], hu_genes['human_mean'],
           c='#377EB8', alpha=0.6, s=25, edgecolors='black', lw=0.3)
texts_b = []
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    sub = hu_genes[hu_genes['pig_gene'] == g]
    if len(sub):
        texts_b.append(ax.annotate(g, (sub['pig_GV_mean'].values[0], sub['human_mean'].values[0]),
                   fontsize=6, fontweight='bold'))
adjust_text(texts_b, ax=ax, force_text=0.3, force_points=0.2,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.4))
ax.set_xlabel('Pig GV mean expression', fontsize=7)
ax.set_ylabel('Human preimplantation mean expr.', fontsize=7)
ax.set_title(f'Pig vs Human (rho={r_hu:.3f}, p={p_hu:.2e})', fontsize=8, fontweight='bold')
ax.tick_params(labelsize=7)

# ── C: Pig vs Mouse scatter ──
ax = axes2[2]; panel_label(ax, 'C')
mo_genes = ortho[ortho['in_mouse']].dropna(subset=['pig_GV_mean','mouse_mean'])
r_mo, p_mo = spearmanr(mo_genes['pig_GV_mean'], mo_genes['mouse_mean'])
ax.scatter(mo_genes['pig_GV_mean'], mo_genes['mouse_mean'],
           c='#FF7F00', alpha=0.6, s=25, edgecolors='black', lw=0.3)
texts_c = []
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    sub = mo_genes[mo_genes['pig_gene'] == g]
    if len(sub):
        texts_c.append(ax.annotate(g, (sub['pig_GV_mean'].values[0], sub['mouse_mean'].values[0]),
                   fontsize=6, fontweight='bold'))
adjust_text(texts_c, ax=ax, force_text=0.3, force_points=0.2,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.4))
ax.set_xlabel('Pig GV mean expression', fontsize=7)
ax.set_ylabel('Mouse preimplantation mean expr.', fontsize=7)
ax.set_title(f'Pig vs Mouse (rho={r_mo:.3f}, p={p_mo:.2e})', fontsize=8, fontweight='bold')
ax.tick_params(labelsize=7)

# ── D: Summary ──
ax = axes2[3]; panel_label(ax, 'D'); ax.axis('off')
ax.set_xlim(0, 10); ax.set_ylim(0, 10)

ax.text(5, 9.6, 'Cross-species summary', ha='center', va='top',
        fontsize=6, fontweight='bold', color='#333333')

# Section 1: gene matching (mini bars)
ax.text(0.4, 8.5, 'Symbol-matched gene pairs (n = 44 F1 genes)', fontsize=6,
        fontweight='bold', color='#377EB8')
match_items = [('Human', in_human, '#4DAF4A'), ('Mouse', in_mouse, '#FF7F00'),
               ('Both species', in_both, '#E41A1C')]
for i, (label, val, c) in enumerate(match_items):
    y = 7.7 - i * 0.8
    ax.text(0.7, y, label, fontsize=6, va='center', color='#333333')
    ax.barh(y, val / pig_only * 3.6, left=4.6, height=0.42, color=c, alpha=0.75, zorder=2)
    ax.text(8.4, y, f'{val}/44', fontsize=6, va='center', fontweight='bold', color=c)

# Section 2: correlation
ax.text(0.4, 5.0, 'Expression correlation (Spearman)', fontsize=6,
        fontweight='bold', color='#4DAF4A')
ax.text(0.7, 4.2, f'Pig–Human:  ρ = {r_hu:.2f}   (nominal P = {p_hu:.2e})', fontsize=6, color='#333333')
ax.text(0.7, 3.5, f'Pig–Mouse:  ρ = {r_mo:.2f}   (nominal P = {p_mo:.2e})', fontsize=6, color='#333333')

# Section 3: take-home
ax.text(5, 1.9, 'F1 correlation is not specific\nrelative to the symbol-matched background',
        ha='center', va='center', fontsize=6, style='italic', color='#555555',
        bbox=dict(boxstyle='round,pad=0.55', facecolor='#F2F2F2', edgecolor='#CCCCCC'))

fig2.tight_layout()
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.svg'), format='svg', bbox_inches='tight')
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.png'), dpi=300, bbox_inches='tight')
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.pdf'), format='pdf', bbox_inches='tight')
    
plt.close(fig2)
print("  -> fig2_crossspecies.svg OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 3: MZT trajectory + OT (3x2 = 6 panels)
# ═══════════════════════════════════════════════════════════════════
print("=== Fig3: MZT trajectory ===")
fig3, axes3 = plt.subplots(3, 2, figsize=(5.8, 6.6))
axes3 = axes3.flatten()

stage_order = ['E0','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10']
stage_num = np.arange(len(stage_order))
stage_colors = {'F1':'#377EB8','F3':'#E41A1C','F4':'#4DAF4A','F6':'#FF7F00'}

# ── A: F3 maternal clearance ──
ax = axes3[0]; panel_label(ax, 'A')
f3_vals = [stg_m.loc[s, 'F3'] if s in stg_m.index else np.nan for s in stage_order]
r_f3, p_f3 = spearmanr(stage_num[~np.isnan(f3_vals)], np.array(f3_vals)[~np.isnan(f3_vals)])
ax.plot(stage_num, f3_vals, '-o', color='#E41A1C', lw=2, markersize=6, label=f'rho={r_f3:.3f}')
ax.fill_between(stage_num, np.array(f3_vals)-np.std(f3_vals)*0.5,
                np.array(f3_vals)+np.std(f3_vals)*0.5, alpha=0.15, color='#E41A1C')
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=6)
ax.set_ylabel('Projected score', fontsize=7)
ax.set_title('F3: stage-associated decline', fontsize=8, fontweight='bold', color='#E41A1C')
ax.tick_params(labelsize=7); ax.legend(fontsize=6)

# ── B: F4 zygotic activation ──
ax = axes3[1]; panel_label(ax, 'B')
f4_vals = [stg_m.loc[s, 'F4'] if s in stg_m.index else np.nan for s in stage_order]
r_f4, p_f4 = spearmanr(stage_num[~np.isnan(f4_vals)], np.array(f4_vals)[~np.isnan(f4_vals)])
ax.plot(stage_num, f4_vals, '-o', color='#4DAF4A', lw=2, markersize=6, label=f'rho={r_f4:.3f}')
ax.fill_between(stage_num, np.array(f4_vals)-np.std(f4_vals)*0.5,
                np.array(f4_vals)+np.std(f4_vals)*0.5, alpha=0.15, color='#4DAF4A')
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=6)
ax.set_ylabel('Projected score', fontsize=7)
ax.set_title('F4: nominal stage-associated rise (descriptive)', fontsize=8, fontweight='bold', color='#4DAF4A')
ax.tick_params(labelsize=7); ax.legend(fontsize=6)

# ── C: F6 zygotic activation ──
ax = axes3[2]; panel_label(ax, 'C')
f6_vals = [stg_m.loc[s, 'F6'] if s in stg_m.index else np.nan for s in stage_order]
r_f6, p_f6 = spearmanr(stage_num[~np.isnan(f6_vals)], np.array(f6_vals)[~np.isnan(f6_vals)])
ax.plot(stage_num, f6_vals, '-o', color='#FF7F00', lw=2, markersize=6, label=f'rho={r_f6:.3f}')
ax.fill_between(stage_num, np.array(f6_vals)-np.std(f6_vals)*0.5,
                np.array(f6_vals)+np.std(f6_vals)*0.5, alpha=0.15, color='#FF7F00')
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=6)
ax.set_ylabel('Projected score', fontsize=7)
ax.set_title('F6: nominal stage-associated rise (descriptive)', fontsize=8, fontweight='bold', color='#FF7F00')
ax.tick_params(labelsize=7); ax.legend(fontsize=6)

# ── D: F1 stable ──
ax = axes3[3]; panel_label(ax, 'D')
f1_vals = [stg_m.loc[s, 'F1'] if s in stg_m.index else np.nan for s in stage_order]
r_f1, p_f1 = spearmanr(stage_num[~np.isnan(f1_vals)], np.array(f1_vals)[~np.isnan(f1_vals)])
ax.plot(stage_num, f1_vals, '-o', color='#377EB8', lw=2, markersize=6, label=f'rho={r_f1:.3f}')
ax.fill_between(stage_num, np.array(f1_vals)-np.std(f1_vals)*0.5,
                np.array(f1_vals)+np.std(f1_vals)*0.5, alpha=0.15, color='#377EB8')
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=6)
ax.set_ylabel('Projected score', fontsize=7)
ax.set_title('F1: modest negative trend (descriptive)', fontsize=8, fontweight='bold', color='#377EB8')
ax.tick_params(labelsize=7); ax.legend(fontsize=6)

# ── E: All factors trajectory overview ──
ax = axes3[4]; panel_label(ax, 'E')
for fc in ['F1','F3','F4','F6']:
    vals = [stg_m.loc[s, fc] if s in stg_m.index else np.nan for s in stage_order]
    vals_arr = np.array(vals)
    vals_norm = (vals_arr - vals_arr.mean()) / vals_arr.std() if vals_arr.std() > 0 else vals_arr
    ax.plot(stage_num, vals_norm, '-o', lw=1.5, markersize=4, color=stage_colors[fc], label=fc, alpha=0.8)
ax.set_xticks(stage_num[::2]); ax.set_xticklabels(stage_order[::2], fontsize=6)
ax.set_ylabel('Normalized score', fontsize=7)
ax.set_title('All factor trajectories (z-normalized)', fontsize=8, fontweight='bold')
ax.axhline(0, color='gray', lw=0.5); ax.legend(fontsize=6); ax.tick_params(labelsize=7)

# ── F: Per-cell projection density (F3) ──
ax = axes3[5]; panel_label(ax, 'F')
stages_in_data = sorted(tproj['stage'].unique(), key=lambda s: int(s[1:]))
stage_idx_map = {s: i for i, s in enumerate(stages_in_data)}
tproj_sorted = tproj.copy()
tproj_sorted['stage_order'] = tproj_sorted['stage'].map(stage_idx_map)
for s in stages_in_data:
    sub = tproj_sorted[tproj_sorted['stage'] == s]
    if len(sub) > 0:
        ax.violinplot(sub['F1'].values, positions=[stage_idx_map[s]],
                     showmeans=True, showmedians=False, widths=0.7)
ax.set_xticks(range(len(stages_in_data)))
ax.set_xticklabels(stages_in_data, fontsize=6, rotation=45, ha='right')
ax.set_ylabel('F1 projected score', fontsize=7)
ax.set_title('Per-cell F1 projection (n=832 E0-E10 cells)', fontsize=8, fontweight='bold')
ax.tick_params(labelsize=7)

fig3.tight_layout()
fig3.savefig(os.path.join(OUT, 'fig3_mzt_trajectory.svg'), format='svg', bbox_inches='tight')
fig3.savefig(os.path.join(OUT, 'fig3_mzt_trajectory.png'), dpi=300, bbox_inches='tight')
fig3.savefig(os.path.join(OUT, 'fig3_mzt_trajectory.pdf'), format='pdf', bbox_inches='tight')
    
plt.close(fig3)
print("  -> fig3_mzt_trajectory.svg OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 4: In silico TF perturbation (2x2 = 4 panels)
# ═══════════════════════════════════════════════════════════════════
print("=== Fig4: Perturbation ===")
fig4, axes4 = plt.subplots(2, 2, figsize=(5.8, 4.8))
axes4 = axes4.flatten()

# ── A: DNMT1 KD F1-F2 shift ──
ax = axes4[0]; panel_label(ax, 'A')
dnmt1 = pert_v[pert_v['tf'] == 'DNMT1']
for label, color, marker in [('100%_KD','#E41A1C','o'),('50%_KD','#FF7F00','s'),
                              ('50%_OE','#377EB8','^'),('100%_OE','#4DAF4A','D')]:
    sub = dnmt1[dnmt1['perturbation'] == label]
    ax.scatter(sub['shift_F1'], sub['shift_F2'], c=color, marker=marker, s=20, alpha=0.7,
              label=label, edgecolors='black', lw=0.3)
ax.axhline(0, color='gray', ls='--', lw=0.5); ax.axvline(0, color='gray', ls='--', lw=0.5)
ax.set_xlabel('Shift in F1', fontsize=7); ax.set_ylabel('Shift in F2', fontsize=7)
ax.set_title('DNMT1 perturbation in MOFA+ space', fontsize=8, fontweight='bold')
ax.legend(fontsize=6); ax.tick_params(labelsize=7)

# ── B: Cosine similarity with PA ──
ax = axes4[1]; panel_label(ax, 'B')
kd_cos = pert_pa[pert_pa['perturbation'].str.contains('KD')]
kd_avg = kd_cos.groupby('tf')['cos_sim'].mean().sort_values()
colors_b = ['#E41A1C' if v > 0.3 else '#FF7F00' if v > 0 else '#999999' for v in kd_avg.values]
ax.barh(range(len(kd_avg)), kd_avg.values, color=colors_b, edgecolor='black', lw=0.5, height=0.6)
ax.set_yticks(range(len(kd_avg))); ax.set_yticklabels(kd_avg.index, fontsize=6)
ax.axvline(0, color='gray', lw=0.8)
for i, v in enumerate(kd_avg.values):
    ax.text(v + 0.02, i, f'{v:+.3f}', fontsize=6, va='center', fontweight='bold')
ax.set_xlabel('Cosine similarity with PA deviation', fontsize=7)
ax.set_title('TF KD alignment with PA deviation', fontsize=8, fontweight='bold')
ax.tick_params(labelsize=7)

# ── C: Per-factor heatmap ──
ax = axes4[2]; panel_label(ax, 'C')
kds = pert_v[pert_v['perturbation'] == '100%_KD']
tfs = sorted(kds['tf'].unique())
hm = np.zeros((len(tfs), 4))
for ti, tf in enumerate(tfs):
    sub = kds[kds['tf'] == tf]
    for fi in range(4):
        hm[ti, fi] = sub[f'shift_{FCOLS[fi]}'].mean()
im = ax.imshow(hm.T, aspect='auto', cmap='RdBu_r', vmin=-0.025, vmax=0.025)
ax.set_xticks(range(len(tfs))); ax.set_xticklabels(tfs, fontsize=6, rotation=45, ha='right')
ax.set_yticks(range(4)); ax.set_yticklabels(FCOLS[:4], fontsize=6)
# 数值标签：用科学记数法，4×9=36 个格子会挤，只标 |shift|>0.005 的（>20% 范围）
for i in range(4):
    for j in range(len(tfs)):
        v = hm[j, i]
        if abs(v) > 0.005:  # 只标显著值
            ax.text(j, i, f'{v*100:.1f}', ha='center', va='center', fontsize=5,
                   color='white' if abs(v) > 0.015 else 'black')
ax.set_title('100% KD: per-factor mean shift ×100', fontsize=8, fontweight='bold')
cbar = fig4.colorbar(im, ax=ax, shrink=0.7, aspect=20)
cbar.ax.tick_params(labelsize=7)

# ── D: Magnitude comparison (KD vs OE, focused y-range) ──
ax = axes4[3]; panel_label(ax, 'D')
tfs_show = ['DNMT1','ATF3','POU5F1','ESRRA','TFAP2C','GATA4','NFE2L3','XBP1']
x_ticks = np.arange(len(tfs_show))
kd_mag = [pert_v[(pert_v['tf']==t)&(pert_v['perturbation']=='100%_KD')]['shift_magnitude'].mean()
          for t in tfs_show]
oe_mag = [pert_v[(pert_v['tf']==t)&(pert_v['perturbation']=='100%_OE')]['shift_magnitude'].mean()
          for t in tfs_show]
ax.bar(x_ticks - 0.15, kd_mag, 0.3, color='#E41A1C', edgecolor='black', lw=0.3, label='KD')
ax.bar(x_ticks + 0.15, oe_mag, 0.3, color='#377EB8', edgecolor='black', lw=0.3, label='OE')
ax.set_xticks(x_ticks); ax.set_xticklabels(tfs_show, fontsize=6, rotation=45, ha='right')
ax.set_ylabel('Mean |shift|', fontsize=7)
ax.set_ylim(min(min(kd_mag), min(oe_mag)) * 0.998, max(max(kd_mag), max(oe_mag)) * 1.002)
ax.set_title('Perturbation magnitude (zoomed)', fontsize=8, fontweight='bold')
ax.legend(fontsize=6); ax.tick_params(labelsize=7)

fig4.tight_layout()
fig4.savefig(os.path.join(OUT, 'fig4_perturbation.svg'), format='svg', bbox_inches='tight')
fig4.savefig(os.path.join(OUT, 'fig4_perturbation.png'), dpi=300, bbox_inches='tight')
fig4.savefig(os.path.join(OUT, 'fig4_perturbation.pdf'), format='pdf', bbox_inches='tight')
    
plt.close(fig4)
print("  -> fig4_perturbation.svg OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 5: PA validation + 3-layer model (3x2 = 6 panels)
# ═══════════════════════════════════════════════════════════════════
print("=== Fig5: PA validation ===")

# Build PA data
pa_all = []
for f in sorted(os.listdir(DESEQ)):
    if not f.startswith('m15_deseq2_') or not (f.endswith('_IVFvsPA.csv') or f.endswith('_IVFvsPA_stageAdj.csv')):
        continue
    if f.endswith('_IVFvsPA_stageAdj.csv'):
        stage = f.replace('m15_deseq2_', '').replace('_IVFvsPA_stageAdj.csv', '')
    else:
        stage = f.replace('m15_deseq2_', '').replace('_IVFvsPA.csv', '')
    df = pd.read_csv(os.path.join(DESEQ, f))
    if 'gene_symbol' not in df.columns:
        continue
    s2i = {}
    for i, row in df.iterrows():
        s = str(row['gene_symbol'])
        if not s.startswith('ENSSSCG'):
            s2i[s] = i
    for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2',
              'IDH2','PKM','EXOSC9','MDH1','GNL3',
              'DPPA5','NANOG']:
        if g in s2i:
            r = df.iloc[s2i[g]]
            pa_all.append({'gene': g, 'stage': stage,
                          'log2FC': r['log2FoldChange'],
                          'padj': r['padj']})
pa_df = pd.DataFrame(pa_all)
overall = pa_df[pa_df['stage'] == 'overall']

fig5 = plt.figure(figsize=(5.8, 6.6))
gs5 = fig5.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

# ── A: Layer 1 bar ──
ax = fig5.add_subplot(gs5[0, 0])
panel_label(ax, 'A')
l1_genes = ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2']
l1_val = []; l1_padj = []; l1_names = []
for g in l1_genes:
    sub = overall[overall['gene'] == g]
    if len(sub) and not np.isnan(sub.iloc[0]['log2FC']):
        l1_val.append(sub.iloc[0]['log2FC']); l1_padj.append(sub.iloc[0]['padj']); l1_names.append(g)
x1 = np.arange(len(l1_names))
c1 = ['#E41A1C' if v > 0.5 else '#FF7F00' if v > -0.5 else '#377EB8' for v in l1_val]
ax.bar(x1, l1_val, color=c1, edgecolor='black', lw=0.5, width=0.6)
ax.axhline(0, color='gray', lw=0.8)
for i, (v, p) in enumerate(zip(l1_val, l1_padj)):
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    yoff = 0.15 if v >= 0 else -0.25
    ax.text(i, v + yoff, f'{v:+.2f} {sig}', ha='center', fontsize=6, fontweight='bold')
ax.set_xticks(x1); ax.set_xticklabels(l1_names, fontsize=6)
ax.set_ylabel('log2 Fold Change (PA / IVF)', fontsize=7)
ax.set_title('F1 gene set (heterogeneous)', fontsize=8, fontweight='bold', color='#E41A1C')
ax.set_ylim(-1.3, 1.85)
ax.axhspan(-1, 1, alpha=0.03, color='gray')

# ── B: Layer 3 bar ──
ax = fig5.add_subplot(gs5[0, 1])
panel_label(ax, 'B')
l3_genes = ['IDH2','PKM','EXOSC9','MDH1','GNL3']
l3_val = []; l3_padj = []; l3_names = []
for g in l3_genes:
    sub = overall[overall['gene'] == g]
    if len(sub) and not np.isnan(sub.iloc[0]['log2FC']):
        l3_val.append(sub.iloc[0]['log2FC']); l3_padj.append(sub.iloc[0]['padj']); l3_names.append(g)
x3 = np.arange(len(l3_names))
c3 = ['#377EB8' if v < -2 else '#999999' for v in l3_val]
ax.bar(x3, l3_val, color=c3, edgecolor='black', lw=0.5, width=0.6)
ax.axhline(0, color='gray', lw=0.8)
for i, (v, p) in enumerate(zip(l3_val, l3_padj)):
    sig = '***' if p < 1e-6 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    ax.text(i, v + 0.5, f'{v:+.2f} {sig}', ha='center', fontsize=6, fontweight='bold',
           color='#377EB8')
ax.set_xticks(x3); ax.set_xticklabels(l3_names, fontsize=6)
ax.set_ylabel('log2 Fold Change (PA / IVF)', fontsize=7)
ax.set_title('Metabolism/translation (not F3)', fontsize=8, fontweight='bold', color='#377EB8')
ax.axhspan(-1, 1, alpha=0.03, color='gray')

# ── C: Per-stage heatmap ──
ax = fig5.add_subplot(gs5[1, 0])
panel_label(ax, 'C')
stg = ['1cell','2cell','4cell','8cell']
gpl = ['DNMT1','ZP3','GDF9','IDH2','PKM','EXOSC9','GNL3','DPPA5','NANOG']
hm = np.full((len(gpl), len(stg)), np.nan)
for gi, g in enumerate(gpl):
    for si, s in enumerate(stg):
        sub = pa_df[(pa_df['gene'] == g) & (pa_df['stage'] == s)]
        if len(sub) and not np.isnan(sub.iloc[0]['log2FC']):
            hm[gi, si] = sub.iloc[0]['log2FC']
hm_masked = np.ma.masked_where(np.isnan(hm), hm)
im = ax.imshow(hm_masked, aspect='auto', cmap='RdBu_r', vmin=-10, vmax=10)
ax.set_xticks(np.arange(len(stg))); ax.set_xticklabels(stg, fontsize=6)
ax.set_yticks(np.arange(len(gpl))); ax.set_yticklabels(gpl, fontsize=6)
for i in range(len(gpl)):
    for j in range(len(stg)):
        v = hm[i, j]
        if not np.isnan(v):
            ax.text(j, i, f'{v:+.1f}', ha='center', va='center', fontsize=6,
                   fontweight='bold', color='white' if abs(v) > 5 else 'black')
ax.set_title('Per-stage log2FC (PA vs IVF)', fontsize=8, fontweight='bold')
cbar = fig5.colorbar(im, ax=ax, shrink=0.5, aspect=15, pad=0.08)
cbar.ax.tick_params(labelsize=5)
cbar.ax.tick_params(labelsize=7)
# ── D: gene-set aggregate summary (read from DESeq2 stage-adjusted) ──
ax = fig5.add_subplot(gs5[1, 1])
panel_label(ax, 'D')
des2 = pd.read_csv('E:/Workbuddy/2026-07-27-11-58-27/data/processed/deseq2/m15_deseq2_overall_IVFvsPA_stageAdj.csv')
sym2lfc = {}
for _, r in des2.iterrows():
    s = str(r.get('gene_symbol',''))
    if s and not s.startswith('ENSSSCG') and pd.notna(r['log2FoldChange']):
        sym2lfc[s] = r['log2FoldChange']
w_df = pd.read_csv('E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa/mofa_multiomics_weights_RNA_v2.csv', index_col=0)
sets_lfc = {}
for fc in ['F1','F3','F4','F6']:
    top = list(w_df[fc].abs().sort_values(ascending=False).head(44).index)
    mapped = [sym2lfc[g] for g in top if g in sym2lfc]
    if mapped: sets_lfc[fc] = float(np.mean(mapped))
labels = ['F1 set (maternal)', 'F3 set (stage-decreasing)', 'F4 set (stage-increasing)', 'F6 set (translation)']
vals = [sets_lfc.get('F1',0), sets_lfc.get('F3',0), sets_lfc.get('F4',0), sets_lfc.get('F6',0)]
clrs = ['#E41A1C','#E41A1C','#4DAF4A','#377EB8']
y = np.arange(len(labels))
ax.barh(y, vals, color=clrs, edgecolor='black', lw=0.5, height=0.6)
ax.axvline(0, color='gray', lw=0.8)
for i, v in enumerate(vals):
    ax.text(v + (0.05 if v >= 0 else -0.05), i, f'{v:+.2f}', va='center', ha='left' if v >= 0 else 'right', fontsize=6, fontweight='bold')
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6)
ax.set_xlabel('Mean log2FC (PA vs IVF)', fontsize=7)
ax.set_title('Gene-set mean log2FC', fontsize=7, fontweight='bold')
ax.tick_params(labelsize=6)


# ── E: Analytical framework schematic (real top genes) ──
ax = fig5.add_subplot(gs5[2, :])
panel_label(ax, 'E'); ax.set_xlim(0, 21); ax.set_ylim(-1, 11); ax.axis('off')

# Left: Normal development (compact, 3 horizontal boxes)
ax.text(4.5, 10.5, 'NORMAL MZT (atlas)', ha='center', fontsize=7, fontweight='bold', color='#333')
gv_c = Circle((1.3, 9.4), 0.35, facecolor='#FFF3CD', edgecolor='#FF7F00', lw=1.5)
ax.add_patch(gv_c)
ax.text(1.3, 9.4, 'GV', ha='center', va='center', fontsize=6, fontweight='bold')

for y, label, ec, fc, genes, note in [
    (8.5, 'F1 (maternal)', '#377EB8', '#E5F0FF', 'RARRES1 . DNMT1 . ZP3', 'r=+0.84'),
    (7.0, 'F3 (decreasing)', '#E41A1C', '#FFE5E5', 'FKBP10 . ACTA2 . PTGES', 'rho=-0.90'),
    (5.5, 'F4/F6 (increasing)', '#4DAF4A', '#E5FFE5', 'DNAJB9 . GDF9 + GMPR2', 'nominal'),
]:
    fb = FancyBboxPatch((0.2, y-0.85), 8.6, 1.7, boxstyle='round,pad=0.08',
                        facecolor=fc, edgecolor=ec, lw=1.2, alpha=0.7)
    ax.add_patch(fb)
    ax.text(4.5, y+0.45, label, fontsize=6, fontweight='bold', color=ec, ha='center')
    ax.text(4.5, y, genes, fontsize=5.5, ha='center', color=ec)
    ax.text(4.5, y-0.45, note, fontsize=5, ha='center', color='#666', style='italic')

# Middle arrow
ax.annotate('', xy=(13.5, 8), xytext=(9.2, 8),
           arrowprops=dict(arrowstyle='->', color='#666', lw=1.5))
ax.text(11.3, 8.3, 'PA vs IVF', ha='center', fontsize=6, color='#E41A1C', fontweight='bold')

# Right: PA embryos (compact, 3 horizontal boxes)
ax.text(17, 10.5, 'PA EMBRYOS', ha='center', fontsize=7, fontweight='bold', color='#E41A1C')
for y, label, ec, fc, mean_lfc in [
    (8.5, 'F1 (heterogeneous)', '#377EB8', '#E5F0FF', 'mean LFC=-0.78'),
    (7.0, 'F3 (negative shift)', '#E41A1C', '#FFE5E5', 'mean LFC=-3.97**'),
    (5.5, 'F4/F6 (no direction)', '#4DAF4A', '#E5FFE5', 'LFC=-0.17/-0.41'),
]:
    fb2 = FancyBboxPatch((13.7, y-0.85), 6.8, 1.7, boxstyle='round,pad=0.08',
                         facecolor=fc, edgecolor=ec, lw=1.2, alpha=0.7)
    ax.add_patch(fb2)
    ax.text(17, y+0.45, label, fontsize=6, fontweight='bold', color=ec, ha='center')
    ax.text(17, y, mean_lfc, fontsize=5.5, ha='center', color='#333')

ax.text(10, 0.5, 'GV oocyte multi-omics (MOFA+) -> Atlas projection -> 42-sample independent validation (SRP301735)',
        ha='center', fontsize=6, style='italic', color='gray')


fig5.savefig(os.path.join(OUT, 'fig5_pa_validation.svg'), format='svg', bbox_inches='tight')
fig5.savefig(os.path.join(OUT, 'fig5_pa_validation.png'), dpi=300, bbox_inches='tight')
fig5.savefig(os.path.join(OUT, 'fig5_pa_validation.pdf'), format='pdf', bbox_inches='tight')
    
plt.close(fig5)
print("  -> fig5_pa_validation.svg OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 6: SCENIC regulon landscape + metabolism (3 panels vertically)
# ═══════════════════════════════════════════════════════════════════
print("=== Fig6: Regulon + Metabolism ===")

# Load SCENIC data
regulon_activity = pd.read_csv(os.path.join(BASE, 'm7_regulon_activity.csv'))

# Clean regulon names from Unnamed:0 column
reg_names = regulon_activity['Unnamed: 0'].str.extract(r'^(\w+)')[0].fillna('Unknown')
# Get numeric columns (E0-E8 stages)
ra_num_cols = [c for c in regulon_activity.columns
               if c.startswith('E') and c[1:].isdigit()]
regulon_stage = regulon_activity[ra_num_cols].copy()
regulon_stage.index = reg_names.values

# Sort stages
ra_stages = sorted(ra_num_cols, key=lambda s: int(s[1:]))

fig6 = plt.figure(figsize=(5.8, 5.8))
gs6 = fig6.add_gridspec(3, 1, height_ratios=[1.2, 1, 1], hspace=0.4)

# ── A: Regulon heatmap ──
ax = fig6.add_subplot(gs6[0])
panel_label(ax, 'A')
# Top 20 regulons by mean activity
top_20 = regulon_stage[ra_stages].mean(axis=1).sort_values(ascending=False).head(15).index.tolist()
hm_data = regulon_stage.loc[top_20, ra_stages].values
im = ax.imshow(hm_data, aspect='auto', cmap='YlOrRd')
ax.set_xticks(range(len(ra_stages)))
ax.set_xticklabels(ra_stages, fontsize=6, rotation=45, ha='right')
ax.set_yticks(range(len(top_20)))
ax.set_yticklabels([str(t)[:18] for t in top_20], fontsize=5)
ax.set_title('SCENIC regulon activity: top 15', fontsize=8, fontweight='bold')
cbar = fig6.colorbar(im, ax=ax, shrink=0.8, aspect=20)
cbar.ax.tick_params(labelsize=7)

# ── B: Key regulon dynamics ──
ax = fig6.add_subplot(gs6[1])
panel_label(ax, 'B')
key_tfs = ['DNMT1','POU5F1','ESRRA','TFAP2C','GATA4','ATF3','NFE2L3','XBP1']
found_tfs = [t for t in key_tfs if t in regulon_stage.index]
if found_tfs:
    colors_tf = plt.cm.tab10(np.linspace(0, 1, len(found_tfs)))
    for ti, tf in enumerate(found_tfs):
        vals = regulon_stage.loc[tf, ra_stages].values
        ax.plot(range(len(ra_stages)), vals, '-o', lw=1.5, markersize=4,
               color=colors_tf[ti], label=tf)
    ax.set_xticks(range(len(ra_stages)))
    ax.set_xticklabels(ra_stages, fontsize=6, rotation=45, ha='right')
    ax.set_ylabel('Regulon activity', fontsize=7)
    ax.set_title('Key regulon temporal dynamics', fontsize=8, fontweight='bold')
    ax.legend(fontsize=6, ncol=4)
    ax.tick_params(labelsize=7)
else:
    ax.text(0.5, 0.5, 'No key regulons found in data', ha='center', va='center', fontsize=6)

# ── C: Metabolism scoring (3 sub-panels with GridSpecFromSubplotSpec) ──
ax_c = fig6.add_subplot(gs6[2])
panel_label(ax_c, 'C')
gs_c = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs6[2], wspace=0.3)

ax_c.axis('off')

stages_m = met['stage'].values
x_met = np.arange(len(stages_m))

# (i) Ratio
ax1 = fig6.add_subplot(gs_c[0])
ax1.bar(x_met, met['ratio'].values, color='#377EB8', edgecolor='black', lw=0.3, width=0.7)
ax1.axhline(1, color='red', ls='--', lw=0.8, alpha=0.5)
ax1.set_xticks(x_met[::2]); ax1.set_xticklabels(stages_m[::2], fontsize=6, rotation=45, ha='right')
ax1.set_ylabel('Glycolysis / OxPhos', fontsize=7)
ax1.set_title('Glycolysis/OxPhos ratio', fontsize=8, fontweight='bold')
ax1.tick_params(labelsize=6)
ax1.set_ylim(0, max(met['ratio'].values) * 1.3)

# (ii) Glycolysis
ax2 = fig6.add_subplot(gs_c[1])
ax2.bar(x_met, met['gly_mean'].values, color='#E41A1C', edgecolor='black', lw=0.3, width=0.7)
ax2.set_xticks(x_met[::2]); ax2.set_xticklabels(stages_m[::2], fontsize=6, rotation=45, ha='right')
ax2.set_ylabel('Module score', fontsize=7)
ax2.set_title('Glycolysis module score', fontsize=8, fontweight='bold')
ax2.tick_params(labelsize=6)

# (iii) OxPhos
ax3 = fig6.add_subplot(gs_c[2])
ax3.bar(x_met, met['ox_mean'].values, color='#4DAF4A', edgecolor='black', lw=0.3, width=0.7)
ax3.set_xticks(x_met[::2]); ax3.set_xticklabels(stages_m[::2], fontsize=6, rotation=45, ha='right')
ax3.set_ylabel('Module score', fontsize=7)
ax3.set_title('OxPhos module score', fontsize=8, fontweight='bold')
ax3.tick_params(labelsize=6)

fig6.savefig(os.path.join(OUT, 'fig6_regulon_metabolism.svg'), format='svg', bbox_inches='tight')
fig6.savefig(os.path.join(OUT, 'fig6_regulon_metabolism.png'), dpi=300, bbox_inches='tight')
fig6.savefig(os.path.join(OUT, 'fig6_regulon_metabolism.pdf'), format='pdf', bbox_inches='tight')
    
plt.close(fig6)
print("  -> fig6_regulon_metabolism.svg OK")

# ═══════════════════════════════════════════════════════════════════
# FIGURE 7: Developmental divergence + volcano (2x2 = 4 panels)
# ═══════════════════════════════════════════════════════════════════
print("=== Fig7: Divergence + Volcano ===")
fig7, axes7 = plt.subplots(2, 2, figsize=(5.8, 5.0))
axes7 = axes7.flatten()

# ── A: Mean divergence per comparison ──
ax = axes7[0]; panel_label(ax, 'A')
ivf_pa = div_df['IVF_vs_PA'].mean()
ivf_vivo = div_df['IVF_vs_vivo'].mean()
pa_vivo = div_df['PA_vs_vivo'].mean()
ax.bar(['IVF-PA','IVF-InVivo','PA-InVivo'], [ivf_pa, ivf_vivo, pa_vivo],
       color=['#377EB8','#4DAF4A','#E41A1C'], edgecolor='black', lw=0.5, width=0.5)
for i, v in enumerate([ivf_pa, ivf_vivo, pa_vivo]):
    ax.text(i, v + 0.005, f'{v:.3f}', ha='center', fontsize=6, fontweight='bold')
ax.set_ylabel('Mean divergence (8,276 genes)', fontsize=7)
ax.set_title('Cross-condition divergence', fontsize=8, fontweight='bold')

# ── B: PA vs in vivo: top divergent genes ──
ax = axes7[1]; panel_label(ax, 'B')
# ID2SYM mapping (same as rebuild_supp_svgs_real.py)
_id2sym = {}
_gtf_path = 'E:/Workbuddy/2026-07-27-11-58-27/data/reference/Sus_scrofa.Sscrofa11.1.113.gtf.gz'
if os.path.exists(_gtf_path):
    with gzip.open(_gtf_path, 'rt') as _f:
        for _line in _f:
            if _line.startswith('#'): continue
            _p = _line.strip().split('\t')
            if len(_p) < 9 or _p[2] != 'gene': continue
            _attrs = dict(re.findall(r'(\w+)\s*"([^"]+)"', _p[8]))
            _gid = _attrs.get('gene_id')
            _sym = _attrs.get('gene_name', '')
            if _gid and _sym and _sym != _gid:
                _id2sym[_gid] = _sym
import gzip  # ensure imported
# Fallback: use DESeq2 mapping for symbols not in GTF
_des = pd.read_csv('E:/Workbuddy/2026-07-27-11-58-27/data/processed/deseq2/m15_deseq2_overall_IVFvsPA_stageAdj.csv')
for _, _r in _des.iterrows():
    _gid = str(_r['gene'])
    _sym = str(_r.get('gene_symbol', ''))
    if _sym and not _sym.startswith('ENSSSCG') and _sym.lower() != 'nan':
        _id2sym[_gid] = _sym

def _to_symbol(g):
    if not g or not str(g).startswith('ENSSSCG'):
        return g
    return _id2sym.get(g, g)

top15 = div_df.sort_values('PA_vs_vivo', ascending=False).head(15)
top15['gene_sym'] = top15['gene'].astype(str).map(_to_symbol)
colors_g = ['#E41A1C' if v > 1.0 else '#FF7F00' if v > 0.5 else '#999' for v in top15['PA_vs_vivo'].values]
ax.barh(range(15), top15['PA_vs_vivo'].values, color=colors_g, edgecolor='black', lw=0.3, height=0.7)
ax.set_yticks(range(15))
ax.set_yticklabels(top15['gene_sym'].values, fontsize=6)
ax.invert_yaxis()
ax.set_xlabel('PA vs in vivo divergence', fontsize=7)
ax.set_title('Top 15 PA-divergent genes', fontsize=8, fontweight='bold')
ax.tick_params(labelsize=6.5)

# ── C: Gene-wise scatter PA vs IVF ──
ax = axes7[2]; panel_label(ax, 'C')
ax.scatter(div_df['IVF_vs_vivo'].values, div_df['PA_vs_vivo'].values,
           alpha=0.15, s=2, c='gray', rasterized=True)
mx = max(div_df['IVF_vs_vivo'].max(), div_df['PA_vs_vivo'].max()) * 1.05
ax.plot([0, mx], [0, mx], 'r--', lw=0.8, alpha=0.5, label='equality')
ax.set_xlabel('IVF vs in vivo divergence', fontsize=7)
ax.set_ylabel('PA vs in vivo divergence', fontsize=7)
ax.set_title('Gene-wise divergence comparison', fontsize=8, fontweight='bold')
ax.legend(fontsize=6); ax.tick_params(labelsize=7)

# ── D: Volcano plot ──
ax = axes7[3]; panel_label(ax, 'D')
# Use m3_de for volcano
sig_up = m3_de[(m3_de['significant'] == True) | (m3_de['padj'] < 0.05)]
sig_dn = m3_de[m3_de['log2FC_PA_vs_IVF'] < 0]
ns = m3_de[(m3_de['padj'] >= 0.05) | (m3_de['significant'] != True)]
ax.scatter(ns['log2FC_PA_vs_IVF'].values, -np.log10(ns['pvalue'].clip(1e-300)),
          c='gray', alpha=0.2, s=2, rasterized=True, label='NS')
sig_data = m3_de[(m3_de['padj'] < 0.05) | (m3_de['significant'] == True)]
if len(sig_data) > 0:
    up_s = sig_data[sig_data['log2FC_PA_vs_IVF'] > 1]
    dn_s = sig_data[sig_data['log2FC_PA_vs_IVF'] < -1]
    ax.scatter(up_s['log2FC_PA_vs_IVF'], -np.log10(up_s['pvalue'].clip(1e-300)),
              c='#E41A1C', alpha=0.6, s=6, label=f'Up in PA (n={len(up_s)})')
    ax.scatter(dn_s['log2FC_PA_vs_IVF'], -np.log10(dn_s['pvalue'].clip(1e-300)),
              c='#377EB8', alpha=0.6, s=6, label=f'Down in PA (n={len(dn_s)})')
ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.5)
ax.axvline(-1, color='gray', ls='--', lw=0.5); ax.axvline(1, color='gray', ls='--', lw=0.5)
ax.set_xlabel('log2FC (PA vs IVF)', fontsize=7)
ax.set_ylabel('-log10(p-value)', fontsize=7)
ax.set_title('IVF vs PA differential expression (6,860 genes)', fontsize=8, fontweight='bold')
ax.legend(fontsize=6); ax.tick_params(labelsize=7)
ax.set_xlim(-6, 6)

fig7.tight_layout()
fig7.savefig(os.path.join(OUT, 'fig7_divergence_volcano.svg'), format='svg', bbox_inches='tight')
fig7.savefig(os.path.join(OUT, 'fig7_divergence_volcano.png'), dpi=300, bbox_inches='tight')
fig7.savefig(os.path.join(OUT, 'fig7_divergence_volcano.pdf'), format='pdf', bbox_inches='tight')
    
plt.close(fig7)
print("  -> fig7_divergence_volcano.svg OK")

# ═══════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("ALL 7 FIGURES GENERATED")
print("="*60)
for fn in sorted(os.listdir(OUT)):
    if fn.endswith('.svg') and fn.startswith('fig'):
        size_kb = os.path.getsize(os.path.join(OUT, fn)) // 1024
        print(f"  {fn:35s} {size_kb:>6d} KB")
print("\nAll SVG with fonttype='none' -- fully editable in Adobe Illustrator.")
print("No inset_axes overlap -- clean GridSpec layout throughout.")
