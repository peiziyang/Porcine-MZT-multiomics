#!/usr/bin/env python
"""Regenerate Fig 2 (cross-species) with corrected F1-top44 gene set + genome-wide background."""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial', 'svg.fonttype': 'none', 'savefig.dpi': 300,
    'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7,
})
import matplotlib.pyplot as plt
from adjustText import adjust_text
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import os

MF = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa'
OUT = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/submission/figures'
os.makedirs(OUT, exist_ok=True)

ortho = pd.read_csv(os.path.join(MF, 'cross_species_orthologs.csv'))
perm = pd.read_csv(os.path.join(MF, 'cross_species_permutation_results.csv'))

# extract key values
f1_hu = perm[(perm['test']=='F1_observed') & (perm['species']=='human')]
f1_mo = perm[(perm['test']=='F1_observed') & (perm['species']=='mouse')]
gw_hu = perm[(perm['test']=='genome_wide') & (perm['species']=='human')]
gw_mo = perm[(perm['test']=='genome_wide') & (perm['species']=='mouse')]

r_hu = float(f1_hu['rho']); p_hu = float(f1_hu['p_value']); n_hu = int(f1_hu['n_genes']); pp_hu = float(f1_hu['expr_matched_perm_p'])
r_mo = float(f1_mo['rho']); p_mo = float(f1_mo['p_value']); n_mo = int(f1_mo['n_genes']); pp_mo = float(f1_mo['expr_matched_perm_p'])
r_gwh = float(gw_hu['rho']); n_gwh = int(gw_hu['n_genes'])
r_gwm = float(gw_mo['rho']); n_gwm = int(gw_mo['n_genes'])

def panel_label(ax, label, x=0.02, y=0.97, fontsize=6):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize,
            fontweight='bold', va='top', ha='left')

fig2, axes2 = plt.subplots(2, 2, figsize=(5.8, 4.8))
axes2 = axes2.flatten()

# A: overlap
ax = axes2[0]; panel_label(ax, 'A')
pig_only = ortho['pig_gene'].nunique()
in_human = int(ortho['in_human'].sum())
in_mouse = int(ortho['in_mouse'].sum())
in_both = int((ortho['in_human'] & ortho['in_mouse']).sum())
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
ax.set_title(f'F1 genes (n={pig_only})', fontsize=8, fontweight='bold')

# B: pig vs human
ax = axes2[1]; panel_label(ax, 'B')
hu_genes = ortho[ortho['in_human']].dropna(subset=['pig_GV_mean', 'human_mean'])
ax.scatter(hu_genes['pig_GV_mean'], hu_genes['human_mean'],
           c='#377EB8', alpha=0.6, s=25, edgecolors='black', lw=0.3)
texts_b = []
for g in ['DNMT1', 'ZP3', 'ZP4', 'RARRES1', 'EIF4G2']:
    sub = hu_genes[hu_genes['pig_gene'] == g]
    if len(sub):
        texts_b.append(ax.annotate(g, (sub['pig_GV_mean'].values[0], sub['human_mean'].values[0]),
                                   fontsize=6, fontweight='bold'))
adjust_text(texts_b, ax=ax, force_text=0.3, force_points=0.2,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.4))
ax.set_xlabel('Pig GV mean expression', fontsize=7)
ax.set_ylabel('Human preimplantation mean expr.', fontsize=7)
ax.set_title(f'F1 pig vs human (rho={r_hu:.3f}, n={n_hu})', fontsize=8, fontweight='bold')
ax.tick_params(labelsize=7)

# C: pig vs mouse
ax = axes2[2]; panel_label(ax, 'C')
mo_genes = ortho[ortho['in_mouse']].dropna(subset=['pig_GV_mean', 'mouse_mean'])
ax.scatter(mo_genes['pig_GV_mean'], mo_genes['mouse_mean'],
           c='#FF7F00', alpha=0.6, s=25, edgecolors='black', lw=0.3)
texts_c = []
for g in ['DNMT1', 'ZP3', 'ZP4', 'RARRES1', 'EIF4G2']:
    sub = mo_genes[mo_genes['pig_gene'] == g]
    if len(sub):
        texts_c.append(ax.annotate(g, (sub['pig_GV_mean'].values[0], sub['mouse_mean'].values[0]),
                                   fontsize=6, fontweight='bold'))
adjust_text(texts_c, ax=ax, force_text=0.3, force_points=0.2,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.4))
ax.set_xlabel('Pig GV mean expression', fontsize=7)
ax.set_ylabel('Mouse preimplantation mean expr.', fontsize=7)
ax.set_title(f'F1 pig vs mouse (rho={r_mo:.3f}, n={n_mo})', fontsize=8, fontweight='bold')
ax.tick_params(labelsize=7)

# D: summary
ax = axes2[3]; panel_label(ax, 'D'); ax.axis('off')
ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.text(5, 9.6, 'Cross-species summary', ha='center', va='top',
        fontsize=6, fontweight='bold', color='#333333')
ax.text(0.4, 8.5, 'Symbol-matched gene pairs (n = 44 F1 genes)', fontsize=6,
        fontweight='bold', color='#377EB8')
match_items = [('Human', in_human, '#4DAF4A'), ('Mouse', in_mouse, '#FF7F00'),
               ('Both species', in_both, '#E41A1C')]
for i, (label, val, c) in enumerate(match_items):
    y = 7.7 - i * 0.8
    ax.text(0.7, y, label, fontsize=6, va='center', color='#333333')
    ax.barh(y, val / pig_only * 3.6, left=4.6, height=0.42, color=c, alpha=0.75, zorder=2)
    ax.text(8.4, y, f'{val}/44', fontsize=6, va='center', fontweight='bold', color=c)

ax.text(0.4, 5.0, 'Expression correlation (Spearman rho)', fontsize=6,
        fontweight='bold', color='#4DAF4A')
ax.text(0.7, 4.25, f'F1 pig-human:  rho = {r_hu:.2f}  (n = {n_hu})', fontsize=6, color='#333333')
ax.text(0.7, 3.6, f'   genome-wide: rho = {r_gwh:.2f}  (n = {n_gwh})', fontsize=6, color='#888888')
ax.text(0.7, 2.85, f'F1 pig-mouse:  rho = {r_mo:.2f}  (n = {n_mo})', fontsize=6, color='#333333')
ax.text(0.7, 2.2, f'   genome-wide: rho = {r_gwm:.2f}  (n = {n_gwm})', fontsize=6, color='#888888')

ax.text(5, 1.15, f'F1 modestly exceeds background;\nexpression-matched perm P = {pp_hu:.3f} (human), {pp_mo:.3f} (mouse)',
        ha='center', va='center', fontsize=5.5, style='italic', color='#555555',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#F2F2F2', edgecolor='#CCCCCC'))

fig2.tight_layout()
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.svg'), format='svg', bbox_inches='tight')
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.png'), dpi=300, bbox_inches='tight')
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.pdf'), format='pdf', bbox_inches='tight')
plt.close(fig2)
print('fig2 regenerated OK')
print(f'F1 human rho={r_hu:.3f} n={n_hu} permP={pp_hu:.3f}')
print(f'F1 mouse rho={r_mo:.3f} n={n_mo} permP={pp_mo:.3f}')
print(f'GW human rho={r_gwh:.3f} n={n_gwh}; GW mouse rho={r_gwm:.3f} n={n_gwm}')
print(f'overlap: both={in_both} hu={in_human} mo={in_mouse} pig={pig_only}')
