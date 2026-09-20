#!/usr/bin/env python
"""Rebuild Fig2 with an infographic Panel D instead of raw text."""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial','svg.fonttype': 'none',
    'savefig.dpi': 300, 'font.size': 9,
})
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Arc
import pandas as pd, numpy as np, os
from scipy.stats import spearmanr

BASE = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
MF   = os.path.join(BASE, 'm4_mofa')
OUT  = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures'

ortho = pd.read_csv(os.path.join(MF, 'cross_species_orthologs.csv'))
pig_only = ortho['pig_gene'].nunique()
in_human = ortho[ortho['in_human']].shape[0]
in_mouse = ortho[ortho['in_mouse']].shape[0]
in_both  = ortho[ortho['in_human'] & ortho['in_mouse']].shape[0]
hu_genes = ortho[ortho['in_human']].dropna(subset=['pig_GV_mean','human_mean'])
mo_genes = ortho[ortho['in_mouse']].dropna(subset=['pig_GV_mean','mouse_mean'])
r_hu, p_hu = spearmanr(hu_genes['pig_GV_mean'], hu_genes['human_mean'])
r_mo, p_mo = spearmanr(mo_genes['pig_GV_mean'], mo_genes['mouse_mean'])

def pl(ax, label):
    ax.text(-0.04, 1.04, label, transform=ax.transAxes, fontsize=13,
            fontweight='bold', va='bottom', ha='left')

fig2 = plt.figure(figsize=(12, 10))
gs = fig2.add_gridspec(2, 2, hspace=0.4, wspace=0.35)

# ═══════════════════════════════════════════════
# A: Species overlap (keep original)
# ═══════════════════════════════════════════════
ax = fig2.add_subplot(gs[0, 0]); pl(ax, 'A')
human_only = in_human - in_both
mouse_only = in_mouse - in_both
pig_excl = pig_only - in_human - in_mouse + in_both
cats = ['Conserved\nin 3 species', 'Pig+Human\nonly', 'Pig+Mouse\nonly', 'Pig only']
vals = [in_both, human_only, mouse_only, pig_excl]
cc = ['#377EB8', '#4DAF4A', '#FF7F00', '#CCCCCC']
ax.bar(cats, vals, color=cc, edgecolor='black', lw=0.5, width=0.6)
for i, v in enumerate(vals):
    ax.text(i, v + 0.5, str(v), ha='center', fontsize=10, fontweight='bold', color=cc[i])
ax.set_ylabel('Number of F1 genes', fontsize=9)
ax.set_title(f'F1 maternal blueprint genes (n={pig_only})', fontsize=11, fontweight='bold')
ax.tick_params(labelsize=8)

# ═══════════════════════════════════════════════
# B: Pig vs Human (keep original)
# ═══════════════════════════════════════════════
ax = fig2.add_subplot(gs[0, 1]); pl(ax, 'B')
ax.scatter(hu_genes['pig_GV_mean'], hu_genes['human_mean'],
           c='#377EB8', alpha=0.6, s=30, edgecolors='black', lw=0.3)
label_offsets = {'DNMT1': (10, 10), 'ZP3': (-30, 10), 'ZP4': (-40, -15),
                  'GDF9': (10, 10), 'RARRES1': (15, -20)}
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    sub = hu_genes[hu_genes['pig_gene'] == g]
    if len(sub):
        x0, y0 = sub['pig_GV_mean'].values[0], sub['human_mean'].values[0]
        off = label_offsets.get(g, (5, 5))
        mx = hu_genes['pig_GV_mean'].max(); my = hu_genes['human_mean'].max()
        ax.annotate(g, (x0, y0), fontsize=8, fontweight='bold',
                   xytext=(x0+off[0]*0.005*mx, y0+off[1]*0.005*my),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
ax.set_xlabel('Pig GV mean expression', fontsize=9)
ax.set_ylabel('Human preimplantation mean expr.', fontsize=9)
ax.set_title(f'Pig vs Human   rho = {r_hu:.3f}   p = {p_hu:.2e}', fontsize=10, fontweight='bold')
ax.tick_params(labelsize=7)

# ═══════════════════════════════════════════════
# C: Pig vs Mouse (keep original)
# ═══════════════════════════════════════════════
ax = fig2.add_subplot(gs[1, 0]); pl(ax, 'C')
ax.scatter(mo_genes['pig_GV_mean'], mo_genes['mouse_mean'],
           c='#FF7F00', alpha=0.6, s=30, edgecolors='black', lw=0.3)
mo_offsets = {'DNMT1': (5, 5), 'ZP3': (5, 5), 'ZP4': (-30, -15),
              'GDF9': (5, 5), 'RARRES1': (-35, 5)}
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    sub = mo_genes[mo_genes['pig_gene'] == g]
    if len(sub):
        x0, y0 = sub['pig_GV_mean'].values[0], sub['mouse_mean'].values[0]
        off = mo_offsets.get(g, (5, 5))
        mx = mo_genes['pig_GV_mean'].max(); my = mo_genes['mouse_mean'].max()
        ax.annotate(g, (x0, y0), fontsize=8, fontweight='bold',
                   xytext=(x0+off[0]*0.005*mx, y0+off[1]*0.005*my),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
ax.set_xlabel('Pig GV mean expression', fontsize=9)
ax.set_ylabel('Mouse preimplantation mean expr.', fontsize=9)
ax.set_title(f'Pig vs Mouse   rho = {r_mo:.3f}   p = {p_mo:.2e}', fontsize=10, fontweight='bold')
ax.tick_params(labelsize=7)

# ═══════════════════════════════════════════════
# D: INFOGRAPHIC (redesigned from scratch)
# ═══════════════════════════════════════════════
ax = fig2.add_subplot(gs[1, 1])
ax.set_xlim(0, 14); ax.set_ylim(0, 12); ax.axis('off')
pl(ax, 'D')

# ── Title ──
ax.text(7, 11.5, 'Cross-Species MZT Conservation', ha='center',
        fontsize=12, fontweight='bold', color='#1a1a2e')

# ── Stat boxes row ──
stat_data = [
    (f'{in_both}/{pig_only}', f'{100*in_both//pig_only}%', '3-Species\nConserved', '#377EB8', '#E5F0FF'),
    (f'{in_human}/{pig_only}', f'{100*in_human//pig_only}%', 'Pig-Human\nOverlap',    '#4DAF4A', '#E5FFE5'),
    (f'{in_mouse}/{pig_only}', f'{100*in_mouse//pig_only}%', 'Pig-Mouse\nOverlap',    '#FF7F00', '#FFF5E5'),
]
for i, (frac, pct, label, fg, bg) in enumerate(stat_data):
    x = 2 + i * 5
    fb = FancyBboxPatch((x-1.8, 8.0), 3.6, 2.8, boxstyle='round,pad=0.15',
                        facecolor=bg, edgecolor=fg, lw=1.5, alpha=0.5)
    ax.add_patch(fb)
    ax.text(x, 10.2, frac, ha='center', fontsize=18, fontweight='bold', color=fg)
    ax.text(x, 9.5, pct, ha='center', fontsize=13, fontweight='bold', color=fg)
    ax.text(x, 8.5, label, ha='center', fontsize=7, color='#555')

# ── Correlation row ──
y_cor = 7.2
fb_cor = FancyBboxPatch((0.8, y_cor-0.1), 12.4, 0.75, boxstyle='round,pad=0.08',
                         facecolor='#FFFDE7', edgecolor='#FFC107', lw=1, alpha=0.6)
ax.add_patch(fb_cor)

# Pig-Human
ax.text(2.5, y_cor+0.45, 'Pig-Human', ha='center', fontsize=8, fontweight='bold', color='#377EB8')
ax.text(2.5, y_cor+0.1, f'rho = {r_hu:.3f}', ha='center', fontsize=9, fontweight='bold', color='#377EB8')
# Pig-Mouse
ax.text(7, y_cor+0.45, 'Pig-Mouse', ha='center', fontsize=8, fontweight='bold', color='#FF7F00')
ax.text(7, y_cor+0.1, f'rho = {r_mo:.3f}', ha='center', fontsize=9, fontweight='bold', color='#FF7F00')
# P values
ax.text(11.5, y_cor+0.45, 'Significance', ha='center', fontsize=8, fontweight='bold', color='#333')
ax.text(11.5, y_cor+0.1, f'p = {p_hu:.1e} / {p_mo:.1e}', ha='center', fontsize=7, color='#555')
# Dividers
ax.axvline(5.5, ymin=0.32, ymax=0.42, color='#DDD', lw=1)
ax.axvline(10, ymin=0.32, ymax=0.42, color='#DDD', lw=1)

# ── Core conserved genes badge row ──
ax.text(7, 5.8, 'Core Conserved Maternal Blueprint Genes', ha='center',
        fontsize=9, fontweight='bold', color='#333')

core_genes = ['DNMT1', 'ZP3', 'ZP4', 'GDF9', 'RARRES1', 'ACTB', 'ANXA1']
gene_colors = ['#E41A1C']*5 + ['#666']*2
cx_start = 2.5
for i, (g, gc) in enumerate(zip(core_genes, gene_colors)):
    xg = cx_start + i * 1.5
    fb_g = FancyBboxPatch((xg-0.55, 4.8), 1.1, 0.6, boxstyle='round,pad=0.05',
                          facecolor='white', edgecolor=gc, lw=1.5, alpha=0.9)
    ax.add_patch(fb_g)
    ax.text(xg, 5.1, g, ha='center', va='center', fontsize=7,
            fontweight='bold', color=gc)

# ── Additional note ──
ax.text(1.5, 3.8, '+DNMT3A, DNMT3B, TET1, EIF4G2, AKIRIN2, ANAPC4, BTG4, ...',
        fontsize=6.5, color='#888', style='italic')

# ── Bottom info (single clean section) ──
fb_info = FancyBboxPatch((0.3, 0.15), 13.4, 3.4, boxstyle='round,pad=0.1',
                          facecolor='#F5F5F5', edgecolor='#CCC', lw=1, alpha=0.9)
ax.add_patch(fb_info)

# Two-row layout
# Row 1: big takeaway
ax.text(7.0, 3.05, '44 maternal blueprint genes - 79% conserved across all three species',
        ha='center', fontsize=9, color='#1a1a2e', fontweight='bold', va='center')
ax.text(7.0, 2.55, 'Pig-Human rho = 0.501 (p=9.8e-04)   |   Pig-Mouse rho = 0.363 (p=3.0e-02)',
        ha='center', fontsize=8, color='#555', fontstyle='italic', va='center')

# Divider
ax.axhline(2.10, xmin=0.04, xmax=0.96, color='#DDD', lw=0.5)

# Row 2: data sources
ax.text(7.0, 1.80, 'DATA SOURCES', ha='center', fontsize=7.5,
        color='#333', fontweight='bold', va='center')
sources = [
    ('Pig:', '#377EB8', '32 GV oocytes, scRNA-seq + WGBS'),
    ('Human+Mouse:', '#FF7F00', 'GSE44183 (Xue et al., 2013; preimplantation qPCR)'),
    ('Matching:', '#555', 'orthologous gene symbols (case-insensitive for mouse)'),
]
y = 1.45
for label, color, body in sources:
    # Label at fixed x=0.6
    ax.text(0.6, y, label, fontsize=7.5, color=color, fontweight='bold', va='center')
    # Body text at x=3.5 to leave clear gap
    ax.text(3.5, y, body, fontsize=7.5, color='#333', va='center')
    y -= 0.30

plt.tight_layout()
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.svg'), format='svg', bbox_inches='tight')
fig2.savefig(os.path.join(OUT, 'fig2_crossspecies.png'), dpi=300, bbox_inches='tight')
plt.close()
print(f'Fig2 rebuilt: Panel D now infographic')
