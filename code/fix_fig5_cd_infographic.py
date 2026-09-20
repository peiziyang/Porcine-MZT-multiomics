#!/usr/bin/env python
"""Redesign Fig5 Panel C+D as infographic panels."""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'Arial', 'svg.fonttype': 'none',
    'savefig.dpi': 300, 'font.size': 9,
})
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from matplotlib.gridspec import GridSpec
import pandas as pd, numpy as np, os

BASE  = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed'
DESEQ = os.path.join(BASE, 'deseq2')
OUT   = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures'
MF    = os.path.join(BASE, 'm4_mofa')

# ── Load PA data ──
pa_all = []
for f in sorted(os.listdir(DESEQ)):
    if not f.startswith('m15_deseq2_'):
        continue
    if not (f.endswith('_IVFvsPA.csv') or f.endswith('_IVFvsPA_stageAdj.csv')):
        continue
    name = f.replace('m15_deseq2_','').replace('_IVFvsPA.csv','').replace('_stageAdj.csv','')
    df = pd.read_csv(os.path.join(DESEQ, f))
    if 'gene_symbol' not in df.columns: continue
    s2i = {}
    for i, row in df.iterrows():
        s = str(row['gene_symbol'])
        if not s.startswith('ENSSSCG'): s2i[s] = i
    for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2',
              'IDH2','PKM','EXOSC9','MDH1','GNL3','DPPA5','NANOG']:
        if g in s2i:
            r = df.iloc[s2i[g]]
            pa_all.append({'gene': g, 'stage': name, 'log2FC': r['log2FoldChange'], 'padj': r['padj']})
pa_df = pd.DataFrame(pa_all)
overall = pa_df[pa_df['stage'].str.startswith('overall')]

def pl(ax, label):
    ax.text(-0.04, 1.04, label, transform=ax.transAxes, fontsize=13,
            fontweight='bold', va='bottom', ha='left')

# ═══════════════════════════════════════════════
# Build complete Fig5
# ═══════════════════════════════════════════════
fig5 = plt.figure(figsize=(15, 18))
gs5 = fig5.add_gridspec(3, 2, hspace=0.45, wspace=0.35,
                        height_ratios=[0.9, 1.2, 1.0])

# ── A: Layer 1 bar (keep as before) ──
ax = fig5.add_subplot(gs5[0, 0]); pl(ax, 'A')
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
    yoff = 0.15 if v >= 0 else -0.35
    ax.text(i, v + yoff, f'{v:+.2f} {sig}', ha='center', fontsize=7.5, fontweight='bold')
ax.set_xticks(x1); ax.set_xticklabels(l1_names, fontsize=8)
ax.set_ylabel('log2 Fold Change (PA / IVF)', fontsize=9)
ax.set_title('Layer 1: Maternal blueprint -- MAINTAINED', fontsize=11, fontweight='bold', color='#E41A1C')
ax.set_ylim(-1.5, 2.5)

# ── B: Layer 3 bar (keep as before) ──
ax = fig5.add_subplot(gs5[0, 1]); pl(ax, 'B')
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
    ax.text(i, 0.4, f'{v:+.2f} {sig}', ha='center', fontsize=7.5, fontweight='bold', color='#377EB8')
ax.set_xticks(x3); ax.set_xticklabels(l3_names, fontsize=8)
ax.set_ylabel('log2 Fold Change (PA / IVF)', fontsize=9)
ax.set_title('Layer 3: Zygotic execution -- ATTENUATED', fontsize=11, fontweight='bold', color='#377EB8')
ax.set_ylim(-10, 1)

# ═══════════════════════════════════════════════
# C: Heatmap REDESIGNED with layer grouping lines + per-gene bar annotation
# ═══════════════════════════════════════════════
ax = fig5.add_subplot(gs5[1, :])
pl(ax, 'C')
stg = ['1cell','2cell','4cell','8cell']
gpl_l1 = ['DNMT1','ZP3','GDF9']       # Layer 1 genes
gpl_l3 = ['IDH2','PKM','EXOSC9','GNL3']  # Layer 3 genes
gpl_zga = ['DPPA5','NANOG']             # ZGA markers
all_genes = gpl_l1 + gpl_l3 + gpl_zga
n_genes = len(all_genes)
n_stages = len(stg)

hm = np.full((n_genes, n_stages), np.nan)
for gi, g in enumerate(all_genes):
    for si, s in enumerate(stg):
        sub = pa_df[(pa_df['gene'] == g) & (pa_df['stage'] == s)]
        if len(sub) and not np.isnan(sub.iloc[0]['log2FC']):
            hm[gi, si] = sub.iloc[0]['log2FC']

hm_masked = np.ma.masked_where(np.isnan(hm), hm)
im = ax.imshow(hm_masked, aspect='equal', cmap='RdBu_r', vmin=-10, vmax=10,
               extent=[-0.5, n_stages-0.5, n_genes-0.5, -0.5])
sep1 = len(gpl_l1) - 0.5   # between L1 and L3
sep2 = sep1 + len(gpl_l3)   # between L3 and ZGA
ax.axhline(sep1, color='#E41A1C', lw=2.5, ls='-', alpha=0.7)
ax.axhline(sep2, color='#4DAF4A', lw=2.5, ls='-', alpha=0.7)

# Layer labels on left side (avoid colorbar)
ax.text(-1.4, (0 + sep1)/2, 'Layer 1', fontsize=8, fontweight='bold',
        color='#E41A1C', va='center', ha='right')
ax.text(-1.4, (sep1 + sep2)/2, 'Layer 3', fontsize=8, fontweight='bold',
        color='#377EB8', va='center', ha='right')
ax.text(-1.4, (sep2 + n_genes - 1)/2, 'ZGA', fontsize=8, fontweight='bold',
        color='#4DAF4A', va='center', ha='right')

# Cell value annotations
for i in range(n_genes):
    for j in range(n_stages):
        v = hm[i, j]
        if not np.isnan(v):
            ax.text(j, i, f'{v:+.1f}', ha='center', va='center', fontsize=7.5,
                   fontweight='bold', color='white' if abs(v) > 5 else 'black')

ax.set_xticks(np.arange(n_stages)); ax.set_xticklabels(stg, fontsize=9)
ax.set_yticks(np.arange(n_genes)); ax.set_yticklabels(all_genes, fontsize=9)
ax.set_xlim(-2.2, n_stages-0.5)
ax.set_title('Per-stage log2FC (PA vs IVF) -- grouped by regulatory layer', fontsize=11, fontweight='bold')
cbar = fig5.colorbar(im, ax=ax, fraction=0.02, pad=0.04, shrink=0.8)
cbar.ax.tick_params(labelsize=8)

# ── Near-heatmap annotation: per-gene bar ──
# Add small horizontal bar chart for overall log2FC next to gene names
gene_stats = []
for g in all_genes:
    sub = overall[overall['gene'] == g]
    if len(sub) and not np.isnan(sub.iloc[0]['log2FC']):
        gene_stats.append((g, sub.iloc[0]['log2FC'], sub.iloc[0]['padj']))
    else:
        gene_stats.append((g, np.nan, np.nan))

# ═══════════════════════════════════════════════
# D: Lollipop chart REDESIGNED
# ═══════════════════════════════════════════════
ax = fig5.add_subplot(gs5[2, :])
pl(ax, 'D')
ax.set_xlim(-12, 6.8); ax.set_ylim(-1, n_genes + 1)

# Gene rows (from Layer 3 up to Layer 1, bottom-up)
ordered_genes = gpl_l1 + gpl_l3 + gpl_zga  # Layer1 top, Layer3 middle, ZGA bottom
layer_colors = (['#E41A1C'] * len(gpl_l1) +
                ['#377EB8'] * len(gpl_l3) +
                ['#4DAF4A'] * len(gpl_zga))
layer_bgs   = (['#FFE5E5'] * len(gpl_l1) +
                ['#E5F0FF'] * len(gpl_l3) +
                ['#E5FFE5'] * len(gpl_zga))

for i, (g, col, bg) in enumerate(zip(ordered_genes, layer_colors, layer_bgs)):
    sub = overall[overall['gene'] == g]
    if len(sub) == 0 or np.isnan(sub.iloc[0]['log2FC']):
        continue
    v = sub.iloc[0]['log2FC']
    p = sub.iloc[0]['padj']
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

    y = n_genes - 1 - i  # flip so Layer1 is at top

    # Background strip
    bg_rect = FancyBboxPatch((-11.5, y-0.35), 17, 0.7, boxstyle='round,pad=0.04',
                             facecolor=bg, edgecolor=col, lw=0.8, alpha=0.4)
    ax.add_patch(bg_rect)

    # Lollipop: stem from 0 to value
    ax.plot([0, v], [y, y], '-', color=col, lw=3, alpha=0.7, solid_capstyle='round')
    # Head
    ax.plot(v, y, 'o', color=col, markersize=14, markeredgecolor='black', markeredgewidth=0.8, zorder=5)

    # Value label inside/next to marker
    label_x = v + 0.4 if v >= 0 else v - 0.8
    ax.text(label_x, y, f'{v:+.2f}', fontsize=8, fontweight='bold', color=col,
            va='center', ha='left' if v >= 0 else 'right')

    # Gene name at left
    ax.text(-11.2, y, g, fontsize=9, fontweight='bold', color='#333', va='center')

    # Significance badge
    if sig != 'ns':
        badge_c = '#CC0000' if p < 0.01 else '#EF6C00' if p < 0.05 else '#888'
        ax.text(3.5, y, sig, fontsize=8, fontweight='bold', color=badge_c, va='center')

# Zero line
ax.axvline(0, color='black', lw=1.2)

# Layer separator lines between groups
for y_line, color in [(len(gpl_zga) - 0.5, '#4DAF4A'),    # between ZGA and L3
                      (len(gpl_zga) + len(gpl_l3) - 0.5, '#E41A1C')]:  # between L3 and L1
    ax.axhline(y_line, color=color, lw=1.5, ls='--', alpha=0.5)

# Layer badges (top-right of each group's first gene)
# Layer 1: top group (y from len(gpl_l3)+len(gpl_zga) to n_genes-1)
# Layer 3: middle group (y from len(gpl_zga) to len(gpl_zga)+len(gpl_l3)-1)
# ZGA: bottom group (y from 0 to len(gpl_zga)-1)
for y_start, y_end, name, col in [
    (len(gpl_zga) + len(gpl_l3), n_genes - 1, 'Layer 1 (Blueprint)', '#E41A1C'),
    (len(gpl_zga), len(gpl_zga) + len(gpl_l3) - 1, 'Layer 3 (Execution)', '#377EB8'),
    (0, len(gpl_zga) - 1, 'ZGA Markers', '#4DAF4A')]:
    y_mid = (y_start + y_end) / 2
    fb = FancyBboxPatch((5.0, y_mid-0.3), 1.5, 0.6, boxstyle='round,pad=0.05',
                       facecolor='white', edgecolor=col, lw=1.5)
    ax.add_patch(fb)
    ax.text(5.75, y_mid, name, fontsize=7, fontweight='bold', color=col,
           ha='center', va='center')

# Legend mini-box
ax.text(-11.2, n_genes + 0.5, 'Gene', fontsize=8, fontweight='bold', color='#666', va='center')
ax.text(0, n_genes + 0.5, 'log2FC (PA vs IVF)', fontsize=8, fontweight='bold', color='#666', ha='center', va='center')
ax.text(4.5, n_genes + 0.5, 'padj', fontsize=8, fontweight='bold', color='#666', ha='center', va='center')

# Significance legend
for xi, (sigl, sigd) in enumerate([('***', 'p < 0.001'), ('**', 'p < 0.01'), ('*', 'p < 0.05'), ('ns', 'not significant')]):
    ax.text(-2 + xi*3.5, n_genes + 0.9, f'{sigl} = {sigd}', fontsize=6, color='#888', va='center')

ax.axis('off')
ax.set_title('Stage-adjusted log2FC in PA vs IVF embryos (SRP301735, 42 samples)', fontsize=11, fontweight='bold')

# ── Save ──
plt.tight_layout()
fig5.savefig(os.path.join(OUT, 'fig5_pa_validation.svg'), format='svg', bbox_inches='tight')
fig5.savefig(os.path.join(OUT, 'fig5_pa_validation.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Fig5 C+D redesigned: grouped heatmap + lollipop chart")
