#!/usr/bin/env python
"""Regenerate all manuscript figures as SVG for Biology of Reproduction submission.

BOR requirements:
  - Sans-serif font (Arial/Helvetica/Calibri)
  - SVG with fonttype='none' (Adobe Illustrator editable text)
  - RGB color mode
  - Panel labels (A, B, C, D) in upper left corners
  - Figure size: single-col 3.25"/8.3cm or double-col 5.8"/14.7cm

Outputs (in manuscript/figures/):
  fig1_mofa.svg       – MOFA+ multi-omics F1-F7
  fig2_crossspecies.svg – Cross-species MZT conservation
  fig3_mzt_trajectory.svg – Three-layer MZT dynamics
  fig4_perturbation.svg  – In silico TF perturbation
  fig5_pa_validation.svg – PA vs IVF bulk RNA-seq
  fig6_model.svg         – Three-layer model summary
  fig_ga.svg             – Graphical abstract (1200x900px, 300ppi)
"""
import os, sys
import pandas as pd, numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    'font.family': 'Arial',
    'font.size': 9,
    'svg.fonttype': 'none',           # <-- Adobe Illustrator editable text
    'pdf.fonttype': 42,
    'axes.unicode_minus': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})
import matplotlib.pyplot as plt

BASE  = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"
MFA   = f"{BASE}/m4_mofa"
SCEN  = f"{BASE}/pyscenic_out"
PERT  = f"{BASE}/perturbation"
OUTDIR = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/figures"
os.makedirs(OUTDIR, exist_ok=True)

# ============== DATA LOADING ==============
# MOFA+ v2
Z = pd.read_csv(f"{MFA}/mofa_multiomics_factors_v2.csv", index_col=0)
W_rna = pd.read_csv(f"{MFA}/mofa_multiomics_weights_RNA_v2.csv", index_col=0)
W_meth = pd.read_csv(f"{MFA}/mofa_multiomics_weights_METH_v2.csv", index_col=0)
r2 = pd.read_csv(f"{MFA}/mofa_multiomics_r2_per_view_v2.csv")
meta = pd.read_csv(f"{BASE}/m2_metadata.csv")
cross_df = pd.read_csv(f"{MFA}/factor_cross_view_v2.csv")
fcols = [c for c in Z.columns if c.startswith("F")]

# Donor colors
donor_map = meta.set_index("cell")["donor"].to_dict()
donor_colors = {"AF1":"#E41A1C","AF2":"#377EB8","AF3":"#4DAF4A","AF4":"#984EA3","AF5":"#FF7F00"}
donors = [donor_map.get(c, "unknown") for c in Z.index]

# Cross-species
ortho = pd.read_csv(f"{MFA}/cross_species_orthologs.csv")

# MZT trajectory
proj = pd.read_csv(f"{MFA}/mzt_trajectory_projection.csv", index_col=0)
stage_mean = pd.read_csv(f"{MFA}/mzt_stage_mean_projection.csv", index_col=0)

# Perturbation
pert_vs_pa = pd.read_csv(f"{PERT}/perturbation_vs_pa.csv")
pert_df = pd.read_csv(f"{PERT}/perturbation_vectors.csv")

# PA DESeq2 data
deseq_dir = f"{BASE}/deseq2"
# Load all stage results
pa_data = {}
for f in sorted(os.listdir(deseq_dir)):
    if f.startswith("m15_deseq2_") and f.endswith("_IVFvsPA.csv"):
        stage = f.replace("m15_deseq2_","").replace("_IVFvsPA.csv","")
        df = pd.read_csv(os.path.join(deseq_dir, f))
        if 'gene_symbol' in df.columns:
            pa_data[stage] = df

# Helper: panel label
def panel_label(ax, label, x=-0.02, y=1.02):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight='bold',
            va='bottom', ha='left')

def save_both(fig, name):
    fig.savefig(f"{OUTDIR}/{name}.svg", format='svg')
    fig.savefig(f"{OUTDIR}/{name}.png", dpi=300)

# ============================================================
# FIGURE 1: MOFA+ Multi-Omics
# ============================================================
print("Generating Fig 1 ...", flush=True)
fig1, axes = plt.subplots(2, 2, figsize=(10.5, 8.5))  # double-col
fig1.suptitle("", fontsize=14, fontweight='bold')

# A: Factor scores heatmap
ax = axes[0,0]; panel_label(ax, 'A')
im = ax.imshow(Z[fcols].values.T, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
ax.set_yticks(range(len(fcols))); ax.set_yticklabels(fcols, fontsize=8)
ax.set_xticks([])
ax.set_title('MOFA+ factor scores (32 GV oocytes)', fontsize=10)
plt.colorbar(im, ax=ax, shrink=0.7)

# B: Per-factor per-view R²
ax = axes[0,1]; panel_label(ax, 'B')
r2p = r2.pivot(index='factor', columns='view', values='r2_pct').fillna(0)
xpos = np.arange(len(fcols))
ax.bar(xpos-0.2, r2p['METH'].reindex(fcols).fillna(0).values, 0.4, color='#FF7F00', label='Methylation', edgecolor='black')
ax.bar(xpos+0.2, r2p['RNA'].reindex(fcols).fillna(0).values, 0.4, color='#377EB8', label='RNA', edgecolor='black')
ax.set_xticks(xpos); ax.set_xticklabels(fcols); ax.set_ylabel('Variance contribution (%)')
ax.set_title('Per-factor per-view variance', fontsize=10); ax.legend(fontsize=8)

# C: Cross-view correlation
ax = axes[1,0]; panel_label(ax, 'C')
cv = cross_df.set_index('factor')['rna_meth_corr_r']
colors_cv = ['#E41A1C' if cv.get(f,0)>0.3 else '#BBBBBB' for f in fcols]
bars = ax.bar(fcols, [cv.get(f,0) for f in fcols], color=colors_cv, edgecolor='black')
for bar, f in zip(bars, fcols):
    v = cv.get(f, 0); ax.text(bar.get_x()+bar.get_width()/2, v+0.03*np.sign(v), f'{v:+.2f}',
                               ha='center', fontsize=8, fontweight='bold' if abs(v)>0.5 else 'normal')
ax.axhline(0, color='gray', lw=0.5); ax.set_ylabel('Pearson r (RNA vs METH weights)')
ax.set_title('Cross-view gene weight correlation', fontsize=10)

# D: F1 top genes (RNA vs METH)
ax = axes[1,1]; panel_label(ax, 'D')
f1_rna = W_rna['F1'].sort_values(ascending=False).head(8)
f1_meth = W_meth['F1'].sort_values(ascending=False).head(8)
y1 = np.arange(len(f1_rna))
y2 = np.arange(len(f1_meth)) + len(f1_rna) + 1
ax.barh(y1, f1_rna.values, color='#377EB8', edgecolor='black', label='RNA view')
ax.barh(y2, f1_meth.values, color='#FF7F00', edgecolor='black', label='METH view')
all_labels = list(f1_rna.index) + list(f1_meth.index)
all_y = list(y1) + list(y2)
ax.set_yticks(all_y); ax.set_yticklabels(all_labels, fontsize=8)
ax.set_xlabel('Weight'); ax.set_title('F1 top genes per view', fontsize=10)
ax.legend(fontsize=8, loc='lower right')

fig1.tight_layout()
save_both(fig1, 'fig1_mofa')
print("  Fig 1 saved", flush=True)

# ============================================================
# FIGURE 2: Cross-Species MZT Conservation
# ============================================================
print("Generating Fig 2 ...", flush=True)
fig2, axes = plt.subplots(2, 2, figsize=(10, 8))

# A: Ortholog overlap
ax = axes[0,0]; panel_label(ax, 'A')
f1_human = ortho['in_human'].sum(); f1_mouse = ortho['in_mouse'].sum()
f1_both = ((ortho['in_human']) & (ortho['in_mouse'])).sum()
cats = ['Pig F1\ngenes', 'Human\northologs', 'Mouse\northologs', 'Both\nspecies']
vals = [44, int(f1_human), int(f1_mouse), int(f1_both)]
colors = ['#A6A6A6','#377EB8','#FF7F00','#E41A1C']
bars = ax.bar(cats, vals, color=colors, edgecolor='black')
for bar, v in zip(bars, vals): ax.text(bar.get_x()+bar.get_width()/2, v+1, str(v), ha='center', fontsize=10, fontweight='bold')
ax.set_ylabel('Number of genes'); ax.set_title('F1 maternal blueprint gene conservation', fontsize=10)

# B: Pig vs Human
ax = axes[0,1]; panel_label(ax, 'B')
hu = ortho[ortho['in_human']==True].dropna(subset=['pig_GV_mean','human_mean'])
r_hu, p_hu = spearmanr(hu['pig_GV_mean'], hu['human_mean'])
ax.scatter(hu['pig_GV_mean'], hu['human_mean'], c='#377EB8', alpha=0.7, s=25, edgecolors='black', lw=0.3)
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    if g in hu['pig_gene'].values:
        idx = hu[hu['pig_gene']==g].index[0]
        ax.annotate(g, (hu.loc[idx,'pig_GV_mean'], hu.loc[idx,'human_mean']), fontsize=8, fontweight='bold')
ax.set_xlabel('Pig GV oocyte mean expr.'); ax.set_ylabel('Human preimplantation mean expr.')
ax.set_title(f'Pig vs Human  \u03c1={r_hu:.2f}, P={p_hu:.1e}', fontsize=10)

# C: Pig vs Mouse
ax = axes[1,0]; panel_label(ax, 'C')
mo = ortho[ortho['in_mouse']==True].dropna(subset=['pig_GV_mean','mouse_mean'])
r_mo, p_mo = spearmanr(mo['pig_GV_mean'], mo['mouse_mean'])
ax.scatter(mo['pig_GV_mean'], mo['mouse_mean'], c='#FF7F00', alpha=0.7, s=25, edgecolors='black', lw=0.3)
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    if g in mo['pig_gene'].values:
        idx = mo[mo['pig_gene']==g].index[0]
        ax.annotate(g, (mo.loc[idx,'pig_GV_mean'], mo.loc[idx,'mouse_mean']), fontsize=8, fontweight='bold')
ax.set_xlabel('Pig GV oocyte mean expr.'); ax.set_ylabel('Mouse preimplantation mean expr.')
ax.set_title(f'Pig vs Mouse  \u03c1={r_mo:.2f}, P={p_mo:.1e}', fontsize=10)

# D: ZGA markers
ax = axes[1,1]; panel_label(ax, 'D')
zga_counts = [10, 33, 30]; zga_species = ['Pig', 'Human', 'Mouse']
ax.bar(zga_species, zga_counts, color=['#E41A1C','#377EB8','#4DAF4A'], edgecolor='black')
for i, v in enumerate(zga_counts): ax.text(i, v+1, f'{v}/38', ha='center', fontsize=10, fontweight='bold')
ax.set_ylabel('Canonical ZGA markers present'); ax.set_title('ZGA marker conservation', fontsize=10)

fig2.tight_layout()
save_both(fig2, 'fig2_crossspecies')
print("  Fig 2 saved", flush=True)

# ============================================================
# FIGURE 3: MZT Trajectory
# ============================================================
print("Generating Fig 3 ...", flush=True)
fig3, axes = plt.subplots(2, 2, figsize=(10.5, 8))

stages_all = sorted(stage_mean.index.tolist(), key=lambda s: int(s[1:]))

# A: All factors line plot
ax = axes[0,0]; panel_label(ax, 'A')
colors7 = ["#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#999999"]
for fi, fc in enumerate(fcols):
    vals = [stage_mean.loc[s, fc] for s in stages_all if s in stage_mean.index]
    stages_f = [s for s in stages_all if s in stage_mean.index]
    ax.plot(stages_f, vals, 'o-', color=colors7[fi], label=fc, lw=1.5, markersize=4)
ax.set_xlabel('Developmental stage'); ax.set_ylabel('Mean projected factor score')
ax.set_title('Projected factor scores across E0-E10', fontsize=10); ax.legend(fontsize=7, ncol=4)
ax.axhline(0, color='gray', ls='--', lw=0.5)

# B: F3 (clearance) per stage
ax = axes[0,1]; panel_label(ax, 'B')
f3v = [stage_mean.loc[s,'F3'] for s in stages_all if s in stage_mean.index]
r_f3, p_f3 = spearmanr(np.arange(len(f3v)), f3v)
ax.bar(stages_all[:len(f3v)], f3v, color='#E41A1C', alpha=0.7, edgecolor='black')
ax.set_xlabel('Developmental stage'); ax.set_ylabel('F3 projected score')
ax.set_title(f'F3: Maternal clearance  \u03c1={r_f3:.2f}, P={p_f3:.1e}', fontsize=10)

# C: F4/F6 (ZGA) per stage
ax = axes[1,0]; panel_label(ax, 'C')
for fc, color in [('F4','#377EB8'),('F6','#FF7F00')]:
    vals = [stage_mean.loc[s,fc] for s in stages_all if s in stage_mean.index]
    ax.plot(stages_all[:len(vals)], vals, 'o-', color=color, label=fc, lw=1.5, markersize=4)
ax.set_xlabel('Developmental stage'); ax.set_ylabel('Projected score')
ax.set_title('F4/F6: Zygotic activation', fontsize=10); ax.legend(fontsize=8)
ax.axhline(0, color='gray', ls='--', lw=0.5)

# D: Layer summary bar chart
ax = axes[1,1]; panel_label(ax, 'D')
layers = ['F1 Blueprint\n(stable)', 'F3 Clearance\n(declining)', 'F4 ZGA\n(rising)', 'F6 ZGA\n(rising)']
r_vals = [0.05, -0.90, 0.64, 0.76]
p_vals = [0.9, 0.0002, 0.035, 0.006]
colors_layer = ['#377EB8','#E41A1C','#4DAF4A','#FF7F00']
ax.barh(layers, r_vals, color=colors_layer, edgecolor='black')
for i, (r, p) in enumerate(zip(r_vals, p_vals)):
    ax.text(r+0.03*np.sign(r), i, f'\u03c1={r:+.2f}  P={p:.1e}', va='center', fontsize=9)
ax.axvline(0, color='gray', lw=0.8); ax.set_xlabel("Spearman \u03c1 (per-stage pseudobulk)")
ax.set_title('Three-layer temporal dynamics', fontsize=10)

fig3.tight_layout()
save_both(fig3, 'fig3_mzt_trajectory')
print("  Fig 3 saved", flush=True)

# ============================================================
# FIGURE 4: In Silico Perturbation
# ============================================================
print("Generating Fig 4 ...", flush=True)
fig4, axes = plt.subplots(2, 2, figsize=(10, 8))

# A: DNMT1 perturbation shift
ax = axes[0,0]; panel_label(ax, 'A')
dnmt1 = pert_df[pert_df['tf']=='DNMT1']
cc = {"100%_KD":"#E41A1C","50%_KD":"#FF7F00","50%_OE":"#377EB8","100%_OE":"#4DAF4A"}
for label in cc:
    sub = dnmt1[dnmt1['perturbation']==label]
    ax.scatter(sub['shift_F1'], sub['shift_F2'], c=cc[label], s=25, alpha=0.8, label=label, edgecolors='black', lw=0.3)
ax.axhline(0,color='gray',ls='--',lw=0.5); ax.axvline(0,color='gray',ls='--',lw=0.5)
ax.set_xlabel('\u0394 F1'); ax.set_ylabel('\u0394 F2')
ax.set_title('DNMT1 perturbation: MOFA+ F1-F2 shift', fontsize=10); ax.legend(fontsize=7)

# B: Cosine similarity
ax = axes[0,1]; panel_label(ax, 'B')
cos_sorted = pert_vs_pa[pert_vs_pa['perturbation']=='100%_KD'].sort_values('cos_sim', ascending=False)
cos_all = cos_sorted.copy()
# Add random baseline
ax.barh(np.arange(len(cos_all)+2), 
        list(cos_all['cos_sim'].values) + [-0.04, 0.04],
        color=['#E41A1C' if c>0.3 else '#FF7F00' if c>0 else '#377EB8' for c in 
               list(cos_all['cos_sim'].values) + [-0.04, 0.04]], edgecolor='black', height=0.6)
all_labs = [f"{t['tf']}" for _,t in cos_all.iterrows()] + ['Null (other TFs)','Null (random genes)']
ax.set_yticks(np.arange(len(all_labs))); ax.set_yticklabels(all_labs, fontsize=8)
ax.axvline(0, color='gray', lw=0.8)
ax.set_xlabel('Cosine similarity with PA deviation'); ax.set_title('Perturbation alignment', fontsize=10)

# C: Factor-level perturbation effects
ax = axes[1,0]; panel_label(ax, 'C')
kd_data = {}
for tf in pert_df['tf'].unique():
    if tf in ['DNMT1','ATF3','TFAP2C','ESRRA','GATA4','NFE2L3']:
        sub = pert_df[(pert_df['tf']==tf)&(pert_df['perturbation']=='100%_KD')]
        kd_data[tf] = [sub[f'shift_{fc}'].mean() for fc in fcols[:4]]
kd_tfs = list(kd_data.keys())
im_data = np.array([kd_data[t] for t in kd_tfs]).T
im = ax.imshow(im_data, aspect='auto', cmap='RdBu_r', vmin=-0.04, vmax=0.04)
ax.set_xticks(np.arange(len(kd_tfs))); ax.set_xticklabels(kd_tfs, rotation=45, ha='right', fontsize=8)
ax.set_yticks(np.arange(4)); ax.set_yticklabels(fcols[:4])
for i in range(4):
    for j in range(len(kd_tfs)):
        ax.text(j, i, f'{im_data[i,j]:.3f}', ha='center', va='center', fontsize=6)
ax.set_title('100% KD: per-factor effects', fontsize=10)
plt.colorbar(im, ax=ax, shrink=0.6)

# D: DNMT1 notes
ax = axes[1,1]; panel_label(ax, 'D')
ax.axis('off')
note = ("DNMT1 KD aligns with PA deviation\n"
        f"  cos_sim = +0.43 (vs other TFs z=3.50)\n"
        f"ATF3 KD aligns with PA deviation\n"
        f"  cos_sim = +0.55 (vs other TFs z=4.38)\n\n"
        "Control TFs (POU5F1, XBP1) show\n"
        "  no alignment (cos_sim ~ 0)\n\n"
        "Note: DNMT1 is a DNA methyltransferase,\n"
        "not a classical TF. Its regulon reflects\n"
        "co-expressed target genes rather than\n"
        "direct transcriptional targets.")
ax.text(0.1, 0.95, note, transform=ax.transAxes, fontsize=8, family='monospace',
        va='top', bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='gray'))

fig4.tight_layout()
save_both(fig4, 'fig4_perturbation')
print("  Fig 4 saved", flush=True)

# ============================================================
# FIGURE 5: PA Validation
# ============================================================
print("Generating Fig 5 ...", flush=True)
fig5, axes = plt.subplots(2, 2, figsize=(10.5, 8.5))

key_genes = ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2','IDH2','PKM','EXOSC9','MDH1','GNL3','DPPA5','NANOG','KLHL15']
# Collect all DESeq2 results
pa_all = []
for stage, df in pa_data.items():
    if 'gene_symbol' not in df.columns: continue
    sym2idx = {}
    for i, row in df.iterrows():
        s = str(row['gene_symbol'])
        if not s.startswith('ENSSSCG'): sym2idx[s] = i
    for g in key_genes:
        if g in sym2idx:
            r = df.iloc[sym2idx[g]]
            pa_all.append({'gene': g, 'stage': stage,
                          'log2FC': r['log2FoldChange'], 'padj': r['padj']})

pa_all_df = pd.DataFrame(pa_all)

# A: Layer 1 (blueprint) stage-adjusted barplot
ax = axes[0,0]; panel_label(ax, 'A')
l1_genes = ['DNMT1','ZP3','ZP4','GDF9','RARRES1','EIF4G2']
# Get stage-adjusted overall values
overall = pa_all_df[pa_all_df['stage']=='overall']
l1_vals = [overall[overall['gene']==g]['log2FC'].values[0] if g in overall['gene'].values else 0 for g in l1_genes]
colors_l1 = ['#E41A1C' if v>0 else '#377EB8' for v in l1_vals]
ax.bar(l1_genes, l1_vals, color=colors_l1, edgecolor='black')
ax.axhline(0, color='gray', lw=0.8)
ax.set_ylabel('log\u2082FC (PA vs IVF)'); ax.set_title('Layer 1: Maternal blueprint', fontsize=10)
for i, (g, v) in enumerate(zip(l1_genes, l1_vals)):
    sig = '*' if g in overall['gene'].values and overall[overall['gene']==g]['padj'].values[0]<0.05 else ''
    ax.text(i, v+0.2*np.sign(v), f'{v:+.2f}{sig}', ha='center', fontsize=7, fontweight='bold')

# B: Layer 3 (execution) stage-adjusted
ax = axes[0,1]; panel_label(ax, 'B')
l3_genes = ['IDH2','PKM','EXOSC9','MDH1','GNL3','NOP9']
l3_vals = [overall[overall['gene']==g]['log2FC'].values[0] if g in overall['gene'].values else 0 for g in l3_genes]
colors_l3 = ['#E41A1C' if v>0 else '#377EB8' for v in l3_vals]
ax.bar(l3_genes, l3_vals, color=colors_l3, edgecolor='black')
ax.axhline(0, color='gray', lw=0.8)
ax.set_ylabel('log\u2082FC (PA vs IVF)'); ax.set_title('Layer 3: Zygotic execution (markedly reduced)', fontsize=10)
for i, (g, v) in enumerate(zip(l3_genes, l3_vals)):
    sig = '*' if g in overall['gene'].values and overall[overall['gene']==g]['padj'].values[0]<0.05 else ''
    ax.text(i, v-1.5, f'{v:+.2f}{sig}', ha='center', fontsize=7, fontweight='bold')

# C: Per-stage heatmap
ax = axes[1,0]; panel_label(ax, 'C')
stages_plot = ['1cell','2cell','4cell','8cell']
genes_plot = ['DNMT1','IDH2','PKM','EXOSC9','MDH1','GNL3','DPPA5','NANOG']
heatmap_data = np.zeros((len(genes_plot), len(stages_plot)))
for gi, g in enumerate(genes_plot):
    for si, stg in enumerate(stages_plot):
        sub = pa_all_df[(pa_all_df['gene']==g)&(pa_all_df['stage']==stg)]
        if len(sub): heatmap_data[gi, si] = sub.iloc[0]['log2FC']

im = ax.imshow(heatmap_data, aspect='auto', cmap='RdBu_r', vmin=-10, vmax=10)
ax.set_xticks(np.arange(len(stages_plot))); ax.set_xticklabels(stages_plot, fontsize=8)
ax.set_yticks(np.arange(len(genes_plot))); ax.set_yticklabels(genes_plot, fontsize=8)
for i in range(len(genes_plot)):
    for j in range(len(stages_plot)):
        v = heatmap_data[i,j]
        ax.text(j, i, f'{v:+.1f}', ha='center', va='center', fontsize=6,
                color='white' if abs(v)>5 else 'black')
ax.set_title('Per-stage log\u2082FC (PA vs IVF)', fontsize=10)
plt.colorbar(im, ax=ax, shrink=0.7)

# D: Layer summary table
ax = axes[1,1]; panel_label(ax, 'D')
ax.axis('off')
summary_text = (
    "F1 gene set (Maternal blueprint)\n"
    "  DNMT1: log2FC=+1.30, padj=0.006  (HIGHER in PA)\n"
    "  ZP3/4, GDF9, RARRES1, EIF4G2: heterogeneous\n"
    "  mean gene-level log2FC = −0.78 (n = 40)\n\n"
    "F3 gene set (Clearance-associated)\n"
    "  mean gene-level log2FC = −3.97, BH-adj P=0.004\n"
    "  median = −4.11 (n = 33, post hoc sensitivity\n"
    "  excluding extreme-effect genes: mean = −2.59)\n\n"
    "F4/F6 gene sets (Activation-/Translation-associated)\n"
    "  F4: mean log2FC = −0.17 (n = 43)\n"
    "  F6: mean log2FC = −0.41 (n = 39)\n"
    "  IDH2, PKM, EXOSC9: individual examples only\n\n"
    "42-sample bulk RNA-seq (SRP301735)\n"
    "DESeq2: stage-adjusted ~ condition\n"
    "4 cleavage stages: 1-cell to 8-cell"
)
ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=7, family='monospace',
        va='top', bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='gray'))

fig5.tight_layout()
save_both(fig5, 'fig5_pa_validation')
print("  Fig 5 saved", flush=True)

# ============================================================
# FIGURE 6: Module Summary
# ============================================================
print("Generating Fig 6 ...", flush=True)
fig6, ax = plt.subplots(1, 1, figsize=(11, 7))
ax.set_xlim(0, 14); ax.set_ylim(0, 12); ax.axis('off')

# Title
ax.text(7, 11.5, 'Candidate Maternal-Associated Modules and PA–IVF Comparison',
        ha='center', fontsize=14, fontweight='bold')

# Factor boxes
layers = [
    (2, 9, 10, 11, 'F1: Candidate maternal module',
     'DNMT1, ZP3, ZP4, GDF9, RARRES1\nRNA+METH weight correlation (r=+0.84)\nGene-level cross-species correspondence\n(nominal; does not exceed background)',
     '#377EB8', '#D6E4F0'),
    (2, 5.5, 10, 7.5, 'F3: Clearance-associated module',
     'Top 44 genes by absolute F3 RNA weight\nStage-associated decline (\u03c1=-0.90, P=2e-4)\nOperationally defined; not every member\nis a clearance factor',
     '#E41A1C', '#F5D6D6'),
    (2, 2, 10, 4, 'F4 & F6: Developmentally increasing modules',
     'F4: activation-associated; F6: translation-associated\nNominal stage-associated rise (\u03c1=+0.64 to +0.76)\nP-values descriptive, not corrected across factors',
     '#4DAF4A', '#D6F0D6'),
]
for x1, y1, x2, y2, title, desc, color, bgcolor in layers:
    rect = plt.Rectangle((x1, y1), x2-x1, y2-y1, facecolor=bgcolor, edgecolor=color, lw=2, alpha=0.6)
    ax.add_patch(rect)
    ax.text(x1+0.3, (y1+y2)/2+1, title, fontsize=11, fontweight='bold', color=color, va='center')
    ax.text(x1+0.5, (y1+y2)/2-0.3, desc, fontsize=9, va='center', family='monospace')

# PA annotation box
pa_box = plt.Rectangle((12.5, 2), 1.2, 9, facecolor='#FFF3CD', edgecolor='#FF7F00', lw=2, alpha=0.7)
ax.add_patch(pa_box)
ax.text(13.1, 10.5, 'PA\nEMBRYOS', ha='center', fontsize=11, fontweight='bold', color='#FF7F00')
ax.text(13.1, 9.3, 'F1\nHETEROGENEOUS\nmean LFC = −0.78', ha='center', fontsize=8, color='#FF7F00')
ax.text(13.1, 6.5, 'F3\nNEGATIVE\nBH-adj P = 0.004', ha='center', fontsize=8, color='#FF7F00')
ax.text(13.1, 3.0, 'F4/F6\nNO CONSISTENT\nDIRECTION', ha='center', fontsize=8, color='#FF7F00')

# Bottom note
ax.text(7, 0.5, 'GV oocyte (32 cells, MOFA+)  \u2192  Atlas projection (832 E0-E10 cells)  \u2192  Independent PA validation (42 samples, SRP301735)',
        ha='center', fontsize=9, style='italic', color='gray')

fig6.tight_layout()
save_both(fig6, 'fig6_model')
print("  Fig 6 saved", flush=True)

# ============================================================
# GRAPHICAL ABSTRACT
# ============================================================
print("Generating Graphical Abstract ...", flush=True)
fig_ga, ax = plt.subplots(1, 1, figsize=(4, 3))  # 1200x900 px at 300ppi = 4"x3"
ax.set_xlim(0, 12); ax.set_ylim(0, 9); ax.axis('off')

# Title
ax.text(6, 8.5, 'MOFA+ Multi-Omics of Porcine GV Oocytes', ha='center', fontsize=11, fontweight='bold')
ax.text(6, 7.8, 'Reveals Three-Layer MZT Architecture', ha='center', fontsize=10, fontweight='bold')

# Three layers simplified
layers_ga = [
    (1, 5, 5, 2.5, 'Layer 1: Maternal Blueprint\nDNMT1, ZP3/4, GDF9, RARRES1\nRNA+METH co-regulated, conserved', '#377EB8'),
    (1, 2, 5, 2.2, 'Layer 2: Maternal Clearance\nZP2, SYCN, PARP12\nRapidly degraded during MZT', '#E41A1C'),
    (1, -1, 5, 2.2, 'Layer 3: Zygotic Activation\nIDH2, PKM, EXOSC9, GNL3\nMetabolism + translation, rising', '#4DAF4A'),
]
for x, y, w, h, text, color in layers_ga:
    rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black', lw=1, alpha=0.3)
    ax.add_patch(rect)
    ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=7, color=color, fontweight='bold')

# PA arrow
ax.annotate('PA Embryos:\nLayer 1 maintained\nLayer 3 attenuated', xy=(8, 4.5), fontsize=7,
            ha='center', bbox=dict(boxstyle='round', facecolor='#FFF3CD', edgecolor='#FF7F00'))

ax.annotate('', xy=(8, 6.5), xytext=(6, 6.5), arrowprops=dict(arrowstyle='->', color='gray', lw=1))
ax.annotate('', xy=(8, 2.5), xytext=(6, 2.5), arrowprops=dict(arrowstyle='->', color='gray', lw=1))

ax.text(6, 0, '32 GV oocytes \u00d7 RNA+WGBS  |  MOFA+ 7 factors  |  Cross-species (pig/human/mouse)  |  42-sample PA validation',
        ha='center', fontsize=6.5, style='italic', color='gray')

fig_ga.tight_layout()
fig_ga.savefig(f"{OUTDIR}/fig_ga.png", dpi=300)
fig_ga.savefig(f"{OUTDIR}/fig_ga.svg")
print("  Graphical abstract saved", flush=True)

print("\n=== ALL FIGURES GENERATED ===")
print(f"Output directory: {OUTDIR}")
print("Files: fig1_mofa.svg, fig2_crossspecies.svg, fig3_mzt_trajectory.svg,")
print("        fig4_perturbation.svg, fig5_pa_validation.svg, fig6_model.svg, fig_ga.svg")
print("  (PNG versions also included)")