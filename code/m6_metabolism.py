"""
M6: 代谢模块评分 — glycolysis vs OxPhos 发育轮廊
================================================
效仿 Malkowska 2022 Nat Commun hexa-species approach:
per-cell module score for glycolysis and oxidative phosphorylation
across developmental stages and conditions.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mygene
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw'
OUT_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'

print("="*60)
print("M6: Metabolic Module Scoring (Glycolysis vs OxPhos)")
print("="*60)

# ============================================================
# 1. Load data
# ============================================================
print("\n[1/5] Loading data...")
counts_vivo = pd.read_csv(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_genes_counts.txt.gz', sep='\t', index_col=0)
si = pd.read_excel(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_SampleInfo.xlsx')
si_v = si[si['Vaild']=='Yes']
cell_to_day = dict(zip(si_v['ID'], si_v['Stage_II']))
cell_to_lane = dict(zip(si_v['ID'], si_v['Lane']))

valid_cells = [c for c in counts_vivo.columns if c in cell_to_day and pd.notna(cell_to_day[c])]
counts_vivo = counts_vivo[valid_cells]
print(f"  GSE168106: {len(valid_cells)} cells, {counts_vivo.shape[0]} genes")

# Also load GSE164812 for IVF/PA comparison
fpkm_paivf = pd.read_csv(f'{DATA_DIR}/GSE164812/GSE164812_gene_FPKM_matrix.txt.gz', sep='\t', index_col=0)
fpkm_paivf = fpkm_paivf.apply(pd.to_numeric, errors='coerce').fillna(0)
print(f"  GSE164812: {fpkm_paivf.shape[1]} cells")

# ============================================================
# 2. Gene symbol mapping
# ============================================================
print("\n[2/5] Mapping pig Ensembl IDs to gene symbols...")
mg = mygene.MyGeneInfo()

# Get all unique genes
all_genes = sorted(set(counts_vivo.index) | set(fpkm_paivf.index))
all_genes_ens = [g for g in all_genes if g.startswith('ENSSSCG')]
print(f"  Total pig Ensembl genes: {len(all_genes_ens)}")

# Batch query (mygene has limits, do in chunks of 1000)
gene_to_symbol = {}
CHUNK = 1000
for i in range(0, len(all_genes_ens), CHUNK):
    chunk = all_genes_ens[i:i+CHUNK]
    try:
        results = mg.querymany(chunk, scopes='ensembl.gene', species='pig', fields='symbol')
        for r in results:
            if 'symbol' in r:
                gene_to_symbol[r['query']] = r['symbol']
    except Exception as e:
        print(f"  Chunk {i//CHUNK}: error {e}")
    if i % 5000 == 0:
        print(f"  Processed {i}/{len(all_genes_ens)}...")

print(f"  Mapped {len(gene_to_symbol)} genes to symbols")
print(f"  Examples: ENSSSCG00000000002 -> {gene_to_symbol.get('ENSSSCG00000000002','?')}")

# ============================================================
# 3. Define glycolysis and OxPhos gene sets
# ============================================================
print("\n[3/5] Defining metabolic gene sets...")

# KEGG Glycolysis / Gluconeogenesis (hsa00010) core enzymes
GLYCOLYSIS_GENES = [
    'HK1','HK2','HK3','GCK',           # hexokinase
    'GPI',                               # glucose-6-phosphate isomerase
    'PFKL','PFKM','PFKP',               # phosphofructokinase
    'ALDOA','ALDOB','ALDOC',            # aldolase
    'TPI1',                              # triosephosphate isomerase
    'GAPDH','GAPDHS',                    # glyceraldehyde-3-phosphate dehydrogenase
    'PGK1','PGK2',                       # phosphoglycerate kinase
    'PGAM1','PGAM2','PGAM4',            # phosphoglycerate mutase
    'ENO1','ENO2','ENO3','ENO4',        # enolase
    'PKLR','PKM',                        # pyruvate kinase
    'LDHA','LDHB','LDHC','LDHD',        # lactate dehydrogenase
    'PDK1','PDK2','PDK3','PDK4',        # pyruvate dehydrogenase kinase
    'SLC2A1','SLC2A2','SLC2A3','SLC2A4', # glucose transporter
    'SLC16A1','SLC16A3',                 # monocarboxylate transporter (MCT1, MCT4)
]

# KEGG Oxidative Phosphorylation (hsa00190) core subunits
OXPHOS_GENES = [
    # Complex I (NADH dehydrogenase)
    'NDUFV1','NDUFV2','NDUFV3','NDUFS1','NDUFS2','NDUFS3',
    'NDUFS4','NDUFS5','NDUFS6','NDUFS7','NDUFS8',
    'NDUFA1','NDUFA2','NDUFA3','NDUFA4','NDUFA5','NDUFA6',
    'NDUFA7','NDUFA8','NDUFA9','NDUFA10','NDUFAB1','NDUFB1',
    # Complex II (succinate dehydrogenase)
    'SDHA','SDHB','SDHC','SDHD',
    # Complex III (cytochrome bc1)
    'UQCRC1','UQCRC2','UQCRFS1','UQCRB','UQCRQ','UQCRH',
    'UQCR10','UQCR11','CYC1',
    # Complex IV (cytochrome c oxidase)
    'COX4I1','COX5A','COX5B','COX6A1','COX6B1','COX6C',
    'COX7A2','COX7B','COX7C','COX8A',
    # Complex V (ATP synthase)
    'ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D','ATP5F1E',
    'ATP5PO','ATP5PF','ATP5PD','ATP5PB','ATP5MC1','ATP5MC2','ATP5MC3',
    'ATP5MG','ATP5MF','ATP5ME',
]

# Map gene symbols to Ensembl IDs (reverse of gene_to_symbol)
symbol_to_ens = {v: k for k, v in gene_to_symbol.items()}

glycolysis_ens = [symbol_to_ens[s] for s in GLYCOLYSIS_GENES if s in symbol_to_ens]
oxphos_ens = [symbol_to_ens[s] for s in OXPHOS_GENES if s in symbol_to_ens]

print(f"  Glycolysis genes found: {len(glycolysis_ens)}/{len(GLYCOLYSIS_GENES)}")
print(f"  OxPhos genes found: {len(oxphos_ens)}/{len(OXPHOS_GENES)}")

# ============================================================
# 4. Compute module scores
# ============================================================
print("\n[4/5] Computing per-cell module scores...")

def module_score(mat, gene_set):
    """Simple mean expression of gene set per cell."""
    available = [g for g in gene_set if g in mat.index]
    if len(available) == 0:
        return pd.Series(0, index=mat.columns)
    sub = mat.loc[available]
    return np.log1p(sub).mean(axis=0)

# GSE168106 (counts)
gly_score_vivo = module_score(counts_vivo, glycolysis_ens)
oxphos_score_vivo = module_score(counts_vivo, oxphos_ens)

# GSE164812 (FPKM, log2+1)
fpkm_log = np.log2(fpkm_paivf + 1)
gly_score_paivf = module_score(fpkm_log, glycolysis_ens)
oxphos_score_paivf = module_score(fpkm_log, oxphos_ens)

# ============================================================
# 5. Build result table + visualization
# ============================================================
print("\n[5/5] Building Fig.6 panels...")

# Combine scores
vivo_df = pd.DataFrame({
    'cell': gly_score_vivo.index,
    'glycolysis': gly_score_vivo.values,
    'oxphos': oxphos_score_vivo.values,
    'stage': [cell_to_day.get(c, '?') for c in gly_score_vivo.index],
    'condition': 'in_vivo'
})

cond_labels = []
for c in fpkm_paivf.columns:
    cond_labels.append(c.split(' ')[0] if ' ' in c else '?')
paivf_df = pd.DataFrame({
    'cell': gly_score_paivf.index,
    'glycolysis': gly_score_paivf.values,
    'oxphos': oxphos_score_paivf.values,
    'stage': 'embryo',
    'condition': cond_labels
})

# Stage order
stage_order_vivo = ['E0','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10',
                    'E11','E12','E13','E14']

# ============================================================
# 6. Visualization
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# A: Glycolysis vs OxPhos scatter (all cells, colored by stage)
ax = axes[0, 0]
for s in stage_order_vivo[:9]:
    mask = vivo_df['stage'] == s
    if mask.sum() > 1:
        ax.scatter(vivo_df.loc[mask, 'glycolysis'].values,
                  vivo_df.loc[mask, 'oxphos'].values,
                  s=8, alpha=0.6, label=s, rasterized=True)
ax.set_xlabel('Glycolysis Score')
ax.set_ylabel('OxPhos Score')
ax.set_title('A. Metabolic State of In Vivo Pig Embryos\n(E0-E8, colored by stage)')
ax.legend(fontsize=7, ncol=2)

# B: Metabolic ratio (glycolysis/OxPhos) across stages
ax = axes[0, 1]
ratio_data = []
for s in stage_order_vivo:
    mask = vivo_df['stage'] == s
    if mask.sum() > 1:
        gly = vivo_df.loc[mask, 'glycolysis'].mean()
        ox = vivo_df.loc[mask, 'oxphos'].mean()
        ratio_data.append({'stage': s, 'gly_mean': gly, 'ox_mean': ox, 'ratio': gly/(ox+1e-8), 'n': mask.sum()})
ratio_df = pd.DataFrame(ratio_data)

x = range(len(ratio_df))
ax.plot(x, ratio_df['gly_mean'], 'o-', color='#e74c3c', linewidth=2, label='Glycolysis')
ax.plot(x, ratio_df['ox_mean'], 's-', color='#3498db', linewidth=2, label='OxPhos')
ax.set_xticks(x)
ax.set_xticklabels(ratio_df['stage'], fontsize=8)
ax.set_xlabel('Developmental Stage')
ax.set_ylabel('Mean Module Score')
ax.set_title('B. Glycolysis vs OxPhos Across Development\n(in vivo, E0-E14)')
ax.legend(fontsize=9)

# C: Glycolysis/OxPhos ratio
ax = axes[0, 2]
ax.bar(x, ratio_df['ratio'], color=['#e74c3c' if r>1 else '#3498db' for r in ratio_df['ratio']], alpha=0.7)
ax.axhline(1, color='grey', linestyle='--', alpha=0.5)
ax.set_xticks(x)
ax.set_xticklabels(ratio_df['stage'], fontsize=8)
ax.set_ylabel('Glycolysis / OxPhos Ratio')
ax.set_title('C. Metabolic Switch: Glycolysis/OxPhos Ratio\n(ratio > 1 = glycolytic, < 1 = oxidative)')

# D: IVF vs PA metabolic comparison
ax = axes[1, 0]
for cond, c, m in [('IVF', '#d62728', 'o'), ('PA', '#1f77b4', '^')]:
    mask = paivf_df['condition'] == cond
    if mask.sum() > 1:
        ax.scatter(paivf_df.loc[mask, 'glycolysis'].values,
                  paivf_df.loc[mask, 'oxphos'].values,
                  s=40, alpha=0.7, c=c, marker=m, edgecolors='white',
                  linewidth=0.5, label=f'{cond} (n={mask.sum()})')
ax.set_xlabel('Glycolysis Score')
ax.set_ylabel('OxPhos Score')
ax.set_title('D. IVF vs PA Metabolic State\n(Early Cleavage Stage)')
ax.legend(fontsize=9)

# E: Violin/box comparison of glycolysis
ax = axes[1, 1]
data_plot = [paivf_df[paivf_df['condition']=='IVF']['glycolysis'].values,
             paivf_df[paivf_df['condition']=='PA']['glycolysis'].values,
             vivo_df[vivo_df['stage'].isin(['E2','E3','E4','E5'])]['glycolysis'].values]
bp = ax.boxplot(data_plot, labels=['IVF','PA','in_vivo'], patch_artist=True)
for box, c in zip(bp['boxes'], ['#d62728','#1f77b4','#2ca02c']):
    box.set_facecolor(c)
    box.set_alpha(0.6)
ax.set_ylabel('Glycolysis Score')
ax.set_title('E. Glycolysis: IVF vs PA vs in_vivo\n(Matched Stages E2-E5)')

# F: Summary panel
ax = axes[1, 2]
ax.axis('off')
summary_text = [
    "M6: Metabolic Module Scoring",
    "=============================",
    f"Glycolysis genes: {len(glycolysis_ens)}",
    f"OxPhos genes: {len(oxphos_ens)}",
    f"Method: per-cell mean log expr",
    f"  of KEGG pathway genes",
    "",
    "Key findings:",
    f"  Gly/Phos ratio range:",
    f"    {ratio_df['ratio'].min():.2f} - {ratio_df['ratio'].max():.2f}",
    f"  Mean ratio (E0-E8):",
    f"    {ratio_df['ratio'].mean():.2f}",
    "",
    "Cross-species comparison:",
    "  Malkowska 2022 (hexa-species)",
    "  reported conserved switch from",
    "  bivalent → glycolytic at late",
    "  blastocyst. Pig data can be",
    "  added as the 7th species.",
    "",
    "[Note: scMetabolism only supports",
    "  human; pig genes mapped via",
    "  mygene Ensembl→symbol query]"
]
for i, line in enumerate(summary_text):
    ax.text(0.02, 0.98-i*0.048, line, fontsize=8, family='monospace',
           transform=ax.transAxes, verticalalignment='top')

plt.suptitle('Fig.6: Metabolic Profiling of Porcine Preimplantation Embryos\n(Glycolysis vs OxPhos Module Scores)',
            fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(f'{OUT_DIR}/m6_metabolism.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: m6_metabolism.png")

# ============================================================
# 7. Save results
# ============================================================
vivo_df.to_csv(f'{OUT_DIR}/m6_metabolism_vivo.csv', index=False)
paivf_df.to_csv(f'{OUT_DIR}/m6_metabolism_paivf.csv', index=False)
ratio_df.to_csv(f'{OUT_DIR}/m6_metabolism_ratio.csv', index=False)

print(f"\nDone! Outputs: {OUT_DIR}/")
print(f"  m6_metabolism.png")
print(f"  m6_metabolism_vivo.csv")
print(f"  m6_metabolism_paivf.csv")
print(f"  m6_metabolism_ratio.csv")
print(f"\n  Metabolic ratio (Gly/OxPhos) by stage:")
for _, row in ratio_df.iterrows():
    bar = '█' * int(row['ratio']*10)
    print(f"    {row['stage']:4s} (n={row['n']:3d}): ratio={row['ratio']:.3f} {bar}")
