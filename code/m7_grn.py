"""
M7: Gene Regulatory Network — Regulon Co-expression + Activity
===============================================================
简版 SCENIC：共表达网络定义 TF regulon → 沿伪时序评估活动度。

用 M5 体内 E0-E8 数据（8276 common genes），spearman 相关建共表达网络，
已知 TF 列表定义每个 TF + top 100 co-expressed = regulon，
AUCell 风格计算 regulon 活动度，沿发育阶段热图展示。
"""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw'
OUT_DIR = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed'

print("="*60)
print("M7: GRN — Co-expression Regulon Analysis")
print("="*60)

# ============================================================
# 1. Load data + Known TF list
# ============================================================
print("\n[1/5] Loading data and TF list...")

# Load previously mapped gene symbols
# From M6 we have gene_to_symbol. Let me rebuild quickly.
import mygene
mg = mygene.MyGeneInfo()

counts_vivo = pd.read_csv(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_genes_counts.txt.gz', sep='\t', index_col=0)
si = pd.read_excel(f'{DATA_DIR}/GSE168106/GSE168106_scRNA-seq_SampleInfo.xlsx')
si_v = si[si['Vaild']=='Yes']
cell_to_day = dict(zip(si_v['ID'], si_v['Stage_II']))
early_days = ['E0','E1','E2','E3','E4','E5','E6','E7','E8']
vivo_cells = [c for c in counts_vivo.columns if c in cell_to_day and pd.notna(cell_to_day[c]) and cell_to_day[c] in early_days]

# Filter: genes with expression in >= 10% cells
expr_rate = (counts_vivo[vivo_cells] > 5).mean(axis=1)
genes_use = counts_vivo.index[expr_rate >= 0.10]
print(f"  Genes with expr in >=10% cells: {len(genes_use)}")

# Map to symbols (just for the filtered set)
print("  Mapping gene symbols...")
genes_ens = [g for g in genes_use if g.startswith('ENSSSCG')]
gene_to_symbol = {}
CHUNK = 800
for i in range(0, len(genes_ens), CHUNK):
    chunk = genes_ens[i:i+CHUNK]
    results = mg.querymany(chunk, scopes='ensembl.gene', species='pig', fields='symbol')
    for r in results:
        if 'symbol' in r:
            gene_to_symbol[r['query']] = r['symbol']
print(f"  Mapped {len(gene_to_symbol)} symbols")

# Known mammalian TFs and key developmental regulators
# Comprehensive list from AnimalTFDB + literature
KNOWN_TFS = {
    # Pluripotency / early embryo
    'POU5F1','NANOG','SOX2','KLF4','MYC','SALL4','ZFP42','PRDM14','TBX3',
    'ESRRB','NR5A2','TFCP2L1','DPPA3','DPPA4','UTF1','ZSCAN4','DUX4',
    # Trophectoderm / extraembryonic
    'CDX2','TEAD4','GATA3','GATA2','TFAP2C','TFAP2A','EOMES','ELF5','ETS2',
    # Primitive endoderm / hypoblast
    'GATA6','SOX17','PDGFRA','HNF4A',
    # Mesendoderm / gastrulation
    'MIXL1','T','EOMES','FOXA2','GSC','MESP1','MESP2','TBX6',
    # Epiblast / ectoderm
    'OTX2','PAX6','SOX1','GBX2','LHX2',
    # General developmental TFs  
    'FOXD3','DNMT1','DNMT3A','DNMT3B','DNMT3L','TET1','TET2','TET3',
    # Cell cycle / ZGA-related
    'E2F1','E2F4','MYBL2','LIN28A','LIN28B','DPPA2','DPPA4',
    # Wnt / FGF / TGFb pathway TFs
    'TCF3','TCF4','TCF7L2','LEF1','SMAD2','SMAD3','SMAD4','SMAD1','SMAD5',
    'SMAD9','STAT3','FOXO1','FOXO3','FOXO4',
    # Nuclear receptors
    'ESR1','ESR2','NR2F2','PPARG','RARG','RXRA','NR3C1',
    # Homeobox
    'HOXA1','HOXB1','HOXC4','HOXD4','MSX1','MSX2','DLX5','SHOX',
    # Forkhead
    'FOXA1','FOXC1','FOXD1','FOXH1','FOXJ1','FOXM1','FOXO1',
    'FOXO3','FOXP1','FOXP2',
    # SOX family
    'SOX1','SOX2','SOX3','SOX4','SOX5','SOX6','SOX7','SOX9',
    'SOX10','SOX11','SOX15','SOX17','SOX18','SOX21',
    # GATA family
    'GATA1','GATA2','GATA3','GATA4','GATA5','GATA6',
    # PAX family
    'PAX1','PAX2','PAX3','PAX5','PAX6','PAX7','PAX8','PAX9',
    # POU family
    'POU1F1','POU2F1','POU2F2','POU3F1','POU3F2','POU4F1','POU5F1',
    'POU5F2','POU6F1',
    # TBX family
    'TBX1','TBX2','TBX3','TBX4','TBX5','TBX6','TBX10','TBX15',
    'TBX18','TBX20','TBX21','TBX22',
    # IRX, LHX, NKX families
    'LHX1','LHX2','LHX3','LHX4','LHX5','LHX6','LHX8','LHX9',
    'NKX2-1','NKX2-5','NKX6-1',
    # bHLH factors
    'HAND1','HAND2','TWIST1','TWIST2','TCF12','TCF21','TAL1',
    'HES1','HES5','HEY1','HEY2','BHLHE40','BHLHE41',
    # Others
    'CUX1','CUX2','HMGA2','HMGB1','HMGB2','ID1','ID2','ID3','ID4',
    'NFYA','NFYB','NFYC','YY1','YY2','CTCF','RAD21','SMC3',
    # ZGA drivers (conserved)
    'ZNF280A','ZNF280C','ZNF609','ZBTB14','ZBTB24',
    # PRRS/swine specific interest
    'IRF1','IRF2','IRF3','IRF7','IRF9','NFKB1','NFKB2','RELA',
    'STAT1','STAT2','STAT6',
    # Chromatin remodelers
    'SMARCA4','SMARCA2','SMARCC1','SMARCD1','CHD1','CHD2','CHD4',
    'CHD7','CHD8','ARID1A','ARID1B','ARID2','BAZ1A','BAZ2A',
    # Methylation-related  
    'UHRF1','UHRF2','MBD3','MECP2','ZFP57','TRIM28','SETDB1',
}

# Map TF symbols to Ensembl IDs
symbol_to_ens = {v: k for k, v in gene_to_symbol.items()}
tf_ens = [symbol_to_ens[s] for s in KNOWN_TFS if s in symbol_to_ens]
print(f"  TFs found in data: {len(tf_ens)}/{len(KNOWN_TFS)}")
print(f"  Example TFs: {[(e, gene_to_symbol[e]) for e in tf_ens[:10]]}")

# ============================================================
# 2. Expression matrix preparation
# ============================================================
print("\n[2/5] Preparing expression matrix...")
# Use genes that are somewhat variable
expr = counts_vivo.loc[genes_ens, vivo_cells]
log_expr = np.log1p(expr)

# Select top variable genes for network (speed)
gene_var = log_expr.var(axis=1)
top_var_genes = gene_var.nlargest(2000).index
# Include all found TFs
genes_for_network = sorted(set(top_var_genes) | set(tf_ens))
genes_for_network = [g for g in genes_for_network if g in log_expr.index]

print(f"  Genes for network: {len(genes_for_network)} (2000 top var + {len(tf_ens)} TFs)")
log_sub = log_expr.loc[genes_for_network]

# ============================================================
# 3. Co-expression network
# ============================================================
print("\n[3/5] Computing co-expression network...")
# For each TF, find top 100 co-expressed target genes
# Use spearman correlation (robust to outliers)

# Compute correlation between TFs and all other genes
tf_set = set(tf_ens)
other_genes = [g for g in genes_for_network if g not in tf_set]

tf_regulons = {}
n_tfs = len(tf_ens)
for idx, tf in enumerate(tf_ens):
    tf_expr = log_sub.loc[tf].values
    corrs = []
    for target in other_genes:
        target_expr = log_sub.loc[target].values
        if len(set(tf_expr)) > 1 and len(set(target_expr)) > 1:
            r, _ = spearmanr(tf_expr, target_expr)
            if not np.isnan(r):
                corrs.append((target, abs(r)))
    # Top 100 by absolute correlation
    corrs.sort(key=lambda x: x[1], reverse=True)
    tf_regulons[tf] = [g for g, r in corrs[:100]]
    if (idx+1) % 50 == 0:
        print(f"  Processed {idx+1}/{n_tfs} TFs...")

print(f"  Built {len(tf_regulons)} regulons (avg {np.mean([len(v) for v in tf_regulons.values()]):.0f} targets)")

# ============================================================
# 4. Regulon activity (AUCell-like)
# ============================================================
print("\n[4/5] Computing regulon activity...")

# For each regulon, compute mean expression of target genes per cell
regulon_activity = {}
for tf, targets in tf_regulons.items():
    avail = [t for t in targets if t in log_sub.index]
    if len(avail) >= 10:
        # Mean z-score of target genes
        sub = log_sub.loc[avail]
        z_scores = (sub - sub.mean(axis=1).values.reshape(-1,1)) / (sub.std(axis=1).values.reshape(-1,1) + 1e-8)
        regulon_activity[tf] = z_scores.mean(axis=0)

print(f"  Regulons with >=10 targets: {len(regulon_activity)}")

# Aggregate by developmental stage
stage_order = ['E0','E1','E2','E3','E4','E5','E6','E7','E8']
regulon_by_stage = {}
for tf, activity in regulon_activity.items():
    by_stage = {}
    for s in stage_order:
        cells_s = [c for c in activity.index if cell_to_day.get(c) == s]
        if cells_s:
            by_stage[s] = activity[cells_s].mean()
    regulon_by_stage[tf] = by_stage

# Convert to matrix
regulon_matrix = pd.DataFrame(regulon_by_stage).T  # regulons × stages
print(f"  Activity matrix: {regulon_matrix.shape}")

# Find dynamic regulons (high variance across stages)
regulon_var = regulon_matrix.var(axis=1)
dynamic_regulons = regulon_var.nlargest(30).index
print(f"  Top 30 most dynamic regulons")

# ============================================================
# 5. Visualization
# ============================================================
print("\n[5/5] Generating Fig.5...")

fig, axes = plt.subplots(2, 3, figsize=(20, 14))

# A: Regulon activity heatmap (top 30 dynamic regulons)
ax = axes[0, 0]
heat_data = regulon_matrix.loc[dynamic_regulons]
# Normalize per row for visualization
heat_norm = heat_data.apply(lambda x: (x - x.min())/(x.max()-x.min()+1e-8), axis=1)
# Sort by peak stage
peak_stage = heat_norm.idxmax(axis=1)
heat_sorted = heat_norm.loc[peak_stage.sort_values().index]

im = ax.imshow(heat_sorted.values, aspect='auto', cmap='RdBu_r', vmin=0, vmax=1)
ax.set_xticks(range(len(stage_order)))
ax.set_xticklabels(stage_order, fontsize=9)
tf_labels = [gene_to_symbol.get(g, g[:15]) for g in heat_sorted.index]
ax.set_yticks(range(len(tf_labels)))
ax.set_yticklabels(tf_labels, fontsize=7)
ax.set_xlabel('Developmental Stage')
ax.set_title('A. Dynamic Regulon Activity Across Development\n(Top 30 TFs by variance, normalized)')
plt.colorbar(im, ax=ax, shrink=0.8, label='Normalized Activity')

# B: Selected key embryo regulon trajectories
ax = axes[0, 1]
key_tf_of_interest = [tf for tf in tf_ens if gene_to_symbol.get(tf, '') in 
    ['POU5F1','NANOG','SOX2','CDX2','GATA3','GATA6','MYC','KLF4',
     'TEAD4','EOMES','DNMT3A','LIN28A','ZSCAN4','TFAP2C','SOX17','T']]
key_tf_found = [tf for tf in key_tf_of_interest if tf in regulon_matrix.index]
for tf in key_tf_found[:10]:
    vals = regulon_matrix.loc[tf, stage_order]
    ax.plot(range(len(stage_order)), vals.values, 'o-', markersize=4, linewidth=2,
           label=gene_to_symbol.get(tf, tf[:12]), alpha=0.8)
ax.set_xticks(range(len(stage_order)))
ax.set_xticklabels(stage_order, fontsize=8)
ax.set_xlabel('Stage')
ax.set_ylabel('Regulon Activity')
ax.set_title('B. Key Embryo TF Regulon Dynamics')
ax.legend(fontsize=6, ncol=2, loc='upper left')

# C: Regulon activity UMAP (colored by top TF activity)
ax = axes[0, 2]
# Use M2 PCA coords if available, else do simple PCA
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
X_pca_input = StandardScaler().fit_transform(log_sub.T.values)
pca_quick = PCA(n_components=2, random_state=42)
X_p2 = pca_quick.fit_transform(X_pca_input)
# Color by NANOG regulon activity if available
nanog_ens = symbol_to_ens.get('NANOG')
if nanog_ens and nanog_ens in regulon_activity:
    nanog_act = regulon_activity[nanog_ens]
    sc = ax.scatter(X_p2[:,0], X_p2[:,1], c=nanog_act.values, cmap='RdBu_r',
                   s=10, alpha=0.7, rasterized=True)
    plt.colorbar(sc, ax=ax, label='NANOG Regulon Activity')
else:
    ax.scatter(X_p2[:,0], X_p2[:,1], c='grey', s=10, alpha=0.5)
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_title('C. NANOG Regulon Activity on PCA Embedding\n(early embryo data)')

# D: Network visualization (top TFs and their connections)
ax = axes[1, 0]
# Show connectivity of top 10 dynamic TFs
top_tfs = dynamic_regulons[:10]
tf_summary = []
for tf in top_tfs:
    sym = gene_to_symbol.get(tf, tf[:12])
    n_targets = len(tf_regulons.get(tf, []))
    tf_summary.append(f'{sym}: {n_targets} targets')

y_pos = np.arange(len(tf_summary))
ax.barh(y_pos, [len(tf_regulons.get(tf, [])) for tf in top_tfs])
ax.set_yticks(y_pos)
ax.set_yticklabels(tf_summary, fontsize=8)
ax.set_xlabel('Number of Co-expressed Targets (|r| ranked top 100)')
ax.set_title('D. Top Dynamic TF Regulon Sizes')

# E: Stage-specific TF activity peak
ax = axes[1, 1]
stage_tf_count = {}
for s in stage_order:
    top_tf_stage = regulon_matrix[s].nlargest(5).index
    stage_tf_count[s] = [gene_to_symbol.get(tf, tf[:10]) for tf in top_tf_stage]

# Simple bubble-like display
for i, s in enumerate(stage_order):
    tfs_at_stage = stage_tf_count[s]
    for j, tf_name in enumerate(tfs_at_stage):
        ax.text(i, 4-j, tf_name[:10], ha='center', fontsize=7,
               bbox=dict(boxstyle='round,pad=0.3', facecolor=plt.cm.tab20(i%20), alpha=0.6))
ax.set_xlim(-0.5, len(stage_order)-0.5)
ax.set_ylim(-0.5, 4.5)
ax.set_xticks(range(len(stage_order)))
ax.set_xticklabels(stage_order, fontsize=9)
ax.set_yticks([])
ax.set_title('E. Top 5 TFs per Developmental Stage\n(by regulon activity)')

# F: Summary
ax = axes[1, 2]
ax.axis('off')
summary = [
    "M7: GRN — Co-expression Regulon",
    "=================================",
    f"TFs identified: {len(tf_ens)}",
    f"  (from {len(KNOWN_TFS)} curated list)",
    f"Genes in network: {len(genes_for_network)}",
    f"Regulons built: {len(tf_regulons)}",
    f"Regulons with >=10 targets: {len(regulon_activity)}",
    f"Most dynamic regulons: {len(dynamic_regulons)}",
    "",
    "Method: Spearman-based co-expression",
    "TF regulon = TF + top 100 abs(r)",
    "Activity = mean z-score of targets",
    "",
    "Analogy: SCENIC Step1 (GRNBoost2)",
    "without Step2 (cisTarget motif)",
    "(awaiting pig motif database)",
    "",
    "Note: Full SCENIC+ with enhancer",
    "prediction requires scATAC-seq data",
    "(e.g., from PRJEB81663, pending)"
]
for i, line in enumerate(summary):
    ax.text(0.02, 0.98-i*0.047, line, fontsize=8, family='monospace',
           transform=ax.transAxes, verticalalignment='top')

plt.suptitle('Fig.5: Gene Regulatory Network of Porcine Preimplantation Development\n(Co-expression Regulon Analysis)',
            fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(f'{OUT_DIR}/m7_grn_regulon.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: m7_grn_regulon.png")

# ============================================================
# 6. Save
# ============================================================
regulon_matrix.index = [f"{gene_to_symbol.get(g, g)} ({g})" for g in regulon_matrix.index]
regulon_matrix.to_csv(f'{OUT_DIR}/m7_regulon_activity.csv')

pd.DataFrame({
    'TF_ensembl': list(tf_regulons.keys()),
    'TF_symbol': [gene_to_symbol.get(t, t) for t in tf_regulons.keys()],
    'n_targets': [len(v) for v in tf_regulons.values()]
}).to_csv(f'{OUT_DIR}/m7_tf_regulons.csv', index=False)

print(f"\nDone! Outputs: {OUT_DIR}/")
print(f"  m7_grn_regulon.png")
print(f"  m7_regulon_activity.csv")
print(f"  m7_tf_regulons.csv")
print(f"\n  Top 10 most dynamic TFs:")
for tf in dynamic_regulons[:10]:
    sym = gene_to_symbol.get(tf, tf)
    var = regulon_var[tf]
    print(f"    {sym}: var={var:.4f}")
