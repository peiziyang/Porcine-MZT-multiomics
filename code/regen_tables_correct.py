#!/usr/bin/env python
"""Regenerate all 6 supplementary tables with CORRECT filename-content mapping."""
import pandas as pd, os

ROOT = r'E:/Workbuddy/2026-07-27-11-58-27'
MOFA = os.path.join(ROOT, 'data/processed/m4_mofa')
DESEQ = os.path.join(ROOT, 'data/processed/deseq2')
OUT = os.path.join(ROOT, 'manuscript/submission/tables')
os.makedirs(OUT, exist_ok=True)

# Clear old (misaligned) files
for f in os.listdir(OUT):
    if f.endswith('.csv'):
        os.remove(os.path.join(OUT, f))

# ---- Table S1: MOFA+ factor summary ----
s1 = pd.read_csv(os.path.join(MOFA, 'factor_summary_v2.csv'))
s1.to_csv(os.path.join(OUT, 'Table_S1_mofa_factor_summary.csv'), index=False)
print(f'S1 MOFA+ factor summary: {len(s1)} rows, cols={list(s1.columns)[:6]}')

# ---- Table S2: GO enrichment ----
s2 = pd.read_csv(os.path.join(MOFA, 'factor_go_enrichment_v2.csv'))
s2.to_csv(os.path.join(OUT, 'Table_S2_go_enrichment.csv'), index=False)
print(f'S2 GO enrichment: {len(s2)} rows')

# ---- Table S3: DESeq2 results ----
s3 = pd.read_csv(os.path.join(DESEQ, 'm15_deseq2_overall_IVFvsPA_stageAdj.csv'))
s3.to_csv(os.path.join(OUT, 'Table_S3_deseq2_results.csv'), index=False)
print(f'S3 DESeq2 results: {len(s3)} rows')

# ---- Table S4: cross-species mapping ----
s4 = pd.read_csv(os.path.join(MOFA, 'cross_species_orthologs.csv'))
s4.to_csv(os.path.join(OUT, 'Table_S4_cross_species_mapping.csv'), index=False)
print(f'S4 cross-species mapping: {len(s4)} rows')

# ---- Table S5: CpG stats for F1 genes ----
cpg = pd.read_csv(os.path.join(MOFA, 'per_gene_cpg_stats.csv'))
w = pd.read_csv(os.path.join(MOFA, 'mofa_multiomics_weights_RNA_v2.csv'), index_col=0)
f1_genes = list(w['F1'].abs().sort_values(ascending=False).head(44).index)
s5 = cpg[cpg['gene_name'].isin(f1_genes)].copy()
s5.to_csv(os.path.join(OUT, 'Table_S5_cpg_stats.csv'), index=False)
print(f'S5 CpG stats (F1 genes): {len(s5)} rows')

# ---- Table S6: gene-set membership ----
deseq = pd.read_csv(os.path.join(DESEQ, 'm15_deseq2_overall_IVFvsPA_stageAdj.csv'))
# Build symbol -> lfc/baseMean/padj map
sym2lfc, sym2bm, sym2padj = {}, {}, {}
for _, r in deseq.iterrows():
    sym = str(r.get('gene_symbol', ''))
    if sym and not sym.startswith('ENSSSCG'):
        sym2lfc[sym] = r.get('log2FoldChange')
        sym2bm[sym] = r.get('baseMean')
        sym2padj[sym] = r.get('padj')

rows = []
for fc in ['F1', 'F3', 'F4', 'F6']:
    genes = list(w[fc].abs().sort_values(ascending=False).head(44).index)
    for rank, g in enumerate(genes):
        rows.append({
            'factor': fc, 'rank': rank + 1, 'gene': g,
            'F1_member': fc == 'F1', 'F3_member': fc == 'F3',
            'F4_member': fc == 'F4', 'F6_member': fc == 'F6',
            'log2FoldChange': sym2lfc.get(g, None),
            'baseMean': sym2bm.get(g, None),
            'padj': sym2padj.get(g, None),
        })
s6 = pd.DataFrame(rows)
s6.to_csv(os.path.join(OUT, 'Table_S6_gene_set_membership.csv'), index=False)
print(f'S6 gene-set membership: {len(s6)} rows')

print('\n=== 验证：文件名 vs 内容 ===')
for fn in sorted(os.listdir(OUT)):
    df = pd.read_csv(os.path.join(OUT, fn), nrows=2)
    cols = list(df.columns)[:5]
    print(f'  {fn} → {cols}')
