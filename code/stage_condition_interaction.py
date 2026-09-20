#!/usr/bin/env python
"""
Stage x condition interaction test for the PA/IVF bulk RNA-seq cohort (SRP301735).

Design: ~ stage + condition + stage:condition
Reports:
  - overall condition main effect (Wald, condition[T.PA] coefficient in full model = PA effect at reference stage 1cell)
  - per-interaction-coefficient Wald tests (stage[T.x]:condition[T.PA])
  - number of genes with significant interaction (padj<0.05) per coefficient
  - top genes with strongest interaction
Saves results to data/processed/deseq2/.
"""
import os, re, gzip, warnings
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
warnings.filterwarnings('ignore')

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
OUTDIR = os.path.join(BASE, "data", "processed", "salmon_out")
GTF = os.path.join(BASE, "data", "raw", "pig_genes.gtf.gz")
METAF = os.path.join(OUTDIR, "sra_metadata.tsv")
RES = os.path.join(BASE, "data", "processed", "deseq2")
os.makedirs(RES, exist_ok=True)

# ---- 1. tx2gene ----
tx2gene_base, gene_symbol = {}, {}
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if 'transcript_id' not in line:
            continue
        tid_m = re.search(r'transcript_id "([^"]+)"', line)
        gid_m = re.search(r'gene_id "([^"]+)"', line)
        gnm_m = re.search(r'gene_name "([^"]+)"', line)
        if tid_m and gid_m:
            base = tid_m.group(1).split('.')[0]
            gid = gid_m.group(1)
            tx2gene_base[base] = gid
            gene_symbol[gid] = gnm_m.group(1) if gnm_m else gid
print(f"tx2gene entries: {len(tx2gene_base)}")

# ---- 2. gene-level counts ----
quant_dirs = [d for d in os.listdir(OUTDIR)
              if os.path.exists(os.path.join(OUTDIR, d, "quant.sf"))]
all_gene_counts = {}
for d in sorted(quant_dirs):
    srr = d.split('_')[0]
    qf = os.path.join(OUTDIR, d, "quant.sf")
    df = pd.read_csv(qf, sep='\t', index_col=0)
    df.index = df.index.str.replace(r'\.\d+$', '', regex=True)
    df['gene'] = df.index.map(tx2gene_base)
    df = df[df['gene'].notna()]
    gc = df.groupby('gene')['NumReads'].sum()
    all_gene_counts[srr] = gc
count_matrix = pd.DataFrame(all_gene_counts).fillna(0.0).round().astype('int64')
count_matrix.index.name = 'gene'
print(f"Count matrix: {count_matrix.shape[0]} genes x {count_matrix.shape[1]} samples")

# ---- 3. metadata ----
meta_raw = pd.read_csv(METAF, sep='\t')
meta = pd.DataFrame(index=meta_raw['run_accession'])
meta['condition'] = meta_raw['sample_title'].str.extract(r'^(IVF|PA)').iloc[:, 0].values
meta['stage'] = meta_raw['sample_title'].str.extract(r'(\d cell)').iloc[:, 0].str.replace(' ', '', regex=False).values
meta['rep'] = meta_raw['sample_title'].str.extract(r'rep(\d+)').astype(int)
count_matrix = count_matrix[meta.index]
print(f"Conditions: {meta['condition'].value_counts().to_dict()}")
print(f"Stages: {meta['stage'].value_counts().to_dict()}")
for (c, s), grp in meta.groupby(['condition', 'stage']):
    print(f"  {c:4s} {s:6s}: n={len(grp)}")

# ---- 4. low-count filter ----
keep = ((count_matrix >= 10).sum(axis=1) >= 3)
count_f = count_matrix[keep]
print(f"Genes after filter: {count_f.shape[0]}")

# ---- 5. Full interaction model ----
print("\nFitting ~ stage + condition + stage:condition ...")
dds = DeseqDataSet(counts=count_f.T.astype(float), metadata=meta[['condition', 'stage']],
                   design="~ stage + condition + stage:condition", quiet=True)
dds.deseq2()
dm = dds.obsm['design_matrix']
print("Design matrix columns:", list(dm.columns))
lfc_cols = list(dm.columns)

# ---- 6. Wald tests ----
# 6a. main effect of condition at reference stage (condition[T.PA])
# 6b. each interaction coefficient
results = {}
# 主效应: condition[T.PA] coefficient
main_col = [c for c in lfc_cols if c.startswith('condition[T.')]
if main_col:
    idx = lfc_cols.index(main_col[0])
    vec = np.zeros(len(lfc_cols)); vec[idx] = 1
    stat = DeseqStats(dds, contrast=vec, quiet=True)
    stat.summary()
    r = stat.results_df.copy()
    r['gene'] = r.index
    r['gene_symbol'] = r['gene'].map(gene_symbol)
    r['coefficient'] = main_col[0]
    results['main_condition'] = r
    n_sig = (r['padj'] < 0.05).sum()
    n_up = ((r['padj'] < 0.05) & (r['log2FoldChange'] > 0)).sum()
    n_down = ((r['padj'] < 0.05) & (r['log2FoldChange'] < 0)).sum()
    print(f"  Main condition effect ({main_col[0]}): {n_sig} genes padj<0.05 (up={n_up}, down={n_down})")

interact_cols = [c for c in lfc_cols if ':' in c]
print(f"  Interaction coefficients: {interact_cols}")
for ic in interact_cols:
    idx = lfc_cols.index(ic)
    vec = np.zeros(len(lfc_cols)); vec[idx] = 1
    stat = DeseqStats(dds, contrast=vec, quiet=True)
    stat.summary()
    r = stat.results_df.copy()
    r['gene'] = r.index
    r['gene_symbol'] = r['gene'].map(gene_symbol)
    r['coefficient'] = ic
    results[ic] = r
    n_sig = (r['padj'] < 0.05).sum()
    print(f"  {ic}: {n_sig} genes padj<0.05")

# ---- 7. Save ----
allres = pd.concat(results.values())
allres.to_csv(os.path.join(RES, "m15_deseq2_stageXcondition_fullModel.csv"))
print(f"\nSaved full results: {len(allres)} rows")

# summary table
summary = []
for name, r in results.items():
    n_sig = (r['padj'] < 0.05).sum()
    n_sig_abs1 = ((r['padj'] < 0.05) & (r['log2FoldChange'].abs() > 1)).sum()
    summary.append({'coefficient': name, 'n_sig_padj005': int(n_sig), 'n_sig_absLFC1': int(n_sig_abs1)})
pd.DataFrame(summary).to_csv(os.path.join(RES, "m15_deseq2_interaction_summary.csv"), index=False)
print("\n=== INTERACTION SUMMARY ===")
print(pd.DataFrame(summary).to_string(index=False))
print("\nDONE")
