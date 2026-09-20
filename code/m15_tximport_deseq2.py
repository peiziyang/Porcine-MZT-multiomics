#!/usr/bin/env python
"""
M15: tximport-style summarization + DESeq2 (IVF vs PA across preimplantation stages)
====================================================================================
Replaces the flawed m10 approach (which collapsed 42 samples into 2 pseudobulk
columns -> no replication -> no real p-values).

Design (real metadata from ENA/GEO SRP301735, not assumed):
  - condition: IVF vs PA
  - stage: 1cell / 2cell / 4cell / 8cell
  - replicate: rep1..repN  (biological replicate unit, from GEO sample titles)

Primary: per-stage DESeq2 (~ condition) -> stage-specific IVF-vs-PA DEGs (Wald).
Secondary: additive ~ stage + condition -> overall IVF-vs-PA main effect (stage-adjusted).

LA scores: donor/replicate-level stats, report effect size + replicate n.
"""
import os, re, gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings
warnings.filterwarnings('ignore')

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
OUTDIR = os.path.join(BASE, "data", "processed", "salmon_out")
GTF = os.path.join(BASE, "data", "raw", "pig_genes.gtf.gz")
METAF = os.path.join(OUTDIR, "sra_metadata.tsv")
RES = os.path.join(BASE, "data", "processed", "deseq2")
os.makedirs(RES, exist_ok=True)

# ----------------------------------------------------------------------------
# 1. Build tx2gene (transcript -> gene_id, gene_symbol). Strip version numbers.
# ----------------------------------------------------------------------------
print("=" * 70)
print("M15: tximport + DESeq2  (IVF vs PA x stage, replicate-level)")
print("=" * 70)
tx2gene_base = {}
gene_symbol = {}
print("  Building tx2gene from GTF ...")
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
print(f"  tx2gene entries: {len(tx2gene_base)} ; genes with symbol: {sum(1 for v in gene_symbol.values() if v)}")

# ----------------------------------------------------------------------------
# 2. Read 42 quant.sf -> gene-level raw counts (sum NumReads per gene)
# ----------------------------------------------------------------------------
quant_dirs = [d for d in os.listdir(OUTDIR)
              if os.path.exists(os.path.join(OUTDIR, d, "quant.sf"))]
print(f"  Found {len(quant_dirs)} quantified samples")

all_gene_counts = {}
matched = 0
for d in sorted(quant_dirs):
    srr = d.split('_')[0]
    qf = os.path.join(OUTDIR, d, "quant.sf")
    df = pd.read_csv(qf, sep='\t', index_col=0)
    # strip transcript version (ENSSSCT...4 -> ENSSSCT...)
    df.index = df.index.str.replace(r'\.\d+$', '', regex=True)
    df['gene'] = df.index.map(tx2gene_base)
    df = df[df['gene'].notna()]
    matched += len(df)
    gc = df.groupby('gene')['NumReads'].sum()
    all_gene_counts[srr] = gc

count_matrix = pd.DataFrame(all_gene_counts).fillna(0.0)
# DESeq2 requires integer counts; salmon NumReads are fractional estimates -> round
count_matrix = count_matrix.round().astype('int64')
count_matrix.index.name = 'gene'
print(f"  Gene-level count matrix: {count_matrix.shape[0]} genes x {count_matrix.shape[1]} samples")
print(f"  Transcript->gene matches: {matched} (transcripts mapped to genes)")

# ----------------------------------------------------------------------------
# 3. Metadata (real labels from ENA sample_title)
# ----------------------------------------------------------------------------
meta_raw = pd.read_csv(METAF, sep='\t')
meta = pd.DataFrame()
meta['SRR'] = meta_raw['run_accession']
meta['condition'] = meta_raw['sample_title'].str.extract(r'^(IVF|PA)')
meta['stage'] = meta_raw['sample_title'].str.extract(r'(\d cell)')
meta['rep'] = meta_raw['sample_title'].str.extract(r'rep(\d+)').astype(int)
# normalize stage: "2 cell" -> "2cell" so it matches the stages[] list and DESeq2 formula
meta['stage'] = meta['stage'].str.replace(' ', '', regex=False)
meta = meta.set_index('SRR')
# reorder columns to match count_matrix
count_matrix = count_matrix[meta.index]
print(f"  Conditions: {meta['condition'].value_counts().to_dict()}")
print(f"  Stages:     {meta['stage'].value_counts().to_dict()}")
print(f"  Replicates per (cond,stage):")
for (c, s), grp in meta.groupby(['condition', 'stage']):
    print(f"    {c:4s} {s:6s}: n={len(grp)}  reps={sorted(grp['rep'].tolist())}")

# filter low-count genes (present in >=3 samples with >=10 counts)
keep = ((count_matrix >= 10).sum(axis=1) >= 3)
count_f = count_matrix[keep]
print(f"  Genes after low-count filter: {count_f.shape[0]} / {count_matrix.shape[0]}")

# ----------------------------------------------------------------------------
# 4. Primary: per-stage DESeq2 (IVF vs PA)
# ----------------------------------------------------------------------------
stages = ['1cell', '2cell', '4cell', '8cell']
stage_results = {}
for stage in stages:
    sub = meta[meta['stage'] == stage]
    n_ivf = (sub['condition'] == 'IVF').sum()
    n_pa = (sub['condition'] == 'PA').sum()
    sub_counts = count_f[sub.index].T.astype(float)   # samples x genes
    sub_meta = sub[['condition']].copy()
    print(f"\n  [{stage}] DESeq2 IVF vs PA : IVF n={n_ivf}, PA n={n_pa}")
    if n_ivf < 2 or n_pa < 2:
        print(f"    SKIP (need >=2 per group)")
        continue
    dds = DeseqDataSet(counts=sub_counts, metadata=sub_meta, design="~ condition")
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "PA", "IVF"])
    stat.summary()
    res = stat.results_df.copy()
    res['gene'] = res.index
    res['gene_symbol'] = res['gene'].map(gene_symbol)
    res['log2FC'] = res['log2FoldChange']
    res['absLFC'] = res['log2FC'].abs()
    res['sig'] = (res['padj'] < 0.05) & (res['log2FC'].abs() > 1)
    n_up = ((res['padj'] < 0.05) & (res['log2FC'] > 1)).sum()   # PA>IVF
    n_down = ((res['padj'] < 0.05) & (res['log2FC'] < -1)).sum()
    print(f"    DEGs (padj<0.05, |log2FC|>1): UP(PA>IVF)={n_up}, DOWN(PA<IVF)={n_down}")
    stage_results[stage] = res

    # save
    out_csv = os.path.join(RES, f"m15_deseq2_{stage}_IVFvsPA.csv")
    res.sort_values('padj').to_csv(out_csv)
    print(f"    saved {out_csv}")

    # volcano
    fig, ax = plt.subplots(figsize=(8, 6))
    cols = np.where(res['sig'] & (res['log2FC'] > 0), '#d62728',
            np.where(res['sig'] & (res['log2FC'] < 0), '#1f77b4', 'grey'))
    ax.scatter(res['log2FC'], -np.log10(res['padj'].clip(lower=1e-300)),
               c=cols, s=6, alpha=0.5, rasterized=True)
    ax.axhline(-np.log10(0.05), color='grey', ls='--', alpha=0.5)
    ax.axvline(1, color='grey', ls=':', alpha=0.3); ax.axvline(-1, color='grey', ls=':', alpha=0.3)
    ax.set_xlabel('log2 FC (PA - IVF)'); ax.set_ylabel('-log10(padj)')
    ax.set_title(f'Volcano: IVF vs PA @ {stage} (n_IVF={n_ivf}, n_PA={n_pa})')
    fig.savefig(os.path.join(RES, f"m15_volcano_{stage}.png"), dpi=130, bbox_inches='tight')
    plt.close(fig)
    # MA
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(np.log10(res['baseMean'] + 1), res['log2FC'], c=cols, s=4, alpha=0.4, rasterized=True)
    ax.axhline(0, color='grey', alpha=0.5); ax.axhline(1, color='grey', ls=':', alpha=0.3); ax.axhline(-1, color='grey', ls=':', alpha=0.3)
    ax.set_xlabel('log10(baseMean+1)'); ax.set_ylabel('log2 FC (PA - IVF)')
    ax.set_title(f'MA: IVF vs PA @ {stage}')
    fig.savefig(os.path.join(RES, f"m15_maplot_{stage}.png"), dpi=130, bbox_inches='tight')
    plt.close(fig)

# ----------------------------------------------------------------------------
# 5. Secondary: additive ~ stage + condition (overall IVF vs PA, stage-adjusted)
# ----------------------------------------------------------------------------
print("\n  [OVERALL] DESeq2 ~ stage + condition (IVF vs PA, stage-adjusted)")
sub_meta_all = meta[['condition', 'stage']].copy()
dds_all = DeseqDataSet(counts=count_f[meta.index].T.astype(float),
                       metadata=sub_meta_all, design="~ stage + condition")
dds_all.deseq2()
stat_all = DeseqStats(dds_all, contrast=["condition", "PA", "IVF"])
stat_all.summary()
res_all = stat_all.results_df.copy()
res_all['gene'] = res_all.index
res_all['gene_symbol'] = res_all['gene'].map(gene_symbol)
res_all['sig'] = (res_all['padj'] < 0.05) & (res_all['log2FoldChange'].abs() > 1)
res_all.sort_values('padj').to_csv(os.path.join(RES, "m15_deseq2_overall_IVFvsPA_stageAdj.csv"))
n_up = ((res_all['padj'] < 0.05) & (res_all['log2FoldChange'] > 1)).sum()
n_down = ((res_all['padj'] < 0.05) & (res_all['log2FoldChange'] < -1)).sum()
print(f"    Overall DEGs (padj<0.05, |log2FC|>1): UP(PA>IVF)={n_up}, DOWN(PA<IVF)={n_down}")

# ----------------------------------------------------------------------------
# 6. Summary + top genes per stage
# ----------------------------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
summary = {"design": "per-stage ~condition (Wald); overall ~stage+condition",
           "n_samples": int(len(meta)),
           "genes_tested_per_stage": {s: int(len(stage_results[s])) for s in stage_results}}
for s in stages:
    if s in stage_results:
        r = stage_results[s]
        n_up = int(((r['padj'] < 0.05) & (r['log2FC'] > 1)).sum())
        n_down = int(((r['padj'] < 0.05) & (r['log2FC'] < -1)).sum())
        print(f"  {s:6s}: UP(PA>IVF)={n_up:4d}  DOWN(PA<IVF)={n_down:4d}")
        top = r[r['sig']].sort_values('padj').head(10)
        for _, row in top.iterrows():
            sym = row['gene_symbol'] if isinstance(row['gene_symbol'], str) else row['gene']
            print(f"            {sym:12s} log2FC={row['log2FC']:+.2f} padj={row['padj']:.2e}")
print(f"\n  Outputs in: {RES}")
print("  - m15_deseq2_<stage>_IVFvsPA.csv  (per-stage DEG tables)")
print("  - m15_deseq2_overall_IVFvsPA_stageAdj.csv")
print("  - m15_volcano_<stage>.png, m15_maplot_<stage>.png")
print("DONE")
