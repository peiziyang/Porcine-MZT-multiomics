#!/usr/bin/env python
"""CpG sensitivity: ≥1×, ≥5×, ≥10× promoter methylation for F1 genes vs background."""
import gzip, os, sys
import numpy as np
import pandas as pd

ROOT = r'E:/Workbuddy/2026-07-27-11-58-27'
GTF  = os.path.join(ROOT, 'data/raw/pig_genes.gtf.gz')
CGMPAT = os.path.join(ROOT, 'data/raw/*/*.CGmap.gz')  # glob

# ----- 1: parse GTF → gene -> (chr, tss) -----
print('Parsing GTF...', flush=True)
gene_tss = {}  # gene_id -> (chr, tss)
with gzip.open(GTF, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] != 'gene':
            continue
        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        # Extract gene_id from attributes
        attrs = {}
        for kv in parts[8].rstrip(';').split('; '):
            if ' ' in kv:
                k, v = kv.split(' ', 1)
                v = v.strip('"')
                attrs[k] = v
        gid = attrs.get('gene_id', '')
        if gid:
            tss = start if strand == '+' else end
            gene_tss[gid] = (chrom, tss)

print(f'  {len(gene_tss)} genes parsed', flush=True)

# ----- 2: parse CGmap files per-threshold -----
import glob
cgmaps = sorted(glob.glob(os.path.join(ROOT, 'data/raw/*/*.CGmap.gz')))
print(f'  {len(cgmaps)} CGmap files', flush=True)

# Map: gene -> list of (per_sample_cpg_count_1x, per_sample_cpg_count_5x, per_sample_cpg_count_10x)
from collections import defaultdict
gene_cpg = defaultdict(lambda: [0, 0, 0])  # [1x, 5x, 10x]

BATCH_SIZE = 5000000
for fi, fn in enumerate(cgmaps):
    sname = os.path.basename(fn).replace('.CGmap.gz', '')
    print(f'  [{fi+1}/{len(cgmaps)}] {sname}', flush=True)
    with gzip.open(fn, 'rt') as fh:
        cpg_hits = defaultdict(lambda: [0, 0, 0])
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            # CGmap: chr, C, pos, context, dinuc, meth_lvl, C_count, CT_count
            context = parts[3]
            if context != 'CG':
                continue
            chrom_cg = parts[0]
            pos = int(parts[2])
            total_reads = int(parts[7])  # CT_count = total
            if total_reads < 1:
                continue
            
            # Find genes whose promoter covers this CpG
            for gid, (gchr, tss) in gene_tss.items():
                if gchr == chrom_cg and abs(pos - tss) <= 2000:
                    cpg_hits[gid][0] += 1  # ≥1x
                    if total_reads >= 5:
                        cpg_hits[gid][1] += 1  # ≥5x
                    if total_reads >= 10:
                        cpg_hits[gid][2] += 1  # ≥10x
    # Accumulate
    for gid, counts in cpg_hits.items():
        gene_cpg[gid][0] += counts[0]
        gene_cpg[gid][1] += counts[1]
        gene_cpg[gid][2] += counts[2]

print(f'  Genes with CpG data: {len(gene_cpg)}', flush=True)

# ----- 3: Build summary dataframe -----
rows = []
for gid, counts in sorted(gene_cpg.items()):
    rows.append({'gene': gid, 'n_cpg_1x': counts[0], 'n_cpg_5x': counts[1], 'n_cpg_10x': counts[2]})

df = pd.DataFrame(rows)
df['fraction_5x_vs_1x'] = np.where(df['n_cpg_1x'] > 0, df['n_cpg_5x'] / df['n_cpg_1x'], 0)
df['fraction_10x_vs_1x'] = np.where(df['n_cpg_1x'] > 0, df['n_cpg_10x'] / df['n_cpg_1x'], 0)

# ----- 4: Compare F1 genes vs background -----
w = pd.read_csv(os.path.join(ROOT, 'data/processed/m4_mofa/mofa_multiomics_weights_RNA_v2.csv'), index_col=0)
f1_genes = list(w['F1'].abs().sort_values(ascending=False).head(44).index)

# Try to match F1 gene symbols to Ensembl IDs
# Use gene_name -> symbol mapping from processed data
gs = pd.read_csv(os.path.join(ROOT, 'data/processed/m4_mofa/per_gene_cpg_stats.csv'))
symbol2ens = dict(zip(gs['gene_name'], gs['gene_id']))

f1_ens = []
for g in f1_genes:
    if g in symbol2ens:
        f1_ens.append(symbol2ens[g])
f1_matched = [e for e in f1_ens if e in set(df['gene'])]

print(f'\nF1 genes: {len(f1_genes)} symbol, {len(f1_ens)} mapped to Ensembl, {len(f1_matched)} in CGmap data', flush=True)

f1_df = df[df['gene'].isin(f1_matched)]
bg_df  = df[~df['gene'].isin(f1_matched)]

print(f'\n===== RESULTS =====')
print(f'{"Threshold":<12} {"F1 mean":>10} {"BG mean":>10} {"F1 median":>10} {"BG median":>10} {"F1>0":>8}')
for label, col in [('≥1x', 'n_cpg_1x'), ('≥5x', 'n_cpg_5x'), ('≥10x', 'n_cpg_10x')]:
    print(f'{label:<12} {f1_df[col].mean():>10.1f} {bg_df[col].mean():>10.1f} {f1_df[col].median():>10.0f} {bg_df[col].median():>10.0f} {sum(f1_df[col]>0):>8}/{len(f1_df)}')

print(f'\n===== FRACTION RETAINED =====')
print(f'≥5x/≥1x: F1={f1_df["fraction_5x_vs_1x"].mean():.2f}, BG={bg_df["fraction_5x_vs_1x"].mean():.2f}')
print(f'≥10x/≥1x: F1={f1_df["fraction_10x_vs_1x"].mean():.2f}, BG={bg_df["fraction_10x_vs_1x"].mean():.2f}')

# Core 5 genes
print(f'\n===== CORE FIVE GENES =====')
core5 = ['DNMT1','ZP3','ZP4','GDF9','RARRES1']
for g in core5:
    eid = symbol2ens.get(g, '')
    if eid and eid in set(df['gene']):
        r = df[df['gene'] == eid].iloc[0]
        print(f'  {g}: 1x={r["n_cpg_1x"]}, 5x={r["n_cpg_5x"]}, 10x={r["n_cpg_10x"]}, frac_5x={r["fraction_5x_vs_1x"]:.2f}')
    else:
        print(f'  {g}: not found in CGmap (symbol2ens: {eid})')

df.to_csv(os.path.join(ROOT, 'data/processed/m4_mofa/cpg_depth_sensitivity.csv'), index=False)
print(f'\nSaved: data/processed/m4_mofa/cpg_depth_sensitivity.csv')
