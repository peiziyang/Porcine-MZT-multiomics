#!/usr/bin/env python
"""CpG sensitivity v2: fast indexed TSS matching with bisect."""
import gzip, os, bisect
import numpy as np
import pandas as pd
from collections import defaultdict

ROOT = r'E:/Workbuddy/2026-07-27-11-58-27'
GTF  = os.path.join(ROOT, 'data/raw/pig_genes.gtf.gz')

# ----- 1: Parse GTF → chromosome-indexed gene TSS -----
print('Parsing GTF...', flush=True)
chrom_genes = defaultdict(list)  # chrom -> [(tss, gene_id)]
with gzip.open(GTF, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
        attrs = {}
        for kv in parts[8].rstrip(';').split('; '):
            if ' ' in kv:
                k, v = kv.split(' ', 1)
                attrs[k] = v.strip('"')
        gid = attrs.get('gene_id', '')
        if gid:
            tss = start if strand == '+' else end
            chrom_genes[chrom].append((tss, gid))

# Sort each chromosome's genes by TSS
for chrom in chrom_genes:
    chrom_genes[chrom].sort()

print(f'  {sum(len(v) for v in chrom_genes.values())} genes across {len(chrom_genes)} chromosomes', flush=True)

# ----- 2: Parse CGmap files with fast lookup -----
# To avoid re-parsing CGmap 3 times, compute all thresholds in one pass
import glob
cgmaps = sorted(glob.glob(os.path.join(ROOT, 'data/raw/*/*.CGmap.gz')))
print(f'  Processing {len(cgmaps)} CGmap files...', flush=True)

# Per-gene counts (summed across samples): [total_1x, total_5x, total_10x, samples_with_data]
gene_stats = defaultdict(lambda: [0, 0, 0, 0])

for fi, fn in enumerate(cgmaps):
    sname = os.path.basename(fn).replace('.CGmap.gz', '')
    # Per-sample counts for this file
    samp_counts = defaultdict(lambda: [0, 0, 0])

    with gzip.open(fn, 'rt') as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            context = parts[3]
            if context != 'CG':
                continue
            chrom, pos = parts[0], int(parts[2])
            total_reads = int(parts[7])
            if total_reads < 1:
                continue

            # Find genes with TSS within 2kb using bisect
            gene_list = chrom_genes.get(chrom, [])
            if not gene_list:
                continue
            tss_list = [g[0] for g in gene_list]
            # Find range: TSS in [pos-2000, pos+2000]
            left = bisect.bisect_left(tss_list, pos - 2000)
            right = bisect.bisect_right(tss_list, pos + 2000)
            for i in range(left, right):
                gid = gene_list[i][1]
                samp_counts[gid][0] += 1  # ≥1x
                if total_reads >= 5:
                    samp_counts[gid][1] += 1  # ≥5x
                if total_reads >= 10:
                    samp_counts[gid][2] += 1  # ≥10x

    # Accumulate to global stats
    for gid, counts in samp_counts.items():
        gene_stats[gid][0] += counts[0]
        gene_stats[gid][1] += counts[1]
        gene_stats[gid][2] += counts[2]
        gene_stats[gid][3] += 1  # samples with data

    if (fi + 1) % 10 == 0 or fi == len(cgmaps) - 1:
        print(f'  [{fi+1}/{len(cgmaps)}] done, {len(gene_stats)} genes with CpG data', flush=True)

# ----- 3: Build dataframe -----
rows = []
for gid, counts in sorted(gene_stats.items()):
    rows.append({
        'gene': gid,
        'n_cpg_1x': counts[0],
        'n_cpg_5x': counts[1],
        'n_cpg_10x': counts[2],
        'n_samples_with_data': counts[3],
    })

df = pd.DataFrame(rows)
df['median_cpg_1x_per_sample'] = df['n_cpg_1x'] / df['n_samples_with_data']
df['median_cpg_5x_per_sample'] = df['n_cpg_5x'] / df['n_samples_with_data']
df['median_cpg_10x_per_sample'] = df['n_cpg_10x'] / df['n_samples_with_data']
df['frac_5x'] = np.where(df['n_cpg_1x'] > 0, df['n_cpg_5x'] / df['n_cpg_1x'], 0)
df['frac_10x'] = np.where(df['n_cpg_1x'] > 0, df['n_cpg_10x'] / df['n_cpg_1x'], 0)

# ----- 4: F1 genes vs background -----
w = pd.read_csv(os.path.join(ROOT, 'data/processed/m4_mofa/mofa_multiomics_weights_RNA_v2.csv'), index_col=0)
f1_genes = list(w['F1'].abs().sort_values(ascending=False).head(44).index)
gs = pd.read_csv(os.path.join(ROOT, 'data/processed/m4_mofa/per_gene_cpg_stats.csv'))
symbol2ens = dict(zip(gs['gene_name'], gs['gene_id']))

f1_ens = [symbol2ens[g] for g in f1_genes if g in symbol2ens]
f1_matched = [e for e in f1_ens if e in set(df['gene'])]

f1_df = df[df['gene'].isin(f1_matched)]
bg_df  = df[~df['gene'].isin(f1_matched)]

print(f'\nF1 genes: {len(f1_genes)} symbol, {len(f1_ens)} Ensembl, {len(f1_matched)} in CGmap', flush=True)

print(f'\n===== CpG COVERAGE SENSITIVITY =====')
for label, col in [('≥1x total CpGs', 'n_cpg_1x'), ('≥5x total CpGs', 'n_cpg_5x'), ('≥10x total CpGs', 'n_cpg_10x'),
                    ('≥1x per-sample', 'median_cpg_1x_per_sample'), ('≥5x per-sample', 'median_cpg_5x_per_sample'),
                    ('≥10x per-sample', 'median_cpg_10x_per_sample')]:
    print(f'  {label:<22s}: F1 mean={f1_df[col].mean():.1f} median={f1_df[col].median():.0f}  BG mean={bg_df[col].mean():.1f} median={bg_df[col].median():.0f}')

print(f'\n  Fraction retained ≥5x/≥1x: F1={f1_df["frac_5x"].mean():.2f}, BG={bg_df["frac_5x"].mean():.2f}')
print(f'  Fraction retained ≥10x/≥1x: F1={f1_df["frac_10x"].mean():.2f}, BG={bg_df["frac_10x"].mean():.2f}')

# Core 5
print(f'\n===== CORE FIVE F1 GENES =====')
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    eid = symbol2ens.get(g, '')
    if eid and eid in set(df['gene']):
        r = df[df['gene'] == eid].iloc[0]
        print(f'  {g:>10s}: ≥1x={r["n_cpg_1x"]:.0f}  ≥5x={r["n_cpg_5x"]:.0f}  ≥10x={r["n_cpg_10x"]:.0f}  frac_5x={r["frac_5x"]:.2f}')
    else:
        print(f'  {g:>10s}: not found')

# F1 genes with weak coverage (≤5 CpGs at 5x)
weak = f1_df[f1_df['n_cpg_5x'] <= 5]
print(f'\n  F1 genes with ≤5 CpGs at ≥5x: {len(weak)}/{len(f1_df)}')
if len(weak) > 0 and len(weak) <= 10:
    for _, r in weak.iterrows():
        eid = r['gene']
        sym = [g for g, e in symbol2ens.items() if e == eid]
        print(f'    {sym[0] if sym else eid}: ≥1x={r["n_cpg_1x"]:.0f}, ≥5x={r["n_cpg_5x"]:.0f}')

df.to_csv(os.path.join(ROOT, 'data/processed/m4_mofa/cpg_depth_sensitivity.csv'), index=False)
print(f'\nSaved: cpg_depth_sensitivity.csv ({len(df)} genes)')
