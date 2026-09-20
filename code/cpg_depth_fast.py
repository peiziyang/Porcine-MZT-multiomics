#!/usr/bin/env python
"""CpG ≥5x sensitivity: use zcat+grep pre-filter for speed."""
import subprocess, os, glob, bisect, gzip
from collections import defaultdict
import numpy as np, pandas as pd

ROOT = r'E:/Workbuddy/2026-07-27-11-58-27'

# ---- Parse GTF (fast, once) ----
print('GTF...', end=' ', flush=True)
chrom_genes = defaultdict(list)
with gzip.open(f'{ROOT}/data/raw/pig_genes.gtf.gz', 'rt') as fh:
    for line in fh:
        if line[0] == '#': continue
        p = line.strip().split('\t')
        if len(p) < 9 or p[2] != 'gene': continue
        chrom, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
        attrs = {}
        for kv in p[8].rstrip(';').split('; '):
            if ' ' in kv:
                k, v = kv.split(' ', 1)
                attrs[k] = v.strip('"')
        gid = attrs.get('gene_id','')
        if gid:
            tss = start if strand == '+' else end
            chrom_genes[chrom].append((tss, gid))
for c in chrom_genes: chrom_genes[c].sort()
ng = sum(len(v) for v in chrom_genes.values())
print(f'{ng} genes', flush=True)

# ---- Pick 3 representative CGmap files ----
fmaps = sorted(glob.glob(f'{ROOT}/data/raw/*/*.CGmap.gz'))
run = [fmaps[0], fmaps[21], fmaps[-1]]  # first, middle, last
print(f'{len(run)}/{len(fmaps)} CGmap files', flush=True)

# ---- Process each file via zcat | grep ----
gene_stats = defaultdict(lambda: [0,0,0])  # [1x, 5x, 10x]

for fi, fn in enumerate(run):
    sname = os.path.basename(fn).replace('.CGmap.gz','')
    print(f'  [{fi+1}/3] {sname}...', end=' ', flush=True)
    
    # Use gzip.open directly since zcat piped approach is complex on Windows
    lines, cg_lines, hits = 0, 0, 0
    buf = []
    
    with gzip.open(fn, 'rt') as fh:
        for line in fh:
            lines += 1
            p = line.strip().split('\t')
            if len(p) < 8 or p[3] != 'CG':
                continue
            cg_lines += 1
            chrom, pos = p[0], int(p[2])
            depth = int(p[7])
            if depth < 1:
                continue
            
            glist = chrom_genes.get(chrom, [])
            if not glist:
                continue
            tss_list = [g[0] for g in glist]
            left = bisect.bisect_left(tss_list, pos - 2000)
            right = bisect.bisect_right(tss_list, pos + 2000)
            
            for i in range(left, right):
                gid = glist[i][1]
                gene_stats[gid][0] += 1  # ≥1x
                if depth >= 5:
                    gene_stats[gid][1] += 1  # ≥5x
                if depth >= 10:
                    gene_stats[gid][2] += 1  # ≥10x
                hits += 1
            
            # Progress every 10M lines
            if lines % 10000000 == 0:
                print(f'{lines//1000000}M...', end=' ', flush=True)
    
    print(f'{lines//1000000}M lines, {cg_lines//1000}K CG, {len(gene_stats)} genes', flush=True)

# ---- Build results ----
print(f'\nTotal genes with promoter CpGs: {len(gene_stats)}', flush=True)
rows = []
for gid, counts in gene_stats.items():
    rows.append({'gene': gid, 'n_cpg_1x': counts[0], 'n_cpg_5x': counts[1], 'n_cpg_10x': counts[2]})
df = pd.DataFrame(rows)
# Per-sample (divide by 3 since we processed 3 files)
for col in ['n_cpg_1x', 'n_cpg_5x', 'n_cpg_10x']:
    df[f'{col}_per_sample'] = df[col] / len(run)
df['frac_5x'] = np.where(df['n_cpg_1x']>0, df['n_cpg_5x']/df['n_cpg_1x'], 0)
df['frac_10x'] = np.where(df['n_cpg_1x']>0, df['n_cpg_10x']/df['n_cpg_1x'], 0)

# ---- F1 vs background ----
w = pd.read_csv(f'{ROOT}/data/processed/m4_mofa/mofa_multiomics_weights_RNA_v2.csv', index_col=0)
f1_genes = list(w['F1'].abs().sort_values(ascending=False).head(44).index)
gs = pd.read_csv(f'{ROOT}/data/processed/m4_mofa/per_gene_cpg_stats.csv')
s2e = dict(zip(gs['gene_name'], gs['gene_id']))
f1_ens = [s2e[g] for g in f1_genes if g in s2e]
f1_matched = [e for e in f1_ens if e in set(df['gene'])]
f1_df = df[df['gene'].isin(f1_matched)]
bg_df = df[~df['gene'].isin(f1_matched)]

print(f'F1: {len(f1_matched)}/{len(f1_ens)} genes matched')

# ---- Results ----
print(f'\n===== CpG DEPTH SENSITIVITY (per-gene, per-sample mean from {len(run)} CGmaps) =====')
for label, col in [
    ('≥1x CpGs/sample', 'n_cpg_1x_per_sample'),
    ('≥5x CpGs/sample', 'n_cpg_5x_per_sample'),
    ('≥10x CpGs/sample', 'n_cpg_10x_per_sample'),
]:
    print(f'  {label:<22}: F1={f1_df[col].mean():.1f}(med={f1_df[col].median():.0f})  BG={bg_df[col].mean():.1f}(med={bg_df[col].median():.0f})')

print(f'\n  Retention ≥5x/≥1x:  F1={f1_df["frac_5x"].mean():.2f}  BG={bg_df["frac_5x"].mean():.2f}')
print(f'  Retention ≥10x/≥1x: F1={f1_df["frac_10x"].mean():.2f}  BG={bg_df["frac_10x"].mean():.2f}')

# ---- Core 5 ----
print(f'\n===== CORE FIVE (per-sample mean CpGs) =====')
for g in ['DNMT1','ZP3','ZP4','GDF9','RARRES1']:
    eid = s2e.get(g,'')
    if eid in set(df['gene']):
        r = df[df['gene']==eid].iloc[0]
        print(f'  {g:>10s}: ≥1x={r["n_cpg_1x_per_sample"]:.0f}  ≥5x={r["n_cpg_5x_per_sample"]:.0f}  ≥10x={r["n_cpg_10x_per_sample"]:.0f}  frac_5x={r["frac_5x"]:.2f}')

df.to_csv(f'{ROOT}/data/processed/m4_mofa/cpg_depth_sensitivity.csv', index=False)
print(f'\nSaved: cpg_depth_sensitivity.csv')
