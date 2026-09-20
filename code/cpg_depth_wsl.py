#!/usr/bin/env python3
"""CpG depth sensitivity — processes pre-filtered CG data from stdin.
Usage: zcat file.CGmap.gz | awk '$4=="CG"' | python3 cpg_depth_wsl.py
Input: tab-separated CG lines from stdin
"""
import sys, bisect
from collections import defaultdict

# Parse GTF (fast, once)
chrom_genes = defaultdict(list)
with open('/mnt/e/Workbuddy/2026-07-27-11-58-27/data/raw/pig_genes.gtf') as fh:
# The GTF needs to be decompressed first
import gzip
with gzip.open('/mnt/e/Workbuddy/2026-07-27-11-58-27/data/raw/pig_genes.gtf.gz', 'rt') as fh:
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

# Process stdin
gene_cpg = defaultdict(lambda: [0,0,0])
n = 0
for line in sys.stdin:
    n += 1
    if n % 1000000 == 0:
        sys.stderr.write(f'\r  {n//1000000}M CG lines processed')
        sys.stderr.flush()
    p = line.strip().split('\t')
    if len(p) < 8: continue
    chrom, pos = p[0], int(p[2])
    depth = int(p[7])
    if depth < 1: continue
    
    glist = chrom_genes.get(chrom, [])
    if not glist: continue
    tss_list = [g[0] for g in glist]
    left = bisect.bisect_left(tss_list, pos - 2000)
    right = bisect.bisect_right(tss_list, pos + 2000)
    for i in range(left, right):
        gid = glist[i][1]
        gene_cpg[gid][0] += 1
        if depth >= 5: gene_cpg[gid][1] += 1
        if depth >= 10: gene_cpg[gid][2] += 1

# Output
sys.stderr.write(f'\n  {len(gene_cpg)} genes with CpG data\n')
for gid, counts in sorted(gene_cpg.items()):
    print(f'{gid}\t{counts[0]}\t{counts[1]}\t{counts[2]}')
