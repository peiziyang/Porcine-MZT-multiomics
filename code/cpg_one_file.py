#!/usr/bin/env python3
"""Process ONE CGmap file for CpG depth sensitivity. WSL only."""
import gzip, bisect, os, sys, glob
from collections import defaultdict

ROOT = '/mnt/e/Workbuddy/2026-07-27-11-58-27'
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

# 1. GTF
print('GTF...', flush=True)
chrom_genes = defaultdict(list)
with gzip.open(ROOT+'/data/raw/pig_genes.gtf.gz','rt') as f:
    for l in f:
        if l[0]=='#': continue
        p=l.strip().split('\t')
        if len(p)<9 or p[2]!='gene': continue
        c,s,e,st=p[0],int(p[3]),int(p[4]),p[6]
        a={kv.split(' ')[0]:kv.split(' ',1)[1].strip('"') for kv in p[8].rstrip(';').split('; ') if ' ' in kv}
        g=a.get('gene_id','')
        if g:
            t=s if st=='+' else e
            chrom_genes[c].append((t,g))
for c in chrom_genes: chrom_genes[c].sort()
print(f'{sum(len(v) for v in chrom_genes.values())} genes', flush=True)

# 2. Process ONE CGmap
fmaps=sorted(glob.glob(ROOT+'/data/raw/*/*.CGmap.gz'))
fn=fmaps[0]
print(f'File: {os.path.basename(fn)}', flush=True)

gene=defaultdict(lambda:[0,0,0])
ncg=0
with gzip.open(fn,'rt') as f:
    for i,ln in enumerate(f):
        if i%5000000==0:
            print(f'  {i//1000000}M lines, {ncg//1000}K CG, {len(gene)} genes', flush=True)
        p=ln.strip().split('\t')
        if len(p)<8 or p[3]!='CG': continue
        ncg+=1
        c,po=p[0],int(p[2]); d=int(p[7])
        if d<1: continue
        gl=chrom_genes.get(c,[])
        if not gl: continue
        tl=[g[0] for g in gl]
        L=bisect.bisect_left(tl,po-2000)
        R=bisect.bisect_right(tl,po+2000)
        for j in range(L,R):
            gid=gl[j][1]
            gene[gid][0]+=1
            if d>=5: gene[gid][1]+=1
            if d>=10: gene[gid][2]+=1

print(f'\nDone: {i} lines, {ncg} CG, {len(gene)} genes with promoter CpGs', flush=True)

# 3. Write output
out=ROOT+'/data/processed/m4_mofa/cpg_depth_raw.tsv'
with open(out,'w') as fh:
    fh.write('gene\tn_cpg_1x\tn_cpg_5x\tn_cpg_10x\n')
    for gid,counts in sorted(gene.items()):
        fh.write(f'{gid}\t{counts[0]}\t{counts[1]}\t{counts[2]}\n')
print(f'Saved: {out}', flush=True)
