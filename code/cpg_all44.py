#!/usr/bin/env python3
"""Process ALL 44 CGmap files for CpG depth sensitivity. WSL only, background."""
import gzip, bisect, os, sys, glob
from collections import defaultdict

ROOT = '/mnt/e/Workbuddy/2026-07-27-11-58-27'

# 1. GTF
print('Parsing GTF...', flush=True)
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

# 2. Process ALL 44 CGmap files
fmaps=sorted(glob.glob(ROOT+'/data/raw/*/*.CGmap.gz'))
print(f'Processing {len(fmaps)} files...', flush=True)

# Global accumulation: gene -> [total_1x, total_5x, total_10x, n_samples]
gene_all = defaultdict(lambda:[0,0,0,0])

for fi,fn in enumerate(fmaps):
    sname=os.path.basename(fn).replace('.CGmap.gz','')
    print(f'[{fi+1:2d}/{len(fmaps)}] {sname}...', end=' ', flush=True)
    
    gene={}
    with gzip.open(fn,'rt') as f:
        for ln in f:
            p=ln.strip().split('\t')
            if len(p)<8 or p[3]!='CG': continue
            c,po=p[0],int(p[2]); d=int(p[7])
            if d<1: continue
            gl=chrom_genes.get(c,[])
            if not gl: continue
            tl=[g[0] for g in gl]
            L=bisect.bisect_left(tl,po-2000)
            R=bisect.bisect_right(tl,po+2000)
            for j in range(L,R):
                gid=gl[j][1]
                if gid not in gene:
                    gene[gid]=[0,0,0]
                gene[gid][0]+=1
                if d>=5: gene[gid][1]+=1
                if d>=10: gene[gid][2]+=1
    
    # Accumulate
    for gid,counts in gene.items():
        gene_all[gid][0]+=counts[0]
        gene_all[gid][1]+=counts[1]
        gene_all[gid][2]+=counts[2]
        gene_all[gid][3]+=1
    print(f'{len(gene)} genes', flush=True)

# 3. Write output
out=ROOT+'/data/processed/m4_mofa/cpg_depth_all44.tsv'
with open(out,'w') as fh:
    fh.write('gene\tn_cpg_1x\tn_cpg_5x\tn_cpg_10x\tn_samples\n')
    for gid,counts in sorted(gene_all.items()):
        fh.write(f'{gid}\t{counts[0]}\t{counts[1]}\t{counts[2]}\t{counts[3]}\n')
print(f'\nSaved: {out} ({len(gene_all)} genes)', flush=True)
