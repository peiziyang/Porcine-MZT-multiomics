#!/usr/bin/env python
"""
CGmap → gene-level methylation matrix (v5).
Python: parses GTF → BED (works).
WSL awk+bedtools: extracts CpG → maps to genes (fast).
"""
import os, sys, gzip, glob, time
import pandas as pd
import numpy as np

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
CGMAP_DIR = os.path.join(BASE, "data/raw/cgmap")
OUT_DIR = os.path.join(BASE, "data/processed/m4_mofa")
GTF_FILE = "E:/Workbuddy/2026-07-27-11-58-27/data/raw/pig_genes.gtf.gz"
os.makedirs(OUT_DIR, exist_ok=True)

def wsl_path(p):
    p = p.replace("\\", "/")
    if p.startswith("E:"):
        return "/mnt/e" + p[2:]
    return p

gene_bed_py = os.path.join(OUT_DIR, "pig_genes.bed").replace("\\", "/")
gene_bed_wsl = wsl_path(gene_bed_py)

# Step 1: GTF → BED (Python)
print("STEP 1: GTF → BED...", flush=True)
with gzip.open(GTF_FILE, "rt") as gtf, open(gene_bed_py, "w") as out:
    for line in gtf:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        gid = None
        for a in parts[8].split(";"):
            a = a.strip()
            if a.startswith("gene_id "):
                gid = a.split('"')[1] if '"' in a else a.split()[1]
                break
        if gid:
            out.write(f"{parts[0]}\t{parts[3]}\t{parts[4]}\t{gid}\n")

n_genes = sum(1 for _ in open(gene_bed_py))
print(f"  {n_genes} genes in BED", flush=True)

# Step 2: Process CGmap files via WSL
cgmap_files = sorted(glob.glob(os.path.join(CGMAP_DIR, "*.CGmap.gz")))
print(f"\nSTEP 2: {len(cgmap_files)} files via WSL awk+bedtools...", flush=True)

gene_methyl = {}

for fpath in cgmap_files:
    fname = os.path.basename(fpath)
    cell_id = "_".join(fname.split("_")[1:]).replace(".CGmap.gz", "")
    fstart = time.time()
    
    # Extract CpG sites: "zcat file | awk '$4==\"CG\" {print $1,$3,$3+1,$6}' > cpg.bed"
    cpg_bed = os.path.join(OUT_DIR, f"_tmp_{cell_id}_cpg.bed").replace("\\", "/")
    cpg_wsl = wsl_path(cpg_bed)
    fp_wsl = wsl_path(fpath)
    
    cmd1 = f'wsl --distribution Ubuntu --exec bash -c "zcat {fp_wsl} | awk \'\\$4==\\"CG\\"{{print \\$1\\"\\t\\"\\$3\\"\\t\\"\\$3+1\\"\\t\\"\\$6}}\' > {cpg_wsl}"'
    os.system(cmd1)
    
    if not os.path.exists(cpg_bed) or os.path.getsize(cpg_bed) == 0:
        print(f"  {cell_id}: no CpG sites found ({time.time()-fstart:.0f}s)", flush=True)
        continue
    
    # bedtools map: for each gene (BED intervals), compute mean methylation of CpGs
    result_bed = os.path.join(OUT_DIR, f"_tmp_{cell_id}_gene_methyl.txt").replace("\\", "/")
    res_wsl = wsl_path(result_bed)
    
    cmd2 = f'wsl --distribution Ubuntu --exec bash -c "bedtools map -a {gene_bed_wsl} -b {cpg_wsl} -c 4 -o mean -null 0 > {res_wsl}"'
    os.system(cmd2)
    
    # Read results
    n_genes_with_methyl = 0
    if os.path.exists(result_bed):
        with open(result_bed) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 5 and parts[4] != "." and parts[4] != "0" and float(parts[4]) > 0:
                    gid = parts[3]
                    methyl = float(parts[4])
                    if gid not in gene_methyl:
                        gene_methyl[gid] = {}
                    gene_methyl[gid][cell_id] = methyl
                    n_genes_with_methyl += 1
    
    # Cleanup
    for f in [cpg_bed, result_bed]:
        if os.path.exists(f):
            os.remove(f)
    
    elapsed = time.time() - fstart
    print(f"  {cell_id}: {n_genes_with_methyl} genes with CpG ({elapsed:.0f}s)", flush=True)

# Step 3: Methylation matrix
print(f"\nSTEP 3: Methylation matrix...", flush=True)
meth_data = {}
for gid, cells in gene_methyl.items():
    for cid, vals in cells.items():
        if cid not in meth_data:
            meth_data[cid] = {}
        meth_data[cid][gid] = vals

meth_df = pd.DataFrame(meth_data).T.fillna(meth_df.median())
print(f"  {meth_df.shape}", flush=True)
meth_df.to_csv(os.path.join(OUT_DIR, "methylation_matrix.csv"))
print(f"  Saved", flush=True)

# Step 4: MOFA+
print(f"\nSTEP 4: MOFA+...", flush=True)
rna = pd.read_csv(os.path.join(BASE, "data/processed/m7_hvg_symbol.csv"), index_col=0).T
common = sorted(set(meth_df.index) & set(rna.index))
print(f"  Common: {len(common)} cells", flush=True)

if len(common) >= 5:
    common_genes = sorted(set(meth_df.columns) & set(rna.columns))
    top = (rna[common_genes].var().sort_values(ascending=False).head(3000).index
           if len(common_genes) > 3000 else common_genes)
    meth_aligned = pd.DataFrame(0.5, index=rna.index, columns=top)
    for g in meth_df.columns:
        if g in meth_aligned.columns:
            meth_aligned[g] = meth_df[g]
    
    pd.concat([rna.loc[common, top], meth_aligned.loc[common, top]], axis=1).to_csv(
        os.path.join(OUT_DIR, "multiomics_views.csv"))
    
    try:
        from mofapy2.run.entry_point import entry_point
        ent = entry_point()
        ent.set_data_options(scale_views=True)
        view_df = pd.concat([
            rna.loc[common, top].add_prefix("RNA_"),
            meth_aligned.loc[common, top].add_prefix("METH_")
        ], axis=1)
        ent.set_data_df(view_df)
        ent.set_model_options(factors=5, likelihoods="gaussian")
        ent.set_train_options(iter=500, convergence_mode="fast", verbose=True, seed=42)
        ent.build()
        ent.run()
        pd.DataFrame(ent.get_factors(), index=common).to_csv(
            os.path.join(OUT_DIR, "mofa_multiomics_factors.csv"))
        print(f"  MOFA+ done!", flush=True)
    except Exception as e:
        print(f"  MOFA+ failed: {type(e).__name__}", flush=True)

print("\nDONE!", flush=True)
