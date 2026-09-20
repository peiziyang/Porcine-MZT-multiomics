#!/usr/bin/env python
"""
CGmap → gene-level methylation matrix (v6).
Bypasses shell quoting issues by writing awk scripts to files.
"""
import os, sys, gzip, glob, time, subprocess
import pandas as pd
import numpy as np

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
CGMAP_DIR = os.path.join(BASE, "data/raw/cgmap")
OUT_DIR = os.path.join(BASE, "data/processed/m4_mofa")
GTF_FILE = "E:/Workbuddy/2026-07-27-11-58-27/data/raw/pig_genes.gtf.gz"
os.makedirs(OUT_DIR, exist_ok=True)

def wsl(cmd):
    """Run command in WSL, return output."""
    result = subprocess.run(
        ["wsl", "--distribution", "Ubuntu", "--exec", "bash", "-c", cmd],
        capture_output=True, text=True, timeout=3600
    )
    if result.stderr and "error" in result.stderr.lower():
        print(f"  WSL stderr: {result.stderr[:200]}", flush=True)
    return result.stdout

# Step 1: GTF → BED (Python)
print("STEP 1: GTF → BED...", flush=True)
gene_bed = os.path.join(OUT_DIR, "pig_genes.bed").replace("\\", "/")
with gzip.open(GTF_FILE, "rt") as gtf, open(gene_bed, "w") as out:
    for line in gtf:
        if line.startswith("#"): continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "gene": continue
        gid = None
        for a in parts[8].split(";"):
            a = a.strip()
            if a.startswith("gene_id "):
                gid = a.split('"')[1] if '"' in a else a.split()[1]
                break
        if gid:
            out.write(f"{parts[0]}\t{parts[3]}\t{parts[4]}\t{gid}\n")
n_genes = sum(1 for _ in open(gene_bed))
print(f"  {n_genes} genes", flush=True)

# Step 2: Process CGmap files
cgmap_files = sorted(glob.glob(os.path.join(CGMAP_DIR, "*.CGmap.gz")))
print(f"\nSTEP 2: {len(cgmap_files)} files...", flush=True)

gene_methyl = {}

for fpath in cgmap_files:
    fname = os.path.basename(fpath)
    cell_id = "_".join(fname.split("_")[1:]).replace(".CGmap.gz", "")
    fstart = time.time()
    
    fp_wsl = "/mnt/e" + fpath[2:].replace("\\", "/")
    gene_bed_wsl = "/mnt/e" + gene_bed[2:].replace("\\", "/")
    cpg_bed_wsl = f"/mnt/e{OUT_DIR[2:]}/_tmp_{cell_id}_cpg.bed"
    res_bed_wsl = f"/mnt/e{OUT_DIR[2:]}/_tmp_{cell_id}_gene_methyl.txt"
    cpg_bed_win = f"{OUT_DIR}/_tmp_{cell_id}_cpg.bed"
    res_bed_win = f"{OUT_DIR}/_tmp_{cell_id}_gene_methyl.txt"
    
    # Write a temp awk script to avoid quoting hell
    awk_script = f"""{cpg_bed_win}.awk"""
    with open(awk_script, "w") as af:
        af.write('$4=="CG"{print $1"\t"$3"\t"$3+1"\t"$6}\n')
    
    awk_wsl = f"/mnt/e{OUT_DIR[2:]}/_tmp_{cell_id}_cpg.bed.awk"
    
    # Extract CpG sites
    wsl(f"zcat {fp_wsl} | awk -f {awk_wsl} > {cpg_bed_wsl}")
    os.remove(awk_script)
    
    if not os.path.exists(cpg_bed_win) or os.path.getsize(cpg_bed_win) == 0:
        print(f"  {cell_id}: no CpG found ({time.time()-fstart:.0f}s)", flush=True)
        continue
    
    # bedtools map
    wsl(f"bedtools map -a {gene_bed_wsl} -b {cpg_bed_wsl} -c 4 -o mean -null 0 > {res_bed_wsl}")
    
    # Read results
    n_genes = 0
    if os.path.exists(res_bed_win):
        with open(res_bed_win) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 5 and parts[4] not in (".", "0", "") and float(parts[4]) > 0:
                    gid = parts[3]
                    methyl = float(parts[4])
                    if gid not in gene_methyl:
                        gene_methyl[gid] = {}
                    gene_methyl[gid][cell_id] = methyl
                    n_genes += 1
    
    # Cleanup
    for f in [cpg_bed_win, res_bed_win]:
        if os.path.exists(f):
            os.remove(f)
    
    print(f"  {cell_id}: {n_genes} genes ({time.time()-fstart:.0f}s)", flush=True)

# Step 3: Methylation matrix
print(f"\nSTEP 3: Methylation matrix...", flush=True)
meth_data = {}
for gid, cells in gene_methyl.items():
    for cid, methyl in cells.items():
        if cid not in meth_data:
            meth_data[cid] = {}
        meth_data[cid][gid] = methyl
meth_df = pd.DataFrame(meth_data).T
meth_df = meth_df.fillna(meth_df.median())
print(f"  {meth_df.shape}", flush=True)
meth_df.to_csv(os.path.join(OUT_DIR, "methylation_matrix.csv"))

# Step 4: MOFA+
print(f"\nSTEP 4: MOFA+...", flush=True)
rna = pd.read_csv(os.path.join(BASE, "data/processed/m7_hvg_symbol.csv"), index_col=0).T
common = sorted(set(meth_df.index) & set(rna.index))
print(f"  Common: {len(common)} cells", flush=True)
if len(common) >= 5:
    cg = sorted(set(meth_df.columns) & set(rna.columns))
    top = (rna[cg].var().sort_values(ascending=False).head(3000).index
           if len(cg) > 3000 else cg)
    ma = pd.DataFrame(0.5, index=rna.index, columns=top)
    for g in meth_df.columns:
        if g in ma.columns:
            ma[g] = meth_df[g]
    pd.concat([rna.loc[common, top], ma.loc[common, top]], axis=1).to_csv(
        os.path.join(OUT_DIR, "multiomics_views.csv"))
    try:
        from mofapy2.run.entry_point import entry_point
        ent = entry_point()
        ent.set_data_options(scale_views=True)
        vd = pd.concat([rna.loc[common, top].add_prefix("RNA_"),
                        ma.loc[common, top].add_prefix("METH_")], axis=1)
        ent.set_data_df(vd)
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
