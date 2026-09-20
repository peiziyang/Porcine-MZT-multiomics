"""
CGmap processing using WSL awk + bedtools for 100x speedup.
1. Extract CpG sites per CGmap file
2. Map to genes and aggregate methylation
"""
import os, sys, glob, subprocess, time
import pandas as pd
import numpy as np

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
CGMAP_DIR = os.path.join(BASE, "data/raw/cgmap")
OUT_DIR = os.path.join(BASE, "data/processed/m4_mofa")
GTF_FILE = "E:/Workbuddy/2026-07-27-11-58-27/data/raw/pig_genes.gtf.gz"
os.makedirs(OUT_DIR, exist_ok=True)

def wsl_path(p):
    """Convert Windows path to WSL /mnt/e/ path."""
    p = p.replace("\\", "/")
    if p.startswith("E:"):
        return "/mnt/e" + p[2:]
    return p

WSL_PREFIX = "wsl --distribution Ubuntu --exec bash -c"

# Step 1: Convert GTF to BED (gene intervals)
print("STEP 1: GTF → BED...", flush=True)
gene_bed = os.path.join(OUT_DIR, "pig_genes.bed").replace("\\", "/")
cmd = f'{WSL_PREFIX} "zcat {wsl_path(GTF_FILE)} | awk \'$3==\\"gene\\" {{split($NF,a,\\\\\\\"\\\\\\\"); for(i in a) if(a[i]~\"gene_id\") gid=a[i+1]; print $1\\\"\\t\\\"$4\\\"\\t\\\"$5\\\"\\t\\\"gid}}\' > {wsl_path(gene_bed)}"'
os.system(cmd)
n_genes = int(os.popen(f'{WSL_PREFIX} "wc -l < {wsl_path(gene_bed)}"').read().strip())
print(f"  {n_genes} genes in BED", flush=True)

# Step 2: For each CGmap file, extract CpG → aggregate methylation per gene
cgmap_files = sorted(glob.glob(os.path.join(CGMAP_DIR, "*.CGmap.gz")))
print(f"\nSTEP 2: Processing {len(cgmap_files)} files via WSL awk...", flush=True)

# {gene_id: {cell_id: [methyl_values]}}
gene_methyl = {}

for fpath in cgmap_files:
    fname = os.path.basename(fpath)
    cell_id = "_".join(fname.split("_")[1:]).replace(".CGmap.gz", "")
    fstart = time.time()
    fp_wsl = wsl_path(fpath)
    
    # Step 2a: Extract CpG sites: chrom, pos, methyl_level
    cpg_file = os.path.join(OUT_DIR, f"_tmp_{cell_id}_cpg.bed").replace("\\", "/")
    cpg_wsl = wsl_path(cpg_file)
    cmd = f'{WSL_PREFIX} "zcat {fp_wsl} | awk \'$4==\\"CG\\" {{print $1\\\"\\t\\\"$3\\\"\\t\\\"$3+1\\\"\\t\\\"$6}}\' > {cpg_wsl}"'
    os.system(cmd)
    
    # Step 2b: Map CpGs to genes using bedtools map
    result_file = os.path.join(OUT_DIR, f"_tmp_{cell_id}_gene_methyl.txt").replace("\\", "/")
    res_wsl = wsl_path(result_file)
    gene_bed_wsl = wsl_path(gene_bed)
    
    bedtools_cmd = f"bedtools map -a {gene_bed_wsl} -b {cpg_wsl} -c 4 -o mean -null 0"
    cmd = f'{WSL_PREFIX} "{bedtools_cmd} > {res_wsl}"'
    os.system(cmd)
    
    # Step 2c: Read results
    n_genes_with_methyl = 0
    if os.path.exists(result_file):
        with open(result_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 5 and parts[4] != "." and parts[4] != "0" and float(parts[4]) > 0:
                    gid = parts[3]
                    methyl = float(parts[4])
                    if gid not in gene_methyl:
                        gene_methyl[gid] = {}
                    gene_methyl[gid][cell_id] = [methyl]  # bedtools mean = single value
                    n_genes_with_methyl += 1
    
    # Clean up temp files
    for tmp in [cpg_file, result_file]:
        if os.path.exists(tmp):
            os.remove(tmp)
    
    elapsed = time.time() - fstart
    print(f"  {cell_id}: {n_genes_with_methyl} genes with CpG ({elapsed:.0f}s)", flush=True)

# Step 3: Methylation matrix
print(f"\nSTEP 3: Building methylation matrix...", flush=True)
meth_data = {}
for gid, cells in gene_methyl.items():
    for cid, vals in cells.items():
        if cid not in meth_data:
            meth_data[cid] = {}
        meth_data[cid][gid] = np.mean(vals)

meth_df = pd.DataFrame(meth_data).T
meth_df = meth_df.fillna(meth_df.median())
print(f"  {meth_df.shape}", flush=True)
meth_df.to_csv(os.path.join(OUT_DIR, "methylation_matrix.csv"))
print(f"  Saved methylation_matrix.csv", flush=True)

# Step 4: MOFA+ (same as before)
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
