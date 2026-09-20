#!/usr/bin/env python
"""PA mechanism validation: test if MOFA+ F1 key genes are dysregulated in PA embryos.

Uses SRP301735 bulk RNA-seq DESeq2 results to check:
  1. DNMT1 and ATF3 expression in IVF vs PA across cleavage stages
  2. MOFA+ F1 top genes as a module — is the module significantly PA-downregulated?
  3. GSEA-like test: Are F1 positive-weight genes enriched among PA-downregulated?
  4. Additional check: do the methylation-driven genes (F1 METH top) behave differently?

Inputs:
  data/processed/deseq2/m15_tximport_deseq2.py output (or raw salmon quant)
  data/processed/m4_mofa/factor_top_genes_v2.csv

If DESeq2 results not available, re-compute from salmon quant + metadata.
"""
import os, sys
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed"

# ====== 1. Load MOFA+ F1 key genes ======
print("Loading MOFA+ F1 key genes ...", flush=True)
top_df = pd.read_csv(f"{BASE}/m4_mofa/factor_top_genes_v2.csv")

# F1 RNA_POS genes: the shared maternal program (DNMT1, ZP3, ZP4, GDF9, RARRES1, etc.)
f1_pos = top_df[(top_df["factor"]=="F1") & (top_df["view"]=="RNA_pos")]
f1_genes = sorted(set(("RARRES1;DNMT1;ZP3;EIF4G2;ZP4").split(";")))
f1_meth_genes = sorted(set(("ZP3;ZP4;DNMT1;GDF9;ZP2").split(";")))

# ATF3 and its key targets
atf3_genes = ["ATF3"]

# Metabolic genes from F4 METH
f4_meth = ["MDH1", "ALDH2", "PGAM1", "HAT1", "S100A10"]

# ZGA markers
zga_markers = ["DPPA5", "DNMT3A", "DNMT3B", "TET1", "POU5F1", "NANOG"]

# PA divergence genes from earlier analysis (KLHL15, MT-ATP6, MT-ND5)
pa_classic = ["KLHL15", "IDH2", "PKM", "NOP9", "GNL3", "EXOSC9"]

all_key = sorted(set(f1_genes + f1_meth_genes + atf3_genes + f4_meth + zga_markers + pa_classic))
print(f"  Key genes: {all_key}", flush=True)

# ====== 2. Load salmon quant data (if available) ======
# Try to find per-stage DESeq2 results first
deseq_dir = f"{BASE}/deseq2"
pa_results = None
salmon_matrix = None

# Check for raw salmon quant (the matrix used for DESeq2)
salmon_files = ["m15_salmon_quant.csv", "m15_salmon_counts.csv", "tximport_counts.csv",
                "salmon_gene_counts.csv"]
for f in salmon_files:
    path = os.path.join(BASE, "deseq2", f)
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        print(f"  Found: {path} ({os.path.getsize(path)//1024} KB)", flush=True)
        salmon_matrix = pd.read_csv(path, index_col=0)
        break

if salmon_matrix is None:
    # Fallback: attempt to compute PA effect from the DESeq2 results directory
    print("\n  No salmon quant matrix found. Checking DESeq2 output files...", flush=True)
    # Try to find any differential expression file
    for root, dirs, files in os.walk(deseq_dir):
        for f in files:
            if f.endswith('.csv') and any(kw in f.lower() for kw in ['deg', 'deseq', 'diff', 'de_']):
                fp = os.path.join(root, f)
                print(f"  Checking: {fp} ({os.path.getsize(fp)//1024} KB)", flush=True)
                try:
                    df = pd.read_csv(fp)
                    if 'gene' in df.columns or 'log2FoldChange' in df.columns:
                        print(f"    Found DE results! cols: {list(df.columns)}", flush=True)
                        pa_results = df
                        break
                except:
                    pass
        if pa_results is not None:
            break

if salmon_matrix is not None:
    print(f"\n  Salmon matrix: {salmon_matrix.shape}", flush=True)
    # Match genes
    salmon_genes = list(salmon_matrix.index) if salmon_matrix.index.name == 'gene' else list(salmon_matrix.columns)
    
    # Find key genes in salmon data
    found = [g for g in all_key if g in salmon_genes]
    print(f"  Key genes in salmon data: {found}", flush=True)
    
    if found:
        # Compute per-sample stage from column names
        # Columns might look like "IVF_2cell_rep1" or "PA_4cell_rep2" etc.
        cols = salmon_matrix.columns if salmon_matrix.index.name != 'gene' else salmon_matrix.index
        print(f"  Sample columns (first 5): {list(salmon_matrix.columns)[:5]}", flush=True)
        
elif pa_results is not None:
    print(f"\n  Found DESeq2 results, using those directly", flush=True)
    if 'gene' in pa_results.columns and 'log2FoldChange' in pa_results.columns:
        for g in all_key:
            if g in pa_results['gene'].values:
                row = pa_results[pa_results['gene'] == g]
                lfc = row['log2FoldChange'].values[0] if 'log2FoldChange' in pa_results.columns else None
                padj = row['padj'].values[0] if 'padj' in pa_results.columns else None
                print(f"    {g}: log2FC={lfc:.3f}, padj={padj}", flush=True)

# ====== 3. Module-level PA test ======
print("\nModule-level PA dysregulation test:")

# Load MOFA+ F1 positive-weight genes
f1_weights = pd.read_csv(f"{BASE}/m4_mofa/mofa_multiomics_weights_RNA_v2.csv", index_col=0)
f1_pos_weights = f1_weights["F1"].sort_values(ascending=False)

# Top 100 F1-positive genes
f1_top100 = set(f1_pos_weights.head(100).index)
print(f"  F1 top-100 positive-weight genes: {len(f1_top100)}", flush=True)

# Try to find these in DE results
if pa_results is not None and 'gene' in pa_results.columns:
    de_genes = set(pa_results['gene'])
    f1_in_de = f1_top100 & de_genes
    print(f"  F1 top100 in DE results: {len(f1_in_de)}", flush=True)
    
    if 'log2FoldChange' in pa_results.columns:
        f1_de = pa_results[pa_results['gene'].isin(f1_in_de)].copy()
        mean_lfc_f1 = f1_de['log2FoldChange'].mean()
        all_lfc = pa_results['log2FoldChange'].dropna()
        print(f"  F1 top100 mean log2FC: {mean_lfc_f1:.3f}")
        print(f"  All DE genes mean log2FC: {all_lfc.mean():.3f}")
        # MWU test
        in_f1 = pa_results[pa_results['gene'].isin(f1_in_de)]['log2FoldChange'].dropna().values
        not_f1 = pa_results[~pa_results['gene'].isin(f1_in_de)]['log2FoldChange'].dropna().values
        stat, p = mannwhitneyu(in_f1, not_f1, alternative="two-sided")
        print(f"  Mann-Whitney: in-F1 vs not-F1, stat={stat:.0f}, p={p:.2e}", flush=True)

# ====== 4. Try DESeq2 from python (if needed) ======
# If we have salmon quant but no DESeq2 results, run DESeq2
if salmon_matrix is not None:
    print("\nAttempting DESeq2 from salmon matrix ...", flush=True)
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        
        # Matrix format check: need genes as rows, samples as columns
        if salmon_matrix.index.name == 'gene' or salmon_matrix.shape[0] < salmon_matrix.shape[1]:
            counts = salmon_matrix
        else:
            counts = salmon_matrix.T
        
        # Build metadata from column names
        samples = list(counts.columns)
        conditions = ["IVF" if "IVF" in s else "PA" if "PA" in s else "Unknown" for s in samples]
        stages = [s.split("_")[1] if "_" in s else "Unknown" for s in samples]
        
        meta_df = pd.DataFrame({"sample": samples, "condition": conditions, "stage": stages})
        print(f"  Conditions: {meta_df['condition'].unique()}", flush=True)
        print(f"  Stages: {meta_df['stage'].value_counts().to_dict()}", flush=True)
        
        # DESeq2 per stage
        results = []
        for stage in meta_df["stage"].unique():
            stg_samples = meta_df[meta_df["stage"]==stage]["sample"].tolist()
            if len(stg_samples) < 4:
                continue
            stg_meta = meta_df[meta_df["stage"]==stage].set_index("sample")
            stg_counts = counts[stg_samples].round().astype(int)
            stg_counts = stg_counts.loc[stg_counts.sum(axis=1) > 10]  # filter low expressing
            
            if len(stg_counts) == 0:
                continue
            
            dds = DeseqDataSet(counts=stg_counts, metadata=stg_meta, design="~condition")
            dds.deseq2()
            stat_res = DeseqStats(dds, contrast=["condition", "PA", "IVF"])
            stat_res.summary()
            res = stat_res.results_df.copy()
            res["gene"] = res.index
            res["stage"] = stage
            
            for g in all_key:
                if g in res.index:
                    r = res.loc[g]
                    lfc = r.get("log2FoldChange", 0)
                    p = r.get("padj", 1)
                    print(f"    {stage} {g}: log2FC={lfc:.3f}, padj={p:.3f}", flush=True)
            
            results.append(res)
        
        if results:
            all_res = pd.concat(results)
            all_res.to_csv(f"{BASE}/deseq2/m23_pa_validation_deseq.csv", index=False)
            print(f"\n  Saved DESeq2 results: {len(all_res)} rows", flush=True)
            
    except ImportError:
        print("  pydeseq2 not available, cannot run DESeq2.", flush=True)
    except Exception as e:
        print(f"  DESeq2 failed: {e}", flush=True)

# ====== 5. PA deviation test from earlier pseudotime ======
print("\nChecking PA deviation from trajectory analysis ...", flush=True)
traj_file = f"{BASE}/m5_gene_divergence_v2.csv"
if os.path.exists(traj_file):
    div = pd.read_csv(traj_file)
    if "gene" in div.columns and "PA_invivo" in div.columns:
        # Check F1 key genes in divergence
        for g in f1_genes:
            if g in div["gene"].values:
                row = div[div["gene"]==g]
                val = row["PA_invivo"].values[0]
                print(f"  {g}: PA-vs-invivo divergence = {val:.3f}", flush=True)

print("\n=== PA VALIDATION SUMMARY ===")
print("Key genes checked:", ", ".join(all_key))
print("\nDONE.")