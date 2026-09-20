#!/usr/bin/env python
"""
FPKM → raw counts back-calculation → DESeq2 (IVF vs PA pseudobulk).
Replaces the original Wilcoxon analysis with DESeq2 on pseudo-counts.
Maintains honest scope (approximation, not gold-standard FASTQ-derived counts).
"""
import os, re
import pandas as pd
import numpy as np

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
FPKM_FILE = os.path.join(BASE, "data/raw/GSE164812/GSE164812_gene_FPKM_matrix.txt.gz")
OUT_DIR = os.path.join(BASE, "data/processed/m3c_deseq2")
TRANSCRIPTOME = os.path.join(BASE, "data/raw/pig_transcriptome.fa.gz")
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("FPKM → DESeq2: back-calculating counts from GSE164812 FPKM matrix")
print("=" * 60)

# 1. Build transcript length lookup from pig FASTA
import gzip
# Parse FASTA to get lengths
seq_lens = {}
with gzip.open(TRANSCRIPTOME, "rt") as f:
    name, seq = None, []
    for line in f:
        if line.startswith(">"):
            if name and seq:
                seq_lens[name] = len("".join(seq))
            name = line.strip().split()[0][1:]
            seq = []
        else:
            seq.append(line.strip())
    if name and seq:
        seq_lens[name] = len("".join(seq))

print(f"\nLoaded {len(seq_lens)} transcript lengths from FASTA")
median_len = np.median(list(seq_lens.values()))
print(f"Median transcript length: {median_len:.0f} bp")

# 2. Load FPKM matrix
df = pd.read_csv(FPKM_FILE, sep="\t", index_col=0)
sample_cols = [c for c in df.columns if c not in ["Biotype", "Position", "GeneName"]]
print(f"\nLoaded FPKM matrix: {df.shape} genes × {len(sample_cols)} samples")
fpkm = df[sample_cols].astype(float)
genes = list(fpkm.index)  # ENSSSCG IDs

# 3. Create gene-to-length mapping
# For each gene ID in FPKM, assign best transcript length
# The FPKM gene IDs are ENSSSCGxxxxx
# Transcript IDs are ENSSSCTxxxxx
# They share the numeric part: ENSSSCG00000037674 → ENSSSCT00000037674? No, different numbers.
# We'll use median length as approximation
gene_lengths = pd.Series(median_len, index=genes)

# 4. Back-calculate counts from FPKM
# FPKM = count × 10^9 / (length × total_counts)
# count = FPKM × length_kb × total_counts / 10^6
# For SMART-seq2 single cell, estimated 2M mapped reads per cell
EST_LIB_SIZE = 2_000_000  # reasonable estimate for single-cell SMART-seq2
print(f"\nUsing estimated library size: {EST_LIB_SIZE/1e6:.1f}M reads/cell")
length_kb = gene_lengths / 1000.0

pseudo_counts = pd.DataFrame(index=genes)
for col in sample_cols:
    fpkm_col = fpkm[col].values
    counts = fpkm_col * length_kb.values * EST_LIB_SIZE / 1e6
    pseudo_counts[col] = np.round(counts).astype(int)

print(f"\nPseudo-count matrix: {pseudo_counts.shape}, total reads: {pseudo_counts.sum().sum()/1e6:.0f}M")

# 5. Parse condition and stage from column names
sample_info = []
for col in sample_cols:
    parts = col.split()
    cond = "IVF" if "IVF" in col else "PA"
    stage_extracted = re.search(r"(\d+)-cell", col)
    stage = stage_extracted.group(1) if stage_extracted else "unknown"
    sample_info.append({"sample": col, "condition": cond, "stage": stage})

meta = pd.DataFrame(sample_info)
print(f"\nCondition counts: {meta['condition'].value_counts().to_dict()}")
print(f"Stage counts: {meta['stage'].value_counts().to_dict()}")

# 6. Pseudobulk by condition (sum counts within condition)
pseudobulk = pd.DataFrame()
for cond in ["IVF", "PA"]:
    cols = meta[meta["condition"] == cond]["sample"].tolist()
    pseudobulk[cond] = pseudo_counts[cols].sum(axis=1)

# Keep only genes with at least 10 total counts across both conditions
pseudobulk = pseudobulk.loc[(pseudobulk.sum(axis=1) >= 10)]
print(f"\nPseudobulk after filtering: {pseudobulk.shape}")

# 7. DESeq2 with pydeseq2 (if available)
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    
    # With only 2 conditions (1 pseudobulk per condition), we need LRT
    # But LRT requires a reduced design
    # Use simplest approach: compute log2FC + approximate p-value
    
    # Actually, with only 2 pseudobulk groups (IVF, PA), we can't do Wald test
    # because we need within-group variance estimates (n≥2 per group)
    # 
    # Alternative: use the individual sample counts with "stage" as covariate
    # This gives us n=42 observations (21 IVF + 21 PA) at the cell level
    # But this reintroduces pseudoreplication
    #
    # Best practical approach:
    # Pseudobulk by stage-condition, giving 8 groups (4 stages × 2 conditions)
    
    print("\nCreating stage-condition pseudobulk groups...")
    # Aggregate by stage-condition
    pb_stage = {}
    for (cond, stage), group in meta.groupby(["condition", "stage"]):
        key = f"{cond}_{stage}"
        cols = group["sample"].tolist()
        pb_stage[key] = pseudo_counts[cols].sum(axis=1)
    
    pb_df = pd.DataFrame(pb_stage)
    pb_df = pb_df.loc[(pb_df.sum(axis=1) >= 10)]
    
    # Create metadata for pseudobulk samples
    pb_meta = pd.DataFrame({"sample": list(pb_stage.keys())})
    pb_meta["condition"] = pb_meta["sample"].apply(lambda x: x.split("_")[0])
    pb_meta["stage"] = pb_meta["sample"].apply(lambda x: x.split("_")[1])
    pb_meta = pb_meta.set_index("sample")
    
    print(f"Pseudobulk-by-stage: {pb_df.shape}, groups: {list(pb_stage.keys())}")
    
    # Run DESeq2
    dds = DeseqDataSet(counts=pb_df.T.astype(int), metadata=pb_meta,
                       design="~ stage + condition", refit_cooks=False)
    dds.deseq2()
    
    # Extract IVF vs PA comparison
    stat_res = DeseqStats(dds, contrast=["condition", "IVF", "PA"])
    stat_res.summary()
    
    results = stat_res.results_df
    print(f"\nDESeq2 results: {results.shape}")
    print(f"Significant (padj < 0.05): {(results['padj'] < 0.05).sum()}")
    
    # Save
    results.to_csv(os.path.join(OUT_DIR, "deseq2_results.csv"))
    print(f"Saved: {os.path.join(OUT_DIR, 'deseq2_results.csv')}")
    
    # Compare with original Wilcoxon results
    print(f"\n=== Comparison with Wilcoxon ===")
    print(f"Top 20 DESeq2 genes (by padj):")
    top = results.sort_values("padj").head(20)
    for gene, row in top.iterrows():
        print(f"  {gene}: log2FC={row['log2FC']:.3f}, padj={row['padj']:.2e}")

except ImportError as e:
    print(f"\npydeseq2 not available ({e})")
    print("Computing log2FC from pseudobulk...")
    pb_df["log2FC_IVF_vs_PA"] = np.log2(pb_df.filter(like="IVF").sum(axis=1) + 1).values - np.log2(pb_df.filter(like="PA").sum(axis=1) + 1).values
    
    # Compare with original Wilcoxon top genes
    print(f"\nTop 20 genes by pseudobulk log2FC:")
    pb_df["absFC"] = pb_df["log2FC_IVF_vs_PA"].abs()
    top = pb_df.sort_values("absFC", ascending=False).head(20)
    for gene, row in top.iterrows():
        print(f"  {gene}: log2FC={row['log2FC_IVF_vs_PA']:.3f}")

print(f"\n{'='*60}")
print("Done. These are pseudo-count derived from FPKM.")
print("Gold-standard DESeq2 from raw FASTQ is pending (SRA salmon pipeline).")
