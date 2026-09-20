#!/usr/bin/env python
"""
Phase R1: CpG-level methylation reanalysis from raw CGmap files.
Streams through 44 CGmap.gz files (~45 GB), maps CpG sites to gene regions
using pig Sscrofa11.1 GTF, and builds per-gene methylation matrices.

Outputs:
  1. methylation_promoter.csv      (TSS +/- 2kb)
  2. methylation_genebody.csv      (gene body)
  3. methylation_combined.csv      (promoter + gene body)
  4. cpg_qc_report.csv             (per-sample QC: CpG count, coverage, etc.)
  5. per_gene_cpg_stats.csv        (per-gene: CpG count per region, mean coverage)

Processes 44 files @ ~29M lines each (~8M CG sites) = ~350M CG sites total.
Target: < 2 hours on Windows Python.
"""
import gzip, os, re, sys, time
from collections import defaultdict
import pandas as pd
import numpy as np

# ── Paths ──
CG_DIR = r'E:/Workbuddy/2026-07-27-11-58-27/data/raw/cgmap'
GTF    = r'E:/Workbuddy/2026-07-27-11-58-27/data/reference/Sus_scrofa.Sscrofa11.1.113.gtf.gz'
OUT    = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa'
os.makedirs(OUT, exist_ok=True)

TSS_WINDOW = 2000  # promoter = TSS +/- 2kb

# ── Step 1: Parse GTF, build gene region index ──
print("=" * 60)
print("STEP 1: Building gene region index from GTF...")
t0 = time.time()

# gene_regions: {chr: [(start, end, strand, gene_id, gene_name), ...]}
gene_regions = defaultdict(list)
gene_info = {}  # gene_id -> {tss, strand, name}

with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        chrom = parts[0]
        start = int(parts[3]) - 1  # 0-based
        end   = int(parts[4])
        strand = parts[6]
        attrs = parts[8]

        gid_match = re.search(r'gene_id "([^"]+)"', attrs)
        gname_match = re.search(r'gene_name "([^"]+)"', attrs)
        if not gid_match:
            continue
        gid = gid_match.group(1)
        gname = gname_match.group(1) if gname_match else gid

        # Determine TSS based on strand
        if strand == '+':
            tss = start
        else:
            tss = end - 1  # 0-based TSS for negative strand

        gene_info[gid] = {
            'tss': tss,
            'strand': strand,
            'gene_body_start': start,
            'gene_body_end': end,
            'gene_name': gname,
            'chrom': chrom,
        }

        # Promoter region: TSS +/- TSS_WINDOW
        prom_start = max(0, tss - TSS_WINDOW)
        prom_end = tss + TSS_WINDOW

        gene_regions[chrom].append((prom_start, prom_end, strand, gid, gname, 'promoter'))
        gene_regions[chrom].append((start, end, strand, gid, gname, 'gene_body'))

# Sort intervals per chromosome for binary search
for chrom in gene_regions:
    gene_regions[chrom].sort(key=lambda x: x[0])

n_genes = len(gene_info)
n_regions = sum(len(v) for v in gene_regions.values())
print(f"  Loaded {n_genes} genes, {n_regions} regions across {len(gene_regions)} chromosomes")
print(f"  Time: {time.time() - t0:.1f}s")

# ── Step 2: Pre-build interval lookup for fast CG→gene mapping ──
# For each chrom, build sorted starts for binary search
chrom_intervals = {}
for chrom, intervals in gene_regions.items():
    starts = np.array([iv[0] for iv in intervals])
    ends   = np.array([iv[1] for iv in intervals])
    chrom_intervals[chrom] = (starts, ends, intervals)

def find_overlapping_genes(chrom, pos):
    """Find all gene regions overlapping a CpG position using binary search."""
    if chrom not in chrom_intervals:
        return []
    starts, ends, intervals = chrom_intervals[chrom]

    # Binary search to find candidate regions
    idx = np.searchsorted(starts, pos, side='right') - 1
    if idx < 0:
        idx = 0

    results = []
    # Check intervals around idx (some may extend past pos)
    for i in range(max(0, idx - 5), min(len(intervals), idx + 10)):
        if starts[i] <= pos < ends[i]:
            results.append(intervals[i])
    return results

# ── Step 3: Process CGmap files ──
print("\n" + "=" * 60)
print("STEP 3: Streaming CGmap files, mapping CpG sites to genes...")

cgmap_files = sorted([
    f for f in os.listdir(CG_DIR)
    if f.endswith('.CGmap.gz') and f.startswith('GSM')
])

# Data structures: per-gene accumulation
# For each sample: {gene_id: {'promoter': {'meth_sum','total_sum','n_cpg'}, 'gene_body': {...}}}
all_samples = {}
sample_qc = {}

for fi, fname in enumerate(cgmap_files):
    sample_id = fname.replace('.CGmap.gz', '')
    fpath = os.path.join(CG_DIR, fname)
    t1 = time.time()

    # Per-sample accumulators
    gene_data = defaultdict(lambda: {
        'promoter_meth': [], 'promoter_cov': [],
        'genebody_meth': [], 'genebody_cov': [],
    })

    total_cpg = 0
    mapped_cpg = 0
    cg_context = 0

    with gzip.open(fpath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            chrom = parts[0]
            pos = int(parts[2])
            context = parts[3]
            meth_level = float(parts[5])
            meth_reads = int(parts[6])
            total_reads = int(parts[7])

            total_cpg += 1

            # Only CG context is CpG methylation
            if context != 'CG':
                continue
            cg_context += 1

            # Only include sites with coverage >= 1
            if total_reads < 1:
                continue

            # Find overlapping gene regions
            hits = find_overlapping_genes(chrom, pos)
            if not hits:
                continue

            mapped_cpg += 1

            for _, _, _, gid, gname, region_type in hits:
                if region_type == 'promoter':
                    gene_data[gid]['promoter_meth'].append(meth_level)
                    gene_data[gid]['promoter_cov'].append(total_reads)
                else:
                    gene_data[gid]['genebody_meth'].append(meth_level)
                    gene_data[gid]['genebody_cov'].append(total_reads)

    # Aggregate per gene
    agg = {}
    for gid, gd in gene_data.items():
        agg[gid] = {
            'promoter_mean_meth': np.mean(gd['promoter_meth']) if gd['promoter_meth'] else np.nan,
            'promoter_n_cpg': len(gd['promoter_meth']),
            'promoter_mean_cov': np.mean(gd['promoter_cov']) if gd['promoter_cov'] else 0,
            'genebody_mean_meth': np.mean(gd['genebody_meth']) if gd['genebody_meth'] else np.nan,
            'genebody_n_cpg': len(gd['genebody_meth']),
            'genebody_mean_cov': np.mean(gd['genebody_cov']) if gd['genebody_cov'] else 0,
        }

    all_samples[sample_id] = agg
    sample_qc[sample_id] = {
        'total_sites': total_cpg,
        'cg_sites': cg_context,
        'cg_with_cov': sum(1 for _ in []) if False else cg_context,  # placeholder
        'mapped_to_genes': mapped_cpg,
        'genes_with_promoter': sum(1 for gid, gd in gene_data.items() if gd['promoter_meth']),
        'genes_with_genebody': sum(1 for gid, gd in gene_data.items() if gd['genebody_meth']),
    }

    elapsed = time.time() - t1
    print(f"  [{fi+1:2d}/{len(cgmap_files)}] {sample_id}  "
          f"CG={cg_context//1000}K  mapped={mapped_cpg//1000}K  "
          f"genes_prom={sample_qc[sample_id]['genes_with_promoter']}  "
          f"genes_body={sample_qc[sample_id]['genes_with_genebody']}  "
          f"({elapsed:.1f}s)")

# ── Step 4: Build methylation matrices ──
print("\n" + "=" * 60)
print("STEP 4: Building methylation matrices...")

# Determine gene universe: genes with promoter data in >= 50% of samples
gene_presence = defaultdict(int)
for sample_id, agg in all_samples.items():
    for gid in agg:
        if not np.isnan(agg[gid]['promoter_mean_meth']):
            gene_presence[gid] += 1

min_samples = max(1, len(all_samples) // 2)
common_genes = sorted([g for g, cnt in gene_presence.items() if cnt >= min_samples])
print(f"  Genes with promoter meth in >= {min_samples} samples: {len(common_genes)}")

# Build matrices
samples = sorted(all_samples.keys())
n_samples = len(samples)
n_genes = len(common_genes)

prom_mat = np.full((n_genes, n_samples), np.nan)
body_mat = np.full((n_genes, n_samples), np.nan)
comb_mat = np.full((n_genes, n_samples), np.nan)
prom_ncpg_mat = np.zeros((n_genes, n_samples), dtype=int)
prom_cov_mat = np.zeros((n_genes, n_samples))

for gi, gid in enumerate(common_genes):
    for si, sid in enumerate(samples):
        if gid in all_samples[sid]:
            d = all_samples[sid][gid]
            if not np.isnan(d['promoter_mean_meth']):
                prom_mat[gi, si] = d['promoter_mean_meth']
                prom_ncpg_mat[gi, si] = d['promoter_n_cpg']
                prom_cov_mat[gi, si] = d['promoter_mean_cov']
            if not np.isnan(d['genebody_mean_meth']):
                body_mat[gi, si] = d['genebody_mean_meth']
            # Combined: average promoter + gene body if both available
            if (not np.isnan(d['promoter_mean_meth']) and
                not np.isnan(d['genebody_mean_meth'])):
                comb_mat[gi, si] = (d['promoter_mean_meth'] + d['genebody_mean_meth']) / 2
            elif not np.isnan(d['promoter_mean_meth']):
                comb_mat[gi, si] = d['promoter_mean_meth']
            elif not np.isnan(d['genebody_mean_meth']):
                comb_mat[gi, si] = d['genebody_mean_meth']

# Build gene name lookup
gene_names = [gene_info[gid]['gene_name'] if gid in gene_info else gid for gid in common_genes]

# ── Step 5: Save matrices ──
prom_df = pd.DataFrame(prom_mat, index=common_genes, columns=samples)
prom_df.index.name = 'gene_id'
prom_df.to_csv(os.path.join(OUT, 'methylation_promoter.csv'))
print(f"  Saved methylation_promoter.csv: {prom_df.shape}")

body_df = pd.DataFrame(body_mat, index=common_genes, columns=samples)
body_df.index.name = 'gene_id'
body_df.to_csv(os.path.join(OUT, 'methylation_genebody.csv'))
print(f"  Saved methylation_genebody.csv: {body_df.shape}")

comb_df = pd.DataFrame(comb_mat, index=common_genes, columns=samples)
comb_df.index.name = 'gene_id'
comb_df.to_csv(os.path.join(OUT, 'methylation_combined.csv'))
print(f"  Saved methylation_combined.csv: {comb_df.shape}")

# ── Step 6: QC report ──
qc_rows = []
for sid in samples:
    qc_rows.append({
        'sample': sid,
        **sample_qc[sid],
        'genes_final': sum(1 for gi in range(n_genes) if not np.isnan(prom_mat[gi, samples.index(sid)])),
        'mean_cpg_per_gene': np.nanmean(prom_ncpg_mat[:, samples.index(sid)]),
        'median_cpg_per_gene': np.nanmedian(prom_ncpg_mat[:, samples.index(sid)]),
    })
qc_df = pd.DataFrame(qc_rows)
qc_df.to_csv(os.path.join(OUT, 'cpg_qc_report.csv'), index=False)
print(f"  Saved cpg_qc_report.csv: {len(qc_df)} samples")

# Per-gene CpG stats
gene_stats = []
for gi, gid in enumerate(common_genes):
    gene_stats.append({
        'gene_id': gid,
        'gene_name': gene_names[gi],
        'mean_promoter_ncpg': np.nanmean(prom_ncpg_mat[gi, :]),
        'median_promoter_ncpg': np.nanmedian(prom_ncpg_mat[gi, :]),
        'mean_promoter_cov': np.nanmean(prom_cov_mat[gi, :]),
        'samples_with_data': np.sum(~np.isnan(prom_mat[gi, :])),
        'mean_methylation': np.nanmean(prom_mat[gi, :]),
    })
gene_stats_df = pd.DataFrame(gene_stats)
gene_stats_df.to_csv(os.path.join(OUT, 'per_gene_cpg_stats.csv'), index=False)
print(f"  Saved per_gene_cpg_stats.csv: {len(gene_stats_df)} genes")

# ── Step 7: Compare with existing methylation matrix ──
print("\n" + "=" * 60)
print("STEP 7: Comparing with existing methylation matrix...")
old_meth = pd.read_csv(os.path.join(OUT, 'methylation_matrix.csv'), index_col=0)
old_meth_in_genes = [g for g in old_meth.index if g in common_genes]
print(f"  Old matrix: {old_meth.shape[0]} genes, {old_meth.shape[1]} samples")
print(f"  Old genes in new gene set: {len(old_meth_in_genes)} / {old_meth.shape[0]}")

# Compare old vs new (combined) for shared genes/samples
shared_genes = sorted(set(old_meth.index) & set(common_genes))
shared_samples = sorted(set(old_meth.columns) & set(samples))
if shared_genes and shared_samples:
    old_sub = old_meth.loc[shared_genes, shared_samples]
    new_sub = comb_df.loc[shared_genes, shared_samples]
    # Correlation per gene
    cors = []
    for g in shared_genes:
        old_vals = old_sub.loc[g].values
        new_vals = new_sub.loc[g].values
        mask = ~(np.isnan(old_vals) | np.isnan(new_vals))
        if mask.sum() >= 5:
            cors.append(np.corrcoef(old_vals[mask], new_vals[mask])[0, 1])
    print(f"  Old vs new methylation correlation (per gene, n={len(cors)}): "
          f"mean={np.mean(cors):.3f}, median={np.median(cors):.3f}")

print(f"\n{'='*60}")
print(f"DONE. Total time: {(time.time() - t0)/60:.1f} min")
