#!/usr/bin/env python
"""Optimized v2: vectorized CpG→gene mapping using numpy instead of per-site Python loop."""
import gzip, os, re, sys, time
from collections import defaultdict
import pandas as pd
import numpy as np

CG_DIR = r'E:/Workbuddy/2026-07-27-11-58-27/data/raw/cgmap'
GTF    = r'E:/Workbuddy/2026-07-27-11-58-27/data/reference/Sus_scrofa.Sscrofa11.1.113.gtf.gz'
OUT    = r'E:/Workbuddy/2026-07-27-11-58-27/data/processed/m4_mofa'
TSS_WINDOW = 2000

os.makedirs(OUT, exist_ok=True)

# ── Step 1: GTF gene regions ──
print("STEP 1: Building gene regions...")
t0 = time.time()
gene_regions = defaultdict(list)
gene_info = {}

with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene': continue
        chrom = parts[0]
        start = int(parts[3]) - 1
        end   = int(parts[4])
        strand = parts[6]
        attrs = parts[8]
        gid_m = re.search(r'gene_id "([^"]+)"', attrs)
        gname_m = re.search(r'gene_name "([^"]+)"', attrs)
        if not gid_m: continue
        gid = gid_m.group(1)
        gname = gname_m.group(1) if gname_m else gid
        tss = start if strand == '+' else end - 1
        gene_info[gid] = {'tss':tss, 'strand':strand, 'gene_body_start':start,
                          'gene_body_end':end, 'gene_name':gname, 'chrom':chrom}
        prom_start = max(0, tss - TSS_WINDOW)
        prom_end = tss + TSS_WINDOW
        gene_regions[chrom].extend([
            (prom_start, prom_end, strand, gid, gname, 0),  # 0=promoter
            (start, end, strand, gid, gname, 1),              # 1=gene_body
        ])

# Sort and build numpy arrays per chromosome
print(f"  Loaded {len(gene_info)} genes, {sum(len(v) for v in gene_regions.values())} regions")
# Flatten
chrom_data = {}
for chrom, intervals in gene_regions.items():
    arr = np.array(sorted(intervals, key=lambda x: x[0]), dtype=object)
    chrom_data[chrom] = {
        'starts': arr[:, 0].astype(np.int64),
        'ends':   arr[:, 1].astype(np.int64),
        'info':   arr[:, 2:],  # strand, gid, gname, type
    }
print(f"  Chrom index: {len(chrom_data)} chrs, {sum(d['starts'].shape[0] for d in chrom_data.values())} regions")
print(f"  Time: {time.time()-t0:.1f}s")

# ── Step 2: Process CGmap files (optimized per-file batch) ──
cgmap_files = sorted([f for f in os.listdir(CG_DIR) if f.endswith('.CGmap.gz') and f.startswith('GSM')])

all_samples = {}
print(f"\nSTEP 2: Processing {len(cgmap_files)} CGmap files...")

for fi, fname in enumerate(cgmap_files):
    sample_id = fname.replace('.CGmap.gz', '')
    fpath = os.path.join(CG_DIR, fname)
    t1 = time.time()

    # Read all CG CpG positions in one pass (int codes, not object arrays)
    chr_codes = []; positions = []; meth_vals = []; covs = []
    chr_to_code = {}
    with gzip.open(fpath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8 or parts[3] != 'CG':
                continue
            total = int(parts[7])
            if total < 1:
                continue
            ch = parts[0]
            if ch not in chr_to_code:
                chr_to_code[ch] = len(chr_to_code)
            chr_codes.append(chr_to_code[ch])
            positions.append(int(parts[2]))
            meth_vals.append(float(parts[5]))
            covs.append(total)
    code_to_chr = {v: k for k, v in chr_to_code.items()}
    chr_arr = np.array(chr_codes, dtype=np.int32)
    pos_all = np.array(positions, dtype=np.int64)
    meth_all = np.array(meth_vals, dtype=np.float32)
    cov_all = np.array(covs, dtype=np.int32)
    n_cg = len(pos_all)
    
    # Per-gene running sums (memory-efficient, no list storage)
    gene_data = defaultdict(lambda: {'prom_sum':0.0, 'prom_n':0, 'prom_cov_sum':0.0,
                                      'body_sum':0.0, 'body_n':0, 'body_cov_sum':0.0})
    
    # Process chromosome by chromosome (vectorized per chrom)
    mapped = 0
    for code in sorted(code_to_chr.keys()):
        chrom = code_to_chr[code]
        if chrom not in chrom_data:
            continue
        cd = chrom_data[chrom]
        starts = cd['starts']
        ends = cd['ends']
        
        # Get CG sites on this chromosome
        mask = chr_arr == code
        cg_pos = pos_all[mask]
        cg_meth = meth_all[mask]
        cg_cov = cov_all[mask]
        if len(cg_pos) == 0:
            continue
        
        # Vectorized search
        idx_right = np.searchsorted(starts, cg_pos, side='right') - 1
        idx_right = np.clip(idx_right, 0, len(starts) - 1)
        in_range = (cg_pos >= starts[idx_right]) & (cg_pos < ends[idx_right])
        mapped += in_range.sum()
        
        # Map to genes
        for j in np.where(in_range)[0]:
            info = cd['info'][idx_right[j]]
            gid = info[1]
            rtype = info[3]
            if rtype == 0:
                gene_data[gid]['prom_sum'] += float(cg_meth[j])
                gene_data[gid]['prom_n'] += 1
                gene_data[gid]['prom_cov_sum'] += int(cg_cov[j])
            else:
                gene_data[gid]['body_sum'] += float(cg_meth[j])
                gene_data[gid]['body_n'] += 1
                gene_data[gid]['body_cov_sum'] += int(cg_cov[j])

    # Aggregate (running sums → means)
    agg = {}
    for gid, gd in gene_data.items():
        agg[gid] = {
            'promoter_mean_meth': gd['prom_sum'] / gd['prom_n'] if gd['prom_n'] > 0 else np.nan,
            'promoter_n_cpg': gd['prom_n'],
            'promoter_mean_cov': gd['prom_cov_sum'] / gd['prom_n'] if gd['prom_n'] > 0 else 0,
            'genebody_mean_meth': gd['body_sum'] / gd['body_n'] if gd['body_n'] > 0 else np.nan,
            'genebody_n_cpg': gd['body_n'],
            'genebody_mean_cov': gd['body_cov_sum'] / gd['body_n'] if gd['body_n'] > 0 else 0,
        }
    all_samples[sample_id] = agg
    
    elapsed = time.time() - t1
    n_prom = sum(1 for g in gene_data if gene_data[g]['prom_n'] > 0)
    print(f"  [{fi+1:2d}/{len(cgmap_files)}] {sample_id}  CG={n_cg//1000}K  mapped={mapped//1000}K  "
          f"genes_p={n_prom}  ({elapsed:.1f}s)")

# ── Step 3: Build methylation matrices ──
print("\nSTEP 3: Building methylation matrices...")
gene_presence = defaultdict(int)
for sid, agg in all_samples.items():
    for gid in agg:
        if not np.isnan(agg[gid]['promoter_mean_meth']):
            gene_presence[gid] += 1

min_samples = max(1, len(all_samples) // 2)
common_genes = sorted([g for g, cnt in gene_presence.items() if cnt >= min_samples])
samples = sorted(all_samples.keys())
n_genes, n_samples = len(common_genes), len(samples)

prom_mat = np.full((n_genes, n_samples), np.nan)
body_mat = np.full((n_genes, n_samples), np.nan)
prom_ncpg_mat = np.zeros((n_genes, n_samples), dtype=int)

for gi, gid in enumerate(common_genes):
    for si, sid in enumerate(samples):
        if gid in all_samples[sid]:
            d = all_samples[sid][gid]
            if not np.isnan(d['promoter_mean_meth']):
                prom_mat[gi, si] = d['promoter_mean_meth']
                prom_ncpg_mat[gi, si] = d['promoter_n_cpg']
            if not np.isnan(d['genebody_mean_meth']):
                body_mat[gi, si] = d['genebody_mean_meth']

gene_names = [gene_info[gid]['gene_name'] if gid in gene_info else gid for gid in common_genes]

# Save
prom_df = pd.DataFrame(prom_mat, index=common_genes, columns=samples)
prom_df.to_csv(os.path.join(OUT, 'methylation_promoter.csv'))
body_df = pd.DataFrame(body_mat, index=common_genes, columns=samples)
body_df.to_csv(os.path.join(OUT, 'methylation_genebody.csv'))

# QC report
qc = []
for sid in samples:
    agg = all_samples[sid]
    gids_with_prom = sum(1 for gid in common_genes if gid in agg and not np.isnan(agg[gid]['promoter_mean_meth']))
    ncpg_vals = [agg[gid]['promoter_n_cpg'] for gid in common_genes if gid in agg]
    qc.append({'sample': sid, 'genes_with_promoter': gids_with_prom,
               'median_cpg_per_gene': np.median(ncpg_vals) if ncpg_vals else 0})
pd.DataFrame(qc).to_csv(os.path.join(OUT, 'cpg_qc_report.csv'), index=False)

# Gene stats
gene_stats = []
for gi, gid in enumerate(common_genes):
    gene_stats.append({'gene_id': gid, 'gene_name': gene_names[gi],
                       'mean_promoter_ncpg': np.nanmean(prom_ncpg_mat[gi, :]),
                       'samples_with_data': np.sum(~np.isnan(prom_mat[gi, :]))})
pd.DataFrame(gene_stats).to_csv(os.path.join(OUT, 'per_gene_cpg_stats.csv'), index=False)

# Compare with old
old_meth = pd.read_csv(os.path.join(OUT, 'methylation_matrix.csv'), index_col=0)
shared_genes = sorted(set(old_meth.index) & set(common_genes))
shared_samples = sorted(set(old_meth.columns) & set(samples))
print(f"  Old matrix: {old_meth.shape}  New: {prom_df.shape}")
print(f"  Shared genes: {len(shared_genes)}  Shared samples: {len(shared_samples)}")
if shared_genes and shared_samples:
    old_sub = old_meth.loc[shared_genes, shared_samples]
    new_sub = prom_df.loc[shared_genes, shared_samples]
    cors = []
    for g in shared_genes:
        ov = old_sub.loc[g].values; nv = new_sub.loc[g].values
        mask = ~(np.isnan(ov) | np.isnan(nv))
        if mask.sum() >= 5:
            cors.append(np.corrcoef(ov[mask], nv[mask])[0,1])
    print(f"  Old vs new methyl corr: median={np.median(cors):.3f} (n={len(cors)} genes)")

print(f"\nDONE. Total: {(time.time()-t0)/60:.1f} min")
