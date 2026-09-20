"""Debug: process first CGmap file only, see where it hangs."""
import gzip, time, sys
from collections import defaultdict

fpath = "E:/Workbuddy/2026-07-27-11-58-27/data/raw/cgmap/GSM7508586_AF1_12.CGmap.gz"
start = time.time()

n_lines = 0
n_cg = 0
chrom_counts = defaultdict(int)

with gzip.open(fpath, "rt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 6:
            n_lines += 1
            continue
        n_lines += 1
        if parts[3] == "CG":
            n_cg += 1
            chrom_counts[parts[0]] += 1
        
        if n_lines % 1000000 == 0:
            elapsed = time.time() - start
            print(f"  {n_lines/1e6:.0f}M lines, {n_cg/1e6:.1f}M CpG ({elapsed:.0f}s)", flush=True)
        
        if n_lines >= 5000000:
            break

elapsed = time.time() - start
print(f"\nScanned {n_lines:,} lines ({n_cg:,} CpG) in {elapsed:.1f}s", flush=True)
print(f"Chromosome distribution:")
for chrom, count in sorted(chrom_counts.items(), key=lambda x: -x[1])[:10]:
    print(f"  {chrom}: {count:,} CpG")
