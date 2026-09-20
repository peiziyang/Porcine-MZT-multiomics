import gzip
GI = r"E:/迅雷下载/Sus_scrofa.gene_info.gz"
G2G = r"E:/迅雷下载/gene2go.gz"

print("=== gene_info 前3行 ===")
with gzip.open(GI, 'rt') as f:
    for i, l in enumerate(f):
        if i < 3:
            print(repr(l.rstrip()[:400]))
        else:
            break

print("\n=== gene_info 行数(猪基因总数) ===")
n = 0
with gzip.open(GI, 'rt') as f:
    for _ in f:
        n += 1
print("rows:", n)

print("\n=== gene2go 前3行 ===")
with gzip.open(G2G, 'rt') as f:
    for i, l in enumerate(f):
        if i < 3:
            print(repr(l.rstrip()[:400]))
        else:
            break

print("\n=== gene2go 中猪(taxon 9823)行数 + Category 取值样例 ===")
pig_rows = 0
cats = set()
qual_not = 0
sample = []
with gzip.open(G2G, 'rt') as f:
    for i, l in enumerate(f):
        if l.startswith('#'):
            continue
        parts = l.rstrip('\n').split('\t')
        if len(parts) < 8:
            continue
        tax = parts[0]
        if tax != '9823':
            continue
        pig_rows += 1
        cats.add(parts[7])
        if 'NOT' in parts[4]:
            qual_not += 1
        if len(sample) < 3:
            sample.append(parts)
print("pig(9823) rows:", pig_rows)
print("Category values:", cats)
print("rows with NOT qualifier:", qual_not)
for s in sample:
    print("  sample:", s)
