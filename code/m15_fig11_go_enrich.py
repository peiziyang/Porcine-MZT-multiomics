"""
M15 Fig 11: Real GO enrichment of stage-conserved DEGs (IVF vs PA, pig cleavage).
------------------------------------------------------------------------------
Offline, no network. Uses user-downloaded NCBI files:
  - E:/迅雷下载/Sus_scrofa.gene_info.gz  (Ensembl -> Entrez/Symbol, pig only)
  - E:/迅雷下载/gene2go.gz              (Entrez -> GO, filter taxon 9823 = pig)

Method:
  * Recompute 321 conserved DEGs (Ensembl-id indexed) from m15 4-stage CSVs,
    using the same rule as Fig 10 but keeping Ensembl IDs.
  * Background universe = all genes detected in m15 (union of 4 stages, Ensembl),
    mapped to Entrez, restricted to those with >=1 GO annotation.
  * Hypergeometric test per GO term (scipy.stats.hypergeom), BH FDR correction.
  * Bar plot of top GO terms per domain (BP/MF/CC), annotated with -log10(q).
All numbers are computed, not assumed.
"""
import os, gzip, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
try:
    from scipy.stats import false_discovery_control
    HAVE_FDC = True
except Exception:
    HAVE_FDC = False

BASE = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed/deseq2'
GI   = r'E:/迅雷下载/Sus_scrofa.gene_info.gz'
G2G  = r'E:/迅雷下载/gene2go.gz'
OUT_PNG  = os.path.join(BASE, 'm15_fig11_go_enrich.png')
OUT_CSV  = os.path.join(BASE, 'm15_fig11_go_enrich.csv')
OUT_MAP  = os.path.join(BASE, 'm15_conserved_deg_ensembl.tsv')
PIG_TAX = '9823'

# ---------------------------------------------------------------- 1. conserved DEG (Ensembl indexed)
stages = ['1cell', '2cell', '4cell', '8cell']
TH_FC, TH_P = 1.0, 0.05
lfc, padj = {}, {}
for s in stages:
    df = pd.read_csv(f'{BASE}/m15_deseq2_{s}_IVFvsPA.csv')
    df = df.sort_values('padj').drop_duplicates('gene', keep='first')   # Ensembl-id indexed
    lfc[s]  = df.set_index('gene')['log2FoldChange']
    padj[s] = df.set_index('gene')['padj']
universe_genes = set()
for s in stages:
    universe_genes |= set(lfc[s].index)
lfc  = pd.DataFrame(lfc)
padj = pd.DataFrame(padj)
lfc  = lfc.loc[lfc.index.intersection(padj.index)]
padj = padj.loc[lfc.index]
sig   = (padj < TH_P) & (lfc.abs() > TH_FC)
n_sig = sig.sum(axis=1)
pos   = (lfc > 0).sum(axis=1)
neg   = (lfc < 0).sum(axis=1)
consistent = (pos >= 3) | (neg >= 3)
mask = (n_sig >= 2) & consistent
conserved = sorted(lfc.index[mask])
print(f'[1] detected genes (universe candidates): {len(universe_genes)}')
print(f'[1] conserved DEGs (Ensembl): {len(conserved)}')

# ---------------------------------------------------------------- 2. parse annotations
print('[2] parsing gene_info (dbXrefs Ensembl -> Entrez, Symbol bridge)...')
# (a) Ensembl->symbol from our own m15 stage CSVs (GTF-consistent)
ensembl2sym = {}
for s in stages:
    df = pd.read_csv(f'{BASE}/m15_deseq2_{s}_IVFvsPA.csv', usecols=['gene', 'gene_symbol'])
    for g, sym in zip(df['gene'].astype(str), df['gene_symbol'].astype(str)):
        if g not in ensembl2sym and sym and sym != 'nan':
            ensembl2sym[g] = sym
# (b) gene_info: Symbol->Entrez + Ensembl(dbXrefs)->Entrez
symbol2entrez = {}
ensembl2entrez_db = {}
with gzip.open(GI, 'rt') as f:
    header = f.readline()
    cols = header.lstrip('#').rstrip('\n').split('\t')
    i_gid = cols.index('GeneID')
    i_sym = cols.index('Symbol')
    i_db  = cols.index('dbXrefs') if 'dbXrefs' in cols else None
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) <= i_gid:
            continue
        gid = p[i_gid]
        sym = p[i_sym] if i_sym < len(p) else ''
        if sym:
            symbol2entrez[sym] = gid
        if i_db is not None and i_db < len(p):
            m = re.search(r'Ensembl:([A-Za-z0-9_.]+)', p[i_db])
            if m:
                en = m.group(1).split('.')[0]
                if en.startswith('ENSSSCG'):
                    ensembl2entrez_db[en] = gid
# (c) merge: prefer direct Ensembl(dbXrefs), fallback to symbol bridge
ensembl2entrez = dict(ensembl2entrez_db)
fb = 0
for g, sym in ensembl2sym.items():
    if g not in ensembl2entrez and sym in symbol2entrez:
        ensembl2entrez[g] = symbol2entrez[sym]
        fb += 1
print(f'    direct Ensembl(dbXrefs)->Entrez: {len(ensembl2entrez_db)}; via symbol bridge: +{fb}')
print(f'    total Ensembl->Entrez: {len(ensembl2entrez)}; Ensembl->symbol: {len(ensembl2sym)}')

print('[2] parsing gene2go (pig subset, Entrez->GO)...')
entrez2go = {}
go_info = {}   # go_id -> (category, term)
with gzip.open(G2G, 'rt') as f:
    header = f.readline()
    cols = header.lstrip('#').rstrip('\n').split('\t') if header.startswith('#') else header.rstrip('\n').split('\t')
    # find by name
    def idx(name, alts):
        for nm in [name] + alts:
            if nm in cols:
                return cols.index(nm)
        return None
    i_tax = idx('tax_id', ['tax_id'])
    i_gid = idx('GeneID', ['GeneID'])
    i_go  = idx('GO_ID', ['GO_ID', 'go_id'])
    i_qual= idx('Qualifier', ['Qualifier'])
    i_term= idx('GO_term', ['GO_term'])
    i_cat = idx('Category', ['Category'])
    for line in f:
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) <= max(i_tax, i_gid, i_go):
            continue
        if p[i_tax] != PIG_TAX:
            continue
        if i_qual is not None and 'NOT' in p[i_qual]:
            continue
        gid = p[i_gid]; go = p[i_go]
        entrez2go.setdefault(gid, set()).add(go)
        if go not in go_info and i_term is not None and i_cat is not None:
            go_info[go] = (p[i_cat], p[i_term])
print(f'    pig genes with GO: {len(entrez2go)}')
print(f'    GO terms observed: {len(go_info)}')

# INVERT: go -> gene set for enrichment (entrez2go above is gene->GO)
go2genes = {}
for _gid, _gos in entrez2go.items():
    for _go in _gos:
        go2genes.setdefault(_go, set()).add(_gid)

# ---------------------------------------------------------------- 3. build foreground/background Entrez sets
def to_entrez(ens_set):
    out = set()
    miss = 0
    for g in ens_set:
        e = ensembl2entrez.get(g)
        if e:
            out.add(e)
        else:
            miss += 1
    return out, miss

bg_entrez, bg_miss = to_entrez(universe_genes)
deg_entrez, deg_miss = to_entrez(set(conserved))
# restrict to those with GO annotation
bg_with_go  = bg_entrez  & set(entrez2go.keys())
deg_with_go = deg_entrez & set(entrez2go.keys())
# ensure deg is a subset of bg so N <= M (otherwise hypergeom.sf -> NaN)
deg_with_go = deg_with_go & bg_with_go
M = len(bg_with_go)
N = len(deg_with_go)
print(f'[3] background Entrez total={len(bg_entrez)} (unmapped={bg_miss}); with GO: M={M}')
print(f'[3] conserved  Entrez total={len(deg_entrez)} (unmapped={deg_miss}); with GO: N={N}')

# ---------------------------------------------------------------- 4. hypergeometric enrichment
rows = []
for go, term_genes in go2genes.items():
    n_i = len(term_genes & bg_with_go)          # background annotated
    k_i = len(term_genes & deg_with_go)         # DEG annotated
    if k_i == 0:
        continue
    # p = P(X >= k_i)
    p = hypergeom.sf(k_i - 1, M, n_i, N)
    cat, term = go_info.get(go, ('?', go))
    rows.append((go, cat, term, n_i, k_i, k_i / N, p))
res = pd.DataFrame(rows, columns=['GO', 'domain', 'term', 'bg_n', 'deg_n', 'frac_in_deg', 'p'])
# coerce to numeric, then drop non-finite p-values (e.g. hypergeom edge cases) before FDR
res['p'] = pd.to_numeric(res['p'], errors='coerce')
res = res[np.isfinite(res['p'])].copy()
if HAVE_FDC:
    res['q'] = false_discovery_control(res['p'].values, method='bh')
else:
    # manual BH
    m = len(res)
    order = res['p'].argsort().values
    q = np.empty(m)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        val = res['p'].iloc[i] * m / (rank + 1)
        prev = min(prev, val)
        q[i] = prev
    res['q'] = q.clip(upper=1.0)
res['log10q'] = -np.log10(res['q'].clip(lower=1e-300))
res = res.sort_values('q').reset_index(drop=True)
res.to_csv(OUT_CSV, index=False)
print(f'[4] GO terms tested: {len(res)}; significant (q<0.05): {(res["q"]<0.05).sum()}')

# save conserved DEG mapping for reproducibility
cmap = pd.DataFrame({
    'gene': conserved,
    'entrez': [ensembl2entrez.get(g, '') for g in conserved],
    'symbol': [ensembl2sym.get(g, '') for g in conserved],
})
cmap['mean_log2FC'] = [float(lfc.loc[g].mean()) for g in conserved]
cmap['n_sig_stages'] = [int(n_sig.loc[g]) for g in conserved]
cmap['direction'] = np.where(cmap['mean_log2FC'] > 0, 'PA-up', 'PA-down')
cmap.to_csv(OUT_MAP, sep='\t', index=False)
print(f'[4] saved conserved DEG map -> {OUT_MAP}')

# ---------------------------------------------------------------- 5. plot Fig 11
domain_label = {'biological_process': 'Biological Process',
                'molecular_function': 'Molecular Function',
                'cellular_component': 'Cellular Component', '?': 'Other'}
fig, axes = plt.subplots(1, 3, figsize=(18, 9))
for ax, dom in zip(axes, ['biological_process', 'molecular_function', 'cellular_component']):
    sub = res[res['domain'] == dom].sort_values('q').head(12)
    if len(sub) == 0:
        ax.set_title(domain_label.get(dom, dom) + ' (none significant)')
        ax.axis('off'); continue
    sub = sub.iloc[::-1]
    labels = [f'{t[:38]}' for t in sub['term']]
    vals = sub['log10q'].values
    y = np.arange(len(sub))
    colors = ['#c0392b' if v >= 1.3 else '#e67e22' for v in vals]
    ax.barh(y, vals, color=colors)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('-log10(FDR q)', fontsize=10)
    ax.set_title(f'{domain_label.get(dom, dom)}\n(n={len(sub)} shown, deg_n label)', fontsize=11)
    for i, (v, k, n) in enumerate(zip(vals, sub['deg_n'], sub['bg_n'])):
        ax.text(v + 0.05, i, f'{int(k)}/{int(n)}', va='center', fontsize=7.5, color='#333')
    ax.axvline(1.3, color='grey', ls='--', alpha=0.5)
fig.suptitle('Figure 11. GO enrichment of stage-conserved DEGs (IVF vs PA, pig cleavage)\n'
             f'background M={M} genes, DEG N={N} genes; hypergeometric + BH', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_PNG, dpi=150, bbox_inches='tight')
print(f'[5] saved {OUT_PNG}')

# print top terms per domain
print('\n=== TOP BP ===')
print(res[res['domain']=='biological_process'].head(8)[['term','deg_n','bg_n','q']].to_string(index=False))
print('\n=== TOP MF ===')
print(res[res['domain']=='molecular_function'].head(8)[['term','deg_n','bg_n','q']].to_string(index=False))
print('\n=== TOP CC ===')
print(res[res['domain']=='cellular_component'].head(8)[['term','deg_n','bg_n','q']].to_string(index=False))
print('\nDONE')
