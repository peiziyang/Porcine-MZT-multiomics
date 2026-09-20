#!/usr/bin/env python3
# M3c_validation LOCAL part: GTF-based symbol mapping + PEG/MEG Fisher overlap
import pandas as pd, numpy as np, re, os, gzip
from scipy.stats import fisher_exact
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/processed"
REF  = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/reference"
OUT  = f"{BASE}/m3c_validation"; os.makedirs(OUT, exist_ok=True)

# ---------- 1. GTF -> ENSSSCG -> pig symbol ----------
print("Parsing GTF for gene_id->gene_name ...")
g2s = {}
with gzip.open(f"{REF}/Sus_scrofa.Sscrofa11.1.113.gtf.gz", "rt") as f:
    for line in f:
        if not line.startswith("#") and "\tgene\t" in line:
            gid = re.search(r'gene_id "([^"]+)"', line)
            gnm = re.search(r'gene_name "([^"]+)"', line)
            if gid and gnm:
                g2s[gid.group(1)] = gnm.group(1)
print(f"  GTF symbols loaded: {len(g2s)}")

# ---------- 2. load & subset ----------
df = pd.read_csv(f"{BASE}/m3_ivf_vs_pa_results.csv")
pa_down = df[(df.log2FC_PA_vs_IVF < -1) & (df.padj < 0.05)].copy()
pa_up   = df[(df.log2FC_PA_vs_IVF >  1) & (df.padj < 0.05)].copy()
print(f"PA-down: {len(pa_down)} | PA-up: {len(pa_up)}")

def sym(g): return g2s.get(g, "")
pa_down["pig_symbol"] = pa_down.gene.map(sym)
pa_up["pig_symbol"]   = pa_up.gene.map(sym)
down_sym = set(s for s in pa_down.pig_symbol if s)
up_sym   = set(s for s in pa_up.pig_symbol if s)
all_sym  = set(s for s in df.gene.map(sym) if s)
print(f"down symbols: {len(down_sym)} | background tested symbols: {len(all_sym)}")
print("  sample down symbols:", sorted(down_sym)[:12])

# ---------- 3. PEG / MEG overlap (Fisher exact) ----------
# Mammalian paternally(paternal)/maternally expressed imprinted genes; pig uses same symbols
PEG = ["IGF2","MEST","PEG3","ZIM2","PEG10","RTL1","DLK1","SNRPN","NDN","MAGEL2",
       "SGCE","GRB10","PLAGL1","INPP5F","GATM","GPR1","NAP1L5","ZDBF2","SDC1",
       "PON2","FCGRT","KCNK9","PHF17","COPG2","SLC38A4","TRIM50","L3MBTL1",
       "MKRN3","NLRP2","CD81","NNAT","OSBPL5","CHMP4B","C2CD4A","C2CD4B","PEG1"]
MEG = ["H19","CDKN1C","PHLDA2","KCNQ1","GNAS","MEG3","MEG8","SNORD116",
       "ZNF597","RBM8A","MAGI2","PARD6G","FAM50B","SLC22A18","TFPI2","PEG13"]
peg_set, meg_set = set(PEG), set(MEG)
bg = len(all_sym)

def fisher(geneset, label):
    ov = sorted(down_sym & geneset)
    a = len(ov); b = len(geneset) - a
    c = len(down_sym) - a; d = bg - a - b - c
    p = fisher_exact([[a,b],[c,d]], alternative="greater")[1] if a>0 else 1.0
    return {"label":label,"overlap":ov,"a":a,"b":b,"c":c,"d":d,"p":p,"bg":bg}

peg_res = fisher(peg_set, "PEG (paternally expressed imprinted)")
meg_res = fisher(meg_set, "MEG (maternally expressed imprinted)")
print(f"\nPEG overlap with PA-down: {peg_res['overlap']}  (Fisher p={peg_res['p']:.3e})")
print(f"MEG overlap with PA-down: {meg_res['overlap']}  (Fisher p={meg_res['p']:.3e})")

# ---------- 4. save ----------
pd.DataFrame({"ensembl":pa_down.gene,"pig_symbol":pa_down.pig_symbol,
              "log2FC":pa_down.log2FC_PA_vs_IVF,"padj":pa_down.padj}
            ).to_csv(f"{OUT}/pa_down_genes_mapped.csv", index=False)
pd.DataFrame([{"set":peg_res["label"],"n_overlap":peg_res["a"],"set_size":len(peg_set),
               "bg_tested":bg,"fisher_p":peg_res["p"],
               "overlap_genes":";".join(peg_res["overlap"])},
              {"set":meg_res["label"],"n_overlap":meg_res["a"],"set_size":len(meg_set),
               "bg_tested":bg,"fisher_p":meg_res["p"],
               "overlap_genes":";".join(meg_res["overlap"])}]
            ).to_csv(f"{OUT}/peg_meg_overlap.csv", index=False)

# ---------- 5. figure: PEG/MEG overlap ----------
fig, ax = plt.subplots(figsize=(6.5,4))
bars = [peg_res["a"], meg_res["a"]]
labels = [f"PEG (paternal)\nn={len(peg_set)}", f"MEG (maternal)\nn={len(meg_set)}"]
cols = ["#d95f02","#1b9e77"]
bars2 = ax.bar(range(2), bars, color=cols)
ax.set_xticks(range(2)); ax.set_xticklabels(labels)
ax.set_ylabel("# PA-down genes overlapping")
ax.set_title("Imprinted-gene overlap with PA-downregulated set")
for i,b in enumerate(bars):
    ax.text(i, b+0.1, f"{b}\np={[peg_res,meg_res][i]['p']:.1e}",
            ha="center", va="bottom", fontsize=9)
ax.set_ylim(0, max(bars)+1)
plt.tight_layout(); plt.savefig(f"{OUT}/peg_overlap.png", dpi=130); plt.close()
print("Wrote:", f"{OUT}/peg_overlap.png", f"{OUT}/peg_meg_overlap.csv",
      f"{OUT}/pa_down_genes_mapped.csv")
print("LOCAL PART DONE")
