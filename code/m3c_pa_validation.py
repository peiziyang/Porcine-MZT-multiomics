#!/usr/bin/env python3
# M3c validation: PA-downregulated genes -> paternal/imprinted pathway enrichment
import pandas as pd, numpy as np, json, time, sys, subprocess
from scipy.stats import fisher_exact

BASE = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/processed"
OUT  = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/processed/m3c_validation"
import os; os.makedirs(OUT, exist_ok=True)

# ---------- 1. load & subset ----------
df = pd.read_csv(f"{BASE}/m3_ivf_vs_pa_results.csv")
pa_down = df[(df.log2FC_PA_vs_IVF < -1) & (df.padj < 0.05)].copy()
pa_up   = df[(df.log2FC_PA_vs_IVF >  1) & (df.padj < 0.05)].copy()
print(f"PA-down genes (log2FC<-1 & padj<0.05): {len(pa_down)}")
print(f"PA-up   genes (log2FC> 1 & padj<0.05): {len(pa_up)}")

# ---------- 2. map ENSSSCG -> pig symbol + human ortholog symbol ----------
import mygene
mg = mygene.MyGeneInfo()
ids = list(pa_down.gene) + list(pa_up.gene)
# pig symbol
pig = mg.querymany(ids, species="pig", scopes="ensembl.gene", fields="symbol",
                   as_dataframe=True, df_index=True, silent=True)
pig_sym = {}
for g in ids:
    if g in pig.index:
        s = pig.loc[g]
        s = s["symbol"] if not isinstance(s, pd.DataFrame) else s["symbol"].iloc[0]
        if isinstance(s, str): pig_sym[g] = s
print(f"mapped pig symbol: {len(pig_sym)}/{len(ids)}")

# human ortholog via symbol lookup (pig symbols mostly == human ortholog symbol)
syms = sorted(set(pig_sym.values()))
hum = mg.querymany(syms, species="human", scopes="symbol", fields="symbol,entrezgene",
                   as_dataframe=True, df_index=True, silent=True)
ens2human = {}
for g, ps in pig_sym.items():
    if ps in hum.index:
        row = hum.loc[ps]
        row = row if not isinstance(row, pd.DataFrame) else row.iloc[0]
        hs = row.get("symbol")
        if isinstance(hs, str): ens2human[g] = hs
print(f"mapped human ortholog: {len(ens2human)}/{len(ids)}")

pa_down_h = [ens2human[g] for g in pa_down.gene if g in ens2human]
pa_up_h   = [ens2human[g] for g in pa_up.gene   if g in ens2human]
print(f"PA-down human orthologs for Enrichr: {len(pa_down_h)}")
print(f"PA-up   human orthologs for Enrichr: {len(pa_up_h)}")

# map Ensembl->symbol back for reporting
ens2pig = pig_sym

# ---------- 3. Enrichr enrichment ----------
def enrichr(genes, bg="GO_Biological_Process_2023", trials=3):
    for _ in range(trials):
        try:
            import requests
            r = requests.post("https://maayanlab.cloud/Enrichr/addList",
                              files={"list": (None, "\n".join(genes)),
                                     "description": (None, "pa_down")}, timeout=60)
            if r.status_code != 200: continue
            uid = r.json()["userListId"]
            time.sleep(1.5)
            q = requests.get(f"https://maayanlab.cloud/Enrichr/enrich",
                             params={"userListId": uid, "backgroundType": bg}, timeout=60)
            if q.status_code != 200: continue
            res = q.json()[bg]
            out = [{"term": x["term"], "pval": x["pvalue"],
                    "adjp": x["adjustedPvalue"], "genes": x["overlappingGenes"],
                    "n_overlap": len(x["overlappingGenes"])} for x in res]
            out = sorted(out, key=lambda x: x["adjp"])[:15]
            return out
        except Exception as e:
            print("  enrichr err", bg, e); time.sleep(3)
    return []

results = {}
for bg in ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]:
    print(f"Enrichr {bg} ...")
    results[bg] = enrichr(pa_down_h, bg)
    for x in results[bg][:5]:
        print(f"   {x['term'][:60]:60s} adjp={x['adjp']:.2e} n={x['n_overlap']}")

# ---------- 4. curated PEG (paternally expressed imprinted genes) overlap ----------
# Mammalian PEG set (conserved in pig); used as hypothesis-generating validation.
PEG = ["IGF2","MEST","PEG3","ZIM2","PEG10","RTL1","DLK1","SNRPN","NDN","MAGEL2",
       "SGCE","GRB10","PLAGL1","INPP5F","GATM","GPR1","NAP1L5","ZDBF2","SDC1",
       "PON2","FCGRT","KCNK9","PHF17","COPG2","SLC38A4","TRIM50","L3MBTL1",
       "MKRN3","NLRP2","CD81","PEG1","NNAT","OSBPL5","CHMP4B","C2CD4A","C2CD4B"]
MEG = ["H19","CDKN1C","PHLDA2","KCNQ1","GNAS","MEG3","MEG8","SNORD116",
       "ZNF597","RBM8A","MAGI2","PARD6G","FAM50B","SLC22A18","TFPI2"]
peg_set, meg_set = set(PEG), set(MEG)

# background = all tested genes mapped to human
all_h = [ens2human[g] for g in df.gene if g in ens2human]
bg_n = len(set(all_h))
down_set = set(pa_down_h)

def fisher(geneset, label):
    overlap = sorted(down_set & geneset)
    a = len(overlap); b = len(geneset) - a
    c = len(down_set) - a; d = bg_n - a - b - c
    if a == 0:
        return {"label": label, "overlap": [], "a":0,"b":b,"c":c,"d":d,"p":1.0}
    _, p = fisher_exact([[a,b],[c,d]], alternative="greater")
    return {"label": label, "overlap": overlap, "a":a,"b":b,"c":c,"d":d,"p":p}

peg_res = fisher(peg_set, "PEG (paternally expressed)")
meg_res = fisher(meg_set, "MEG (maternally expressed)")
print(f"\nPEG overlap with PA-down: {peg_res['overlap']}  (Fisher p={peg_res['p']:.3e})")
print(f"MEG overlap with PA-down: {meg_res['overlap']}  (Fisher p={meg_res['p']:.3e})")

# ---------- 5. save outputs ----------
# enriched terms CSV
rows = []
for bg, lst in results.items():
    for x in lst:
        rows.append({"database": bg, "term": x["term"], "adj_p": x["adjp"],
                     "raw_p": x["pval"], "n_overlap": x["n_overlap"],
                     "overlap_genes": ";".join(x["genes"])})
pd.DataFrame(rows).to_csv(f"{OUT}/enrichr_pa_down.csv", index=False)

# PEG/MEG summary
pd.DataFrame([{"set": peg_res["label"], "n_overlap": peg_res["a"],
               "set_size": len(peg_set), "bg_tested": bg_n,
               "fisher_p": peg_res["p"], "overlap_genes": ";".join(peg_res["overlap"])},
              {"set": meg_res["label"], "n_overlap": meg_res["a"],
               "set_size": len(meg_set), "bg_tested": bg_n,
               "fisher_p": meg_res["p"], "overlap_genes": ";".join(meg_res["overlap"])}]
            ).to_csv(f"{OUT}/peg_meg_overlap.csv", index=False)

# also save mapped gene lists
pd.DataFrame({"ensembl": pa_down.gene,
              "pig_symbol": [ens2pig.get(g,"") for g in pa_down.gene],
              "human_ortholog": [ens2human.get(g,"") for g in pa_down.gene],
              "log2FC": pa_down.log2FC_PA_vs_IVF,
              "padj": pa_down.padj}).to_csv(f"{OUT}/pa_down_genes_mapped.csv", index=False)

# ---------- 6. figure: top GO/KEGG/Reactome terms ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=(18,7))
titles = {"GO_Biological_Process_2023":"GO Biological Process",
          "KEGG_2021_Human":"KEGG", "Reactome_2022":"Reactome"}
for ax, (bg, lbl) in zip(axes, titles.items()):
    d = results[bg][:10][::-1]
    if not d:
        ax.set_title(f"{lbl}\n(no result)"); continue
    names = [x["term"].split("(")[0][:38] for x in d]
    vals = [-np.log10(max(x["adjp"],1e-300)) for x in d]
    ax.barh(range(len(d)), vals, color="#2c7fb8")
    ax.set_yticks(range(len(d))); ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("-log10 adj p"); ax.set_title(lbl, fontsize=11)
    for i,v in enumerate(vals):
        ax.text(v+0.1, i, f"{v:.1f}", va="center", fontsize=7)
plt.tight_layout()
plt.savefig(f"{OUT}/enrichr_pa_down_terms.png", dpi=130, bbox_inches="tight")
plt.close()

# figure: PEG/MEG overlap bar
fig, ax = plt.subplots(figsize=(6,4))
for i,(res,col) in enumerate([(peg_res,"#d95f02"),(meg_res,"#1b9e77")]):
    ax.bar(i, res["a"], color=col,
           label=f"{res['label']}\noverlap={res['a']}/{len(peg_set) if i==0 else len(meg_set)}\np={res['p']:.1e}")
ax.set_xticks([0,1]); ax.set_xticklabels(["PEG (paternal)","MEG (maternal)"])
ax.set_ylabel("# PA-down genes in set"); ax.set_title("Imprinted-gene overlap with PA-down set")
ax.legend(fontsize=8, loc="upper right")
plt.tight_layout(); plt.savefig(f"{OUT}/peg_overlap.png", dpi=130); plt.close()

print("\nDONE. outputs in", OUT)
print("pa_down mapped file rows:", len(pa_down))
