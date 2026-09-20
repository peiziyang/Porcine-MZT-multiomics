#!/usr/bin/env python3
# M3c validation PART 3: Reactome pathway enrichment via mygene (pig has reactome, no GO/kegg)
import pandas as pd, numpy as np, os
from scipy.stats import fisher_exact
import mygene
BASE = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/processed"
OUT  = f"{BASE}/m3c_validation"; os.makedirs(OUT, exist_ok=True)
df = pd.read_csv(f"{BASE}/m3_ivf_vs_pa_results.csv")
pa_down = df[(df.log2FC_PA_vs_IVF < -1) & (df.padj < 0.05)].copy()
fg = list(pa_down.gene); bg_all = list(df.gene.unique())
print(f"fg={len(fg)} bg={len(bg_all)}", flush=True)
mg = mygene.MyGeneInfo()

def get_path(ids):
    out = {}
    for i in range(0, len(ids), 400):
        try:
            r = mg.querymany(ids[i:i+400], species="pig", scopes="ensembl.gene",
                             fields="pathway.reactome", as_dataframe=True, df_index=True, silent=True)
            for g in ids[i:i+400]:
                if g in r.index:
                    row = r.loc[g]; row = row if not isinstance(row, pd.DataFrame) else row.iloc[0]
                    out[g] = row
        except Exception as e:
            print("err", e, flush=True)
    return out

print("fg pathway...", flush=True); fg_a = get_path(fg)
print("bg pathway...", flush=True); bg_a = get_path(bg_all)

def reactome_terms(anno):
    t2g = {}
    for g, row in anno.items():
        lst = row.get("pathway.reactome") if hasattr(row, "get") else None
        if isinstance(lst, list):
            for t in lst:
                if isinstance(t, dict) and "name" in t and "id" in t:
                    t2g.setdefault((t["id"], t["name"]), set()).add(g)
    return t2g

bg_n = len(bg_a)
bg_terms = reactome_terms(bg_a)
fg_terms = reactome_terms(fg_a)
rows = []
for (tid, tname), bset in bg_terms.items():
    fset = fg_terms.get((tid, tname), set())
    a = len(fset & set(fg))
    if a < 2: continue
    b = len(bset) - a; c = len(fg) - a; d = bg_n - a - b - c
    _, p = fisher_exact([[a,b],[c,d]], alternative="greater")
    rows.append({"term_id":tid,"term":tname,"fg_overlap":a,"bg_total":len(bset),
                 "raw_p":p,"genes":";".join(sorted(fset)[:20])})
res = pd.DataFrame(rows)
if len(res):
    res["adj_p"] = np.minimum(1, res["raw_p"].values * len(res) / (np.argsort(np.argsort(res["raw_p"].values))+1))
    res = res.sort_values("adj_p").reset_index(drop=True)
    res.to_csv(f"{OUT}/reactome_enrichment.csv", index=False)
    print("Top Reactome pathways in PA-down set:", flush=True)
    for _,r in res.head(15).iterrows():
        print(f"  {r['term'][:55]:55s} adjp={r['adj_p']:.2e} n={r['fg_overlap']} bg={r['bg_total']}", flush=True)
    # figure
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    top = res.head(12)[::-1]
    plt.figure(figsize=(9,6))
    plt.barh(range(len(top)), -np.log10(top.adj_p.clip(lower=1e-300)), color="#2c7fb8")
    plt.yticks(range(len(top)), [t[:48] for t in top.term], fontsize=8)
    plt.xlabel("-log10 adj p"); plt.title("Reactome enrichment: PA-downregulated genes (IVF vs PA)")
    for i,v in enumerate(-np.log10(top.adj_p.clip(lower=1e-300))):
        plt.text(v+0.05, i, f"{v:.1f}", va="center", fontsize=7)
    plt.tight_layout(); plt.savefig(f"{OUT}/reactome_enrichment.png", dpi=130, bbox_inches="tight"); plt.close()
else:
    print("NO Reactome terms", flush=True)
print("PART3 DONE", flush=True)
