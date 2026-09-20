#!/usr/bin/env python3
# M3c validation PART 2: mygene-based GO/KEGG enrichment + Ensembl-level PEG/MEG overlap
import pandas as pd, numpy as np, os, sys
from scipy.stats import fisher_exact
import mygene, requests
BASE = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/processed"
OUT  = f"{BASE}/m3c_validation"; os.makedirs(OUT, exist_ok=True)
df = pd.read_csv(f"{BASE}/m3_ivf_vs_pa_results.csv")
pa_down = df[(df.log2FC_PA_vs_IVF < -1) & (df.padj < 0.05)].copy()
pa_up   = df[(df.log2FC_PA_vs_IVF >  1) & (df.padj < 0.05)].copy()
fg = list(pa_down.gene)
bg_all = list(df.gene.unique())
print(f"fg={len(fg)} bg_tested={len(bg_all)}", flush=True)

mg = mygene.MyGeneInfo()

# ---------- (A) GO/KEGG enrichment (local hypergeometric) ----------
def get_anno(ids, fields="go.BP,go.MF,go.CC,pathway.kegg"):
    out = {}
    for i in range(0, len(ids), 400):
        chunk = ids[i:i+400]
        try:
            r = mg.querymany(chunk, species="pig", scopes="ensembl.gene",
                             fields=fields, as_dataframe=True, df_index=True, silent=True)
            for g in chunk:
                if g in r.index:
                    row = r.loc[g]
                    row = row if not isinstance(row, pd.DataFrame) else row.iloc[0]
                    out[g] = row
        except Exception as e:
            print("anno err", e, flush=True)
    return out

print("annotating foreground...", flush=True)
fg_anno = get_anno(fg)
print("annotating background (chunked)...", flush=True)
bg_anno = get_anno(bg_all)

def collect_terms(anno, kind):
    # kind: 'go.BP' or 'kegg'
    term2genes = {}
    for g, row in anno.items():
        if kind.startswith("go"):
            go = row.get("go") if hasattr(row, "get") else None
            if isinstance(go, dict):
                bp = go.get({"go.BP":"BP","go.MF":"MF","go.CC":"CC"}[kind])
                if isinstance(bp, list):
                    for t in bp:
                        if isinstance(t, dict) and "term" in t and "id" in t:
                            term2genes.setdefault((t["id"], t["term"]), set()).add(g)
        else:
            pw = row.get("pathway") if hasattr(row, "get") else None
            if isinstance(pw, dict):
                kegg = pw.get("kegg")
                if isinstance(kegg, list):
                    for t in kegg:
                        if isinstance(t, dict) and "name" in t and "id" in t:
                            term2genes.setdefault((t["id"], t["name"]), set()).add(g)
    return term2genes

bg_n = len(bg_anno)
rows = []
for kind,label in [("go.BP","GO_BP"),("go.MF","GO_MF"),("go.CC","GO_CC"),("kegg","KEGG")]:
    t2g = collect_terms(bg_anno, kind)
    fg_terms = collect_terms(fg_anno, kind)
    for (tid, tname), bset in t2g.items():
        fset = fg_terms.get((tid, tname), set())
        a = len(fset & set(fg)); 
        if a < 2: continue
        b = len(bset) - a
        c = len(fg) - a
        d = bg_n - a - b - c
        if a == 0: continue
        _, p = fisher_exact([[a,b],[c,d]], alternative="greater")
        rows.append({"type":label,"term_id":tid,"term":tname,
                     "fg_overlap":a,"bg_total":len(bset),"raw_p":p,
                     "genes":";".join(sorted(fset)[:15])})
res = pd.DataFrame(rows)
if len(res):
    # BH correction per type
    res["adj_p"] = res.groupby("type")["raw_p"].transform(
        lambda x: np.minimum(1, x.values * len(x) / np.argsort(np.argsort(x.values)+0.5)))
    res = res.sort_values("adj_p").reset_index(drop=True)
    res.to_csv(f"{OUT}/go_kegg_enrichment.csv", index=False)
    print("Top GO/KEGG:", flush=True)
    for _,r in res.head(12).iterrows():
        print(f"  [{r['type']}] {r['term'][:50]} adjp={r['adj_p']:.2e} n={r['fg_overlap']}", flush=True)
else:
    print("NO enrichment terms found", flush=True)

# ---------- (B) PEG/MEG overlap on Ensembl IDs ----------
PEG = ["IGF2","MEST","PEG3","ZIM2","PEG10","RTL1","DLK1","SNRPN","NDN","MAGEL2",
       "SGCE","GRB10","PLAGL1","INPP5F","GATM","GPR1","NAP1L5","ZDBF2","SDC1",
       "PON2","FCGRT","KCNK9","PHF17","COPG2","SLC38A4","TRIM50","L3MBTL1",
       "MKRN3","NLRP2","CD81","NNAT","OSBPL5","CHMP4B","C2CD4A","C2CD4B","PEG1",
       "ZRSR1","APOBEC3B","DDC","ANTXR1","JKAMP","NDNL2","PWAR1","PWAR2","GPR1",
       "SUSD2","H13","YBEY","C2CD4C","INPP5F","BLCAP","MAGI2","SLC22A18","TFPI2"]
MEG = ["H19","CDKN1C","PHLDA2","KCNQ1","GNAS","MEG3","MEG8","SNORD116",
       "ZNF597","RBM8A","PARD6G","FAM50B","PEG13","RTL1","DLK1","IGF2","GRB10",
       "ZAC1","MEST","GATM","NAP1L5","NNAT","PEG10","SGCE","MAGEL2","NDN","SNRPN"]
def to_ensembl(syms):
    s = list(set(syms))
    out = {}
    for i in range(0,len(s),200):
        try:
            r = mg.querymany(s[i:i+200], species="pig", scopes="symbol",
                             fields="ensembl.gene", as_dataframe=True, df_index=True, silent=True)
            for sy in s[i:i+200]:
                if sy in r.index:
                    row=r.loc[sy]; row=row if not isinstance(row,pd.DataFrame) else row.iloc[0]
                    eg=row.get("ensembl")
                    if isinstance(eg,list): eg=[e.get("gene") for e in eg if isinstance(e,dict)]
                    elif isinstance(eg,dict): eg=[eg.get("gene")]
                    else: eg=[]
                    if eg: out[sy]=eg
        except Exception as e:
            print("peg map err",e,flush=True)
    return out
print("mapping PEG/MEG to pig Ensembl...", flush=True)
peg_e = to_ensembl(PEG); meg_e = to_ensembl(MEG)
peg_set=set(x for v in peg_e.values() for x in v)
meg_set=set(x for v in meg_e.values() for x in v)
fg_set=set(fg)
def ftest(s,label):
    ov=sorted(fg_set & s); a=len(ov); b=len(s)-a; c=len(fg)-a; d=bg_n-a-b-c
    p=fisher_exact([[a,b],[c,d]],alternative="greater")[1] if a>0 else 1.0
    return {"set":label,"overlap":ov,"a":a,"set_size":len(s),"bg":bg_n,"p":p}
pr=ftest(peg_set,"PEG"); mr=ftest(meg_set,"MEG")
print(f"PEG Ensembl overlap with PA-down: {pr['overlap']} p={pr['p']:.3e}", flush=True)
print(f"MEG Ensembl overlap with PA-down: {mr['overlap']} p={mr['p']:.3e}", flush=True)
pd.DataFrame([{"set":pr["set"],"n_overlap":pr["a"],"set_size":pr["set_size"],
               "bg":pr["bg"],"fisher_p":pr["p"],"overlap_genes":";".join(pr["overlap"])},
              {"set":mr["set"],"n_overlap":mr["a"],"set_size":mr["set_size"],
               "bg":mr["bg"],"fisher_p":mr["p"],"overlap_genes":";".join(mr["overlap"])}]
            ).to_csv(f"{OUT}/peg_meg_overlap_enriched.csv", index=False)
print("PART2 DONE", flush=True)
