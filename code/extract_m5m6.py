import csv, json, statistics, subprocess, sys
from collections import defaultdict

BASE = r"E:\Workbuddy\2026-07-27-11-58-27\data\processed"
out = {}

# ---------- M5: top50 divergent genes ----------
p = f"{BASE}/m5_top50_divergent_genes.csv"
rows = []
with open(p, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        d = {k: float(v) for k, v in row.items() if k != "gene"}
        d["gene"] = row["gene"]
        rows.append(d)
cols = ["IVF_vs_PA", "IVF_vs_vivo", "PA_vs_vivo"]
means = {c: statistics.mean(x[c] for x in rows) for c in cols}
out["m5_top50_means"] = {c: round(means[c], 3) for c in cols}
# which comparison has the largest divergence per gene?
def argmax(x):
    return max(cols, key=lambda c: x[c])
cnt = defaultdict(int)
for x in rows:
    cnt[argmax(x)] += 1
out["m5_top50_argmax_counts"] = dict(cnt)
out["m5_top50_n"] = len(rows)
# top gene by PA_vs_vivo
top_pavivo = sorted(rows, key=lambda x: -x["PA_vs_vivo"])[:5]
out["m5_top_pavivo_genes"] = [(x["gene"], round(x["PA_vs_vivo"],3)) for x in top_pavivo]

# ---------- M5 v2: global pattern ----------
p2 = f"{BASE}/m5_gene_divergence_v2.csv"
allrows = []
with open(p2, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        try:
            allrows.append({k: float(v) for k, v in row.items() if k != "gene"})
        except:
            pass
n_all = len(allrows)
pa_gt_iv = sum(1 for x in allrows if x["PA_vs_vivo"] > x["IVF_vs_vivo"])
pa_gt_pa = sum(1 for x in allrows if x["PA_vs_vivo"] > x["IVF_vs_PA"])
out["m5_v2_n_genes"] = n_all
out["m5_v2_pa_gt_ivffrac"] = round(pa_gt_iv / n_all, 4)
out["m5_v2_pa_gt_ivfpairfrac"] = round(pa_gt_pa / n_all, 4)
out["m5_v2_totaldiv_mean"] = round(statistics.mean(x["total_div"] for x in allrows), 4)

# ---------- M6: IVF vs PA per-cell metabolism ----------
p6 = f"{BASE}/m6_metabolism_paivf.csv"
ivf, pa = [], []
with open(p6, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        if row["condition"] not in ("IVF", "PA"):
            continue
        try:
            g = float(row["glycolysis"]); o = float(row["oxphos"])
        except:
            continue
        (ivf if row["condition"] == "IVF" else pa).append((g, o))
def mm(lst):
    g = statistics.mean(x[0] for x in lst); o = statistics.mean(x[1] for x in lst)
    return round(g,3), round(o,3), round(g/o,3), len(lst)
out["m6_ivf"] = mm(ivf)
out["m6_pa"] = mm(pa)

# ---------- M6: per-stage ratio ----------
p6r = f"{BASE}/m6_metabolism_ratio.csv"
rrows = []
with open(p6r, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        rrows.append((row["stage"], float(row["gly_mean"]), float(row["ox_mean"]), float(row["ratio"]), int(row["n"])))
ratios = [x[3] for x in rrows]
out["m6_ratio_all_lt_06"] = all(r < 0.6 for r in ratios)
out["m6_ratio_min"] = (min(rrows, key=lambda x: x[3])[0], round(min(ratios),3))
out["m6_ratio_max"] = (max(rrows, key=lambda x: x[3])[0], round(max(ratios),3))
out["m6_ratio_mean"] = round(statistics.mean(ratios),3)
out["m6_ratio_range"] = (round(min(ratios),3), round(max(ratios),3))

print(json.dumps(out, indent=2))

# ---------- map top genes to symbols via mygene (pig) ----------
top_ids = [g for g, _ in out["m5_top_pavivo_genes"]]
try:
    import mygene
    mg = mygene.MyGeneInfo()
    res = mg.querymany(top_ids, species="pig", scopes="ensembl.gene", fields="symbol", size=1)
    sym = {}
    for x in res:
        sym[x.get("query", "")] = x.get("symbol", "?")
    print("GENE_SYMBOLS:", json.dumps(sym, indent=2))
except Exception as e:
    print("MYGENE_FAIL:", repr(e))
