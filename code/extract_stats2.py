import csv, os, statistics, json, ast
from collections import defaultdict, Counter

BASE = r"E:\Workbuddy\2026-07-27-11-58-27\data\processed"

# ---- regulons_ctx.csv: 3-row header ----
# row0: ['','',Enrichment x8]; row1: ['','',AUC,NES,...,TargetGenes,RankAtMax]; row2: [TF,MotifID,'','',''...]
rp = os.path.join(BASE, "pyscenic_out", "regulons_ctx.csv")
sizes = defaultdict(int)
nrows = 0
with open(rp, newline="") as f:
    r = csv.reader(f)
    next(r); next(r); next(r)   # skip 3 header rows
    for row in r:
        if not row or len(row) < 9:
            continue
        nrows += 1
        tf = row[0]
        try:
            tg = ast.literal_eval(row[8]) if row[8].strip() else []
            cnt = len(tg) if isinstance(tg, (list, tuple)) else 0
        except Exception:
            cnt = 0
        sizes[tf] += cnt
out = {}
out["m7_regulon_rows"] = nrows
out["m7_n_tf"] = len(sizes)
out["m7_target_min"] = min(sizes.values())
out["m7_target_max"] = max(sizes.values())
out["m7_target_mean"] = round(statistics.mean(sizes.values()), 1)
out["m7_top_tf_by_targets"] = sorted(sizes.items(), key=lambda x: -x[1])[:8]

# ---- aucell + stage_label ----
amap = os.path.join(BASE, "pyscenic_out", "aucell_scores_full.csv")
smap = os.path.join(BASE, "m7_cell_stage_map.csv")
stage = {}
with open(smap, newline="") as f:
    r = csv.reader(f)
    h = next(r)
    ci = h.index("cell_id")
    si = h.index("stage_label")
    for row in r:
        if not row:
            continue
        stage[row[ci].strip()] = row[si].strip()
out["m7_stage_label_counts"] = dict(Counter(stage.values()))

rows = []
with open(amap, newline="") as f:
    r = csv.reader(f)
    h = next(r)
    regcols = h[1:]
    for row in r:
        if not row:
            continue
        rows.append((row[0].strip(), [float(x) for x in row[1:]]))

mean_auc = [statistics.mean(v[i] for _, v in rows) for i in range(len(regcols))]
top_idx = sorted(range(len(regcols)), key=lambda i: -mean_auc[i])[:10]
out["m7_top_regulons_meanAUC"] = [(regcols[i], round(mean_auc[i], 4)) for i in top_idx]

# stage-resolved profile for top 10 regulons
profile = {}
for i in top_idx:
    name = regcols[i]
    d = defaultdict(list)
    for cell, vals in rows:
        st = stage.get(cell, "NA")
        d[st].append(vals[i])
    profile[name] = {st: round(statistics.mean(lst), 4) for st, lst in d.items() if lst}
out["m7_top10_stage_profile"] = profile

# also dataset-level for sanity
dmap = {}
with open(smap, newline="") as f:
    r = csv.reader(f); h = next(r); ci=h.index("cell_id"); di=h.index("dataset")
    for row in r:
        if row: dmap[row[ci].strip()] = row[di].strip()
dset_prof = {}
for i in top_idx[:6]:
    name = regcols[i]
    d = defaultdict(list)
    for cell, vals in rows:
        d[dmap.get(cell,"NA")].append(vals[i])
    dset_prof[name] = {ds: round(statistics.mean(lst),4) for ds,lst in d.items() if lst}
out["m7_top6_dataset_profile"] = dset_prof

print(json.dumps(out, ensure_ascii=False, indent=2))
