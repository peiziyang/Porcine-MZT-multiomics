import csv, os, statistics, json
from collections import defaultdict, Counter

BASE = r"E:\Workbuddy\2026-07-27-11-58-27\data\processed"
out = {}

# ---------- M7: regulons_ctx.csv ----------
rp = os.path.join(BASE, "pyscenic_out", "regulons_ctx.csv")
reg = defaultdict(list)  # tf -> list of (sign, n_targets, auc)
with open(rp, newline="") as f:
    r = csv.reader(f)
    header = next(r)
    idx_tf = header.index("TF") if "TF" in header else 0
    # find TargetGenes column
    try:
        idx_tg = header.index("TargetGenes")
    except ValueError:
        idx_tg = len(header) - 1
    nrows = 0
    for row in r:
        if not row or len(row) <= idx_tf:
            continue
        nrows += 1
        tf = row[idx_tf]
        tg = row[idx_tg] if idx_tg < len(row) else ""
        nt = len([x for x in tg.split(";") if x.strip()]) if tg else 0
        reg[tf].append(nt)
out["m7_regulon_rows"] = nrows
out["m7_n_tf"] = len(reg)
out["m7_regulon_target_sizes"] = {tf: sizes for tf, sizes in reg.items()}
out["m7_regulons_with_neg"] = sorted([tf for tf, s in reg.items() if len(s) > 1])
out["m7_total_regulon_units"] = sum(len(s) for s in reg.values())

# ---------- M7: aucell + stage ----------
amap = os.path.join(BASE, "pyscenic_out", "aucell_scores_full.csv")
smap = os.path.join(BASE, "m7_cell_stage_map.csv")
# stage map
stage = {}
with open(smap, newline="") as f:
    r = csv.reader(f)
    h = next(r)
    # find cell id col and stage col
    ci = 0
    si = h.index("stage") if "stage" in [x.lower() for x in h] else (1 if len(h) > 1 else 0)
    for row in r:
        if not row:
            continue
        stage[row[ci].strip()] = row[si].strip() if si < len(row) else ""
out["m7_stage_counts"] = dict(Counter(stage.values()))

# aucell
rows = []
with open(amap, newline="") as f:
    r = csv.reader(f)
    h = next(r)
    regcols = h[1:]
    for row in r:
        if not row:
            continue
        cell = row[0].strip()
        vals = [float(x) for x in row[1:]]
        rows.append((cell, vals))
out["m7_aucell_cells"] = len(rows)
out["m7_aucell_regulons"] = len(regcols)

# mean AUC per regulon
mean_auc = [statistics.mean(v[i] for _, v in rows) for i in range(len(regcols))]
top_idx = sorted(range(len(regcols)), key=lambda i: -mean_auc[i])[:12]
out["m7_top_regulons_meanAUC"] = [(regcols[i], round(mean_auc[i], 4)) for i in top_idx]

# stage-resolved mean AUC for top regulons (group stages)
stage_order = ["GV", "MII", "1-cell", "2-cell", "4-cell", "8-cell", "Morula", "Blastocyst",
               "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11",
               "E12", "E13", "E14", "in_vitro", "pgEpiSC"]
# normalize stage labels
def norm(s):
    s = s.strip()
    for key in stage_order:
        if s == key or s.startswith(key):
            return key
    return s
stage_vals = defaultdict(list)
for cell, vals in rows:
    st = norm(stage.get(cell, "NA"))
    stage_vals[st].append(vals)
# pick top 8 regulons, show stage profile
top8 = top_idx[:8]
profile = {}
for i in top8:
    name = regcols[i]
    prof = {}
    for st, lst in stage_vals.items():
        if lst:
            prof[st] = round(statistics.mean(v[i] for v in lst), 4)
    profile[name] = prof
out["m7_top8_stage_profile"] = profile

# ---------- M3c ----------
m3 = os.path.join(BASE, "m3c_validation")
# reactome
rm = os.path.join(m3, "reactome_enrichment.csv")
if os.path.exists(rm):
    with open(rm, newline="") as f:
        r = csv.reader(f)
        h = next(r)
        terms = []
        for row in r:
            if not row:
                continue
            terms.append(row)
    out["m3c_reactome_terms"] = terms[:8]
# peg/meg overlap
pm = os.path.join(m3, "peg_meg_overlap.csv")
if os.path.exists(pm):
    with open(pm, newline="") as f:
        r = csv.reader(f)
        h = next(r)
        out["m3c_peg_meg_header"] = h
        out["m3c_peg_meg_rows"] = [row for row in r]
# pa down mapped
pd = os.path.join(m3, "pa_down_genes_mapped.csv")
if os.path.exists(pd):
    with open(pd, newline="") as f:
        rr = list(csv.reader(f))
    out["m3c_n_pa_down"] = len(rr) - 1 if len(rr) > 1 else 0

# ---------- M4 ----------
m4t = os.path.join(BASE, "m4_oocyte_types.csv")
if os.path.exists(m4t):
    with open(m4t, newline="") as f:
        rr = list(csv.reader(f))
    out["m4_oocyte_types_header"] = rr[0] if rr else []
    out["m4_oocyte_types_n"] = len(rr) - 1
    # cluster counts if a column exists
    if rr and len(rr[0]) > 1:
        col = rr[0][-1]
        out["m4_oocyte_type_counts"] = dict(Counter(row[-1] for row in rr[1:]))
# deg summary (count rows)
m4d = os.path.join(BASE, "m4_oocyte_type_degs.csv")
if os.path.exists(m4d):
    with open(m4d, newline="") as f:
        rr = list(csv.reader(f))
    out["m4_deg_n"] = len(rr) - 1 if len(rr) > 1 else 0

# ---------- M5 ----------
m5 = os.path.join(BASE, "m5_top50_divergent_genes.csv")
if os.path.exists(m5):
    with open(m5, newline="") as f:
        rr = list(csv.reader(f))
    out["m5_top50_header"] = rr[0] if rr else []
    out["m5_top50_genes"] = [row[0] for row in rr[1:11]]  # first 10

# ---------- M6 ----------
m6 = os.path.join(BASE, "m6_metabolism_ratio.csv")
if os.path.exists(m6):
    with open(m6, newline="") as f:
        rr = list(csv.reader(f))
    out["m6_ratio_header"] = rr[0] if rr else []
    out["m6_ratio_rows"] = rr[1:12]

print(json.dumps(out, ensure_ascii=False, indent=2))
