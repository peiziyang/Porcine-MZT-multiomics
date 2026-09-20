#!/usr/bin/env python3
"""
M7 full pySCENIC upgrade (part 2): ctx-enriched regulons + AUCell -> Fig.5 panels.
Reads:
  pyscenic_out/aucell_scores_full.csv  (AUCell matrix; orientation auto-detected)
  m7_cell_stage_map.csv                (cell_id, stage_label, stage_order)
Final convention after load: auc is cells (rows) x regulons (columns).
Writes Fig.5 upgrade panels with true developmental-stage resolution
(GV -> E0..E14 -> in_vitro -> pgEpiSC), superseding the co-expression-only version.
"""
import pandas as pd, numpy as np, os, sys, traceback
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/mnt/e/Workbuddy/2026-07-27-11-58-27'
PROC = f'{BASE}/data/processed'
OUT  = f'{PROC}/pyscenic_out'
LOG  = f'{OUT}/ctx_full.log'

def log(m):
    try:
        with open(LOG, 'a') as f:
            f.write(m + "\n")
    except Exception:
        pass
    print(m)

try:
    # ---- stage map ----
    stage = pd.read_csv(f'{PROC}/m7_cell_stage_map.csv')
    cell_label = dict(zip(stage['cell_id'], stage['stage_label']))
    cell_order = dict(zip(stage['cell_id'], stage['stage_order'].fillna(999)))
    known_cells = set(cell_label.keys())

    # ---- AUCell matrix: force to cells (rows) x regulons (columns) ----
    raw = pd.read_csv(f'{OUT}/aucell_scores_full.csv', index_col=0)
    col_is_cells = len(set(raw.columns) & known_cells) > 0.5 * len(raw.columns)
    if col_is_cells:
        auc = raw.T                      # regulons x cells -> cells x regulons
        log("aucell orientation: transposed (regulons x cells -> cells x regulons)")
    else:
        auc = raw                        # already cells x regulons
        log("aucell orientation: cells x regulons (as-is)")
    auc = auc.loc[[c for c in auc.index if c in known_cells]]   # keep known cells (rows)
    log(f"AUCell matrix used: {auc.shape[0]} cells x {auc.shape[1]} regulons")

    # ---- order cells by developmental stage ----
    cells_sorted = sorted(auc.index, key=lambda c: (cell_order.get(c, 999), c))
    auc_ord = auc.loc[cells_sorted]                              # cells x regulons
    labels_ord = [cell_label.get(c, 'NA') for c in cells_sorted]

    # Panel A: top 25 most variable enriched regulons, cells ordered by stage
    # variance of each regulon across cells = auc.var(axis=0)
    row_var = auc.var(axis=0).sort_values(ascending=False)
    top = row_var.head(25).index.tolist()
    M = auc_ord[top].T.values        # regulons (y) x cells (x)
    fig, ax = plt.subplots(figsize=(13, 8))
    im = ax.imshow(M, aspect='auto', cmap='viridis')
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top, fontsize=8)
    step = max(1, len(labels_ord) // 24)
    xt = list(range(0, len(labels_ord), step))
    ax.set_xticks(xt); ax.set_xticklabels([labels_ord[i] for i in xt], rotation=90, fontsize=7)
    ax.set_title('Top 25 most variable enriched regulons (AUCell; cells ordered by stage)')
    plt.colorbar(im, ax=ax, label='AUCell')
    plt.tight_layout(); fig.savefig(f'{OUT}/fig5b_regulon_heatmap_stage.png', dpi=150); plt.close()

    # Panel B: mean enriched regulon activity by developmental stage
    stage_order_list = ['GV'] + [f'E{i}' for i in range(0, 15)] + ['in_vitro', 'pgEpiSC']
    mean_by_stage = pd.DataFrame(index=top)
    for sg in stage_order_list:
        cs = [c for c in cells_sorted if cell_label.get(c) == sg]
        if cs:
            mean_by_stage[sg] = auc.loc[cs, top].mean(axis=0).values   # per-regulon mean
    mean_by_stage = mean_by_stage.fillna(0)
    fig, ax = plt.subplots(figsize=(11, 8))
    im = ax.imshow(mean_by_stage.values, aspect='auto', cmap='RdBu_r')
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top, fontsize=8)
    ax.set_xticks(range(len(stage_order_list))); ax.set_xticklabels(stage_order_list, rotation=90, fontsize=8)
    ax.set_title('Mean enriched regulon activity by developmental stage')
    plt.colorbar(im, ax=ax, label='AUCell')
    plt.tight_layout(); fig.savefig(f'{OUT}/fig5b_regulon_by_stage.png', dpi=150); plt.close()

    # Panel C: key pluripotency / lineage TF regulons (base-name match)
    interesting_bases = ['NANOG', 'POU5F1', 'SOX2', 'CDX2', 'TEAD4', 'GATA3', 'ELF5', 'OCT4']
    present = list(dict.fromkeys(
        r for r in auc.columns if any(r.split('_')[0] == b for b in interesting_bases)))
    if present:
        fig, axes = plt.subplots(1, len(present), figsize=(3 * len(present), 4))
        if len(present) == 1:
            axes = [axes]
        for ax, r in zip(axes, present):
            ax.scatter(range(len(auc_ord[r].values)), auc_ord[r].values, s=3)
            ax.set_title(r, fontsize=9); ax.set_xlabel('cells (stage-ordered)'); ax.set_ylabel('AUCell')
        plt.tight_layout(); fig.savefig(f'{OUT}/fig5b_key_tf_regulons.png', dpi=150); plt.close()
        log(f"key TF panel regulons: {present}")

    log("PLOT DONE")
except Exception:
    log("PLOT ERROR:\n" + traceback.format_exc())
    sys.exit(1)
