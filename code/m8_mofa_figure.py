#!/usr/bin/env python
"""Fig.8: MOFA+ multi-group factor map of the porcine preimplantation atlas.
Reads m8_mofa_*.csv and renders a 2x2 composite PNG.
Panels: A variance explained per factor; B factor trajectory over canonical
stages (GV→E14); C top genes on the dominant factor; D heatmap of top genes
x top-3 factors.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/"
OUT = BASE + "m8_mofa_figure.png"

r2 = pd.read_csv(BASE + "m8_mofa_r2.csv")
fac = pd.read_csv(BASE + "m8_mofa_factors.csv")
load = pd.read_csv(BASE + "m8_mofa_loadings.csv", index_col=0)
top = pd.read_csv(BASE + "m8_mofa_topgenes.csv")
fbst = pd.read_csv(BASE + "m8_mofa_factor_by_stage.csv")

K = r2.shape[0]
factors = [f"F{k+1}" for k in range(K)]
r2v = r2["R2"].values
total_model = r2v.sum()

# ---- canonical stage ordering (exclude in_vitro / pgEpiSC from trajectory) ----
canonical = ["GV"] + [f"E{i}" for i in range(15)]  # GV, E0..E14
st = fbst.copy().set_index("stage_label")
stages = [s for s in canonical if s in st.index]
xs = np.arange(len(stages))

# correlations with stage_order, computed only on canonical stages
so_map = {s: i - 1 for i, s in enumerate(stages)}  # GV=-1, E0=0, ...
stage_vals = np.array([so_map[s] for s in stages])
corr = np.array([np.corrcoef(st.loc[stages, f].values, stage_vals)[0, 1] for f in factors])
# if correlation is nan (zero variance), set 0
corr = np.where(np.isfinite(corr), corr, 0.0)
traj_ks = list(np.argsort(-np.abs(corr))[:3])

fig, ax = plt.subplots(2, 2, figsize=(13, 10.5))
plt.subplots_adjust(hspace=0.34, wspace=0.30)

# ---- A: variance explained per factor ----
ax[0, 0].barh(range(K)[::-1], r2v * 100, color="#4C72B0")
for i, v in enumerate(r2v * 100):
    ax[0, 0].text(v + 0.15, (K - 1 - i), f"{v:.1f}%", va="center", fontsize=8)
ax[0, 0].set_yticks(range(K)[::-1]); ax[0, 0].set_yticklabels(factors, fontsize=9)
ax[0, 0].set_xlabel("Variance explained (%)")
ax[0, 0].set_title("A. MOFA+ variance explained per latent factor", fontsize=11, fontweight="bold")
ax[0, 0].set_xlim(0, max(r2v * 100) * 1.22)
ax[0, 0].text(0.98, 0.04, f"model explains {total_model*100:.1f}% of total variance;\nF1 alone {r2v[0]*100:.1f}%",
              transform=ax[0, 0].transAxes, ha="right", fontsize=7.5, color="#555")

# ---- B: factor trajectory along canonical developmental stage ----
for k in traj_ks:
    ax[0, 1].plot(xs, st.loc[stages, factors[k]].values, marker="o", label=f"{factors[k]} (r={corr[k]:.2f})")
ax[0, 1].set_xticks(xs); ax[0, 1].set_xticklabels(stages, rotation=90, fontsize=8)
ax[0, 1].set_xlabel("Developmental stage")
ax[0, 1].set_ylabel("Mean factor score")
ax[0, 1].set_title("B. Latent factors track GV → E14 progression", fontsize=11, fontweight="bold")
ax[0, 1].legend(fontsize=8, loc="best")

# ---- C: top genes on the dominant (highest-R2) factor ----
dom = int(np.argmax(r2v))
tg = top[top["factor"] == factors[dom]].head(15).iloc[::-1]
cols = ["#d62728" if w > 0 else "#1f77b4" for w in tg["weight"]]
ax[1, 0].barh(range(len(tg)), tg["weight"].values, color=cols)
ax[1, 0].set_yticks(range(len(tg))); ax[1, 0].set_yticklabels(tg["gene"].values, fontsize=8)
ax[1, 0].axvline(0, color="k", lw=0.6); ax[1, 0].set_xlabel("MOFA+ weight")
ax[1, 0].set_title(f"C. Top genes on {factors[dom]} (dominant axis, {r2v[dom]*100:.1f}%)", fontsize=10.5, fontweight="bold")

# ---- D: heatmap top genes x top-3 factors ----
hk = traj_ks[:3]
genes_union = []
for k in hk:
    genes_union += list(top[top["factor"] == factors[k]].head(10)["gene"])
genes_union = list(dict.fromkeys(genes_union))[:24]
H = load.loc[genes_union, [factors[k] for k in hk]].values
im = ax[1, 1].imshow(H, cmap="RdBu_r", aspect="auto", vmin=-np.abs(H).max(), vmax=np.abs(H).max())
ax[1, 1].set_xticks(range(len(hk))); ax[1, 1].set_xticklabels([factors[k] for k in hk], fontsize=9)
ax[1, 1].set_yticks(range(len(genes_union))); ax[1, 1].set_yticklabels(genes_union, fontsize=7)
ax[1, 1].set_title("D. Gene x factor weight structure (top genes, top-3 factors)", fontsize=10.5, fontweight="bold")
fig.colorbar(im, ax=ax[1, 1], shrink=0.8, label="weight")

fig.suptitle("Fig. 8. MOFA+ latent-factor map of the porcine preimplantation reference atlas",
             fontsize=13, fontweight="bold")
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print("saved", OUT)
print(f"dominant factor = {factors[dom]} (R2={r2v[dom]*100:.2f}%)")
print("top-3 stage-correlated factors:", [factors[k] for k in traj_ks], "corrs:", np.round(corr[traj_ks], 2))
print(f"total model variance explained = {total_model*100:.2f}%")
print("F1 top genes:", list(top[top['factor']=='F1'].head(10)['gene']))
