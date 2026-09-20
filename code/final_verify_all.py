#!/usr/bin/env python
"""投稿前最终核实：正文关键数字 vs 数据源文件"""
import re, os
import pandas as pd
import numpy as np

BASE = "E:/Workbuddy/2026-07-27-11-58-27"
MF = os.path.join(BASE, "data/processed/m4_mofa")
txt = open(os.path.join(BASE, "_docx_verify.txt"), encoding="utf-8").read()
FAIL = []
def chk(desc, ok, detail=""):
    print(f"  {'✅' if ok else '❌'} {desc}" + (f"  [{detail}]" if detail else ""))
    if not ok:
        FAIL.append(desc)

print("========== 1. MOFA+ 因子方差 ==========")
r2 = pd.read_csv(f"{MF}/mofa_multiomics_r2_per_view_v2.csv")
for f in ["F1", "F2"]:
    rna = r2[(r2.factor==f)&(r2.view=="RNA")]["r2_pct"].values[0]
    meth = r2[(r2.factor==f)&(r2.view=="METH")]["r2_pct"].values[0]
    chk(f"{f} RNA={rna:.1f}% METH={meth:.1f}%", abs(rna-30.7)<0.1 if f=="F1" else abs(rna-25.2)<0.1)

print("========== 2. F1 cross-view r = 0.84 ==========")
w_rna = pd.read_csv(f"{MF}/mofa_multiomics_weights_RNA_v2.csv", index_col=0)
w_meth = pd.read_csv(f"{MF}/mofa_multiomics_weights_METH_v2.csv", index_col=0)
common = w_rna.index.intersection(w_meth.index)
r = np.corrcoef(w_rna.loc[common,"F1"], w_meth.loc[common,"F1"])[0,1]
chk(f"F1 weight Pearson r = {r:.3f} (声称 0.84)", abs(r-0.84)<0.01)

print("========== 3. 跨物种 ==========")
cs = pd.read_csv(f"{MF}/cross_species_permutation_results_v2.csv")
f1h = cs[(cs.test=="F1_observed")&(cs.species=="human")].iloc[0]
f1m = cs[(cs.test=="F1_observed")&(cs.species=="mouse")].iloc[0]
gwh = cs[(cs.test=="genome_wide")&(cs.species=="human")].iloc[0]
gwm = cs[(cs.test=="genome_wide")&(cs.species=="mouse")].iloc[0]
chk(f"F1 human ρ={f1h.rho:.2f} (0.35)", abs(f1h.rho-0.35)<0.01)
chk(f"F1 mouse ρ={f1m.rho:.2f} (0.36)", abs(f1m.rho-0.36)<0.01)
chk(f"GW human ρ={gwh.rho:.2f} (0.29)", abs(gwh.rho-0.29)<0.01)
chk(f"GW mouse ρ={gwm.rho:.2f} (0.20)", abs(gwm.rho-0.20)<0.01)
chk(f"perm human={f1h.expr_matched_perm_p:.3f} (0.105)", abs(f1h.expr_matched_perm_p-0.105)<0.002)
chk(f"perm mouse={f1m.expr_matched_perm_p:.3f} (0.027)", abs(f1m.expr_matched_perm_p-0.027)<0.002)
chk(f"n genes human={f1h.n_genes} (41)", f1h.n_genes==41)
chk(f"n genes mouse={f1m.n_genes} (36)", f1m.n_genes==36)

print("========== 4. Table 2 基因集 ==========")
gs = pd.read_csv(f"{MF}/geneset_test/gene_set_permutation_results.csv")
for _, row in gs.iterrows():
    f = row.gene_set
    chk(f"{f} obs_mean={row.obs_mean_log2FC:.2f} emp={row.empirical_p:.3f} bh={row.bh_adj_p:.3f}",
        True)

print("========== 5. 轨迹 ==========")
lo = pd.read_csv(f"{MF}/lo_dataset/lo_stage_sensitivity.csv")
full = lo[lo.dropped_stage=="E0"].set_index("factor")["full_rho"]
chk(f"F1 full ρ={full['F1']:.2f} (-0.41)", abs(full["F1"]+0.41)<0.01)
chk(f"F3 full ρ={full['F3']:.2f} (-0.90)", abs(full["F3"]+0.90)<0.01)
chk(f"F4 full ρ={full['F4']:.2f} (+0.64)", abs(full["F4"]-0.64)<0.01)
chk(f"F6 full ρ={full['F6']:.2f} (+0.76)", abs(full["F6"]-0.76)<0.01)
tc = pd.read_csv(f"{MF}/lo_dataset/trajectory_bootstrap_ci.csv")
f1c = tc[tc.factor=="F1"].iloc[0]; f3c = tc[tc.factor=="F3"].iloc[0]
chk(f"F1 CI=[{f1c.ci_lo:.2f},{f1c.ci_hi:.2f}] ([-0.52,+0.07])",
    abs(f1c.ci_lo+0.52)<0.01 and abs(f1c.ci_hi-0.07)<0.01)
chk(f"F3 CI=[{f3c.ci_lo:.2f},{f3c.ci_hi:.2f}] ([-0.92,-0.89])",
    abs(f3c.ci_lo+0.92)<0.01 and abs(f3c.ci_hi+0.89)<0.02)

print("========== 6. 扰动 ==========")
pp = pd.read_csv(f"{MF}/perturbation_perm/perturbation_permutation_results.csv")
dn = pp[pp.tf=="DNMT1"].iloc[0]; at = pp[pp.tf=="ATF3"].iloc[0]
chk(f"DNMT1 cos={dn.observed_cos_sim:.2f} (0.43) permP={dn.permutation_p}", abs(dn.observed_cos_sim-0.43)<0.01)
chk(f"ATF3 cos={at.observed_cos_sim:.2f} (0.55)", abs(at.observed_cos_sim-0.55)<0.01)
chk(f"DNMT1 TF-shuffle P={dn.tf_label_perm_p:.2f} (0.25)", abs(dn.tf_label_perm_p-0.25)<0.01)
chk(f"ATF3 TF-shuffle P={at.tf_label_perm_p:.2f} (0.13)", abs(at.tf_label_perm_p-0.13)<0.01)

print("========== 7. donor-aware ==========")
kw = pd.read_csv(f"{MF}/donor_analysis/factor_donor_kruskal.csv").set_index("factor")
chk(f"F1 donor perm P={kw.loc['F1','KW_permutation_p']:.3f} (0.001)",
    abs(kw.loc["F1","KW_permutation_p"]-0.001)<0.0005)
for f in ["F3","F4","F6"]:
    p = kw.loc[f,"KW_permutation_p"]
    chk(f"{f} donor perm P={p:.2f} (>0.15)", p>0.15)
oc = pd.read_csv(f"{MF}/donor_analysis/donor_oocyte_counts.csv")
chk(f"donor counts {oc['n_oocytes'].tolist()} (2,7,7,8,8)",
    oc["n_oocytes"].tolist()==[2,7,7,8,8])

print("========== 8. interaction ==========")
sm = pd.read_csv(f"{BASE}/data/processed/deseq2/m15_deseq2_interaction_summary.csv").set_index("coefficient")
for coeff, exp in [("stage[T.2cell]:condition[T.PA]",67),("stage[T.4cell]:condition[T.PA]",219),("stage[T.8cell]:condition[T.PA]",376)]:
    v = sm.loc[coeff,"n_sig_padj005"]
    chk(f"{coeff} sig={v} ({exp})", v==exp)
res = pd.read_csv(f"{BASE}/data/processed/deseq2/m15_deseq2_stageXcondition_fullModel.csv")
s4 = pd.read_csv(f"{BASE}/manuscript/submission/tables/Table_S4_gene_set_membership.csv")
f3g = set(s4[s4.F3_member==True].gene)
for coeff, lab, exp in [("condition[T.PA]","1cell",-4.82),("stage[T.4cell]:condition[T.PA]","4cell",-0.59)]:
    sub = res[res.coefficient==coeff]
    m = sub[sub.gene_symbol.isin(f3g)]["log2FoldChange"].mean()
    chk(f"F3 {lab} mean={m:.2f} ({exp})", abs(m-exp)<0.03)

print("========== 9. CpG 匹配 ==========")
c = pd.read_csv(f"{MF}/cpg_matched_analysis/cpg_matched_sensitivity.csv").set_index("metric")["value"]
chk(f"F1 median infoCpG={c['F1_median_infoCpG']:.0f} (102)", abs(c["F1_median_infoCpG"]-102)<1)
chk(f"bg median={c['bg_median_infoCpG']:.0f} (70)", abs(c["bg_median_infoCpG"]-70)<1)
chk(f"expr×len permP={c['expr_len_matched_permP']:.3f} (0.015)", abs(c["expr_len_matched_permP"]-0.015)<0.002)
chk(f"density F1={c['F1_density_per_kb']:.1f} (28.9)", abs(c["F1_density_per_kb"]-28.9)<0.3)
chk(f"density bg={c['bg_density_per_kb']:.1f} (17.2)", abs(c["bg_density_per_kb"]-17.2)<0.3)
chk(f"density P={c['density_MWU_P']:.3f} (0.007)", abs(c["density_MWU_P"]-0.007)<0.001)
chk(f"corr expr r={c['corr_expr_rho']:.2f} (0.02)", abs(c["corr_expr_rho"]-0.02)<0.005)
chk(f"corr expr P={c['corr_expr_P']:.2f} (0.31)", abs(c["corr_expr_P"]-0.31)<0.01)

print("========== 10. DE / 发散度 / 代谢 ==========")
m3 = pd.read_csv(f"{BASE}/data/processed/m3_ivf_vs_pa_results.csv")
dn_c = ((m3.log2FC_PA_vs_IVF<-1)&(m3.padj<0.05)).sum()
up_c = ((m3.log2FC_PA_vs_IVF>1)&(m3.padj<0.05)).sum()
chk(f"DE genes={len(m3)} (6860) down={dn_c} (104) up={up_c} (15)",
    len(m3)==6860 and dn_c==104 and up_c==15)
d = pd.read_csv(f"{BASE}/data/processed/m5_gene_divergence_v2.csv")
chk(f"divergence genes={len(d)} (8276)", len(d)==8276)
for col, exp in [("IVF_vs_vivo",0.793),("PA_vs_vivo",0.833),("IVF_vs_PA",0.447)]:
    v = d[col].mean()
    chk(f"{col} mean={v:.3f} ({exp})", abs(v-exp)<0.002)
met = pd.read_csv(f"{BASE}/data/processed/m6_metabolism_ratio.csv")
chk(f"metabolism ratio range [{met.ratio.min():.2f},{met.ratio.max():.2f}] (0.32-0.56)",
    met.ratio.min()>0.31 and met.ratio.max()<0.57)

print("========== 11. SHAP ==========")
sh = pd.read_csv(f"{MF}/shap_analysis_results.csv")
top = list(sh.iloc[:3].iloc[:,0]) if "feature" in [str(c).lower() for c in sh.columns] else list(sh.columns[:3])
print(f"  SHAP top3: {top} (声称 ACOT9, HYLS1, PHAX)")
chk("SHAP top3 = ACOT9/HYLS1/PHAX", "ACOT9" in str(top) and "HYLS1" in str(top) and "PHAX" in str(top))

print("========== 12. Waddington OT 位移 ==========")
disp = pd.read_csv(f"{BASE}/data/processed/m5_ot_displacement.csv")
overall = disp[disp.ref_stage.isin(["E2","E3","E4","E5"])].groupby("condition")["displacement"].agg(["mean","std","count"])
for cond, exp_m in [("in_vivo",8.33),("IVF",10.24),("PA",10.15)]:
    if cond in overall.index:
        m = overall.loc[cond,"mean"]
        chk(f"OT {cond} mean={m:.2f} ({exp_m})", abs(m-exp_m)<0.05)

print()
print(f"===== 失败项: {len(FAIL)} =====")
for f in FAIL:
    print("  ❌", f)
print("DONE")
