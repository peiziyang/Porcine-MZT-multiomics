import re

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/"
SRC = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/Results_Discussion_Fig1-6.md"
OUT = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/Results_Discussion_Fig1-6_with_figures.md"

FIG_IMAGES = {
 1: [
   (BASE + "m2_umap_harmony_by_stage_group.png",
    "Fig. 1. Harmony-batch-corrected UMAP of the integrated porcine preimplantation atlas, coloured by developmental stage (GV → E14 → pgEpiSC). The full oocyte-to-post-implantation trajectory is preserved."),
   (BASE + "m2_umap_harmony_by_source.png",
    "Fig. 1 (panel B). Same embedding coloured by reproductive source (IVF / PA / in vivo / GV); cells from the same developmental window co-cluster, showing stage — not batch — is the dominant axis of variation."),
   (BASE + "m2_harmony_comparison.png",
    "Fig. 1 (panel E). Batch-removal diagnostics: dataset variance 7.7% → 2.1%, stage variance 18.3% → 15.5%, stage:batch ratio 2.4 → 7.4."),
 ],
 2: [
   (BASE + "m5_trajectory_comparison.png",
    "Fig. 2A–B. Cross-condition pseudotime divergence. PA embryos diverge most from the in vivo trajectory (PA_vs_vivo mean = 3.24) versus IVF_vs_vivo (2.74) and IVF_vs_PA (0.90)."),
   (BASE + "m3_ivf_vs_pa_volcano.png",
    "Fig. 2D. IVF vs PA volcano (Wilcoxon on FPKM, BH). 186 significant genes; 104 down- vs 15 up-regulated in PA."),
   (BASE + "m3c_validation/reactome_enrichment.png",
    "Fig. 2 inset. Reactome enrichment of the 104 PA-down genes (Rho GTPase, mitotic centromere, SUMOylation)."),
   (BASE + "m3c_validation/peg_overlap.png",
    "Fig. 2 inset. Zero overlap of PA-down genes with paternally/maternally expressed imprint loci (Fisher p = 1.0)."),
 ],
 3: [
   (BASE + "m3_volcano_early_vs_late.png",
    "Fig. 3A. Pseudobulk volcano: Early (E0–E3) vs Late (E11–E14). 11,621 DEGs at padj < 0.05, |log2FC| > 1."),
   (BASE + "m3_maplot_early_vs_late.png",
    "Fig. 3B. MA plot of the same contrast; large effect sizes span ~half the transcriptome."),
 ],
 4: [
   (BASE + "m4_multiomics_overview.png",
    "Fig. 4. Multi-omics re-analysis of GSE234116 GV oocytes: reproduced Type I/II classification, DEG profile, and referenced DMR hypermethylation in Type II."),
 ],
 5: [
   (BASE + "pyscenic_out/fig5b_regulon_heatmap_stage.png",
    "Fig. 5A–B. pySCENIC regulon heatmap (46 motif-validated regulons × cells, ordered by stage) and AUCell activity."),
   (BASE + "pyscenic_out/fig5b_regulon_by_stage.png",
    "Fig. 5C. Stage-resolved mean AUC for key regulons (ESRRA, MXD3, POU5F1/OCT4, TFAP2C, GATA4, BACH1, ATF3)."),
   (BASE + "pyscenic_out/fig5b_key_tf_regulons.png",
    "Fig. 5D. Key TF regulon activity trajectories across development."),
   (BASE + "m9_lineage_rss_figure.png",
    "Fig. 5E. Lineage-level regulon specificity in E6–E14 blastocyst cells (cluster marker profiles + lineage-enriched regulons)."),
   (BASE + "m9_tf_network.png",
    "Fig. 5F. Core TF regulatory network (6 TFs × 20 targets, 120 co-expression edges)."),
 ],
 6: [
   (BASE + "m6_metabolism.png",
    "Fig. 6. Glycolysis/OxPhos module-score ratio across the preimplantation time course and in IVF vs PA. Ratio < 0.6 throughout; OxPhos-dominant."),
 ],
 7: [
   (BASE + "m5_waddington_ot.png",
    "Fig. 7. Waddington-OT lineage reconstruction on the in-vivo GSE168106 backbone (E0–E8): UMAP by fate entropy, entropy/probability trajectories, per-stage IVF/PA displacement, and UMAP by max terminal-fate probability."),
 ],
 8: [
   (BASE + "m8_mofa_figure.png",
    "Fig. 8. MOFA+ latent-factor map: variance explained per factor, stage-correlated factor trajectories, top F1 genes, and gene x factor weight heatmap."),
 ],
}

lines = open(SRC, encoding='utf-8').read().split('\n')
out = []
cur_fig = None
inserted = set()
for line in lines:
    m = re.match(r'^## Figure (\d+) ', line)
    if m:
        cur_fig = int(m.group(1))
    if line.strip() == '---' and cur_fig is not None and cur_fig not in inserted:
        for path, cap in FIG_IMAGES.get(cur_fig, []):
            out.append("![%s](%s)" % (cap, path))
        out.append("")  # spacer
        inserted.add(cur_fig)
        cur_fig = None
    out.append(line)

open(OUT, 'w', encoding='utf-8').write('\n'.join(out))
print("Wrote", OUT, "with figures for:", sorted(inserted))
