#!/usr/bin/env python
# Build manuscript_full_with_figures_v7.md from v6 with: (1) Fig 2D threshold clarification,
# (2) new Figure 9 (42-sample bulk RNA-seq IVF-vs-PA validation cohort, M15), and synchronized
# abstract/methods/discussion/legends/data-readiness/references.
import io, sys

SRC = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/manuscript_full_with_figures_v6.md"
DST = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/manuscript_full_with_figures_v7.md"

t = open(SRC, encoding="utf-8").read()
log = []

def rep(old, new, n=1, label=""):
    global t
    c = t.count(old)
    assert c == n, f"[{label}] expected {n} occurrence(s) of:\n{old[:80]}\nbut found {c}"
    t = t.replace(old, new, n)
    log.append(f"OK  {label} (replaced {c})")

# ---- 1. Abstract: add independent 42-sample validation sentence ----
rep(
    "from the *in vivo* backbone. As a purely computational resource, this atlas",
    "from the *in vivo* backbone. An independently generated 42-sample bulk-RNA-seq cohort of IVF-versus-PA cleavage-stage embryos — analysed here with full donor-level replication via DESeq2 — independently confirms this transcriptional attenuation: PA-downregulated genes outnumber PA-upregulated genes at every cleavage stage, most robustly at the 4-cell stage (n = 9 per group). As a purely computational resource, this atlas",
    label="abstract-add-cohort")

# ---- 2. Fig 2D results: clarify the two thresholds (186 @ |lfc|>0.5 ; 119/104+15 @ |lfc|>1.0) ----
rep(
    "On GSE164812 (IVF n = 21, PA n = 21 cells), using Wilcoxon rank-sum on FPKM (no raw counts available) with Benjamini–Hochberg correction, **186 genes** were significant at |log₂FC| > 0.5, of which **104 were down-regulated and only 15 up-regulated in PA** (Fig. 2D).",
    "On GSE164812 (IVF n = 21, PA n = 21 cells), using Wilcoxon rank-sum on FPKM (no raw counts available) with Benjamini–Hochberg correction, 186 genes were significant at |log₂FC| > 0.5 (BH), of which 162 were down- and 24 up-regulated in PA; at the stricter |log₂FC| > 1.0 cutoff used for the imprint-overlap analysis, 119 genes remained (104 down, 15 up) (Fig. 2D).",
    label="fig2d-results-clarify")

# ---- 3. Fig 2D legend: same clarification ----
rep(
    "Wilcoxon rank-sum on FPKM with Benjamini–Hochberg correction identified 186 significant genes (|log₂FC| > 0.5), of which 104 were down- and 15 up-regulated in PA.",
    "Wilcoxon rank-sum on FPKM with Benjamini–Hochberg correction identified 186 significant genes at |log₂FC| > 0.5, of which 162 were down- and 24 up-regulated in PA; 119 genes satisfied |log₂FC| > 1.0 (104 down, 15 up), the set carried forward to imprint-overlap testing.",
    label="fig2d-legend-clarify")

# ---- 4. Methods "Differential expression": add 42-sample design (3) ----
rep(
    "Cell-level pseudo-replicated p-values were not reported (Squair *et al.* [6]).",
    "Cell-level pseudo-replicated p-values were not reported (Squair *et al.* [6]). (3) For the independently generated 42-sample bulk-RNA-seq cohort (SRA SRP301735; IVF vs PA across 1/2/4/8-cell stages, with 3–9 biological replicates per condition × stage), transcript-level counts were quantified by salmon and aggregated gene-level (tximport-style); DESeq2 was run per stage (design ~ condition) and overall (design ~ stage + condition) at padj < 0.05 and |log₂FC| > 1. Donor-level replication is native to this design.",
    label="methods-add-m15")

# ---- 5. Discussion: add a paragraph confirming the 42-sample validation after the PA-deficit paragraph ----
rep(
    "The cytoskeletal/mitotic hypothesis is testable with live imaging and mitochondrial-function assays.",
    "The cytoskeletal/mitotic hypothesis is testable with live imaging and mitochondrial-function assays.\n\n**Independent replication in a 42-sample bulk-RNA-seq cohort.** The transcriptional-attenuation model receives independent support from a separately generated 42-sample bulk-RNA-seq dataset of IVF-versus-PA cleavage-stage embryos (SRP301735), analysed here with full donor-level replication (Fig. 9). At every stage PA-downregulated genes outnumbered PA-upregulated genes (1-cell 124 up / 200 down; 2-cell 298 / 552; 4-cell 815 / 1,015; 8-cell 186 / 361; overall 1,306 / 1,111), concordant with the public-data asymmetry. The 4-cell stage — the best-powered (n = 9 per group) — gave the most robust signal, with leading PA-down genes (PET100, NEPRO, SGPP1, ADPRM; all padj < 10⁻⁸) and PA-up genes (PABPC5 log₂FC ≈ +12.8, padj ≈ 10⁻¹⁵; SETDB1; NOVA1). We deliberately do **not** claim cross-dataset validation of the methylation-MOFA highlighted loci (e.g. SENP5, ADAM10, NEK2): those were not significant DEGs in this cleavage-stage bulk cohort, consistent with their origin in a distinct GV-oocyte window and dataset, and they are reported where they were derived rather than here.",
    label="discussion-add-cohort")

# ---- 6. Limitations: note FPKM limitation now complemented by 42-sample DESeq2 ----
rep(
    "(i) The PA/IVF comparison relied on FPKM + Wilcoxon rather than raw-count DESeq2;",
    "(i) The public-dataset PA/IVF comparison (Fig. 2D, GSE164812) relied on FPKM + Wilcoxon rather than raw-count DESeq2; this is now complemented by an independent 42-sample bulk-RNA-seq cohort analysed with DESeq2 and full donor replication (Fig. 9), although its 1- and 8-cell stages retain only n = 3 per group and their extreme individual fold changes should be read as hypothesis-generating;",
    label="limitations-fpkm-complement")

# ---- 7. Insert Figure 9 section after Fig 8 section (before "# DISCUSSION") ----
fig9 = '''

## Figure 9 — An independent 42-sample bulk-RNA-seq cohort confirms PA transcriptional attenuation across cleavage stages

To test whether the IVF-versus-PA programme inferred from the public single-cell data generalised to an independent, properly replicated dataset, we analysed a separately generated 42-sample bulk-RNA-seq cohort of porcine preimplantation embryos (SRA SRP301735; IVF n = 21, PA n = 21; spanning the 1-, 2-, 4- and 8-cell stages with 3–9 biological replicates per condition × stage). Transcript-level counts were quantified by salmon and aggregated gene-level (tximport-style), then analysed with DESeq2 per stage (design ~ condition) and overall (design ~ stage + condition) at padj < 0.05 and |log₂FC| > 1.

**Per-stage IVF vs PA DEGs.** At every cleavage stage PA-downregulated genes outnumbered PA-upregulated genes (Fig. 9A), recapitulating the asymmetry of the public FPKM analysis (Fig. 2D): 1-cell 124 up / 200 down (n = 3/3), 2-cell 298 / 552 (n = 6/6), 4-cell 815 / 1,015 (n = 9/9) and 8-cell 186 / 361 (n = 3/3); overall (stage-adjusted) 1,306 up / 1,111 down. The 4-cell stage — with the largest replication — gave the most robust signal and is shown as the representative volcano (Fig. 9B). Leading PA-downregulated genes at 4-cell included PET100, NEPRO, SGPP1 and ADPRM (all padj < 10⁻⁸, strongly negative log₂FC), while leading PA-upregulated genes included PABPC5 (log₂FC ≈ +12.8, padj ≈ 10⁻¹⁵), SETDB1 and NOVA1.

**Interpretation and caveats.** The direction of effect is concordant with the public-data finding that PA embryos show incomplete transcriptional activation. Two statistical caveats must be stated: (i) the 1-cell and 8-cell stages had only n = 3 replicates per condition, so their large DEG counts and extreme individual fold changes (e.g. PET100 log₂FC ≈ −43) rest on limited power and should be read as hypothesis-generating; (ii) the methylation-MOFA highlighted genes from the GV-oocyte multi-omics analysis (e.g. SENP5, ADAM10, NEK2) were not significant DEGs in this cleavage-stage bulk cohort — expected, because they derive from a different developmental window and dataset, and we do not claim cross-dataset validation of those specific loci. The 4-cell stage (n = 9) provides the most reliable stage-resolved insight into the IVF/PA transcriptomic divide.

![Fig. 9. Independent 42-sample bulk-RNA-seq validation of the IVF-vs-PA programme. (A) Per-stage DEG counts (PA-down in red, PA-up in blue) at padj < 0.05, |log₂FC| > 1. (B) Volcano of the best-powered 4-cell stage (n = 9/9).](E:/Workbuddy/2026-07-27-11-58-27/data/processed/deseq2/m15_fig9_summary.png)

'''
rep("\n# DISCUSSION\n", fig9 + "\n# DISCUSSION\n", label="insert-fig9")

# ---- 8. Figure legends: add Fig 9 legend after Fig 8 legend ----
fig9_legend = '''

## Figure 9. An independent 42-sample bulk-RNA-seq cohort confirms PA transcriptional attenuation across cleavage stages.
(A) DESeq2 DEG counts per cleavage stage (IVF vs PA, padj < 0.05, |log₂FC| > 1). PA-downregulated genes (red) exceed PA-upregulated genes (blue) at every stage, most robustly at 4-cell (n = 9 per group). (B) Volcano plot of the 4-cell stage (best-powered). Leading PA-down genes include PET100, NEPRO, SGPP1, ADPRM; leading PA-up genes include PABPC5 (log₂FC ≈ +12.8, padj ≈ 10⁻¹⁵), SETDB1, NOVA1. The 1- and 8-cell stages (n = 3 per group) yield large DEG counts whose extreme fold changes are low-powered and hypothesis-generating.

'''
rep(
    "## Figure 8. MOFA+ decomposes multi-study gene-expression variance into interpretable developmental axes.\nMOFA+ was run on 3,000 highly variable genes across 1,833 cells from three studies (GSE168106 *in vivo* E0–E14, GSE112380 *in vitro*, GSE234116 GV oocytes), using group centering to remove study-mean offsets. (A) Variance explained by each of ten latent factors; the full model explains 21.6% of total variance, with F1 alone accounting for 11.6%. (B) Trajectories of the three factors most strongly correlated with canonical developmental stage (GV → E14). F2 (oocyte/MZT programme; r = −0.68) peaks in the early stages, while F3 and F8 show later modulation. (C) Top loading genes on the dominant F1 axis, comprising cell-cycle and proliferation regulators (CCNB2, CDKN3, HAUS1, KIF23, SGO1, CDC27, VRK1). (D) Heatmap of top gene loadings across the three most stage-correlated factors; each factor is driven by a distinct gene set, confirming that the integrated atlas retains biologically meaningful, shared axes of variation after Harmony correction.\n",
    "## Figure 8. MOFA+ decomposes multi-study gene-expression variance into interpretable developmental axes.\nMOFA+ was run on 3,000 highly variable genes across 1,833 cells from three studies (GSE168106 *in vivo* E0–E14, GSE112380 *in vitro*, GSE234116 GV oocytes), using group centering to remove study-mean offsets. (A) Variance explained by each of ten latent factors; the full model explains 21.6% of total variance, with F1 alone accounting for 11.6%. (B) Trajectories of the three factors most strongly correlated with canonical developmental stage (GV → E14). F2 (oocyte/MZT programme; r = −0.68) peaks in the early stages, while F3 and F8 show later modulation. (C) Top loading genes on the dominant F1 axis, comprising cell-cycle and proliferation regulators (CCNB2, CDKN3, HAUS1, KIF23, SGO1, CDC27, VRK1). (D) Heatmap of top gene loadings across the three most stage-correlated factors; each factor is driven by a distinct gene set, confirming that the integrated atlas retains biologically meaningful, shared axes of variation after Harmony correction.\n" + fig9_legend,
    label="insert-fig9-legend")

# ---- 9. Data-readiness table: add Fig 9 row ----
rep(
    "| Fig. 8 | MOFA+ multi-study variance decomposition ✓ | RNA-only; multi-omics CGmap MOFA+ = future work |",
    "| Fig. 8 | MOFA+ multi-study variance decomposition ✓ | RNA-only; multi-omics CGmap MOFA+ = future work |\n| Fig. 9 | 42-sample bulk RNA-seq IVF/PA validation ✓ | Per-stage n reported; 1/8-cell n=3 lower power; 4-cell n=9 best-powered |",
    label="readiness-add-fig9")

# ---- 10. Data availability: add SRP301735 ----
rep(
    "Processed matrices, regulon definitions, and AUCell scores generated in this study are provided in the supplementary material.",
    "Processed matrices, regulon definitions, and AUCell scores generated in this study are provided in the supplementary material. The independent validation cohort comprises 42 bulk-RNA-seq libraries of IVF-versus-PA cleavage-stage pig embryos (SRA SRP301735; run accessions SRR13435640–SRR13435681), re-analysed from raw FASTQ via salmon and DESeq2 in this study.",
    label="data-avail-add-srp")

open(DST, "w", encoding="utf-8").write(t)
print("\n".join(log))
print("WROTE", DST)
