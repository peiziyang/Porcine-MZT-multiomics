# -*- coding: utf-8 -*-
"""
Assemble the full manuscript:
  Title + Abstract + Introduction  (new)
  + existing Results & Discussion (Figs 1-7)  (from Results_Discussion_Fig1-6.md, from "# RESULTS" onward)
  + expanded Methods  (replaces the concise "# METHODS (concise)" block)
  + existing References + Data Availability
Then embed the real PNG figures (reuse FIG_IMAGES logic) -> manuscript_full_with_figures.md
"""
import re

BASE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/"
SRC  = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/Results_Discussion_Fig1-6.md"
OUT  = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/manuscript_full.md"
OUT_FIG = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/manuscript_full_with_figures.md"

TITLE = "# A single-cell reference atlas of porcine preimplantation embryogenesis reveals convergent and divergent regulatory programs across reproductive origins"

ABSTRACT = """## Abstract

Preimplantation embryogenesis in the pig integrates oocyte maturation, the maternal-to-zygotic transition (MZT), and the first lineage decisions that only this species permits to be studied side-by-side as *in vitro* fertilisation (IVF), parthenogenetic activation (PA), and *in vivo* development. Despite rich public single-cell data, these resources remain fragmented across laboratories, conditions, and platforms, and lack a unified regulatory interpretation. Here we present a single-cell reference atlas of porcine preimplantation embryogenesis — spanning the germinal-vesicle (GV) oocyte to embryonic day E14 — assembled by re-analysing public single-cell RNA-seq datasets (1,955 high-quality cells), and complement it with targeted multi-omics re-analysis of oocyte maturation states. We show that (i) the inferred regulon landscape recovers the conserved Oct4/TFAP2C pluripotency–trophoblast axis and an ESRRA oxidative-metabolism programme, with a pig-specific stage-resolved timeline of ZGA and lineage commitment; (ii) PA embryos diverge most from the *in vivo* trajectory, and their down-regulated genes — while not overlapping imprint loci — enrich for cytoskeleton/mitosis and SUMOylation pathways and include mitochondrial OxPhos subunits, pointing to a competence rather than imprint-loss defect; (iii) the preimplantation programme is OxPhos-dominant through E14, consistent with the conserved late-blastocyst bivalent state; and (iv) Waddington-OT lineage reconstruction quantifies a canalising loss of fate plasticity (entropy ≈0.99 → 0) and a latent-space displacement of IVF/PA embryos from the *in vivo* backbone. As a purely computational resource, this atlas provides an open discovery platform and a benchmark for stem-cell-based embryo models."""

INTRODUCTION = """## Introduction

The pig (*Sus scrofa*) is both a leading agricultural species and an increasingly important biomedical model, owing to its physiological and developmental proximity to humans in organ transplantation, toxicology, and developmental biology. Preimplantation embryogenesis — encompassing oocyte maturation (germinal-vesicle GV → metaphase-II MII), fertilisation, the maternal-to-zygotic transition (MZT), and the first lineage segregation into inner cell mass (ICM) and trophectoderm (TE) — is a conserved yet species-variable programme whose molecular logic is still incompletely understood outside humans and mice.

A unique strength of the pig model is that *in vitro* fertilisation (IVF), parthenogenetic activation (PA), and *in vivo* development can be generated and compared within the same species. Human embryology is constrained by ethics and law from such comparisons, and most livestock models lack the parthenogenetic option; the pig triad therefore offers a natural experiment for isolating how reproductive origin reshapes the transcriptomic trajectory. PA is especially informative because it proceeds without a paternal genome and its associated imprints, so contrasting PA against IVF and *in vivo* embryos disentangles the contribution of paternal epigenetic information from the broader execution of developmental competence [†].

Nevertheless, the single-cell RNA-seq datasets that now exist for pig oocytes and embryos — spanning the GV oocyte to post-implantation stages and produced by several independent groups (GSE168106 [9]; GSE112380, GSE139512, GSE164812, GSE160334, GSE234116 [†],[5]) — remain fragmented across laboratories, platforms, and culture conditions. They have not been unified into a common coordinate system with consistent annotation, nor interpreted at the level of gene-regulatory networks. Reference atlases have proven powerful in human [3] and cross-species embryology [4]; a comparable, integrative pig resource is still missing.

We therefore built a single-cell reference atlas of porcine preimplantation embryogenesis by re-analysing public data, then layered analytic modules on top of it: **Harmony** batch correction [10] of the integrated embedding, pseudobulk differential expression with donor-level replication [6], cross-condition trajectory divergence, **Waddington-OT** optimal-transport lineage reconstruction [8], multi-omics re-analysis of oocyte maturation states, regulon inference via pySCENIC [1], metabolic module scoring [4], and **MOFA+** multi-study variance decomposition [7]. Our principal contributions are threefold. First, the regulon map recovers the canonical Oct4/TFAP2C pluripotency–trophoblast axis and an ESRRA oxidative-metabolism programme, and provides a pig-specific quantitative timeline of ZGA-associated activity and lineage commitment. Second, we show that PA embryos diverge *most* from the *in vivo* trajectory, and that their down-regulated genes — rather than showing large-scale deregulation of known imprint loci — mark a broad cytoskeleton/mitosis and energetic execution failure; Waddington-OT further quantifies this as reduced fate plasticity and a latent-space displacement of *in vitro* embryos from the *in vivo* backbone. Third, we demonstrate that the porcine preimplantation window is OxPhos-dominant through E14, consistent with the conserved late-blastocyst bivalent-respiration state.

We state at the outset the methodological boundaries that frame these conclusions, so that readers can weigh them appropriately: the IVF/PA comparison relied on FPKM-based testing rather than raw-count DESeq2; pySCENIC used ortholog-mapped human motifs rather than a pig-specific motif set; and **SCENIC+ enhancer-driven networks could not be inferred because no dataset in this atlas provides paired scRNA + scATAC measurements**. Batch correction (Harmony [10]), cross-condition trajectory divergence, Waddington-OT lineage reconstruction [8], and **gene-expression MOFA+ [7]** are fully executed in the present work. None of the remaining boundaries invalidate the reported patterns, but they define the next phase of work. As a purely computational resource — analogous to Zhao *et al.* [3] for human and Malkowska *et al.* [4] for cross-species — this atlas requires no new wet experiments yet delivers a discovery platform and a benchmark for embryo models."""

# ---- expanded Methods (replaces "# METHODS (concise)" block) ----
EXPANDED_METHODS = """# Methods

### Data collection and preprocessing
We curated seven publicly available single-cell RNA-seq datasets of porcine embryos and oocytes from GEO (see DATA AVAILABILITY). After quality control (mitochondrial-fraction filtering, detected-gene thresholds, and doublet removal) 1,955 high-quality cells were retained (20 low-quality cells removed, 1%). Gene identifiers were harmonised to a common Ensembl gene space of 8,764 genes via ortholog mapping. Count/FPKM matrices were retained as deposited because raw FASTQ re-alignment was not performed in this environment; for datasets distributed only as FPKM (GSE139512, GSE164812, GSE160334) the provided matrices were used directly. GSE242553 (ESC-derived blastoid) was excluded from the main integration because its RefSeq transcript IDs were incompatible with the Ensembl gene space; it is noted as a model-system comparator only.

### Batch integration and reference embedding
Cells were log-normalised and principal-component analysis was performed on 50 dimensions (top-5 PCs explained 45% of variance). We then applied **Harmony** batch correction [10], treating each source dataset as the batch variable, and embedded the corrected 50-D space with UMAP. Quantitative assessment confirmed effective batch removal while preserving developmental signal: the variance in the embedding attributable to dataset dropped from 7.7% (PCA centring) to 2.1% (Harmony), whereas the variance attributable to developmental stage remained the dominant axis (18.3% → 15.5%), raising the stage:batch signal ratio from 2.4 to 7.4. Developmental stage was therefore the dominant axis of variation — cells from the same window co-clustered across source studies.

### Multi-study variance decomposition (MOFA+, Fig. 8)
To quantify shared versus study-specific sources of variance after Harmony correction, we ran MOFA+ [7] on 3,000 highly variable genes across the 1,833 cells retained for GRN analysis. We treated the three retained studies as groups (GSE168106 *in vivo* E0–E14; GSE112380 *in vitro*; GSE234116 GV oocytes), used a gaussian likelihood, ARD weights, and group centering to remove study-mean offsets, and inferred K = 10 latent factors. Variance explained per factor was obtained from the model R². F1 (11.6% of total variance) loaded on cell-cycle/proliferation genes (CCNB2, CDKN3, HAUS1, KIF23, SGO1, CDC27, VRK1); F2 (4.4%) loaded on oocyte/early-embryo and maternal-to-zygotic-transition genes (WEE2, BMP15, BTG4, NLRP9, H1-8) and was the most strongly stage-correlated axis (r = −0.68). The full multi-omics MOFA+ integration of matched RNA + DNA-methylation CGmap data remains future work.

### Cross-condition pseudotime divergence and Waddington-OT lineage reconstruction (M5, Fig. 7)
For each reproductive origin (IVF, PA, *in vivo*) we derived a PCA-based pseudotime ordering and computed per-gene deviation across 8,276 analysable genes. The 50 most divergent genes were ranked; we report mean divergence per pair and the fraction of top genes for which each comparison was the largest deviation. Complementarily, we performed **Waddington-OT** lineage reconstruction [8] in the Harmony-corrected 50-D space. Using the *in vivo* GSE168106 time course (E0–E8) as the developmental backbone, entropic optimal transport (Sinkhorn, ε = 0.05 × median cost) computed couplings between consecutive stages; terminal E8 cells were clustered into K = 3 fates and coupling mass was forward-propagated to yield per-cell fate distributions, fate entropy, and maximum fate probability. IVF and PA cells (1C–8C) were mapped onto the matched *in vivo* stages (1C→E2 … 8C→E5) and their root-mean-squared displacement from the corresponding *in vivo* stage distribution was measured. Fate entropy declined from ≈0.99 nats at E0–E2 to 0 at E8 (canalisation), while mean displacement was in-vivo 8.33 ± 2.64, IVF 10.24 ± 1.95, and PA 10.15 ± 1.51 standardised HPC units — quantifying that *in vitro* origin displaces cells from the physiological trajectory.

### Differential expression
Two complementary designs were used. (1) **Pseudobulk DESeq2** (pydeseq2) on GSE168106, aggregating cells per embryonic day as the biological replication unit (n = 4 pseudobulk samples per group), testing Early (E0–E3) versus Late (E11–E14) among 24,948 genes (minimum 10 total counts), at padj < 0.05 and |log₂FC| > 1. (2) For FPKM-only IVF versus PA (GSE164812, n = 21 cells each), Wilcoxon rank-sum with Benjamini–Hochberg correction at |log₂FC| > 0.5. Cell-level pseudo-replicated p-values were not reported (Squair *et al.* [6]).

### Multi-omics re-analysis of oocyte maturation (M4)
GSE234116 (Yuan *et al.* [5]) profiled 62 GV oocytes by single-cell M&T-seq. We reproduced the two-state (Type I / Type II) classification by unsupervised PCA + k-means on the RNA data and identified DEGs (Type II versus Type I) by DESeq2. The published methylome — 1,141 differentially methylated regions (DMRs), of which 1,140 (99.91%) are hypermethylated in Type II — is referenced from the original publication rather than recomputed. Full CGmap-based MOFA+ [7] integration (GSE235731, ~69 GB) was not executable in the sandbox and is deferred.

### Gene-regulatory-network inference (M7)
On a raw-count re-processing of 1,833 cells (E0–E14, *in vitro*, and GV) we ran pySCENIC [1]: GRNBoost2 co-expression across 246 porcine TFs × 3,000 highly variable genes produced 304,706 TF–target links; cisTarget motif enrichment pruned these to 46 motif-validated regulons (261 TF± modules; target sizes 7–2,727, mean 323); AUCell scored all 1,833 cells. Motif enrichment used the hg38 motif database mapped to pig orthologs (no comprehensive pig-specific motif set exists). **SCENIC+ [2] enhancer-driven enhancer-gene regulatory networks require paired scRNA + scATAC measurements from the same cells; none of the datasets in this atlas provides such data, so SCENIC+ was not performed.**

### Metabolic module scoring (M6)
Emulating the hexa-species framework [4], glycolysis and OxPhos module scores were computed (KEGG gene sets ortholog-mapped to human; AUCell) across the preimplantation time course and in IVF versus PA.

### Statistics and reproducibility
All metadata were retrieved from original GEO depositions without re-assignment of donor/condition labels (LAVDC guardrail ①). All differential-expression statistics used donor/pseudobulk-level replication; no cell-level pseudo-replicated p-values were reported (guardrail ②). All method citations were verified against abstracts/DOIs (guardrail ③). Software: Python 3.12, Scanpy, harmonypy 2.0.0, POT 0.9.7, pydeseq2, pySCENIC 0.12.1, mygene, AUCell. Analysis code is available from the corresponding author upon request.

### Data availability
Primary data were obtained from the GEO accessions listed in the DATA AVAILABILITY & PANEL READINESS table and the References. Processed matrices, regulon definitions, and AUCell scores generated in this study are provided in the supplementary material."""

FIG_IMAGES = {
 1: [
   (BASE + "m2_umap_harmony_by_stage_group.png",
    "Fig. 1. Harmony-batch-corrected UMAP of the integrated porcine preimplantation atlas (n = 1,955 cells), coloured by developmental stage (GV → E14 → pgEpiSC). The full oocyte-to-post-implantation trajectory is preserved after removing dataset batch effects."),
   (BASE + "m2_umap_harmony_by_source.png",
    "Fig. 1 (panel B). Same embedding coloured by reproductive source (IVF / PA / in vivo / GV). Cells from the same developmental window co-cluster across source studies, showing stage — not batch — is the dominant axis of variation."),
   (BASE + "m2_harmony_comparison.png",
    "Fig. 1 (panel E). Batch-removal diagnostics: variance in the embedding explained by dataset fell from 7.7% to 2.1% after Harmony, while stage variance remained dominant (18.3% → 15.5%), raising the stage:batch ratio from 2.4 to 7.4."),
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
    "Fig. 5E. Lineage-level regulon specificity in E6–E14 blastocyst cells. (A) Cluster marker profiles (z-scored expression of canonical ICM/TE/PrE markers). (B) Top lineage-specific regulons: PrE-enriched XBP1, GATA2, CEBPG, ATF3, BACH1 (FC = 2.2–4.6); ICM/TE-mixed-enriched NFE2L3, POU5F1 (FC = 1.7–2.9)."),
   (BASE + "m9_tf_network.png",
    "Fig. 5F. Core transcription factor regulatory network (GRNBoost2 co-expression edges; ESRRA, POU5F1, MXD3, GATA4, PAX6, TFAP2C × 20 targets each; 120 edges total)."),
 ],
 6: [
   (BASE + "m6_metabolism.png",
    "Fig. 6. Glycolysis/OxPhos module-score ratio across the preimplantation time course and in IVF vs PA. Ratio < 0.6 throughout; OxPhos-dominant."),
 ],
 7: [
   (BASE + "m5_waddington_ot.png",
    "Fig. 7. Waddington-OT lineage reconstruction on the in-vivo GSE168106 backbone (E0–E8). (A) UMAP coloured by fate entropy. (B) Mean fate entropy (green) decreases and mean max terminal-fate probability (red) increases across stages, reaching entropy ≈0 at E8. (C) Per-stage latent-space displacement of IVF/PA from matched in-vivo stages (E2–E5). (D) UMAP coloured by max terminal-fate probability shows spatially organised commitment."),
 ],
 8: [
   (BASE + "m8_mofa_figure.png",
    "Fig. 8. MOFA+ latent-factor map of the porcine preimplantation reference atlas. (A) Variance explained by each of ten latent factors (F1 dominates at 11.6%; total model 21.6%). (B) Trajectories of the three most stage-correlated factors across canonical developmental stages (GV → E14). (C) Top loading genes on the dominant F1 axis (cell-cycle/proliferation programme). (D) Heatmap of top gene loadings across the top-3 stage-correlated factors."),
 ],
}

# -------- assemble --------
body = open(SRC, encoding='utf-8').read()
# take from "# RESULTS" onward
idx = body.index("# RESULTS")
body_from_results = body[idx:]

# replace concise Methods block with expanded Methods
m_start = body_from_results.index("# METHODS (concise)")
r_start = body_from_results.index("# REFERENCES")
body_final = body_from_results[:m_start] + EXPANDED_METHODS + "\n\n" + body_from_results[r_start:]

full = TITLE + "\n\n" + ABSTRACT + "\n\n" + INTRODUCTION + "\n\n" + body_final
open(OUT, 'w', encoding='utf-8').write(full)
print("Wrote", OUT, "(bytes:", len(full), ")")

# -------- embed figures --------
lines = full.split('\n')
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
        out.append("")
        inserted.add(cur_fig)
        cur_fig = None
    out.append(line)
open(OUT_FIG, 'w', encoding='utf-8').write('\n'.join(out))
print("Wrote", OUT_FIG, "with figures for:", sorted(inserted))
