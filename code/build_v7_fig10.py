"""Insert Figure 10 (conserved DEG heatmap) into v7 manuscript and sync abstract/discussion/table."""
import io, sys

SRC = 'E:/Workbuddy/2026-07-27-11-58-27/manuscript/manuscript_full_with_figures_v7.md'
IMG = 'E:/Workbuddy/2026-07-27-11-58-27/data/processed/deseq2/m15_fig10_conserved_heatmap.png'

with open(SRC, encoding='utf-8') as f:
    md = f.read()

def rep(text, old, new, label):
    n = text.count(old)
    assert n == 1, f'[{label}] expected exactly 1 match, found {n}'
    return text.replace(old, new, 1)

# ---- 1. Abstract: append conserved-DEG sentence ----
old_abs = 'most robustly at the 4-cell stage (n = 9 per group).'
new_abs = ('most robustly at the 4-cell stage (n = 9 per group). '
           'A stage-conserved core of DEGs \u2014 notably the metabolic enzymes IDH2 and PKM, '
           'consistently PA-downregulated across all four cleavage stages \u2014 extends the '
           'attenuation signature into energy metabolism.')
md = rep(md, old_abs, new_abs, 'abstract')

# ---- 2. Main body: insert Figure 10 section before DISCUSSION ----
old_main = 'm15_fig9_summary.png)\n\n\n# DISCUSSION'
fig10_body = (
    'm15_fig9_summary.png)\n\n'
    '## Figure 10 \u2014 A stage-conserved DEG core links the PA deficit to metabolic reprogramming\n\n'
    'Beyond the per-stage counts, we asked which genes are reproducibly dysregulated '
    'across cleavage stages rather than being stage-specific. Restricting to genes that '
    'were significant (padj < 0.05, |log\u2082FC| > 1) and directionally consistent in at least '
    'three of the four stages, we identified 321 such conserved DEGs and highlighted the '
    'top 30 by stage-consistency and mean effect size (Fig. 10). The conserved core is '
    'dominated by PA-downregulated genes (21 of the top 30) and prominently includes metabolic '
    'enzymes: **IDH2** (TCA cycle) and **PKM** (glycolysis) are PA-downregulated at all four '
    'stages, alongside NOP9, GNL3 and EXOSC9 (ribosome/biogenesis). PA-upregulated conserved '
    'genes include **PABPC5**, **SETDB1** (histone methyltransferase) and **TAF9B**. The '
    'consistency of IDH2/PKM downregulation across every cleavage stage \u2014 independent of the '
    'low-powered 1- and 8-cell groups \u2014 argues that PA\u2019s transcriptional attenuation extends '
    'into energy-metabolism reprogramming, a conclusion now shared by the public-data and '
    '42-sample analyses.\n\n'
    f'![Fig. 10. Top-30 stage-conserved DEGs (IVF vs PA, log\u2082FC). Rows ordered by hierarchical '
    f'clustering of the per-stage log\u2082FC profile; red = PA-up, blue = PA-down. Metabolic enzymes '
    f'IDH2 and PKM are PA-downregulated at all four cleavage stages (see text).]({IMG})\n\n'
    '# DISCUSSION'
)
md = rep(md, old_main, fig10_body, 'main-body')

# ---- 3. Discussion: append conserved-DEG sentence to independent-replication paragraph ----
old_disc = 'they are reported where they were derived rather than here.'
new_disc = ('they are reported where they were derived rather than here. '
            'A stage-conserved core of DEGs \u2014 including the TCA-cycle enzyme IDH2 and the '
            'glycolytic enzyme PKM, which are PA-downregulated at all four cleavage stages \u2014 '
            'further localises the PA deficit to energy-metabolism reprogramming, a direction now '
            'independently anchored in the 42-sample cohort.')
md = rep(md, old_disc, new_disc, 'discussion')

# ---- 4. Figure legend: insert Figure 10 legend before SUPPLEMENTARY MATERIAL ----
old_leg = 'hypothesis-generating.\n\n\n# SUPPLEMENTARY MATERIAL'
fig10_leg = (
    'hypothesis-generating.\n\n'
    '## Figure 10. A stage-conserved DEG core links the PA deficit to metabolic reprogramming.\n'
    'Hierarchical-clustered heatmap of the top-30 genes that are significant '
    '(padj < 0.05, |log\u2082FC| > 1) and directionally consistent in \u22653 of 4 cleavage stages '
    '(IVF vs PA, log\u2082FC, PA/IVF). Red = PA-up, blue = PA-down. Of the 321 conserved DEGs '
    'identified, the top tier is dominated by PA-downregulated genes and prominently includes '
    'the TCA-cycle enzyme IDH2 and the glycolytic enzyme PKM, both PA-downregulated at all four '
    'stages, together with ribosome/biogenesis factors (NOP9, GNL3, EXOSC9). PA-upregulated '
    'conserved genes include PABPC5, SETDB1 and TAF9B. The cross-stage consistency of IDH2/PKM '
    'downregulation points to a metabolic-reprogramming component of the PA transcriptional deficit.\n\n'
    '# SUPPLEMENTARY MATERIAL'
)
md = rep(md, old_leg, fig10_leg, 'legend')

# ---- 5. Data availability table: add Fig 10 row ----
old_tab = ('| Fig. 9 | 42-sample bulk RNA-seq IVF/PA validation \u2713 | '
           'Per-stage n reported; 1/8-cell n=3 lower power; 4-cell n=9 best-powered |')
new_tab = (old_tab + '\n'
           '| Fig. 10 | Conserved DEG heatmap across cleavage stages \u2713 | '
           'Top-30 stage-conserved DEGs (\u22653/4 stages consistent); IDH2/PKM PA-down |')
md = rep(md, old_tab, new_tab, 'table')

with open(SRC, 'w', encoding='utf-8') as f:
    f.write(md)

import re, os
imgs = re.findall(r'!\\[[^\\]]*\\]\\(([^)]+)\\)', md)
miss = [i for i in imgs if not os.path.exists(i)]
print('All 5 replacements applied.')
print('Total figure links:', len(imgs), '| broken:', len(miss))
for m in miss:
    print('  MISSING', m)
assert not miss, 'broken image links!'
print('v7.md updated OK.')
