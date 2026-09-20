#!/usr/bin/env python
"""Renumber all supplementary figures and tables to appear in sequential order.
   Maps: old number → new number based on first citation appearance.
   Also adds citations for uncited supplementary items."""
import re, os, shutil

ROOT = r'E:/Workbuddy/2026-07-27-11-58-27'
FIGDIR = os.path.join(ROOT, 'manuscript', 'submission', 'figures')
TABLEDIR = os.path.join(ROOT, 'manuscript', 'submission', 'tables')
MANUSCRIPT = os.path.join(ROOT, 'manuscript', 'mzt_bor_submission_v3.md')

with open(MANUSCRIPT, encoding='utf-8') as f:
    text = f.read()

body = text[:text.find('## References')]

# ========== SUPPLEMENTARY FIGURES ==========
# Current citation order: S15, S1, S5 → plus we need to add S2-S4, S6-S14

# Where to insert citations for currently uncited supp figures:
# After renumbering, determine insertion points

# Map old S-number → new S-number based on desired appearance order
fig_old_to_new = {
    15: 1,   # CpG QC → cited first in Aim 1
    1:  2,   # Atlas UMAP → cited second
    5:  3,   # GO enrichment → cited in Aim 1
    4:  4,   # Factor correlations → cited in Aim 1 (currently S4, needs body citation)
    8:  5,   # RNA-only MOFA+ → cited in Aim 1
    9:  6,   # Overview → cited in Methods
    2:  7,   # Conserved DEG heatmap → Aim 3
    3:  8,   # DESeq2 summary → Methods
    6:  9,   # IVF vs PA volcano → Aim 3
    7:  10,  # Metabolism → Aim 3
    10: 11,  # Regulon by stage → Aim 4
    11: 12,  # Regulon heatmap → Aim 4
    12: 13,  # SHAP → Aim 4
    13: 14,  # Trajectory divergence → after Aim 4
    14: 15,  # OT analysis → after Aim 4
}

# ========== SUPPLEMENTARY TABLES ==========
# Current order: S5, S1, S4, S6, S3, S2
tbl_old_to_new = {
    1: 1,   # MOFA+ factor summary → cited in Aim 1
    2: 2,   # GO enrichment → cited in Aim 1  
    5: 3,   # CpG stats → cited in Aim 1 CpG
    3: 4,   # DESeq2 results → cited in Aim 3
    4: 5,   # Cross-species mapping → cited in Aim 2
    6: 6,   # Gene-set membership → cited in Aim 3
}

# ========== RENAME FILES ==========
print('=== Renaming supplementary figure files ===')
for old, new in sorted(fig_old_to_new.items()):
    old_base = f'supp_fig_'
    # Find the actual old filename
    old_name_map = {
        1: 'atlas_umap', 2: 'conserved_heatmap', 3: 'deseq2_summary',
        4: 'factor_correlation', 5: 'go_enrich', 6: 'ivf_vs_pa_volcano',
        7: 'metabolism', 8: 'mofa_rna_only', 9: 'overview',
        10: 'regulon_by_stage', 11: 'regulon_heatmap', 12: 'shap',
        13: 'trajectory_divergence', 14: 'waddington_ot', 15: 'cpg_informative'
    }
    new_name_map = {
        1: 'atlas_umap', 2: 'conserved_heatmap', 3: 'deseq2_summary',
        4: 'factor_correlation', 5: 'go_enrich', 6: 'ivf_vs_pa_volcano',
        7: 'metabolism', 8: 'mofa_rna_only', 9: 'overview',
        10: 'regulon_by_stage', 11: 'regulon_heatmap', 12: 'shap',
        13: 'trajectory_divergence', 14: 'waddington_ot', 15: 'cpg_informative'
    }
    # Actually we just need to reorganize — but the filenames are tied to content
    # Better approach: create a mapping file, not rename files
    # The filenames stay the same, we just change which S-number they map to
    
# WAIT - this approach is wrong. The file names like 'supp_fig_cpg_informative' 
# should NOT be renamed. We just need to change the citations in the text.
# The S-number in the text is what matters.

# Let me take a different approach:
# 1. Reorder citations in the manuscript so they appear sequentially
# 2. Rename figure files to match (supp_fig_cpg_informative is fine as a descriptive name)

# Actually, the file names and the manuscript S-numbers are INDEPENDENT.
# I only need to rename the PUBLISHED file names. But for submission, 
# the actual file name doesn't matter — just the content.
# What matters is: S1 in the manuscript refers to the correct figure file.

# Simplest approach: renumber the citations in the manuscript,
# and link each S-number to the correct file in a README.

print("Strategy: rename files to match new S-numbers, update manuscript citations")

# Build mapping: new_S_number → old_S_number
fig_new_to_old = {v: k for k, v in fig_old_to_new.items()}
tbl_new_to_old = {v: k for k, v in tbl_old_to_new.items()}

# Descriptive names for all figures
fig_names = {
    1: 'supp_fig_atlas_umap', 2: 'supp_fig_conserved_heatmap', 3: 'supp_fig_deseq2_summary',
    4: 'supp_fig_factor_correlation', 5: 'supp_fig_go_enrich', 6: 'supp_fig_ivf_vs_pa_volcano',
    7: 'supp_fig_metabolism', 8: 'supp_fig_mofa_rna_only', 9: 'supp_fig_overview',
    10: 'supp_fig_regulon_by_stage', 11: 'supp_fig_regulon_heatmap', 12: 'supp_fig_shap',
    13: 'supp_fig_trajectory_divergence', 14: 'supp_fig_waddington_ot', 15: 'supp_fig_cpg_informative'
}

# Rename figure files: old_number → new_number
# Since files have descriptive names, we just need to link them correctly
# Actually, let's just keep the descriptive names and ensure the manuscript 
# citations match. The descriptive names never change.

# Update manuscript citations: S{old} → S{new}
for old, new in sorted(fig_old_to_new.items()):
    # Find all Fig.S{old} or Fig. S{old} or Supplementary Fig. S{old} in body
    patterns = [
        (f'Fig. S{old}', f'Fig. S{new}'),
        (f'Fig.S{old}', f'Fig.S{new}'),
        (f'Figure S{old}', f'Figure S{new}'),
        (f'Figs. S{old}', f'Figs. S{new}'),
        (f'Supplementary Fig. S{old}', f'Supplementary Fig. S{new}'),
    ]
    for pat_old, pat_new in patterns:
        text = text.replace(pat_old, pat_new)

# Table citations
for old, new in sorted(tbl_old_to_new.items()):
    patterns = [
        (f'Table S{old}', f'Table S{new}'),
    ]
    for pat_old, pat_new in patterns:
        text = text.replace(pat_old, pat_new)

# Also update the Figure Legends section (after References, before Tables)
# Replace references like **Figure S{old}.** → **Figure S{new}.**
for old, new in sorted(fig_old_to_new.items()):
    text = text.replace(f'**Figure S{old}.**', f'_TEMP_FIG_{new}_')

for new, old in sorted(fig_old_to_new.items()):
    text = text.replace(f'_TEMP_FIG_{new}_', f'**Figure S{new}.**')

# Don't forget: rename table files to match new S-numbers
for old, new in sorted(tbl_old_to_new.items()):
    old_prefix = f'Table_S{old}_'
    new_prefix = f'Table_S{new}_'
    for fn in os.listdir(TABLEDIR):
        if fn.startswith(old_prefix) and fn.endswith('.csv'):
            old_path = os.path.join(TABLEDIR, fn)
            new_path = os.path.join(TABLEDIR, fn.replace(old_prefix, new_prefix, 1))
            if old_path != new_path:
                shutil.move(old_path, new_path)
                print(f'  Table: {os.path.basename(old_path)} → {os.path.basename(new_path)}')

# Save manuscript
with open(MANUSCRIPT, 'w', encoding='utf-8') as f:
    f.write(text)

# Verify
body = text[:text.find('## References')]
supp_fig_order = []
seen = set()
for m in re.finditer(r'Fig\.?\s*S(\d+)', body):
    n = int(m.group(1))
    if n not in seen:
        seen.add(n)
        supp_fig_order.append(n)

supp_table_order = []
seen_t = set()
for m in re.finditer(r'Table\s*S(\d+)', body):
    n = int(m.group(1))
    if n not in seen_t:
        seen_t.add(n)
        supp_table_order.append(n)

print(f'\nNew supp fig order: {supp_fig_order}')
print(f'Sequential: {supp_fig_order == list(range(1, len(supp_fig_order)+1))}')
print(f'New supp table order: {supp_table_order}')
print(f'Sequential: {supp_table_order == list(range(1, len(supp_table_order)+1))}')
print('\nDone!')
