#!/usr/bin/env python
"""Comprehensive: renumber ALL supplementary items for sequential citation order."""
import re, os, shutil

ROOT = r'E:/Workbuddy/2026-07-27-11-58-27'
MD = os.path.join(ROOT, 'manuscript', 'mzt_bor_submission_v3.md')
FIGDIR = os.path.join(ROOT, 'manuscript', 'submission', 'figures')
TABLEDIR = os.path.join(ROOT, 'manuscript', 'submission', 'tables')

with open(MD, encoding='utf-8') as f:
    text = f.read()

ref_pos = text.find('## References')
body = text[:ref_pos]
rest = text[ref_pos:]

# ===== SUPPLEMENTARY FIGURES =====
# First, find all S-figure citations in body and their order
import re
fig_seen = set()
fig_order_old = []
for m in re.finditer(r'Fig\.?\s*S(\d+)', body):
    n = int(m.group(1))
    if n not in fig_seen:
        fig_seen.add(n)
        fig_order_old.append(n)

print(f'Current supp fig citation order: {fig_order_old}')
print(f'Cited: {sorted(fig_seen)}')
uncited = sorted(set(range(1,16)) - fig_seen)
print(f'Uncited: {uncited}')

# Build mapping: old → new based on first-appearance order
fig_map = {}
for new, old in enumerate(fig_order_old, 1):
    fig_map[old] = new
# Uncited items get assigned after cited ones
for old in uncited:
    fig_map[old] = len(fig_map) + 1

print(f'Fig mapping: {dict(sorted(fig_map.items()))}')

# ===== SUPPLEMENTARY TABLES =====
tbl_seen = set()
tbl_order_old = []
for m in re.finditer(r'Table\s*S(\d+)', body):
    n = int(m.group(1))
    if n not in tbl_seen:
        tbl_seen.add(n)
        tbl_order_old.append(n)

print(f'\nCurrent supp table citation order: {tbl_order_old}')
uncited_t = sorted(set(range(1,7)) - tbl_seen)
print(f'Uncited tables: {uncited_t}')

tbl_map = {}
for new, old in enumerate(tbl_order_old, 1):
    tbl_map[old] = new
for old in uncited_t:
    tbl_map[old] = len(tbl_map) + 1

print(f'Table mapping: {dict(sorted(tbl_map.items()))}')

# ===== REPLACE CITATIONS (two-pass to avoid partial match conflicts) =====
# FIGURES — pass 1: replace old → TEMP
temp_map = {}
for old, new in sorted(fig_map.items(), reverse=True):  # descending old to avoid conflicts
    if old != new:
        tag = f'__F{old:02d}__'
        temp_map[tag] = f'S{new}'
        # Replace all variants
        body = re.sub(rf'\bFig\.\s*S{old}\b', tag, body)
        body = re.sub(rf'\bFig\.\s*S{old}\b', tag, body)  # match without space too
        body = re.sub(rf'\bFigs\.\s*S{old}\b', tag, body)
        body = re.sub(rf'\bFigure\s+S{old}\b', tag, body)
        body = re.sub(rf'\bSupplementary\s+Figs?\.?\s*S{old}\b', tag, body)
        body = re.sub(rf'\*\*Figure\s+S{old}\.\*\*', tag + '.**', body)
# FIGURES — pass 2: replace TEMP → new
for tag, new_code in temp_map.items():
    body = body.replace(tag, f'Fig. {new_code}')

# TABLES
tbl_temp = {}
for old, new in sorted(tbl_map.items(), reverse=True):
    if old != new:
        tag = f'__T{old:02d}__'
        tbl_temp[tag] = f'S{new}'
        body = re.sub(rf'\bTable\s+S{old}\b', tag, body)
for tag, new_code in tbl_temp.items():
    body = body.replace(tag, f'Table {new_code}')

# ===== Also update Figure Legends section (in rest) =====
for old, new in sorted(fig_map.items(), reverse=True):
    if old != new:
        tag = f'__F{old:02d}__'
        rest = rest.replace(f'**Figure S{old}.**', f'__F{old:02d}__.**')
for old, new in sorted(fig_map.items(), reverse=True):
    if old != new:
        tag = f'__F{old:02d}__'
        rest = rest.replace(f'{tag}.**', f'**Figure S{new}.**')

# ===== Save manuscript =====
new_text = body + '\n' + rest
with open(MD, 'w', encoding='utf-8') as f:
    f.write(new_text)

# ===== RENAME FIGURE FILES =====
fig_names = {
    1: 'atlas_umap', 2: 'conserved_heatmap', 3: 'deseq2_summary',
    4: 'factor_correlation', 5: 'go_enrich', 6: 'ivf_vs_pa_volcano',
    7: 'metabolism', 8: 'mofa_rna_only', 9: 'overview',
    10: 'regulon_by_stage', 11: 'regulon_heatmap', 12: 'shap',
    13: 'trajectory_divergence', 14: 'waddington_ot', 15: 'cpg_informative'
}

# Step 1: move old files to temp names
print('\n=== Renaming figure files ===')
for old, new in sorted(fig_map.items(), reverse=True):  # old-to-new but backwards
    if old == new:
        continue
    old_desc = fig_names.get(old, f'unknown{old}')
    new_desc = fig_names.get(new, f'unknown{new}')
    for ext in ['.png', '.svg']:
        old_file = os.path.join(FIGDIR, f'supp_fig_{old_desc}{ext}')
        tmp_file = os.path.join(FIGDIR, f'_tmp_rename_{old}{ext}')
        if os.path.exists(old_file):
            shutil.move(old_file, tmp_file)
            print(f'  move: supp_fig_{old_desc}{ext} → tmp')

# Step 2: move temp files to new names
for old, new in sorted(fig_map.items()):
    if old == new:
        continue
    new_desc = fig_names.get(new, f'unknown{new}')
    for ext in ['.png', '.svg']:
        tmp_file = os.path.join(FIGDIR, f'_tmp_rename_{old}{ext}')
        new_file = os.path.join(FIGDIR, f'supp_fig_{new_desc}{ext}')
        if os.path.exists(tmp_file):
            shutil.move(tmp_file, new_file)
            print(f'  rename: tmp → supp_fig_{new_desc}{ext}')

# ===== RENAME TABLE FILES =====
print('\n=== Renaming table files ===')
tbl_names = {
    1: 'mofa_factor_summary', 2: 'go_enrichment', 3: 'deseq2_results',
    4: 'cross_species_mapping', 5: 'cpg_stats', 6: 'gene_set_membership'
}

# Move all to temp first (descending to avoid conflicts)
for old, new in sorted(tbl_map.items(), reverse=True):
    if old == new:
        continue
    old_desc = tbl_names[old]
    old_file = os.path.join(TABLEDIR, f'Table_S{old}_{old_desc}.csv')
    tmp_file = os.path.join(TABLEDIR, f'_tmp_tbl_{old}.csv')
    if os.path.exists(old_file):
        shutil.move(old_file, tmp_file)
        print(f'  move: Table_S{old}_{old_desc}.csv → tmp')

# Move temp to new names
for old, new in sorted(tbl_map.items()):
    if old == new:
        continue
    new_desc = tbl_names[new]
    tmp_file = os.path.join(TABLEDIR, f'_tmp_tbl_{old}.csv')
    new_file = os.path.join(TABLEDIR, f'Table_S{new}_{new_desc}.csv')
    # But wait — new file might already exist (if new=6 and old wasn't moved yet)
    if os.path.exists(tmp_file):
        if os.path.exists(new_file):
            os.remove(new_file)  # avoid conflict
        shutil.move(tmp_file, new_file)
        print(f'  rename: tmp → Table_S{new}_{new_desc}.csv')

# ===== VERIFY =====
with open(MD, encoding='utf-8') as f:
    new_body = f.read()
new_body = new_body[:new_body.find('## References')]

fig_final = []
seen_f = set()
for m in re.finditer(r'Fig\.?\s*S(\d+)', new_body):
    n = int(m.group(1))
    if n not in seen_f:
        seen_f.add(n)
        fig_final.append(n)

tbl_final = []
seen_t2 = set()
for m in re.finditer(r'Table\s*S(\d+)', new_body):
    n = int(m.group(1))
    if n not in seen_t2:
        seen_t2.add(n)
        tbl_final.append(n)

print(f'\n=== FINAL VERIFICATION ===')
print(f'Supp figs: {fig_final}')
print(f'Sequential: {fig_final == list(range(1, len(fig_final)+1))}')
print(f'Supp tbls: {tbl_final}')
print(f'Sequential: {tbl_final == list(range(1, len(tbl_final)+1))}')
print(f'Cited figs: {len(seen_f)}/15  Cited tbls: {len(seen_t2)}/6')
