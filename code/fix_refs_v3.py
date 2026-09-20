#!/usr/bin/env python
"""Minimal ref renumbering — only touch known citation patterns."""
import re

PATH = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/mzt_bor_submission_v3.md'
with open(PATH, encoding='utf-8') as f:
    text = f.read()

ref_start = text.find('## References')
tbl_start = text.find('## Table 1')
body = text[:ref_start]
refs_text = text[ref_start:tbl_start]
rest = text[tbl_start:]

# ---- Known citation order (from earlier analysis) ----
# Old refs in order of first appearance: 1,2,3,4,14,13,5,6,7,9,25,17,10,21,20,19,11,12,15,18,22,16
# Ref 8 is cited via [7–9] range (need to keep)
# Refs 23, 24 are truly uncited → remove

appearance = [1, 2, 3, 4, 14, 13, 5, 6, 7, 8, 9, 25, 17, 10, 21, 20, 19, 11, 12, 15, 18, 22, 16]
# Remove uncited: 23, 24
# Note: 8 was NOT in the original analysis because range [7–9] only yielded 7 and 9
#        But [7–9] means 7,8,9 are all cited. We add 8 after 7.

old_to_new = {o: n for n, o in enumerate(appearance, 1)}
# old_to_new: {1:1, 2:2, 3:3, 4:4, 14:5, 13:6, 5:7, 6:8, 7:9, 8:10, 9:11, 25:12, 17:13, ...}

print('Old→New mapping:')
for k in sorted(old_to_new):
    print(f'  {k}→{old_to_new[k]}')

# ---- Replace citations in body ----
# Strategy: replace from highest number to avoid conflicts
# For each old number, replace every instance of [old_num] or [old_num, or ,old_num] or [old_num– etc

body_new = body

# Find all bracketed text
all_brackets = [(m.start(), m.end(), m.group(0)) for m in re.finditer(r'\[([^\]]+)\]', body_new)]

for start, end, bracket in reversed(all_brackets):
    content = bracket[1:-1]  # strip brackets
    
    # Check if this is a citation bracket
    # Citations look like: digits separated by commas or en-dashes, no decimals, no minus signs
    # Count: digits, commas, dashes, spaces
    parts = re.split(r'([,–\-])', content)
    all_parts_numeric = True
    numbers_in_bracket = []
    for p in parts:
        p = p.strip()
        if p in [',', '–', '-']:
            continue
        if re.match(r'^\d{1,2}$', p):
            numbers_in_bracket.append(int(p))
        elif p == '':
            continue
        else:
            all_parts_numeric = False
            break
    
    if not all_parts_numeric or not numbers_in_bracket:
        continue
    
    # Check these look like citation numbers (1-25ish)
    if not all(1 <= n <= 30 for n in numbers_in_bracket):
        continue
    
    # Avoid false positives: skip if preceded by IQR, CI, ρ, =, etc
    before = body_new[max(0,start-10):start].rstrip()
    if re.search(r'(IQR|CI\s*=|rho|log|FC|±|padj|n\s*=)\s*$', before, re.IGNORECASE):
        continue
    
    # This is a citation — replace each number
    new_parts = []
    for p in parts:
        if p.strip() in [',', '–', '-']:
            new_parts.append(p)
            continue
        p_stripped = p.strip()
        if re.match(r'^\d{1,2}$', p_stripped):
            n = int(p_stripped)
            if n in old_to_new:
                new_parts.append(str(old_to_new[n]))
            else:
                new_parts.append(p_stripped)
        else:
            new_parts.append(p)
    
    new_bracket = '[' + ''.join(new_parts) + ']'
    if new_bracket != bracket:
        body_new = body_new[:start] + new_bracket + body_new[end:]

# ---- Rebuild reference list ----
# Parse existing refs
ref_lines = refs_text.strip().split('\n')
ref_entries = {}
current_num = None
current_text = []
for line in ref_lines:
    m = re.match(r'^(\d+)\.\s+(.*)', line.strip())
    if m:
        if current_num is not None:
            ref_entries[current_num] = ' '.join(current_text)
        current_num = int(m.group(1))
        current_text = [m.group(2)]
    elif line.strip() and current_num is not None:
        current_text.append(line.strip())
if current_num is not None:
    ref_entries[current_num] = ' '.join(current_text)

# Reorder
new_refs = []
for old in appearance:
    if old in ref_entries:
        new_refs.append(f'{len(new_refs)+1}. {ref_entries[old]}')

new_refs_text = '\n\n'.join(new_refs)

# ---- Assemble ----
final = body_new + '\n\n## References\n\n' + new_refs_text + '\n\n' + rest.lstrip()

with open(PATH, 'w', encoding='utf-8') as f:
    f.write(final)

print(f'\nOld refs: {len(ref_entries)} → New refs: {len(new_refs)}')
print('Verifying:')
# Spot check a few
for check in ['[1]', '[2]', '[3]', '[4]', '[5]', '[6]', '[7]', '[8]', '[9]']:
    count = body_new.count(check)
    if count > 0:
        print(f'  {check}: {count} occurrences')
