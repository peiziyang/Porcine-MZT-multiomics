#!/usr/bin/env python
"""Fix: reference renumbering, remove uncited refs, check figure/table order."""
import re

PATH = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/mzt_bor_submission_v3.md'
with open(PATH, encoding='utf-8') as f:
    text = f.read()

# Split
ref_start = text.find('## References')
tbl1_start = text.find('## Table 1')
body = text[:ref_start]
refs_section = text[ref_start:tbl1_start]
rest = text[tbl1_start:]

# Find first-appearance order of all cited refs
seen = set()
order = []
for m in re.finditer(r'\[([^\]]+)\]', body):
    content = m.group(1)
    # Skip non-citation brackets: CIs, IQRs, etc
    if '−' in content and not content.replace('−','').replace('.','').replace(',','').replace('-','').replace(' ','').replace('+','').isdigit():
        nums = re.findall(r'\d+', content)
        for n in nums:
            ni = int(n)
            if ni <= 30 and ni not in seen:  # assume ref numbers ≤ 30
                seen.add(ni)
                order.append(ni)

print(f'Refs in appearance order: {order}')
print(f'Uncited: {set(range(1,26)) - set(order)}')

# Build mapping: old_num -> new_num
old_to_new = {}
for new_num, old_num in enumerate(order, 1):
    old_to_new[old_num] = new_num

# Build the replacement function
def replace_ref(match):
    content = match.group(1)
    # Handle ranges like [7–9]
    parts = re.split(r'([,\-–])', content)
    new_parts = []
    for p in parts:
        nums = re.findall(r'\d+', p)
        for n in nums:
            ni = int(n)
            if ni in old_to_new:
                p = p.replace(n, str(old_to_new[ni]), 1) if n in p else p
        # Special: replace each number individually
        new_p = p
        for n in sorted(re.findall(r'\d+', p), key=lambda x: -len(x)):
            ni = int(n)
            if ni in old_to_new:
                new_p = new_p.replace(n, str(old_to_new[ni]))
        new_parts.append(new_p)
    new_content = ''.join(new_parts)
    return '[' + new_content + ']'

# Replace citations in body
body_new = body
for m in re.finditer(r'\[([^\]]+)\]', body):
    old = m.group(0)
    if re.search(r'\d{4}', old):  # years or large numbers - skip
        continue
    new = replace_ref(m)
    if old != new:
        body_new = body_new.replace(old, new, 1)

# Reorder references
ref_lines = refs_section.strip().split('\n')
ref_entries = {}
current_num = None
current_text = []
for line in ref_lines:
    line = line.strip()
    m = re.match(r'^(\d+)\.\s+(.*)', line)
    if m:
        if current_num is not None:
            ref_entries[current_num] = '\n'.join(current_text)
        current_num = int(m.group(1))
        current_text = [m.group(2)]
    elif line and current_num is not None:
        current_text.append(line)
if current_num is not None:
    ref_entries[current_num] = '\n'.join(current_text)

# Remove uncited refs and reorder
new_refs = []
for old_num in order:
    if old_num in ref_entries:
        entry = ref_entries[old_num]
        # Update embedded citation numbers in the entry text if needed
        new_refs.append(f'{len(new_refs)+1}. {entry}')
    else:
        print(f'WARNING: Ref {old_num} not found in reference list!')

# Build final text
new_refs_text = '\n\n'.join(new_refs)
final = body_new + '\n\n' + new_refs_text + '\n\n' + rest.lstrip()

with open(PATH, 'w', encoding='utf-8') as f:
    f.write(final)

print(f'\nOld refs: {len(ref_entries)}, New refs: {len(new_refs)}')
print(f'Mapping: {dict(sorted(old_to_new.items()))}')
print('Done!')
