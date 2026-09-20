#!/usr/bin/env python
"""Robust: renumber all refs in order of first appearance, remove truly uncited ones."""
import re

PATH = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/mzt_bor_submission_v3.md'
with open(PATH, encoding='utf-8') as f:
    text = f.read()

ref_start = text.find('## References')
tbl1_start = text.find('## Table 1')
body = text[:ref_start]
refs_text = text[ref_start:tbl1_start]
rest = text[tbl1_start:]

# ---- Phase 1: Determine which refs are cited ----
# Find all citation groups: [digit(s)] that are citation markers.
# Patterns: [1], [1,2], [3,4,14], [7–9], [7-9], [1,3–5]
# Avoid: [−0.58, −0.15] (has minus), [0.4, 0.8] (has decimal)

all_cited_old = set()
cit_pattern = re.compile(r'\[((?:\d{1,2}(?:[,–-]\d{1,2})*))\]')
for m in cit_pattern.finditer(body):
    before = body[max(0,m.start()-8):m.start()]
    # Skip if preceded by CI, IQR, =, or has minus/decimal nearby
    if re.search(r'(CI|IQR|rho|log|±|[a-z])\s*$', before, re.IGNORECASE):
        continue
    raw = m.group(0)
    nums = [int(x) for x in re.findall(r'\d+', m.group(1))]
    if all(1 <= n <= 30 for n in nums):
        for n in nums:
            all_cited_old.add(n)

# Also check ranges — include implicit middle numbers
# [7–9] cites 7,8,9 but only 7 and 9 appear explicitly
for m in re.finditer(r'\[(\d{1,2})[–-](\d{1,2})\]', body):
    a, b = int(m.group(1)), int(m.group(2))
    if a <= 30 and b <= 30:
        for n in range(a, b+1):
            all_cited_old.add(n)

print(f'Cited refs (old numbers): {sorted(all_cited_old)}')
uncited = set(range(1, 26)) - all_cited_old
print(f'Uncited (will be removed): {sorted(uncited)}')

# ---- Phase 2: Build first-appearance order ----
seen = set()
order = []
for m in cit_pattern.finditer(body):
    before = body[max(0,m.start()-8):m.start()]
    if re.search(r'(CI|IQR|rho|log|±|[a-z])\s*$', before, re.IGNORECASE):
        continue
    nums = [int(x) for x in re.findall(r'\d+', m.group(1))]
    if all(1 <= n <= 30 for n in nums):
        for n in nums:
            if n not in seen:
                seen.add(n)
                order.append(n)

# For ranges, insert middle numbers at the range's position
# Find [7–9] style ranges and insert 8 after 7
expanded_order = []
for n in order:
    expanded_order.append(n)
    # Check if n is start of a range in the body
    for m in re.finditer(rf'\[{n}[–-](\d{{1,2}})\]', body):
        end = int(m.group(1))
        for mid in range(n+1, end+1):
            if mid not in expanded_order:
                expanded_order.append(mid)

print(f'Appearance order (expanded): {expanded_order}')

# ---- Phase 3: Build mapping ----
old_to_new = {}
new_to_old = {}
for new_num, old_num in enumerate(expanded_order, 1):
    if old_num in all_cited_old:  # skip uncited ones
        old_to_new[old_num] = new_num
        new_to_old[new_num] = old_num

print(f'Mapping: {dict(sorted(old_to_new.items()))}')

# ---- Phase 4: Replace all citation markers in body ----
def replace_citation(match):
    content = match.group(1)
    # Split on commas and dashes
    tokens = re.split(r'([,–-])', content)
    result = []
    for tok in tokens:
        nums = re.findall(r'\d+', tok)
        for n_str in sorted(nums, key=lambda x: -len(x)):
            n = int(n_str)
            if n in old_to_new:
                tok = tok.replace(n_str, str(old_to_new[n]))
        result.append(tok)
    return '[' + ''.join(result) + ']'

body_new = body
# Replace from right to left to preserve positions
replacements = []
for m in cit_pattern.finditer(body):
    before = body[max(0,m.start()-8):m.start()]
    if re.search(r'(CI|IQR|rho|log|±|[a-z])\s*$', before, re.IGNORECASE):
        continue
    nums = [int(x) for x in re.findall(r'\d+', m.group(1))]
    if all(1 <= n <= 30 for n in nums):
        replacements.append((m.start(), m.end(), m.group(0)))

for start, end, old_str in reversed(replacements):
    # Find old_str in body_new starting from start (with offset due to prior edits)
    pos = body_new.find(old_str, max(0, start - 100), start + 100)
    if pos >= 0:
        new_str = cit_pattern.sub(replace_citation, old_str)
        body_new = body_new[:pos] + new_str + body_new[pos+len(old_str):]

# ---- Phase 5: Reorder references ----
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
    elif line.strip():
        current_text.append(line.strip())
if current_num is not None:
    ref_entries[current_num] = ' '.join(current_text)

# Build new ref list in order
new_refs = []
for new_num in range(1, len(expanded_order) + 1):
    old_num = new_to_old.get(new_num)
    if old_num and old_num in ref_entries:
        entry = ref_entries[old_num]
        # Fix any embedded old ref numbers in the entry text
        for o_n in sorted(old_to_new.keys(), reverse=True):
            if o_n != old_num:  # don't replace the ref's own number
                entry = entry.replace(f'[{o_n}]', f'[{old_to_new[o_n]}]')
        new_refs.append(f'{new_num}. {entry}')

new_refs_text = '\n\n'.join(new_refs)

# ---- Phase 6: Assemble final text ----
final = body_new + '\n\n' + new_refs_text + '\n\n' + rest

with open(PATH, 'w', encoding='utf-8') as f:
    f.write(final)

print(f'\nOld: {len(ref_entries)} refs, New: {len(new_refs)} refs')
print('Saved!')
