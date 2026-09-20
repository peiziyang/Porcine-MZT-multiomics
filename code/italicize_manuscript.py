#!/usr/bin/env python
"""给稿件加斜体标注：in vivo/in vitro + 基因名（大写符号）。避免重复处理已斜体内容。"""
import re

MD = r"E:/Workbuddy/2026-07-27-11-58-27/manuscript/mzt_bor_submission_v3.md"

with open(MD, encoding='utf-8') as f:
    text = f.read()

original = text

def italicize_phrase(text, phrase):
    """把短语（如 'in vivo'）加斜体，避免重复。"""
    pat = re.compile(r'(?<!\*)\b' + re.escape(phrase) + r'\b(?!\*)')
    return pat.sub(lambda m: '*' + m.group(0) + '*', text)

# 1. 拉丁语斜体
for phrase in ['in vivo', 'in vitro']:
    text = italicize_phrase(text, phrase)

# 2. 基因名斜体（大写符号，精确词边界）
genes = [
    'DNMT1', 'ZP3', 'ZP4', 'GDF9', 'RARRES1',
    'ANXA1', 'WWTR1', 'RHOD', 'CDC42EP1',
    'MDH1', 'ALDH2', 'PGAM1',
    'IDH2', 'PKM', 'EXOSC9', 'GNL3',
    'E2F1', 'TET3', 'UHRF1', 'FOXO1', 'POU5F1', 'TFAP2C', 'GATA4', 'ATF3',
    'EIF4G2',
]
for g in genes:
    pat = re.compile(r'(?<!\*)\b' + re.escape(g) + r'\b(?!\*)')
    text = pat.sub(lambda m: '*' + m.group(0) + '*', text)

# 3. 参考文献里的小鼠基因名 Dnmt1（ref 28）
text = text.replace('the Dnmt1 gene', 'the *Dnmt1* gene')

# 统计变化
added = text.count('*') - original.count('*')
print(f"斜体星号增量: {added} (应≈ (in vivo/in vitro 次数 + 基因名次数) × 2)")

# 验证：检查几个关键位置
checks = [
    'in vivo, by *in vitro* fertilization',
    '*in vivo* development',
    '*DNMT1*, *ZP3*',
    '*DNMT1* (maintenance DNA methyltransferase)',
]
print("\n=== 验证 ===")
for c in checks:
    print(f'  {"✓" if c in text else "✗"} {c}')

# 检查是否误伤了已有斜体（物种名/期刊名不应重复加星）
print(f'\n  *Sus scrofa* 仍正确: {"*Sus scrofa*" in text}')
print(f'  **Development** 期刊名: {"*Development*" in text}')
# 检查没有 *** 三连星（说明没重复）
print(f'  三连星 *** 数量（应为0）: {text.count("***")}')
# 检查没有 ** in vivo** 这种错误
print(f'  ** in vivo 误匹配: {text.count("**in vivo")}')

with open(MD, 'w', encoding='utf-8') as f:
    f.write(text)
print("\n已写回文件")
