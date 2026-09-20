#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""补 5 处缺引用 + 全文重编号 + 重建 30 条参考文献列表。"""
import re, shutil

md = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/mzt_bor_submission_v3.md'

with open(md, encoding='utf-8') as f:
    text = f.read()

shutil.copy(md, md + '.bak3.md')

# ===== Step 1: 插入 5 处新引用（占位符 [@@N@@]，避免被重编号正则误伤）=====
insertions = [
    ('and the timely initiation of zygotic transcription.',
     'and the timely initiation of zygotic transcription [@@3@@].'),
    ('during oocyte growth and maturation.',
     'during oocyte growth and maturation [@@4@@].'),
    ('due to similarities in physiology and preimplantation kinetics.',
     'due to similarities in physiology and preimplantation kinetics [@@11@@].'),
    ('consistent with known roles of DNA methylation maintenance machinery in oocyte maturation.',
     'consistent with known roles of DNA methylation maintenance machinery in oocyte maturation [@@27@@].'),
    ('RNA degradation kinetics (e.g., SLAM-seq) during the oocyte-to-embryo transition',
     'RNA degradation kinetics (e.g., SLAM-seq [@@28@@]) during the oocyte-to-embryo transition'),
]
for old, new in insertions:
    if old not in text:
        print(f'  WARN 未找到: {old[:60]}')
    text = text.replace(old, new, 1)
print('Step 1 插入完成')

# ===== Step 2: 重编号所有旧引用标记 =====
old_to_new = {
    1:1, 2:2, 3:5, 4:6, 5:7, 6:8, 7:9, 8:10,
    9:12, 10:13, 11:14, 12:15, 13:16, 14:17, 15:18,
    16:19, 17:20, 18:21, 19:22, 20:23, 21:24, 22:25,
    23:26, 24:29, 25:30,
}

def renumber(m):
    inner = m.group(1)
    nums = [int(x) for x in re.findall(r'\d+', inner)]
    new = [old_to_new.get(n, n) for n in nums]
    if '\u2013' in inner or '-' in inner:
        return f'[{new[0]}\u2013{new[-1]}]'
    return '[' + ','.join(str(n) for n in new) + ']'

text = re.sub(r'\[(\d{1,2}(?:[,\u2013\-]\d{1,2})*)\]', renumber, text)
print('Step 2 重编号完成')

# ===== Step 3: 占位符 → 最终编号 =====
for ph in ['3', '4', '11', '27', '28']:
    text = text.replace(f'[@@{ph}@@]', f'[{ph}]')
print('Step 3 占位符替换完成')

# ===== Step 4: 重建参考文献列表 =====
new_refs = """## References

1. Tadros W, Lipshitz HD. The maternal-to-zygotic transition: a play in two acts. *Development* 2009; 136:3033\u20133042.

2. Lee MT, Bonneau AR, Giraldez AJ. Zygotic genome activation during the maternal-to-zygotic transition. *Annu Rev Cell Dev Biol* 2014; 30:581\u2013613.

3. Eckersley-Maslin MA, Alda-Catalinas C, Reik W. Dynamics of the epigenetic landscape during the maternal-to-zygotic transition. *Nat Rev Mol Cell Biol* 2018; 19:436\u2013450.

4. Zhang K, Smith GW. Maternal control of early embryogenesis in mammals. *Reprod Fertil Dev* 2015; 27:880\u2013896.

5. Surani MA, Barton SC, Norris ML. Development of reconstituted mouse eggs suggests imprinting of the genome during gametogenesis. *Nature* 1984; 308:548\u2013550.

6. McGrath J, Solter D. Completion of mouse embryogenesis requires both the maternal and paternal genomes. *Cell* 1984; 37:179\u2013183.

7. Kono T, Obata Y, Wu Q, et al. Birth of parthenogenetic mice that can develop to adulthood. *Nature* 2004; 428:860\u2013864.

8. Jukam D, Shariati SAM, Skotheim JM. Zygotic genome activation in vertebrates. *Dev Cell* 2017; 42:316\u2013332.

9. Xue Z, Huang K, Cai C, et al. Genetic programs in human and mouse early embryos revealed by single-cell RNA sequencing. *Nature* 2013; 500:593\u2013597.

10. Yan L, Yang M, Guo H, et al. Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. *Nat Struct Mol Biol* 2013; 20:1131\u20131139.

11. Jalali BM, Wasielak-Politowska M. From zygote to blastocyst\u2014molecular aspects of porcine early embryonic development. *Cells* 2026; 15:15.

12. Du ZQ, Liang H, Liu XM, Liu YH, Wang C, Yang CX. Single cell RNA-seq reveals genes vital to in vitro fertilized embryos and parthenotes in pigs. *Sci Rep* 2021; 11:18725.

13. Kong Q, Yang X, Zhang H, et al. Lineage specification and pluripotency revealed by transcriptome analysis from oocyte to blastocyst in pig. *FASEB J* 2020; 34:691\u2013705.

14. Zhi M, Zhang J, Tang Q, et al. Generation and characterization of stable pig pregastrulation epiblast stem cell lines. *Cell Res* 2022; 32:383\u2013400.

15. Alberio R, Sang F, et al. Pluripotency and X chromosome dynamics revealed in pig pre-gastrulating embryos by single cell analysis. *Nat Commun* 2019; 10:623.

16. Yuan X, Chen N, Feng Y, et al. Single-cell multi-omics profiling reveals key regulatory mechanisms that poise germinal vesicle oocytes for maturation in pigs. *Cell Mol Life Sci* 2023; 80:222.

17. Argelaguet R, Arnol D, Bredikhin D, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biol* 2020; 21:111.

18. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis in Python. *Genome Biol* 2018; 19:15.

19. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods* 2019; 16:1289\u20131296.

20. Becht E, McInnes L, Healy J, et al. Dimensionality reduction for visualizing single-cell data using UMAP. *Nat Biotechnol* 2018; 36:38\u201344.

21. Chen P-Y, Chen Z-Y, et al. An integrated single-cell atlas of porcine preimplantation development. *bioRxiv* 2024. Preprint.

22. Xue Z, Huang K, Cai C, et al. Human and mouse preimplantation single-cell RNA-seq. GEO Accession GSE44183. 2013.

23. Aibar S, Gonz\u00e1lez-Blas CB, Moerman T, et al. SCENIC: single-cell regulatory network inference and clustering. *Nat Methods* 2017; 14:1083\u20131086.

24. Van de Sande B, Flerin C, Davie K, et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nat Protoc* 2020; 15:2247\u20132276.

25. Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon: fast and bias-aware quantification of transcript expression. *Nat Methods* 2017; 14:417\u2013419.

26. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 2014; 15:550.

27. Howell CY, Bestor TH, Ding F, et al. Genomic imprinting disrupted by a maternal effect mutation in the Dnmt1 gene. *Cell* 2001; 104:829\u2013838.

28. Herzog VA, Reichholf B, Neumann T, et al. Thiol-linked alkylation of RNA to assess expression dynamics. *Nat Methods* 2017; 14:1198\u20131204.

29. Flamary R, Courty N, Gramfort A, et al. POT: Python Optimal Transport. *J Mach Learn Res* 2021; 22:1\u20138.

30. Schiebinger G, Shu J, Tabaka M, et al. Optimal-transport analysis of single-cell gene expression identifies developmental trajectories in reprogramming. *Cell* 2019; 176:928\u2013943.

"""

ref_start = text.find('## References')
ref_end = text.find('## Table 1')
if ref_start < 0 or ref_end < 0:
    print(f'ERROR: References 定位失败 ref_start={ref_start} ref_end={ref_end}')
else:
    text = text[:ref_start] + new_refs + '\n' + text[ref_end:]

with open(md, 'w', encoding='utf-8') as f:
    f.write(text)

print('Step 4 参考文献重建完成，共 30 条')
