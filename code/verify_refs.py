#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""用 NCBI E-utilities 逐条核验 32 条参考文献的标题/期刊/卷期页真实性。"""
import json, time, urllib.parse, urllib.request

# 32 条文献：编号 -> (搜索标题, 稿件写的期刊年份卷期页)
refs = {
 1: ("The maternal-to-zygotic transition a play in two acts", "Development 2009;136:3033-3042"),
 2: ("Zygotic genome activation during the maternal-to-zygotic transition", "Annu Rev Cell Dev Biol 2014;30:581-613"),
 3: ("Dynamics of the epigenetic landscape during the maternal-to-zygotic transition", "Nat Rev Mol Cell Biol 2018;19:436-450"),
 4: ("Maternal control of early embryogenesis in mammals", "Reprod Fertil Dev 2015;27:880-896"),
 5: ("Development of reconstituted mouse eggs suggests imprinting of the genome during gametogenesis", "Nature 1984;308:548-550"),
 6: ("Completion of mouse embryogenesis requires both the maternal and paternal genomes", "Cell 1984;37:179-183"),
 7: ("Birth of parthenogenetic mice that can develop to adulthood", "Nature 2004;428:860-864"),
 8: ("Zygotic genome activation in vertebrates", "Dev Cell 2017;42:316-332"),
 9: ("Genetic programs in human and mouse early embryos revealed by single-cell RNA sequencing", "Nature 2013;500:593-597"),
 10: ("Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells", "Nat Struct Mol Biol 2013;20:1131-1139"),
 11: ("From zygote to blastocyst molecular aspects of porcine early embryonic development", "Cells 2026;15:15"),
 12: ("Single cell RNA-seq reveals genes vital to in vitro fertilized embryos and parthenotes in pigs", "Sci Rep 2021;11:14393"),
 13: ("Lineage specification and pluripotency revealed by transcriptome analysis from oocyte to blastocyst in pig", "FASEB J 2020;34:691-705"),
 14: ("Generation and characterization of stable pig pregastrulation epiblast stem cell lines", "Cell Res 2022;32:383-400"),
 15: ("Pluripotency and X chromosome dynamics revealed in pig pre-gastrulating embryos by single cell analysis", "Nat Commun 2019;10:500"),
 16: ("Single-cell multi-omics profiling reveals key regulatory mechanisms that poise germinal vesicle oocytes for maturation in pigs", "Cell Mol Life Sci 2023;80:222"),
 17: ("Deciphering differences in DNA methylation and transcriptome profiles of oocytes from pigs with high and low developmental competence", "Environ Epigenet 2025;11:dvaf018"),
 18: ("MOFA+ a statistical framework for comprehensive integration of multi-modal single-cell data", "Genome Biol 2020;21:111"),
 19: ("Single-cell RNA-seq reveals mRNAs and lncRNAs important for oocytes in vitro matured in pigs", "Reprod Domest Anim 2021;56:642-657"),
 20: ("SCANPY large-scale single-cell gene expression data analysis", "Genome Biol 2018;19:15"),
 21: ("Fast sensitive and accurate integration of single-cell data with Harmony", "Nat Methods 2019;16:1289-1296"),
 22: ("Dimensionality reduction for visualizing single-cell data using UMAP", "Nat Biotechnol 2018;36:38-44"),
 24: ("SCENIC single-cell regulatory network inference and clustering", "Nat Methods 2017;14:1083-1086"),
 25: ("A scalable SCENIC workflow for single-cell gene regulatory network analysis", "Nat Protoc 2020;15:2247-2276"),
 26: ("Salmon fast and bias-aware quantification of transcript expression", "Nat Methods 2017;14:417-419"),
 27: ("Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2", "Genome Biol 2014;15:550"),
 28: ("Genomic imprinting disrupted by a maternal effect mutation in the Dnmt1 gene", "Cell 2001;104:829-838"),
 29: ("Single-Cell RNA Sequencing Reveals Differences in Chromatin Remodeling and Energy Metabolism among In Vivo-Developed In Vitro-Fertilized and Parthenogenetically Activated Embryos", "Animals 2024;14:465"),
 30: ("Thiol-linked alkylation of RNA to assess expression dynamics", "Nat Methods 2017;14:1198-1204"),
 31: ("POT Python Optimal Transport", "J Mach Learn Res 2021;22:1-8"),
 32: ("Optimal-transport analysis of single-cell gene expression identifies developmental trajectories in reprogramming", "Cell 2019;176:928-943"),
}

def esearch(title):
    q = urllib.parse.quote(title)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={q}[Title]&retmode=json&retmax=3"
    with urllib.request.urlopen(url, timeout=20) as r:
        d = json.load(r)
    return d['esearchresult'].get('idlist', [])

def esummary(pmid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    with urllib.request.urlopen(url, timeout=20) as r:
        d = json.load(r)
    return d['result'][pmid]

results = []
for num in sorted(refs.keys()):
    title, manuscript = refs[num]
    try:
        pmids = esearch(title)
        time.sleep(0.4)
        if not pmids:
            results.append((num, title, '✗ 未找到 PMID', manuscript, ''))
            continue
        pmid = pmids[0]
        s = esummary(pmid)
        time.sleep(0.4)
        journal = s.get('fulljournalname') or s.get('source', '')
        year = s.get('pubdate', '')[:4]
        vol = s.get('volume', '')
        pages = s.get('pages', '')
        actual = f"{journal} {year};{vol}:{pages}"
        results.append((num, title, f'✓ PMID {pmid}', manuscript, actual))
    except Exception as e:
        results.append((num, title, f'⚠ 错误 {e}', manuscript, ''))

print('=== 核验结果 ===')
for num, title, status, manuscript, actual in results:
    print(f'[{num}] {status}')
    print(f'    稿件: {manuscript}')
    if actual:
        print(f'    PubMed: {actual}')
    print()
