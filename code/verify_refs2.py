#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""精确核验存疑条目（第一作者[Author] AND 标题[Title]）。"""
import json, time, urllib.parse, urllib.request

# 编号 -> (搜索 query, 稿件卷期页)
refs = {
 2: ('Lee MT[Author] AND "Zygotic genome activation during the maternal-to-zygotic transition"[Title]', "Annu Rev Cell Dev Biol 2014;30:581-613"),
 8: ('Jukam D[Author] AND "Zygotic genome activation in vertebrates"[Title]', "Dev Cell 2017;42:316-332"),
 9: ('Xue Z[Author] AND "Genetic programs in human and mouse early embryos"[Title]', "Nature 2013;500:593-597"),
 10: ('Yan L[Author] AND "Single-cell RNA-Seq profiling of human preimplantation embryos"[Title]', "Nat Struct Mol Biol 2013;20:1131-1139"),
 15: ('Ramos-Ibeas[Author] AND "Pluripotency and X chromosome dynamics"[Title]', "Nat Commun 2019;10:500"),
 22: ('Becht E[Author] AND "Dimensionality reduction for visualizing single-cell data using UMAP"[Title]', "Nat Biotechnol 2018;36:38-44"),
 24: ('Aibar S[Author] AND "SCENIC"[Title]', "Nat Methods 2017;14:1083-1086"),
 32: ('Schiebinger G[Author] AND "Optimal-transport analysis of single-cell gene expression"[Title]', "Cell 2019;176:928-943"),
 3: ('Eckersley-Maslin[Author] AND "Dynamics of the epigenetic landscape"[Title]', "Nat Rev Mol Cell Biol 2018;19:436-450"),
 4: ('Zhang K[Author] AND "Maternal control of early embryogenesis"[Title]', "Reprod Fertil Dev 2015;27:880-896"),
}

def esearch(q):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={urllib.parse.quote(q)}&retmode=json&retmax=3"
    with urllib.request.urlopen(url, timeout=20) as r:
        return json.load(r)['esearchresult'].get('idlist', [])

def esummary(pmid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    with urllib.request.urlopen(url, timeout=20) as r:
        return json.load(r)['result'][pmid]

for num in sorted(refs.keys()):
    q, manuscript = refs[num]
    try:
        pmids = esearch(q)
        time.sleep(0.4)
        if not pmids:
            print(f'[{num}] ✗ 未找到 PMID | 稿件: {manuscript}')
            continue
        pmid = pmids[0]
        s = esummary(pmid)
        time.sleep(0.4)
        journal = s.get('source', '')
        year = s.get('pubdate', '')[:4]
        vol = s.get('volume', '')
        pages = s.get('pages', '')
        title = s.get('title', '')[:70]
        print(f'[{num}] PMID {pmid} | {title}')
        print(f'    稿件: {manuscript}')
        print(f'    PubMed: {journal} {year};{vol}:{pages}')
        print()
    except Exception as e:
        print(f'[{num}] ⚠ {e} | 稿件: {manuscript}')
        print()
