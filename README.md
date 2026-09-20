# Porcine-MZT-multiomics

Analysis code and processed data accompanying the manuscript:

> Liu S, Liu Z, Pei Y. *Multi-omics profiling identifies candidate maternal-associated modules in porcine preimplantation embryos.* BMC Genomics (submitted).

## Overview

This repository provides the reproducible analysis code and intermediate/processed
result tables for a multi-omics reanalysis of porcine preimplantation development.
We integrated matched single-oocyte transcriptomic and DNA-methylation profiles
with multi-omics factor analysis (MOFA+), projected the inferred modules onto a
single-cell developmental atlas, and cross-compared them with an independent
parthenogenetic (PA) / *in vitro* fertilization (IVF) bulk RNA-seq cohort.

## Data availability

All **raw** data are publicly available and were downloaded from the following
repositories (no new sequencing was generated in this study):

| Dataset | Repository | Accession |
| --- | --- | --- |
| Matched scRNA + methylation, 32 porcine GV oocytes (5 donors) | GEO | **GSE235729** |
| Independent PA/IVF bulk RNA-seq cohort (42 samples) | SRA | **SRP301735** |
| 1,955-cell porcine preimplantation developmental atlas | GEO | *see manuscript Methods* |

Download the raw files and place them under `data/raw/` before running the
pipeline. Sequence read alignment/quantification was performed with **salmon**
and SRA conversion with **sratoolkit** (`fasterq-dump`).

## Repository structure

```
Porcine-MZT-multiomics/
├── code/            # analysis scripts (Python + bash)
├── processed/       # intermediate / processed result tables (CSV, TXT)
├── requirements.txt # Python dependencies
├── LICENSE          # MIT
├── CITATION.cff     # citation metadata (Zenodo)
└── README.md        # this file
```

## Requirements

- Python 3.8+ (packages in `requirements.txt`)
- External command-line tools: `salmon`, `sratoolkit` (`fasterq-dump`)
- Reference genome: *Sus scrofa* Sscrofa11.1 (Ensembl)

Install Python dependencies with:

```bash
pip install -r requirements.txt
```

## Usage

1. Download raw data from GEO/SRA into `data/raw/`.
2. **Edit the hard-coded local paths** in the scripts to match your environment:
   - `code/cgmap_to_gene_v2.py` — variables `CG_DIR`, `GTF`, `OUT` (top of file)
   - `code/run_pipeline.sh` — variables `BIN`, `SRADIR`, `TMP`, `LOGDIR`, `OUTDIR`
3. Run the individual scripts in `code/` as needed (see descriptions below).

## Scripts (`code/`)

The `code/` directory contains the full analysis pipeline (145 scripts),
organized by analysis module. The main modules are:

| Module | Purpose |
| --- | --- |
| `m2_*.py` | Data integration — Harmony batch correction and atlas projection |
| `m3_*.py` | PA/IVF bulk RNA-seq differential expression and validation |
| `m4_*.py` | MOFA+ multi-omics integration and promoter CpG analysis |
| `m5_*.py` | Developmental trajectory and Waddington optimal transport |
| `m6_*.py` | Transcriptome-based metabolic scoring |
| `m7_*.py` / `m9_*.py` | Gene-regulatory network (SCENIC) and TF/lineage analysis |
| `m8_*.py` | RNA-only MOFA+ robustness analysis |
| `m10_*.py` / `m11_*.py` | CGmap processing and CGmap-to-MOFA integration |
| `m12–m16` | Methylation matrix, R², MOFA+ comparison, DESeq2, GO enrichment |
| `m19_*.py` | In-silico regulon perturbation |
| `m20_*.py` | SHAP donor classifier |
| `m21–m24` | Cross-species, trajectory, PA validation, F1 stability |
| `regenerate_*.py` | Final publication-figure regeneration |
| `*.sh` | SRA→fastq→salmon batch quantification (WSL2) |

> Note: the MOFA+ factor training step produces the tables in `processed/`
> (e.g. `m4_mofa/mofa_multiomics_factors_v2.csv`,
> `m4_mofa/mofa_multiomics_weights_RNA_v2.csv`). MOFA+ can be run with
> `mofapy2` (Python) or the `MOFA2` Bioconductor R package; the trained outputs
> are provided here so the downstream analyses are fully reproducible.

## License

MIT — see [LICENSE](LICENSE).

## How to cite

Liu S, Liu Z, Pei Y. (2026). Multi-omics profiling identifies candidate
maternal-associated modules in porcine preimplantation embryos.
*BMC Genomics*. [journal DOI pending].

Code archive: [Zenodo DOI — to be added after GitHub release].
