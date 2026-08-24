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

| Script | Purpose |
| --- | --- |
| `cgmap_to_gene_v2.py` | Vectorized CpG→gene mapping using numpy (promoter statistics). |
| `generate_s15_cpg_qc.py` | Generates CpG-quality-control report (Supplementary Table 15). |
| `md_to_docx.py` | Converts markdown tables into the .docx supplementary file. |
| `salmon_batch_wsl2.sh` | Batch salmon quantification of RNA-seq samples. |
| `rerun_failed.sh` | Re-runs the 2 failed PA 2-cell samples. |
| `run_pipeline.sh` | Orchestrates SRA→fastq→salmon re-run for the 2 failed samples. |

> Note: the MOFA+ factor training step produces the tables in `processed/`
> (e.g. `mofa_factors.csv`, `mofa_weights_RNA.csv`). MOFA+ can be run with
> `mofapy2` (Python) or the `MOFA2` Bioconductor R package; the trained outputs
> are provided here so the downstream analyses are fully reproducible.

## License

MIT — see [LICENSE](LICENSE).

## How to cite

Liu S, Liu Z, Pei Y. (2026). Multi-omics profiling identifies candidate
maternal-associated modules in porcine preimplantation embryos.
*BMC Genomics*. [journal DOI pending].

Code archive: [Zenodo DOI — to be added after GitHub release].
