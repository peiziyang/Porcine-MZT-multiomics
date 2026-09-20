#!/usr/bin/env python
"""Package all manuscript deliverables into a single upload-ready zip."""
import zipfile, os, shutil

ROOT = r'E:/Workbuddy/2026-07-27-11-58-27'
ZIP = os.path.join(ROOT, 'manuscript', 'pig_mzt_deposit.zip')
os.makedirs(os.path.dirname(ZIP), exist_ok=True)

FILES = {
    # Core processed data
    'data/processed/m4_mofa/mofa_multiomics_weights_RNA_v2.csv': 'processed/mofa_weights_RNA.csv',
    'data/processed/m4_mofa/mofa_multiomics_weights_METH_v2.csv': 'processed/mofa_weights_METH.csv',
    'data/processed/m4_mofa/mofa_multiomics_factors_v2.csv': 'processed/mofa_factors.csv',
    'data/processed/m4_mofa/factor_top_genes_v2.csv': 'processed/factor_top_genes.csv',
    'data/processed/m4_mofa/factor_cross_view_v2.csv': 'processed/factor_cross_view.csv',
    'data/processed/m4_mofa/cross_species_permutation_results.csv': 'processed/cross_species_permutation.csv',
    'data/processed/m4_mofa/perturbation_perm/perturbation_permutation_results.csv': 'processed/perturbation_permutation.csv',
    'data/processed/m4_mofa/per_gene_cpg_stats.csv': 'processed/per_gene_cpg_stats.csv',
    'data/processed/m4_mofa/cpg_qc_report.csv': 'processed/cpg_qc_report.csv',
    'data/processed/deseq2/m15_deseq2_overall_IVFvsPA_stageAdj.csv': 'processed/deseq2_stage_adjusted.csv',
    # Key scripts
    'scripts/cgmap_to_gene_v2.py': 'code/cgmap_to_gene_v2.py',
    'scripts/generate_s15_cpg_qc.py': 'code/generate_s15_cpg_qc.py',
    'scripts/md_to_docx.py': 'code/md_to_docx.py',
    'scripts/salmon_batch_wsl2.sh': 'code/salmon_batch_wsl2.sh',
    'scripts/run_pipeline.sh': 'code/run_pipeline.sh',
    'scripts/rerun_failed.sh': 'code/rerun_failed.sh',
}

missing = []
count = 0
with zipfile.ZipFile(ZIP, 'w', zipfile.ZIP_DEFLATED) as zf:
    # README
    zf.writestr('README.txt', """Processed Data and Code for:
"Integrated Transcriptomic and Methylation Analysis Identifies Candidate Maternal-Associated Programs during Porcine Preimplantation Development"

Contents:
  processed/ — MOFA+ results, DESeq2 output, CpG statistics, cross-species and perturbation permutations
  code/      — Python scripts for CpG reanalysis, figure generation, and salmon batch quantification

Primary raw data accession numbers: GSE235729, GSE44183, SRP301735, GSE168106

Contact: peiya@??? (or as appropriate)
""")
    
    for src, dst in FILES.items():
        sp = os.path.join(ROOT, src)
        if os.path.exists(sp):
            zf.write(sp, dst)
            count += 1
            print(f'  + {dst}')
        else:
            missing.append(src)

    # Also write a gene_set_list.txt
    zf.writestr('processed/gene_set_list.txt', 
"""F1/F3/F4/F6 gene sets: top 44 genes by absolute RNA weight for each MOFA+ factor, 
as defined in factor_top_genes.csv (column: factor, rows sorted by absolute weight descending).

F1 core five genes: DNMT1, ZP3, ZP4, GDF9, RARRES1
""")

size_mb = os.path.getsize(ZIP) / 1024 / 1024
print(f'\nZip: {ZIP} ({size_mb:.1f} MB)')
print(f'Files: {count} added, {len(missing)} missing')
if missing:
    print(f'Missing: {missing}')
