# Sunflower Aquaporin Gene Family: Genome-Wide Characterization and Multi-Stress Transcriptomic Profiling

This repository holds the analysis code, result tables and figures behind:

**Delen, Y., Palali Delen, S., Haliloglu, K., Türkoğlu, A., and Alipour, H.** (2026). Genome-wide characterization of the aquaporin gene family in sunflower (*Helianthus annuus* L.) and transcriptomic profiling under multiple abiotic and biotic stresses.

## Overview

We identified 87 aquaporin genes (99 protein isoforms) in the sunflower reference genome HanXRQr2.0 and profiled their expression across six public RNA-seq datasets (249 libraries), covering 16 stress contrasts and seven stress types (cold, heat, drought, salt, flooding, *Sclerotinia sclerotiorum* and *Orobanche cumana*).

The repository contains the scripts that were run on the cluster, the result tables they produced, the Python scripts that draw the figures and assemble the supplementary workbook, and the figures and tables used in the paper. Raw FASTQ files, BAM files, count matrices and other bulky intermediates are not tracked (see `.gitignore`); they can be regenerated from the SRA run lists in `03_rnaseq/accession_lists/`.

## Repository structure

```
00_setup/                  Environment setup and dependency installation
01_references/             Reference genome download and preparation
02_gene_family/            HMM search, BLAST, domain verification, gene tables, duplication and Ka/Ks
03_rnaseq/                 RNA-seq pipeline (run lists, QC, alignment, counting)
04_expression/             DESeq2, aquaporin expression tables, WGCNA, GO enrichment
05_phylogenetics/          MAFFT alignment, IQ-TREE phylogeny, MEME motifs
06_cis_elements/           Promoter cis-element scanning
07_synteny/                MCScanX synteny analysis
08_figures/                Figure scripts and the supplementary-table builder (Python)
config.sh                  Shared configuration (paths, module names)
results/                   Pipeline outputs read by the figure scripts (same layout as on the cluster)
tables_and_figures/
    figures/               Main figures (Fig 1 to 5, PDF and PNG)
    supplementary_figures/ Supplementary figures (Fig S1 to S7, PDF and PNG)
    tables/                Main tables (Table 1 to 3) and Table S12 as TSV
    supplementary_tables/  Supplementary_Tables.xlsx (Tables S1 to S11 and S13)
```

## Running the pipeline

The shell scripts were written for the SLURM scheduler on the HCC Swan cluster (University of Nebraska-Lincoln). Paths, module versions and sample lists live in `config.sh` (project-wide) and `03_rnaseq/config.sh` (RNA-seq specific); edit those first to match your own system.

Only `00_setup/01_create_dirs.sh` is a plain shell script; everything else is submitted with `sbatch`. The stages are meant to be run in the order below, and each expects the previous stage to have finished. Scripts with a `_v2` suffix supersede the script of the same number: the earlier script was kept because its outputs are inputs of the later one or are referenced in the comparison tables under `results/`.

```bash
# 1. Environment and reference data
bash   00_setup/01_create_dirs.sh
sbatch 00_setup/02_install_r_packages.sh
sbatch 00_setup/03_setup_conda.sh
sbatch 01_references/01_download_genome.sh

# 2. Aquaporin gene family identification
sbatch 02_gene_family/01_hmm_search.sh
sbatch 02_gene_family/02_blast_validation.sh
sbatch 02_gene_family/03_domain_verification_v2.sh   # DeepTMHMM transmembrane helices
sbatch 02_gene_family/04_characterization.sh
sbatch 02_gene_family/05_gene_structure.sh
sbatch 02_gene_family/06_chromosomal_location.sh

# 3. Phylogeny, subfamily assignment, gene tables, duplication, motifs
sbatch 05_phylogenetics/01_msa_and_tree_v2.sh        # reference aquaporins from HMM scans of four proteomes
sbatch 02_gene_family/08_classify_by_tree_v2.sh      # subfamily by tree placement (classify_by_tree.py)
sbatch 02_gene_family/07_build_gene_tables_v2.sh     # gene table, LOC mapping, promoters, Table S1
sbatch 02_gene_family/09_duplication_kaks_v2.sh      # paralog pairs, duplication class, Ka/Ks
sbatch 05_phylogenetics/02_motif_analysis_v2.sh      # MEME protein motifs

# 4. Cis-elements and synteny
sbatch 06_cis_elements/01_cis_element_analysis.sh
sbatch 07_synteny/01_synteny_analysis.sh
sbatch 07_synteny/03_reparse_synteny_v2.sh           # pairs restricted to the aquaporin inventory

# 5. RNA-seq processing
sbatch 03_rnaseq/02_get_srr_ids.sh     # run this before 01_download_sra.sh
sbatch 03_rnaseq/01_download_sra.sh    # array job, one task per BioProject
sbatch 03_rnaseq/03_qc_trim.sh         # array job, one task per run
sbatch 03_rnaseq/04_multiqc.sh
sbatch 03_rnaseq/05_hisat2_align.sh    # array job, one task per run
sbatch 03_rnaseq/06_featurecounts.sh
sbatch 03_rnaseq/07_alignment_stats.sh
# PRJNA492303 (flooding, 96 runs): trimming and alignment of the remaining runs,
# per-BAM counting and a merged count matrix
sbatch 03_rnaseq/03_qc_trim_492303_v2.sh
sbatch 03_rnaseq/05_hisat2_492303_v2.sh
sbatch 03_rnaseq/06_featurecounts_492303_v2.sh
sbatch 03_rnaseq/06b_merge_counts_492303_v2.sh

# 6. Expression analysis
sbatch 04_expression/01_deseq2_analysis.sh           # DESeq2 per BioProject
sbatch 04_expression/06_deseq2_492303_full_v2.sh     # flooding on all 96 runs, ~ genotype + tissue + age + condition
sbatch 04_expression/02_aquaporin_expression_v2.sh   # aquaporin normalized counts and heatmaps (outputs not tracked here)
sbatch 04_expression/03_wgcna_analysis_v2.sh         # WGCNA on PRJNA869183
sbatch 04_expression/04_aquaporin_deg_table_v2.sh    # aquaporin log2FC/padj matrix (Table S8)
sbatch 04_expression/08_go_all_v2.sh                 # GO enrichment, 16 contrasts, one common universe
sbatch 04_expression/09_aquaporin_tpm_v2.sh          # aquaporin TPM in all 249 libraries, control means by tissue

# 7. Figures and supplementary workbook (local machine, Python 3)
cd 08_figures
python fig_tree.py                # Fig 1
python fig_synteny.py             # Fig 2
python fig_expression_set.py      # Fig 3, Fig S3 to S6
python fig_volcano.py             # Fig 4
python cis_vs_expression.py       # Fig 5 and the cis-element tables
python fig_gene_structure.py      # Fig S1
python fig_chromosome_map.py      # Fig S2
python fig_tissue_baseline.py     # Fig S7
python build_supplementary_tables.py   # Supplementary_Tables.xlsx
```

A few things worth knowing before starting:

- `02_get_srr_ids.sh` has to run before `01_download_sra.sh`, because the download job reads the accession lists that the first script fetches from NCBI. The lists used for the paper are tracked in `03_rnaseq/accession_lists/`.
- The array sizes in `03_qc_trim.sh` and `05_hisat2_align.sh` are upper bounds; tasks past the end of `all_samples.txt` simply exit. The download, trimming and alignment jobs skip runs whose output already exists, so re-submitting them is safe.
- `03_domain_verification_v2.sh` predicts transmembrane helices with DeepTMHMM through the BioLib client; the prediction runs on the BioLib servers, so the job needs network access (or run the `biolib` command given in the script header once and resubmit).
- `04_go_kegg_enrichment.sh` is the first GO run. The GO table in the paper (Table S11) comes from `08_go_all_v2.sh`, which tests all 16 contrasts against the same background (all genes with a DESeq2 result in any contrast, 49,914 genes).
- The figure scripts read `results/` by default; set `AQP_RESULTS` to point them at another copy of the pipeline outputs. They write into `tables_and_figures/`.

## Figures

### Main figures

| Figure | Description |
|---|---|
| Fig 1 | Maximum-likelihood phylogeny of 305 aquaporin proteins from sunflower (99), Arabidopsis (59), rice (36), tomato (48) and lettuce (63) |
| Fig 2 | Collinear aquaporin gene pairs from MCScanX: intra-genomic pairs and pairs with lettuce and Arabidopsis |
| Fig 3 | Up- and downregulated aquaporin DEGs in each of the 16 stress contrasts |
| Fig 4 | Volcano plots of the cold, drought 14 d, salt and *Sclerotinia* contrasts with aquaporin genes highlighted |
| Fig 5 | Promoter cis-element content versus stress responsiveness (Fisher tests and Spearman correlations) |

### Supplementary figures

| Figure | Description |
|---|---|
| Fig S1 | Exon-intron structure and MEME motifs of the 87 aquaporin genes |
| Fig S2 | Chromosomal distribution of the 87 aquaporin genes on the 17 chromosomes |
| Fig S3 | log2 fold-change heatmap of the 87 genes across the 16 contrasts |
| Fig S4 | Overlap of the aquaporin DEG sets of the seven stress types |
| Fig S5 | The 30 aquaporin genes significant in five or more contrasts |
| Fig S6 | Module-trait correlations of the WGCNA modules of PRJNA869183 that contain aquaporin genes |
| Fig S7 | Baseline expression of the 87 aquaporin genes in the control libraries, by BioProject and tissue |

## Tables

Main tables (`tables_and_figures/tables/`): Table 1 (RNA-seq datasets and DESeq2 designs), Table 2 (stress-responsive cis-elements in the aquaporin promoters), Table 3 (aquaporin genes responding in the largest number of contrasts) and Table S12 (genome-wide DEG counts per contrast). The `cis_vs_expression_*.tsv` files are the outputs of `08_figures/cis_vs_expression.py`.

`Supplementary_Tables.xlsx` (Tables S1 to S11): S1 gene inventory, S2 excluded candidates, S3 DESeq2 design rationale, S4 LOC-to-TAIR mapping, S5 physicochemical properties, S6 MEME motifs, S7a duplication pairs, S7b synteny pairs, S7c Ka/Ks, S8 aquaporin DEG matrix, S9 WGCNA module membership, S10 module-trait correlations, S11 GO enrichment, S13 baseline TPM of the control libraries by BioProject and tissue (Table S12 is in the Additional file 2 of the paper and is provided here as `TableS12_genome_wide_degs.tsv`).

## Data availability

All RNA-seq data are public in NCBI SRA:

| BioProject | Stress type | Runs used |
|---|---|---|
| PRJNA869183 | Cold, heat, drought, salt, rehydration (seedling leaves) | 57 |
| PRJNA797473 | Drought, 0/7/14/21 d (leaves) | 24 |
| PRJNA1041959 | Drought, PEG 72 h (leaves and roots) | 24 |
| PRJNA492303 | Flooding (leaves and roots, two genotypes, three ages) | 96 |
| PRJNA908908 | *Sclerotinia sclerotiorum* (leaves) | 12 |
| PRJNA850121 | *Orobanche cumana*, stages A to E (roots) | 36 |

Reference genome: HanXRQr2.0 (GCF_002127325.2) from NCBI RefSeq.

## Software

| Software | Version | Purpose |
|---|---|---|
| HMMER | 3.4 | PF00230 (MIP) domain search |
| BLAST+ | 2.17 | Candidate validation, orthologs, self-BLAST for duplications |
| DeepTMHMM | BioLib | Transmembrane helix prediction |
| MAFFT | 7.526 | Multiple sequence alignment |
| trimAl | conda | Alignment trimming |
| IQ-TREE | 2.2.2.7 | Phylogenetic inference (ModelFinder, UFBoot, SH-aLRT) |
| MEME Suite | 5.5 | Protein motif discovery |
| MCScanX | bioconda | Collinearity detection |
| fastp | 0.23 | Read trimming |
| FastQC | 0.12 | Read quality control |
| MultiQC | py37/1.8 module | QC aggregation |
| HISAT2 | 2.2 | Read alignment |
| SAMtools | 1.23 | BAM processing |
| Subread/featureCounts | 2.1 | Read counting |
| R | 4.1.3 | Statistical computing (conda environment aqp_env) |
| DESeq2 | 1.34.0 | Differential expression |
| clusterProfiler | 4.2.0 | GO enrichment (org.At.tair.db) |
| WGCNA | 1.71 | Co-expression networks |
| ComplexHeatmap | 2.10.0 | Expression heatmaps |
| Python | 3.10 | Pipeline helpers; figures need matplotlib, numpy, scipy, biopython and openpyxl |

Numbered versions are the cluster modules that were loaded (see `config.sh` and the `module load` lines of each script) or, for R and its packages, the versions installed in the conda environment `conda_envs/aqp_env` that the scripts in `04_expression/` activate (R 4.1.3 with Bioconductor 3.14 packages). The full package list of that environment is in `00_setup/aqp_env.yml`; it can be recreated with `conda env create -p conda_envs/aqp_env -f 00_setup/aqp_env.yml`. The annotation package org.At.tair.db is not part of the environment and was installed into the project R library (`R_libs`) with BiocManager (`00_setup/02_install_r_packages.sh`). The remaining command-line tools were installed at their then-current release through conda in `00_setup/03_setup_conda.sh`; the environment created there (`aquaporin_env`) is the one used by the synteny scripts.

## Revision history

The first version of this repository described 61 aquaporin genes. The transmembrane-helix screen of that version (a hydropathy sliding window) removed annotated aquaporins with a complete MIP domain; replacing it with DeepTMHMM (`03_domain_verification_v2.sh`) gives the 87-gene inventory used in the paper. The reference aquaporins for the phylogeny are now identified by HMM scans of the four reference proteomes rather than from a fixed list of UniProt accessions, subfamilies are assigned from tree placement, and the duplication, Ka/Ks and synteny tables were rebuilt for the new inventory. The flooding dataset PRJNA492303 was re-analysed with all 96 runs and the design `~ genotype + tissue + age + condition` (the first version used 46 runs), GO enrichment was re-run for all 16 contrasts against one common gene universe, and the figures are now produced by the scripts in `08_figures/`. In `02_gene_family/02_blast_validation.sh` the UniProt accession listed for AtTIP3;2 was Q9ZV07, which is AtPIP2;6 (already in the list); it was corrected to O22588, and the script now stops when any query sequence fails to download instead of only warning below ten sequences. The query set only adds candidates to the HMM hits, so this correction does not change the gene inventory, which is fixed by the PF00230 domain and the DeepTMHMM criteria.

## Citation

If you use this code or the derived tables, please cite the paper above. The individual tools should be cited separately; the versions used are listed in the software table.

## License

This project is licensed under the MIT License.

## Contact

Yavuz Delen - ydelen2@unl.edu
