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
02_gene_family/            HMM search, BLAST, domain verification, protein properties, NPA motifs and ar/R filter,
                           gene tables, duplication and Ka/Ks
03_rnaseq/                 RNA-seq pipeline (run lists, QC, alignment, counting)
04_expression/             DESeq2, aquaporin expression tables, WGCNA, GO enrichment
05_phylogenetics/          MAFFT alignment, IQ-TREE phylogeny, MEME motifs
06_cis_elements/           Promoter cis-element scanning
07_synteny/                MCScanX synteny analysis
08_figures/                Figure scripts and the supplementary-table builder (Python)
config.sh                  Shared configuration (paths, module names)
results/                   Pipeline outputs read by the figure scripts (same layout as on the cluster)
tables_and_figures/
    figures/               Main figures (Fig 1 to 4, PDF and PNG)
    supplementary_figures/ Supplementary figures (Fig S1 to S9, PDF and PNG)
    tables/                Main tables (Table 1 to 3) and Table S12 as TSV
    supplementary_tables/  Supplementary_Tables.xlsx (Tables S1 to S11 and S13)
```

## Running the pipeline

The shell scripts were written for the SLURM scheduler on the HCC Swan cluster (University of Nebraska-Lincoln). Paths, module versions and sample lists live in `config.sh` (project-wide) and `03_rnaseq/config.sh` (RNA-seq specific); edit those first to match your own system.

Only `00_setup/01_create_dirs.sh` is a plain shell script; everything else is submitted with `sbatch`. The stages are meant to be run in the order below, and each expects the previous stage to have finished. Scripts with a `_v2` or `_v3` suffix replace the script of the same number. A script without suffix is kept where a later step reads its outputs; the list below the commands says which of its outputs the paper uses.

```bash
# 1. Environment and reference data
bash   00_setup/01_create_dirs.sh
sbatch 00_setup/02_install_r_packages.sh
sbatch 00_setup/03_setup_conda.sh
sbatch 01_references/01_download_genome.sh

# 2. Aquaporin gene family identification
sbatch 02_gene_family/01_hmm_search.sh
sbatch 02_gene_family/02_blast_validation.sh         # writes the candidate FASTA of the HMM hits
sbatch 02_gene_family/03_domain_verification_v2.sh   # DeepTMHMM transmembrane helices
sbatch 02_gene_family/05_gene_structure.sh
sbatch 02_gene_family/06_chromosomal_location.sh
sbatch 02_gene_family/11_tblastn_unannotated_v3.sh   # aquaporin-like loci of the genome without a gene model

# 3. Phylogeny, subfamily assignment, gene tables, duplication, motifs
sbatch 05_phylogenetics/01_msa_and_tree_v2.sh        # reference aquaporins from HMM scans of four proteomes
sbatch 02_gene_family/02_blast_validation_v3.sh      # BLASTP of the Arabidopsis and rice reference aquaporins against the proteome
sbatch 02_gene_family/08_classify_by_tree_v2.sh      # subfamily by tree placement (classify_by_tree.py)
sbatch 02_gene_family/07_build_gene_tables_v2.sh     # gene table, LOC mapping, promoters, Table S1
sbatch 02_gene_family/09_duplication_kaks_v3.sh      # paralog pairs (all-against-all BLASTP without a subject limit), Ka/Ks
sbatch 05_phylogenetics/02_motif_analysis_v2.sh      # MEME protein motifs

# 4. Cis-elements and synteny
sbatch 06_cis_elements/01_cis_element_analysis.sh
sbatch 07_synteny/01_synteny_analysis.sh                   # gene position files for MCScanX
sbatch 07_synteny/05_synteny_one_protein_per_gene_v3.sh    # BLASTP with one protein per gene, MCScanX (three comparisons)

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
sbatch 04_expression/01_deseq2_analysis.sh           # DESeq2 per BioProject, variance-stabilized counts
sbatch 04_expression/06_deseq2_492303_full_v2.sh     # flooding on all 96 runs, ~ genotype + tissue + age + condition
sbatch 04_expression/10_deseq2_designs_v3.sh         # the 16 contrasts with their designs, interaction tests, contrasts within genotype/tissue, joint FDR
sbatch 04_expression/02_aquaporin_expression_v2.sh   # aquaporin normalized counts and heatmaps (outputs not tracked here)
sbatch 04_expression/03_wgcna_analysis_v3.sh         # WGCNA on all 57 libraries of PRJNA869183
sbatch 04_expression/04_aquaporin_deg_table_v2.sh    # aquaporin log2FC/padj matrix (Table S8)
sbatch 02_gene_family/13_blast_arabidopsis_best_hit.sh   # best Arabidopsis hit of every sunflower protein (source of Table S4)
sbatch 04_expression/08_go_all_v2.sh                 # GO enrichment, 16 contrasts, one common universe
sbatch 04_expression/09_aquaporin_tpm_v2.sh          # aquaporin TPM in all 249 libraries, control means by tissue

# 7. Tables from the cluster outputs, figures and supplementary workbook (local machine, Python 3)
python 02_gene_family/04_protein_properties_v3.py           # protein properties (Table S5)
python 02_gene_family/12_npa_arR_from_alignment_v3.py       # NPA motifs and ar/R filter read from the alignment (Table S5)
python 07_synteny/06_reparse_synteny_v3.py                  # aquaporin pairs in the collinear blocks (Table S7b)
python 02_gene_family/10_duplication_synteny_class_v3.py    # duplication classes with synteny support (Tables S7a, S7c)
cd 08_figures
python fig_tree.py                # Fig 1 (fan layout)
python fig_tree.py "" "" radial   # Fig S1 (radial layout)
python fig_synteny.py             # Fig 2
python fig_expression_set.py      # Fig 3, Fig S4, S5, S6 and S8
python fig_volcano.py             # Fig 4
python fig_gene_structure.py      # Fig S2
python fig_chromosome_map.py      # Fig S3
python cis_vs_expression.py       # Fig S7 and the cis-element tables (Table 2)
python fig_tissue_baseline.py     # Fig S9
python build_supplementary_tables.py   # Supplementary_Tables.xlsx
```

Which outputs the paper uses, and what the steps need:

- `02_get_srr_ids.sh` has to run before `01_download_sra.sh`, because the download job reads the accession lists that the first script fetches from NCBI. The lists used for the paper are tracked in `03_rnaseq/accession_lists/`.
- The array sizes in `03_qc_trim.sh` and `05_hisat2_align.sh` are upper bounds; tasks past the end of `all_samples.txt` simply exit. The download, trimming and alignment jobs skip runs whose output already exists, so re-submitting them is safe.
- `03_domain_verification_v2.sh` predicts transmembrane helices with DeepTMHMM through the BioLib client; the prediction runs on the BioLib servers, so the job needs network access (or run the `biolib` command given in the script header once and resubmit).
- `02_blast_validation.sh` writes `candidate_aquaporins.faa`, the sequences of the HMM hits that the later steps read; its BLASTP search is not used. The BLASTP search of the paper is `02_blast_validation_v3.sh`, with the Arabidopsis and rice reference aquaporins of the phylogenetic analysis as queries. The gene inventory is defined by the PF00230 domain and the DeepTMHMM criterion.
- `06_chromosomal_location.sh` provides the gene coordinates; its BLASTP search is not used. The paralog search of the paper is the one inside `09_duplication_kaks_v3.sh`.
- `07_synteny/01_synteny_analysis.sh` downloads the annotations and writes the gene position files (`07_synteny/data/`); its BLASTP searches and collinearity files are not used. The collinearity files of the paper come from `05_synteny_one_protein_per_gene_v3.sh` and are in `results/07_synteny/mcscanx/`.
- Table 2 of the paper and the promoter tests come from `08_figures/cis_vs_expression.py`, which scans the promoters with the consensus sequences given in the Methods and counts every site once. `06_cis_elements/01_cis_element_analysis.sh` uses other consensus definitions; its outputs are not used in the paper.
- The LOC-to-TAIR table used for the GO enrichment (`results/02_gene_family/results/arabidopsis_orthologs_v2.tsv`, Table S4) comes from the output of `13_blast_arabidopsis_best_hit.sh`: a gene is kept when one of its proteins has a hit with identity above 30% and query coverage above 50%, and the hit with the highest bit score gives the TAIR locus; hits to organelle-encoded proteins and to proteins without a nuclear TAIR locus identifier are not mapped. This step is not part of the script.
- The DESeq2 designs and contrasts of the paper are written out in `04_expression/10_deseq2_designs_v3.sh` (two-group fits for PRJNA869183, additive models for the other datasets). It refits every contrast and compares the result with the stored tables under `results/04_expression/results/deseq2/` (`deseq2_designs_v3/reproduction_report.tsv`: same DEG numbers for all 16 contrasts). `01_deseq2_analysis.sh` writes the variance-stabilized matrix that the WGCNA step reads; its designs and DEG tables are not used.
- Two stored intermediate tables should be read with care. In `results/02_gene_family/03_domain_verification_v2/failed_verification.tsv` only the rows with `has_mip = yes` are excluded candidates (the six of Table S2a); the rows whose identifier has the form `ref|XP_...|` come from the candidate list of `02_blast_validation.sh` and were not tested. In `results/02_gene_family/05_gene_structure/aqp_gene_structure_summary.tsv` there is one row per mRNA, and every row of a gene carries the protein of its last-listed isoform; the gene tables take exons and coding sequences from `aqp_genes.gff3`.
- `11_tblastn_unannotated_v3.sh` merges hits within 10 kb on one strand into a locus. One partial aquaporin copy (NC_035447.2:62,074,918-62,075,636, minus strand) lies 3.7 kb from the gene LOC110911759 and is therefore part of the locus of that gene in `tblastn_loci.tsv`.
- `04_go_kegg_enrichment.sh` is not used for the paper. The GO table of the paper (Table S11) comes from `08_go_all_v2.sh`, which tests all 16 contrasts against the same background (the 49,914 genes with an adjusted P value in at least one contrast).
- The figure scripts read `results/` by default; set `AQP_RESULTS` to point them at another copy of the pipeline outputs. They write into `tables_and_figures/`.

## Figures

### Main figures

| Figure | Description |
|---|---|
| Fig 1 | Maximum-likelihood phylogeny of 305 aquaporin proteins from sunflower (99), Arabidopsis (59), rice (36), tomato (48) and lettuce (63), fan layout |
| Fig 2 | Collinear aquaporin gene pairs from MCScanX: intra-genomic pairs and pairs with lettuce and Arabidopsis |
| Fig 3 | Up- and downregulated aquaporin DEGs in each of the 16 stress contrasts |
| Fig 4 | Volcano plots of the cold, drought 14 d, salt and *Sclerotinia* contrasts with aquaporin genes highlighted |

### Supplementary figures

| Figure | Description |
|---|---|
| Fig S1 | The tree of Fig 1 without rooting, radial layout |
| Fig S2 | Exon-intron structure and MEME motifs of the 87 aquaporin genes |
| Fig S3 | Chromosomal distribution of the 87 aquaporin genes on the 17 chromosomes |
| Fig S4 | log2 fold-change heatmap of the 87 genes across the 16 contrasts |
| Fig S5 | Overlap of the aquaporin DEG sets of the seven stress types |
| Fig S6 | The 30 aquaporin genes significant in five or more contrasts |
| Fig S7 | Promoter cis-element content versus stress responsiveness (Fisher tests and Spearman correlations) |
| Fig S8 | Module-trait correlations of the WGCNA modules of PRJNA869183 that contain aquaporin genes |
| Fig S9 | Baseline expression of the 87 aquaporin genes in the control libraries, by BioProject and tissue |

## Tables

Main tables (`tables_and_figures/tables/`): Table 1 (RNA-seq datasets and DESeq2 designs), Table 2 (stress-responsive cis-elements in the aquaporin promoters; `cis_vs_expression_table2.tsv`), Table 3 (aquaporin genes responding in the largest number of contrasts) and Table S12 (genome-wide DEG counts per contrast). The `cis_vs_expression_*.tsv` files are the outputs of `08_figures/cis_vs_expression.py`.

`Supplementary_Tables.xlsx` (Tables S1 to S11 and S13): S1 gene inventory, S2a excluded candidates, S2b loci found by tblastn in the genome sequence, S3a DESeq2 design rationale, S3b contrasts and replicates, S3c interaction tests, S3d contrasts within genotype, tissue or age, S3e corrections for multiple testing, S4 LOC-to-TAIR mapping, S5 physicochemical properties, NPA motifs and ar/R filter, S6 MEME motifs, S7a duplication pairs, S7b synteny pairs, S7c Ka/Ks, S8 aquaporin DEG matrix, S9 WGCNA module membership, S10 module-trait correlations, S11 GO enrichment, S13 baseline TPM of the control libraries by BioProject and tissue (Table S12 is in the Additional file 2 of the paper and is provided here as `TableS12_genome_wide_degs.tsv`).

## Data availability

All RNA-seq data are public in NCBI SRA:

| BioProject | Stress type | Runs used |
|---|---|---|
| PRJNA869183 | Cold, heat, drought, salt, rehydration (seedling leaves) | 57 |
| PRJNA797473 | Drought, 0/7/14/21 d (leaves) | 24 |
| PRJNA1041959 | Drought, irrigation withheld, 72 h (leaves and roots) | 24 |
| PRJNA492303 | Flooding, days 1 and 8, and day 6 of recovery (leaves and roots, two genotypes) | 96 |
| PRJNA908908 | *Sclerotinia sclerotiorum* (leaves) | 12 |
| PRJNA850121 | *Orobanche cumana*, 4, 12, 24, 48 and 72 h after inoculation (seedling roots) | 36 |

Group labels in the result files: `drought_PEG_72h` and `PEG_72h` (PRJNA1041959) denote the samples taken 72 h after the soil water content had fallen to 30% (irrigation withheld; the experiment used no PEG). `parasitized_A` to `parasitized_E` and `Orobanche_A` to `Orobanche_E` (PRJNA850121) denote the roots sampled 4, 12, 24, 48 and 72 h after inoculation, each compared with the roots sampled at 0 h. In PRJNA492303 the ages `d16`, `d24` and `d30` are the three sampling times (day 1 and day 8 of flooding, day 6 of recovery). The supplementary tables and the figures give these descriptions instead of the labels.

Reference genome: HanXRQr2.0 (GCF_002127325.2) from NCBI RefSeq.

## Software

| Software | Version | Purpose |
|---|---|---|
| HMMER | 3.4 | PF00230 (MIP) domain search |
| BLAST+ | 2.17 | Candidate validation, orthologs, self-BLAST for duplications |
| DeepTMHMM | BioLib | Transmembrane helix prediction |
| MAFFT | 7.526 | Multiple sequence alignment |
| trimAl | 1.5.1 | Alignment trimming |
| IQ-TREE | 2.2.2.7 | Phylogenetic inference (ModelFinder, UFBoot, SH-aLRT) |
| MEME Suite | 5.5 | Protein motif discovery |
| MCScanX | bioconda | Collinearity detection (default parameters, minimum block size 5) |
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

## Citation

If you use this code or the derived tables, please cite the paper above. The individual tools should be cited separately; the versions used are listed in the software table.

## License

This project is licensed under the MIT License.

## Contact

Yavuz Delen - ydelen2@unl.edu
