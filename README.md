# Sunflower Aquaporin Gene Family: Genome-Wide Characterization and Multi-Stress Transcriptomic Profiling

This repository holds the analysis code, supplementary data and figures behind:

**Delen, Y., Delen, S.P., and Haliloglu, K.** (2026). Genome-Wide Characterization of the Aquaporin Gene Family in Sunflower (*Helianthus annuus* L.) and Transcriptomic Profiling Under Multiple Abiotic and Biotic Stress Conditions.

## Overview

We identified 61 aquaporin genes (71 protein isoforms) in the sunflower genome (HanXRQr2.0) and profiled their expression across six public RNA-seq datasets, covering 16 stress contrasts and seven stress types (cold, heat, drought, salt, flooding, *Sclerotinia*, *Orobanche*).

Everything needed to reproduce the analysis is here: the scripts that were actually run, the result tables they produced, and the figures used in the paper. Raw FASTQ files, BAM files and other bulky intermediates are deliberately not tracked (see `.gitignore`), since they can be regenerated from the SRA accessions listed below.

## Repository Structure

```
├── 00_setup/                  # Environment setup and dependency installation
├── 01_references/             # Reference genome download and preparation
├── 02_gene_family/            # HMM search, BLAST, gene identification
├── 03_rnaseq/                 # RNA-seq pipeline (QC, alignment, counting)
├── 04_expression/             # DESeq2, WGCNA, GO enrichment
├── 05_phylogenetics/          # MAFFT alignment, IQ-TREE phylogenetic analysis
├── 06_cis_elements/           # Promoter cis-element scanning
├── 07_synteny/                # MCScanX synteny analysis
├── config.sh                  # Shared configuration (paths, parameters)
├── results/                       # Pipeline output data
│   ├── 02_gene_family/        #   Gene family identification results
│   ├── 04_expression/         #   WGCNA module membership and trait tables
│   └── 05_phylogenetics/      #   IQ-TREE consensus and ML tree files
└── tables_and_figures/        # All publication tables and figures
    ├── figures/               #   Main figures (Fig 1-5, PDF + PNG)
    ├── tables/                #   Main tables (Table 1-4, TSV)
    ├── supplementary_figures/ #   Supplementary figures (Fig S1-S4, PDF + PNG)
    └── supplementary_tables/  #   Supplementary tables (S1-S10, TSV + Excel)
```

## Running the pipeline

The scripts were written for the SLURM scheduler on the HCC Swan cluster. Paths, module versions and sample lists live in `config.sh` (project-wide) and `03_rnaseq/config.sh` (RNA-seq specific), so start by editing those to match your own system.

Only `00_setup/01_create_dirs.sh` is a plain shell script; everything else is submitted with `sbatch`. The stages are meant to be run in order, and each one expects the previous stage to have finished.

```bash
# 1. Environment and reference data
bash   00_setup/01_create_dirs.sh
sbatch 00_setup/02_install_r_packages.sh
sbatch 00_setup/03_setup_conda.sh
sbatch 01_references/01_download_genome.sh

# 2. Aquaporin gene family identification
sbatch 02_gene_family/01_hmm_search.sh
sbatch 02_gene_family/02_blast_validation.sh
sbatch 02_gene_family/03_domain_verification.sh
sbatch 02_gene_family/04_characterization.sh
sbatch 02_gene_family/05_gene_structure.sh
sbatch 02_gene_family/06_chromosomal_location.sh

# 3. RNA-seq processing
sbatch 03_rnaseq/02_get_srr_ids.sh     # run this before 01_download_sra.sh
sbatch 03_rnaseq/01_download_sra.sh    # array job, one task per BioProject
sbatch 03_rnaseq/03_qc_trim.sh         # array job, one task per sample
sbatch 03_rnaseq/04_multiqc.sh
sbatch 03_rnaseq/05_hisat2_align.sh    # array job, one task per sample
sbatch 03_rnaseq/06_featurecounts.sh
sbatch 03_rnaseq/07_alignment_stats.sh

# 4. Expression analysis
sbatch 04_expression/01_deseq2_analysis.sh
sbatch 04_expression/02_aquaporin_expression.sh
sbatch 04_expression/03_wgcna_analysis.sh
sbatch 04_expression/04_go_kegg_enrichment.sh

# 5. Phylogenetics, motifs, cis-elements and synteny
sbatch 05_phylogenetics/01_msa_and_tree.sh
sbatch 05_phylogenetics/02_motif_analysis.sh
sbatch 06_cis_elements/01_cis_element_analysis.sh
sbatch 07_synteny/01_synteny_analysis.sh
sbatch 07_synteny/03_reparse_synteny.sh
sbatch 07_synteny/04_synteny_figure.sh
```

A few things worth knowing before starting:

- `02_get_srr_ids.sh` has to run before `01_download_sra.sh`, because the download job reads the accession lists that the first script fetches from NCBI. The numbering is a leftover from an earlier version of the pipeline.
- The array sizes in `03_qc_trim.sh` and `05_hisat2_align.sh` are upper bounds; tasks past the end of `all_samples.txt` simply exit.
- Downloading the six BioProjects takes by far the longest and is the step most likely to need restarting. The download and trimming jobs skip samples whose output already exists, so re-submitting them is safe.

## Figures

### Main Figures (5)

| Figure | Description |
|---|---|
| Fig 1 | Chromosomal distribution of 61 aquaporin genes across 15 sunflower chromosomes |
| Fig 2 | Synteny relationships: intra-genomic and interspecies conservation |
| Fig 3 | Summary bar chart of upregulated/downregulated aquaporin DEGs per stress contrast |
| Fig 4 | Volcano plots for selected stress contrasts |
| Fig 5 | WGCNA module-trait correlation heatmap |

### Supplementary Figures (4)

| Figure | Description |
|---|---|
| Fig S1 | Maximum-likelihood phylogenetic tree of 228 aquaporin proteins (5 species) |
| Fig S2 | Expression heatmap of 61 aquaporin genes across 16 stress contrasts |
| Fig S3 | UpSet plot showing overlap of aquaporin DEGs across seven stress types |
| Fig S4 | Dot plot of broadly responsive aquaporin genes across 16 stress contrasts |

### Supplementary Tables (12 sheets in `Supplementary_Tables.xlsx`)

S1 (gene list), S2 (physicochemical), S3a (duplication), S3b (synteny), S3c (Ka/Ks), S4 (DEG matrix), S5 (WGCNA membership), S6 (module-trait), S7 (GO enrichment), S8 (MEME motifs), S9 (ortholog mapping), S10 (DESeq2 design rationale).

## Data Availability

All RNA-seq data are publicly available from NCBI SRA:

| BioProject | Stress Type | Samples |
|---|---|---|
| PRJNA869183 | Cold, Heat, Drought, Salt, Rehydration | 57 |
| PRJNA797473 | Drought (progressive, 7/14/21d) | 24 |
| PRJNA1041959 | Drought (PEG 72h) | 24 |
| PRJNA492303 | Flooding | 46 |
| PRJNA908908 | *Sclerotinia sclerotiorum* | 12 |
| PRJNA850121 | *Orobanche cumana* (stages A-E) | 36 |

Reference genome: HanXRQr2.0 (GCF_002127325.2) from NCBI.

## Software Dependencies

| Software | Version | Purpose |
|---|---|---|
| HMMER | 3.4 | Domain search (PF00230) |
| BLAST+ | 2.15 | Reciprocal validation |
| MAFFT | 7.526 | Multiple sequence alignment |
| trimAl | bioconda | Alignment trimming |
| IQ-TREE | 2.2 | Phylogenetic inference |
| HISAT2 | 2.2 | Read alignment |
| SAMtools | 1.23 | BAM processing |
| Subread/featureCounts | 2.1 | Read quantification |
| DESeq2 | 1.38 | Differential expression |
| WGCNA | Bioconductor | Co-expression networks |
| MCScanX | bioconda | Synteny detection |
| MEME Suite | 5.5 | Motif discovery |
| FastQC | 0.12 | Read quality control |
| fastp | 0.23 | Read trimming |
| MultiQC | py37/1.8 module | QC aggregation |
| clusterProfiler | 4.6 | GO enrichment |
| R | 4.1 | Statistical computing |
| Python | 3.10 | Data processing |
| BioPython | conda | Sequence analysis |

Numbered versions are the cluster modules that were loaded. The rest were installed at their then-current release: R packages through BiocManager in `00_setup/02_install_r_packages.sh`, Python and command-line tools through conda in `00_setup/03_setup_conda.sh`. trimAl is the exception; `05_phylogenetics/01_msa_and_tree.sh` installs it into its own conda environment the first time it runs.

## Computational Environment

All computationally intensive analyses were performed on the University of Nebraska-Lincoln Holland Computing Center (HCC) Swan high-performance computing cluster using the SLURM job scheduler. Conda environments were used for dependency management.

## Citation

If you use this code or the derived tables, please cite the paper above. The individual tools should be cited separately; the versions used are listed in the dependency table.

## License

This project is licensed under the MIT License.

## Contact

Yavuz Delen - ydelen2@unl.edu
