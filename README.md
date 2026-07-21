# Sunflower Aquaporin Gene Family: Genome-Wide Characterization and Multi-Stress Transcriptomic Profiling

This repository contains all analysis scripts, supplementary data, and figures used in:

**Delen, Y., Delen, S.P., and Haliloglu, K.** (2026). Genome-Wide Characterization of the Aquaporin Gene Family in Sunflower (*Helianthus annuus* L.) and Transcriptomic Profiling Under Multiple Abiotic and Biotic Stress Conditions.

## Overview

We identified 61 aquaporin genes (71 protein isoforms) in the sunflower genome (HanXRQr2.0) and performed multi-stress transcriptomic profiling using six public RNA-seq datasets covering 16 stress contrasts across seven stress types (cold, heat, drought, salt, flooding, *Sclerotinia*, *Orobanche*).

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
├── hpc/                       # HPC output data and results
│   ├── 02_gene_family/        #   Gene family identification results
│   ├── 04_expression/         #   WGCNA module membership and trait tables
│   └── 05_phylogenetics/      #   IQ-TREE consensus and ML tree files
└── tables_and_figures/        # All publication tables and figures
    ├── figures/               #   Main figures (Fig 1-5, PDF + PNG)
    ├── tables/                #   Main tables (Table 1-4, TSV)
    ├── supplementary_figures/ #   Supplementary figures (Fig S1-S4, PDF + PNG)
    └── supplementary_tables/  #   Supplementary tables (S1-S10, TSV + Excel)
```

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
| IQ-TREE | 2 | Phylogenetic inference |
| HISAT2 | 2.2 | Read alignment |
| SAMtools | 1.23 | BAM processing |
| Subread/featureCounts | 2.1 | Read quantification |
| DESeq2 | 1.38 | Differential expression |
| WGCNA | -- | Co-expression networks |
| MCScanX | -- | Synteny detection |
| MEME Suite | 5.5 | Motif discovery |
| FastQC | 0.12 | Read quality control |
| fastp | 0.23 | Read trimming |
| MultiQC | -- | QC aggregation |
| clusterProfiler | 4.6 | GO enrichment |
| R | 4.1 | Statistical computing |
| Python | 3.12 | Data processing |
| BioPython | -- | Sequence analysis |

## Computational Environment

All computationally intensive analyses were performed on the University of Nebraska-Lincoln Holland Computing Center (HCC) Swan high-performance computing cluster using the SLURM job scheduler. Conda environments were used for dependency management.

## License

This project is licensed under the MIT License.

## Contact

Yavuz Delen - ydelen2@unl.edu
