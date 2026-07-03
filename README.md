# Sunflower Aquaporin Gene Family: Genome-Wide Characterization and Multi-Stress Transcriptomic Profiling

This repository contains all analysis scripts and pipelines used in:

**Delen, Y.** (2026). Genome-Wide Characterization of the Aquaporin Gene Family in Sunflower (*Helianthus annuus* L.) and Transcriptomic Profiling Across Eight Abiotic and Biotic Stress Conditions.

## Overview

We identified 61 aquaporin genes (71 protein isoforms) in the sunflower genome (HanXRQr2.0) and performed multi-stress transcriptomic profiling using six public RNA-seq datasets covering 16 stress contrasts across 8 stress types.

## Repository Structure

```
├── 00_setup/              # Environment setup and dependency installation
├── 01_references/         # Reference genome download and preparation
├── 02_gene_family/        # HMM search, BLAST, gene identification
├── 03_rnaseq/             # RNA-seq pipeline (QC, alignment, counting)
├── 04_expression/         # DESeq2 differential expression analysis
├── 05_phylogenetics/      # MAFFT, IQ-TREE phylogenetic analysis
├── 06_cis_elements/       # Promoter cis-element scanning
├── 07_synteny/            # MCScanX synteny analysis
├── config.sh              # Shared configuration (paths, parameters)
└── manuscript/            # Manuscript source (Markdown)
```

## Data Availability

All RNA-seq data are publicly available from NCBI SRA:

| BioProject | Stress Type | Samples |
|---|---|---|
| PRJNA869183 | Cold, Heat, Drought, Salt, Rehydration | 57 |
| PRJNA797473 | Drought (progressive) | 24 |
| PRJNA1041959 | Drought (PEG) | 24 |
| PRJNA492303 | Flooding | 46 |
| PRJNA908908 | *Sclerotinia sclerotiorum* | 12 |
| PRJNA850121 | *Orobanche cumana* | 36 |

Reference genome: HanXRQr2.0 (GCF_002127325.2) from NCBI.

## Software Dependencies

| Software | Version | Citation |
|---|---|---|
| HMMER | 3.4 | Eddy, 2011 |
| BLAST+ | 2.15 | Altschul et al., 1990 |
| MAFFT | 7.526 | Katoh and Standley, 2013 |
| IQ-TREE | 2 | Minh et al., 2020 |
| HISAT2 | 2.2 | Kim et al., 2019 |
| SAMtools | 1.23 | Danecek et al., 2021 |
| Subread/featureCounts | 2.1 | Liao et al., 2014 |
| DESeq2 | 1.38 | Love et al., 2014 |
| WGCNA | - | Langfelder and Horvath, 2008 |
| MCScanX | - | Wang et al., 2012 |
| MEME Suite | 5.5 | Bailey et al., 2015 |
| FastQC | 0.12 | Andrews, 2010 |
| fastp | 0.23 | Chen et al., 2018 |
| MultiQC | - | Ewels et al., 2016 |
| R | 4.1 | R Core Team, 2021 |
| Python | 3.12 | - |
| BioPython | - | Cock et al., 2009 |

## Computational Environment

All analyses were performed on the University of Nebraska-Lincoln SWAN HPC cluster using SLURM job scheduler. Conda environments were used for dependency management.

## License

This project is licensed under the MIT License.

## Contact

Yavuz Delen - ydelen2@unl.edu
