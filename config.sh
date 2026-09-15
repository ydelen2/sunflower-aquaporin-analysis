#!/usr/bin/env bash
# =============================================================================
# config.sh - Central configuration for sunflower aquaporin gene family project
# SWAN HPC, University of Nebraska-Lincoln
# Project: /work/dweikat/ydelen2/aquaporin_study
# =============================================================================

# ---------------------------------------------------------------------------
# Strict mode
# ---------------------------------------------------------------------------
set -euo pipefail

# ---------------------------------------------------------------------------
# Project directory structure
# ---------------------------------------------------------------------------
export PROJ_DIR="/work/dweikat/ydelen2/aquaporin_study"
export SETUP_DIR="${PROJ_DIR}/00_setup"
export REF_DIR="${PROJ_DIR}/01_references"
export GENE_FAM_DIR="${PROJ_DIR}/02_gene_family"
export RNASEQ_DIR="${PROJ_DIR}/03_rnaseq"
export FASTQ_DIR="${RNASEQ_DIR}/fastq"
export QC_DIR="${RNASEQ_DIR}/qc"
export TRIMMED_DIR="${RNASEQ_DIR}/trimmed"
export ALIGNED_DIR="${RNASEQ_DIR}/aligned"
export COUNTS_DIR="${RNASEQ_DIR}/counts"
export EXPR_DIR="${PROJ_DIR}/04_expression"
export PHYLO_DIR="${PROJ_DIR}/05_phylogenetics"
export CIS_DIR="${PROJ_DIR}/06_cis_elements"
export SYNTENY_DIR="${PROJ_DIR}/07_synteny"
export STRUCT_DIR="${PROJ_DIR}/08_protein_structure"
export FIG_DIR="${PROJ_DIR}/09_figures"
export LOG_DIR="${PROJ_DIR}/logs"

# Local R library
export R_LIBS_USER="${PROJ_DIR}/R_libs"

# Conda environment
export CONDA_ENV_NAME="aquaporin_env"
export CONDA_ENV_PREFIX="${PROJ_DIR}/conda_envs/${CONDA_ENV_NAME}"

# ---------------------------------------------------------------------------
# Reference genome paths - HanXRQr2.0 (INRAE / NCBI GCF_002127325.2)
# ---------------------------------------------------------------------------
export GENOME_FASTA="${REF_DIR}/sunflower/GCF_002127325.2_HanXRQr2.0_genomic.fna"
export GENOME_GFF="${REF_DIR}/sunflower/GCF_002127325.2_HanXRQr2.0_genomic.gff"
export GENOME_GTF="${REF_DIR}/sunflower/GCF_002127325.2_HanXRQr2.0_genomic.gtf"
export GENOME_PROT="${REF_DIR}/sunflower/GCF_002127325.2_HanXRQr2.0_protein.faa"
export GENOME_CDS="${REF_DIR}/sunflower/GCF_002127325.2_HanXRQr2.0_cds_from_genomic.fna"
export HISAT2_INDEX="${REF_DIR}/sunflower/hisat2_index/HanXRQr2"

# Comparative proteomes
export ATHA_PROT="${REF_DIR}/arabidopsis/Athaliana_TAIR10_proteome.faa"
export OSAT_PROT="${REF_DIR}/rice/Osativa_IRGSP1_proteome.faa"
export SLYC_PROT="${REF_DIR}/tomato/Slycopersicum_ITAG4_proteome.faa"
export LSAT_PROT="${REF_DIR}/lettuce/Lsativa_v8_proteome.faa"

# ---------------------------------------------------------------------------
# SRA run lists
# ---------------------------------------------------------------------------
# The SRR accessions of each BioProject are fetched from NCBI by
# 03_rnaseq/02_get_srr_ids.sh and stored in 03_rnaseq/accession_lists/
# (<BioProject>_srr.txt, one run per line). The lists used for the paper are
# tracked in that directory.

# Combined list of all BioProject IDs
declare -a BIOPROJECTS=(
    PRJNA869183
    PRJNA492303
    PRJNA1041959
    PRJNA797473
    PRJNA908908
    PRJNA850121
)

# ---------------------------------------------------------------------------
# Module load commands - SWAN HPC available modules
# ---------------------------------------------------------------------------
load_hisat2()       { module load hisat2/2.2;       }
load_star()         { module load star/2.7.9a;      }
load_subread()      { module load subread/2.1;      }
load_blast()        { module load blast/2.17;       }
load_mafft()        { module load mafft/7.526;      }
load_iqtree()       { module load iqtree/3.1;       }
load_hmmer()        { module load hmmer/3.4;        }
load_meme()         { module load MEME/5.5;         }
load_fastqc()       { module load fastqc/0.12;      }
load_fastp()        { module load fastp/0.23;       }
load_samtools()     { module load samtools/1.23;    }
load_r()            { module load R/4.1;            }
load_trimmomatic()  { module load trimmomatic/0.39; }
load_bedtools()     { module load bedtools/2.31;    }
load_alphafold3()   { module load alphafold3/3.0;   }
load_miniforge()    { module load miniforge/24.5;   }
load_multiqc()      { module load py37/1.8;         }  # multiqc available under py37
load_gffread()      { module load gffread/0.12;     }
load_muscle()       { module load muscle/5.1;       }
load_java()         { module load java/19;          }

# ---------------------------------------------------------------------------
# Common SLURM settings
# ---------------------------------------------------------------------------
export SLURM_PARTITION="batch"
export SLURM_DEFAULT_TIME="24:00:00"
export SLURM_DEFAULT_MEM="16G"
export SLURM_DEFAULT_CPUS=8
export SLURM_MAIL_TYPE="END,FAIL"
export SLURM_MAIL_USER="ydelen2@huskers.unl.edu"

# ---------------------------------------------------------------------------
# Color codes for logging
# ---------------------------------------------------------------------------
export RED='\033[0;31m'
export GREEN='\033[0;32m'
export YELLOW='\033[1;33m'
export BLUE='\033[0;34m'
export MAGENTA='\033[0;35m'
export CYAN='\033[0;36m'
export WHITE='\033[1;37m'
export NC='\033[0m'  # No Color

# ---------------------------------------------------------------------------
# Logging functions
# ---------------------------------------------------------------------------
timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

log_info() {
    echo -e "${GREEN}[INFO  $(timestamp)]${NC} $*"
}

log_warn() {
    echo -e "${YELLOW}[WARN  $(timestamp)]${NC} $*"
}

log_error() {
    echo -e "${RED}[ERROR $(timestamp)]${NC} $*" >&2
}

log_step() {
    echo -e "${CYAN}[STEP  $(timestamp)]${NC} $*"
}

log_done() {
    echo -e "${MAGENTA}[DONE  $(timestamp)]${NC} $*"
}

# Trap handler for clean error reporting
trap_error() {
    log_error "Script failed at line $1 (exit code $2)"
    exit "$2"
}
trap 'trap_error ${LINENO} $?' ERR

# ---------------------------------------------------------------------------
# Utility: check if a command/module is available
# ---------------------------------------------------------------------------
check_cmd() {
    command -v "$1" &>/dev/null || {
        log_error "Required command not found: $1"
        return 1
    }
}

# ---------------------------------------------------------------------------
# Utility: safe download with retry
# ---------------------------------------------------------------------------
safe_download() {
    local url="$1"
    local outfile="$2"
    local max_retries="${3:-3}"
    local attempt=1

    while [[ $attempt -le $max_retries ]]; do
        log_info "Downloading (attempt ${attempt}/${max_retries}): ${url}"
        if wget --no-verbose --tries=3 --timeout=60 -O "${outfile}" "${url}"; then
            log_info "Download complete: ${outfile}"
            return 0
        fi
        log_warn "Download attempt ${attempt} failed"
        ((attempt++))
        sleep 10
    done

    log_error "Failed to download after ${max_retries} attempts: ${url}"
    return 1
}
