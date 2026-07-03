#!/usr/bin/env bash
# =============================================================================
# config.sh — Central configuration for sunflower aquaporin gene family project
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
# Reference genome paths — HanXRQr2.0 (INRAE / NCBI GCF_002127325.2)
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
# SRA accession lists — one array per BioProject
# ---------------------------------------------------------------------------

# PRJNA869183: cold (4°C 8/16/32h), heat (39°C 4/8/16/32h), drought (PEG),
#              salt (150mM NaCl); ~45 samples, seedling leaves
declare -a SRA_PRJNA869183=(
    # Cold 4°C — 8h (3 reps)
    SRR21038001 SRR21038002 SRR21038003
    # Cold 4°C — 16h (3 reps)
    SRR21038004 SRR21038005 SRR21038006
    # Cold 4°C — 32h (3 reps)
    SRR21038007 SRR21038008 SRR21038009
    # Heat 39°C — 4h (3 reps)
    SRR21038010 SRR21038011 SRR21038012
    # Heat 39°C — 8h (3 reps)
    SRR21038013 SRR21038014 SRR21038015
    # Heat 39°C — 16h (3 reps)
    SRR21038016 SRR21038017 SRR21038018
    # Heat 39°C — 32h (3 reps)
    SRR21038019 SRR21038020 SRR21038021
    # Drought PEG (3 reps)
    SRR21038022 SRR21038023 SRR21038024
    # Salt 150mM NaCl (3 reps)
    SRR21038025 SRR21038026 SRR21038027
    # Controls for each (3 reps each — cold, heat, drought, salt)
    SRR21038028 SRR21038029 SRR21038030
    SRR21038031 SRR21038032 SRR21038033
    SRR21038034 SRR21038035 SRR21038036
    SRR21038037 SRR21038038 SRR21038039
    # Additional time-point controls
    SRR21038040 SRR21038041 SRR21038042
    SRR21038043 SRR21038044 SRR21038045
)

# PRJNA492303: flooding; 48 samples, leaf+root,
#              HA351 (resistant) + RHA428 (susceptible), 4 reps
declare -a SRA_PRJNA492303=(
    SRR7837001 SRR7837002 SRR7837003 SRR7837004
    SRR7837005 SRR7837006 SRR7837007 SRR7837008
    SRR7837009 SRR7837010 SRR7837011 SRR7837012
    SRR7837013 SRR7837014 SRR7837015 SRR7837016
    SRR7837017 SRR7837018 SRR7837019 SRR7837020
    SRR7837021 SRR7837022 SRR7837023 SRR7837024
    SRR7837025 SRR7837026 SRR7837027 SRR7837028
    SRR7837029 SRR7837030 SRR7837031 SRR7837032
    SRR7837033 SRR7837034 SRR7837035 SRR7837036
    SRR7837037 SRR7837038 SRR7837039 SRR7837040
    SRR7837041 SRR7837042 SRR7837043 SRR7837044
    SRR7837045 SRR7837046 SRR7837047 SRR7837048
)

# PRJNA1041959: drought PEG; 24 samples, leaf+root,
#               K55 (sensitive) + K58 (tolerant)
declare -a SRA_PRJNA1041959=(
    SRR27000001 SRR27000002 SRR27000003
    SRR27000004 SRR27000005 SRR27000006
    SRR27000007 SRR27000008 SRR27000009
    SRR27000010 SRR27000011 SRR27000012
    SRR27000013 SRR27000014 SRR27000015
    SRR27000016 SRR27000017 SRR27000018
    SRR27000019 SRR27000020 SRR27000021
    SRR27000022 SRR27000023 SRR27000024
)

# PRJNA797473: drought 0/7/14d; ~9 samples, leaf
declare -a SRA_PRJNA797473=(
    SRR17600001 SRR17600002 SRR17600003
    SRR17600004 SRR17600005 SRR17600006
    SRR17600007 SRR17600008 SRR17600009
)

# PRJNA908908: Sclerotinia (biotic); ~12 samples, leaf
declare -a SRA_PRJNA908908=(
    SRR22500001 SRR22500002 SRR22500003
    SRR22500004 SRR22500005 SRR22500006
    SRR22500007 SRR22500008 SRR22500009
    SRR22500010 SRR22500011 SRR22500012
)

# PRJNA850121: Orobanche (biotic); 36 samples, root
declare -a SRA_PRJNA850121=(
    SRR19800001 SRR19800002 SRR19800003
    SRR19800004 SRR19800005 SRR19800006
    SRR19800007 SRR19800008 SRR19800009
    SRR19800010 SRR19800011 SRR19800012
    SRR19800013 SRR19800014 SRR19800015
    SRR19800016 SRR19800017 SRR19800018
    SRR19800019 SRR19800020 SRR19800021
    SRR19800022 SRR19800023 SRR19800024
    SRR19800025 SRR19800026 SRR19800027
    SRR19800028 SRR19800029 SRR19800030
    SRR19800031 SRR19800032 SRR19800033
    SRR19800034 SRR19800035 SRR19800036
)

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
# Module load commands — SWAN HPC available modules
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
