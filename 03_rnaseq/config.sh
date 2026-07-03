#!/bin/bash
# =============================================================================
# config.sh - Shared configuration for sunflower aquaporin RNA-seq pipeline
# Project: Aquaporin gene family analysis in Helianthus annuus (HanXRQr2.0)
# Cluster: SWAN HPC, University of Nebraska-Lincoln
# =============================================================================

# -- Base paths ---------------------------------------------------------------
export PROJECT_DIR="/work/dweikat/ydelen2/aquaporin_study"
export RNASEQ_DIR="${PROJECT_DIR}/03_rnaseq"
export REF_DIR="${PROJECT_DIR}/01_references"
export GENOME_DIR="${REF_DIR}/genome"
export HISAT2_INDEX="${REF_DIR}/hisat2_index/HanXRQr2"
export GTF_FILE="${GENOME_DIR}/HanXRQr2.0.gtf"

# -- RNA-seq subdirectories ---------------------------------------------------
export FASTQ_DIR="${RNASEQ_DIR}/fastq"
export TRIMMED_DIR="${RNASEQ_DIR}/trimmed"
export FASTQC_RAW_DIR="${RNASEQ_DIR}/fastqc_raw"
export FASTQC_TRIM_DIR="${RNASEQ_DIR}/fastqc_trimmed"
export FASTP_DIR="${RNASEQ_DIR}/fastp_reports"
export BAM_DIR="${RNASEQ_DIR}/bam"
export COUNTS_DIR="${RNASEQ_DIR}/counts"
export MULTIQC_DIR="${RNASEQ_DIR}/multiqc"
export STATS_DIR="${RNASEQ_DIR}/stats"
export LOG_DIR="${RNASEQ_DIR}/logs"
export ACCESSION_DIR="${RNASEQ_DIR}/accession_lists"
export SCRIPTS_DIR="${RNASEQ_DIR}/scripts"

# -- BioProjects --------------------------------------------------------------
# Format: BIOPROJECT_ID:SHORT_NAME:DESCRIPTION
export BIOPROJECTS=(
    "PRJNA869183:abiotic_seedling:cold_heat_drought_salt_stress_seedling_leaves"
    "PRJNA492303:flooding:flooding_leaf_root_2genotypes_4reps"
    "PRJNA1041959:drought_lr:drought_leaf_root_2genotypes"
    "PRJNA797473:drought_leaf:drought_leaf_only"
    "PRJNA908908:sclerotinia:sclerotinia_biotic_leaf"
    "PRJNA850121:orobanche:orobanche_biotic_root"
)

# Approximate sample counts (for SLURM array sizing)
declare -A SAMPLE_COUNTS
SAMPLE_COUNTS[PRJNA869183]=45
SAMPLE_COUNTS[PRJNA492303]=48
SAMPLE_COUNTS[PRJNA1041959]=24
SAMPLE_COUNTS[PRJNA797473]=9
SAMPLE_COUNTS[PRJNA908908]=12
SAMPLE_COUNTS[PRJNA850121]=36

export TOTAL_SAMPLES=174

# -- Module versions ----------------------------------------------------------
export MOD_HISAT2="hisat2/2.2"
export MOD_STAR="star/2.7.9a"
export MOD_SUBREAD="subread/2.1"
export MOD_FASTQC="fastqc/0.12"
export MOD_FASTP="fastp/0.23"
export MOD_SAMTOOLS="samtools/1.23"
export MOD_MULTIQC="py37/1.8"
export MOD_TRIMMOMATIC="trimmomatic/0.39"
export MOD_MINIFORGE="miniforge/24.5"
export MOD_ENTREZ="entrez-direct/16.2"
export MOD_SRATOOLKIT="sratoolkit/3.0"

# -- Pipeline parameters ------------------------------------------------------
export FASTP_QUAL=20
export FASTP_MINLEN=50
export FASTP_THREADS=8
export HISAT2_THREADS=8
export FEATURECOUNTS_THREADS=8

# Alignment QC thresholds
export MIN_MAPPING_RATE=70
export MIN_TOTAL_READS=5000000

# -- Helper functions ---------------------------------------------------------

# Create all output directories
init_dirs() {
    local dirs=(
        "$FASTQ_DIR" "$TRIMMED_DIR" "$FASTQC_RAW_DIR" "$FASTQC_TRIM_DIR"
        "$FASTP_DIR" "$BAM_DIR" "$COUNTS_DIR" "$MULTIQC_DIR" "$STATS_DIR"
        "$LOG_DIR" "$ACCESSION_DIR" "$SCRIPTS_DIR"
    )
    for d in "${dirs[@]}"; do
        mkdir -p "$d"
    done
}

# Parse BioProject entry -> sets BP_ID, BP_NAME, BP_DESC
parse_bioproject() {
    local entry="$1"
    BP_ID=$(echo "$entry" | cut -d: -f1)
    BP_NAME=$(echo "$entry" | cut -d: -f2)
    BP_DESC=$(echo "$entry" | cut -d: -f3)
}

# Build master sample list from all accession files
# Output: one SRR ID per line to RNASEQ_DIR/all_samples.txt
build_master_sample_list() {
    local outfile="${RNASEQ_DIR}/all_samples.txt"
    > "$outfile"
    for entry in "${BIOPROJECTS[@]}"; do
        parse_bioproject "$entry"
        local acc_file="${ACCESSION_DIR}/${BP_ID}_srr.txt"
        if [[ -f "$acc_file" ]]; then
            cat "$acc_file" >> "$outfile"
        else
            echo "WARNING: Missing accession file: $acc_file" >&2
        fi
    done
    echo "Master sample list: $(wc -l < "$outfile") samples -> $outfile"
}

# Get sample ID from master list by 1-based index
get_sample_by_index() {
    local idx=$1
    local list="${RNASEQ_DIR}/all_samples.txt"
    if [[ ! -f "$list" ]]; then
        echo "ERROR: Master sample list not found: $list" >&2
        return 1
    fi
    sed -n "${idx}p" "$list"
}

# Find which BioProject a given SRR belongs to
find_bioproject_for_srr() {
    local srr="$1"
    for entry in "${BIOPROJECTS[@]}"; do
        parse_bioproject "$entry"
        local acc_file="${ACCESSION_DIR}/${BP_ID}_srr.txt"
        if [[ -f "$acc_file" ]] && grep -qw "$srr" "$acc_file"; then
            echo "$BP_ID"
            return 0
        fi
    done
    echo "UNKNOWN"
    return 1
}

# Logging with timestamp
log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Error handler
die() {
    log_msg "ERROR: $*" >&2
    exit 1
}
