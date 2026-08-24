#!/usr/bin/env bash
#SBATCH --job-name=dl_genome
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/download_genome_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/download_genome_%j.err
# =============================================================================
# 01_download_genome.sh - Download reference genomes and build HISAT2 index
#
# Downloads:
#   1. Sunflower HanXRQr2.0 (NCBI GCF_002127325.2)
#   2. Arabidopsis TAIR10 proteome
#   3. Rice IRGSP-1.0 proteome
#   4. Tomato ITAG4.0 proteome
#   5. Lettuce v8 proteome
# Then builds HISAT2 index for sunflower.
#
# Submit: sbatch 01_download_genome.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

NCPUS="${SLURM_CPUS_PER_TASK:-16}"

# ---------------------------------------------------------------------------
# NCBI FTP base URLs
# ---------------------------------------------------------------------------
NCBI_FTP="https://ftp.ncbi.nlm.nih.gov/genomes/all"

# Sunflower HanXRQr2.0 - GCF_002127325.2
SUNFLOWER_BASE="${NCBI_FTP}/GCF/002/127/325/GCF_002127325.2_HanXRQr2.0-SUNRISE"
SUNFLOWER_FASTA="${SUNFLOWER_BASE}/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.fna.gz"
SUNFLOWER_GFF="${SUNFLOWER_BASE}/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.gff.gz"
SUNFLOWER_PROT="${SUNFLOWER_BASE}/GCF_002127325.2_HanXRQr2.0-SUNRISE_protein.faa.gz"
SUNFLOWER_CDS="${SUNFLOWER_BASE}/GCF_002127325.2_HanXRQr2.0-SUNRISE_cds_from_genomic.fna.gz"
SUNFLOWER_MD5="${SUNFLOWER_BASE}/md5checksums.txt"

# Arabidopsis thaliana TAIR10 - GCF_000001735.4
ATHA_BASE="${NCBI_FTP}/GCF/000/001/735/GCF_000001735.4_TAIR10.1"
ATHA_PROT_URL="${ATHA_BASE}/GCF_000001735.4_TAIR10.1_protein.faa.gz"

# Rice Oryza sativa IRGSP-1.0 - GCF_001433935.1
OSAT_BASE="${NCBI_FTP}/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0"
OSAT_PROT_URL="${OSAT_BASE}/GCF_001433935.1_IRGSP-1.0_protein.faa.gz"

# Tomato Solanum lycopersicum SL3.0 - GCF_000188115.5
SLYC_BASE="${NCBI_FTP}/GCF/000/188/115/GCF_000188115.5_SL3.1"
SLYC_PROT_URL="${SLYC_BASE}/GCF_000188115.5_SL3.1_protein.faa.gz"

# Lettuce Lactuca sativa v8 - GCF_002870075.4
LSAT_BASE="${NCBI_FTP}/GCF/002/870/075/GCF_002870075.4_Lsat_Salinas_v11"
LSAT_PROT_URL="${LSAT_BASE}/GCF_002870075.4_Lsat_Salinas_v11_protein.faa.gz"

# ============================================================================
# SECTION 1: Download sunflower genome
# ============================================================================
log_step "=== Downloading sunflower HanXRQr2.0 genome ==="

SF_DIR="${REF_DIR}/sunflower"
mkdir -p "${SF_DIR}"

# Download md5 checksums
safe_download "${SUNFLOWER_MD5}" "${SF_DIR}/md5checksums.txt"

# Genome FASTA
if [[ ! -f "${GENOME_FASTA}" ]]; then
    safe_download "${SUNFLOWER_FASTA}" "${SF_DIR}/genomic.fna.gz"
    log_info "Decompressing genome FASTA"
    gunzip -c "${SF_DIR}/genomic.fna.gz" > "${GENOME_FASTA}"
    log_info "Genome FASTA: $(wc -l < "${GENOME_FASTA}") lines"
else
    log_info "Genome FASTA already exists, skipping"
fi

# GFF annotation
if [[ ! -f "${GENOME_GFF}" ]]; then
    safe_download "${SUNFLOWER_GFF}" "${SF_DIR}/genomic.gff.gz"
    gunzip -c "${SF_DIR}/genomic.gff.gz" > "${GENOME_GFF}"
    log_info "GFF annotation: $(grep -c '^[^#]' "${GENOME_GFF}") feature lines"
else
    log_info "GFF already exists, skipping"
fi

# Protein sequences
if [[ ! -f "${GENOME_PROT}" ]]; then
    safe_download "${SUNFLOWER_PROT}" "${SF_DIR}/protein.faa.gz"
    gunzip -c "${SF_DIR}/protein.faa.gz" > "${GENOME_PROT}"
    log_info "Protein sequences: $(grep -c '^>' "${GENOME_PROT}") entries"
else
    log_info "Protein FASTA already exists, skipping"
fi

# CDS sequences
if [[ ! -f "${GENOME_CDS}" ]]; then
    safe_download "${SUNFLOWER_CDS}" "${SF_DIR}/cds_from_genomic.fna.gz"
    gunzip -c "${SF_DIR}/cds_from_genomic.fna.gz" > "${GENOME_CDS}"
    log_info "CDS sequences: $(grep -c '^>' "${GENOME_CDS}") entries"
else
    log_info "CDS FASTA already exists, skipping"
fi

# Verify md5 checksums for downloaded .gz files
log_step "Verifying MD5 checksums for sunflower downloads"
cd "${SF_DIR}"
for gz_file in genomic.fna.gz genomic.gff.gz protein.faa.gz cds_from_genomic.fna.gz; do
    if [[ -f "${gz_file}" ]]; then
        expected=$(grep "${gz_file}" md5checksums.txt | awk '{print $1}' | head -1)
        if [[ -n "${expected}" ]]; then
            actual=$(md5sum "${gz_file}" | awk '{print $1}')
            if [[ "${expected}" == "${actual}" ]]; then
                log_info "MD5 OK: ${gz_file}"
            else
                log_error "MD5 MISMATCH: ${gz_file} (expected ${expected}, got ${actual})"
                exit 1
            fi
        else
            log_warn "No MD5 found in checksums file for ${gz_file}"
        fi
    fi
done

# Convert GFF to GTF using gffread
log_step "Converting GFF to GTF"
if [[ ! -f "${GENOME_GTF}" ]]; then
    load_gffread
    gffread "${GENOME_GFF}" -T -o "${GENOME_GTF}"
    log_info "GTF created: $(wc -l < "${GENOME_GTF}") lines"
else
    log_info "GTF already exists, skipping"
fi

# ============================================================================
# SECTION 2: Download comparative proteomes
# ============================================================================
log_step "=== Downloading comparative proteomes ==="

# Arabidopsis thaliana
ATHA_DIR="${REF_DIR}/arabidopsis"
mkdir -p "${ATHA_DIR}"
if [[ ! -f "${ATHA_PROT}" ]]; then
    safe_download "${ATHA_PROT_URL}" "${ATHA_DIR}/protein.faa.gz"
    gunzip -c "${ATHA_DIR}/protein.faa.gz" > "${ATHA_PROT}"
    log_info "Arabidopsis proteome: $(grep -c '^>' "${ATHA_PROT}") proteins"
else
    log_info "Arabidopsis proteome already exists"
fi

# Rice (Oryza sativa)
OSAT_DIR="${REF_DIR}/rice"
mkdir -p "${OSAT_DIR}"
if [[ ! -f "${OSAT_PROT}" ]]; then
    safe_download "${OSAT_PROT_URL}" "${OSAT_DIR}/protein.faa.gz"
    gunzip -c "${OSAT_DIR}/protein.faa.gz" > "${OSAT_PROT}"
    log_info "Rice proteome: $(grep -c '^>' "${OSAT_PROT}") proteins"
else
    log_info "Rice proteome already exists"
fi

# Tomato (Solanum lycopersicum)
SLYC_DIR="${REF_DIR}/tomato"
mkdir -p "${SLYC_DIR}"
if [[ ! -f "${SLYC_PROT}" ]]; then
    safe_download "${SLYC_PROT_URL}" "${SLYC_DIR}/protein.faa.gz"
    gunzip -c "${SLYC_DIR}/protein.faa.gz" > "${SLYC_PROT}"
    log_info "Tomato proteome: $(grep -c '^>' "${SLYC_PROT}") proteins"
else
    log_info "Tomato proteome already exists"
fi

# Lettuce (Lactuca sativa)
LSAT_DIR="${REF_DIR}/lettuce"
mkdir -p "${LSAT_DIR}"
if [[ ! -f "${LSAT_PROT}" ]]; then
    safe_download "${LSAT_PROT_URL}" "${LSAT_DIR}/protein.faa.gz"
    gunzip -c "${LSAT_DIR}/protein.faa.gz" > "${LSAT_PROT}"
    log_info "Lettuce proteome: $(grep -c '^>' "${LSAT_PROT}") proteins"
else
    log_info "Lettuce proteome already exists"
fi

# ============================================================================
# SECTION 3: Build HISAT2 index
# ============================================================================
log_step "=== Building HISAT2 index for sunflower genome ==="

HISAT2_DIR="${REF_DIR}/sunflower/hisat2_index"
mkdir -p "${HISAT2_DIR}"

if [[ ! -f "${HISAT2_INDEX}.1.ht2" ]]; then
    load_hisat2

    # Extract splice sites and exons from GTF for splice-aware alignment
    log_info "Extracting splice sites from GTF"
    hisat2_extract_splice_sites.py "${GENOME_GTF}" > "${HISAT2_DIR}/splice_sites.tsv"
    log_info "Splice sites: $(wc -l < "${HISAT2_DIR}/splice_sites.tsv")"

    log_info "Extracting exons from GTF"
    hisat2_extract_exons.py "${GENOME_GTF}" > "${HISAT2_DIR}/exons.tsv"
    log_info "Exons: $(wc -l < "${HISAT2_DIR}/exons.tsv")"

    # Build index with splice site and exon information
    log_info "Building HISAT2 index (this will take a while for ~3.6 Gb genome)"
    hisat2-build \
        --seed 42 \
        -p "${NCPUS}" \
        --ss "${HISAT2_DIR}/splice_sites.tsv" \
        --exon "${HISAT2_DIR}/exons.tsv" \
        "${GENOME_FASTA}" \
        "${HISAT2_INDEX}"

    log_info "HISAT2 index files:"
    ls -lh "${HISAT2_INDEX}"*.ht2
else
    log_info "HISAT2 index already exists, skipping"
fi

# ============================================================================
# SECTION 4: Generate genome summary statistics
# ============================================================================
log_step "=== Generating genome summary statistics ==="

load_samtools

# Create FASTA index
if [[ ! -f "${GENOME_FASTA}.fai" ]]; then
    samtools faidx "${GENOME_FASTA}"
fi

# Genome size summary
log_info "Genome assembly statistics:"
awk '{sum+=$2; n++} END {
    printf "  Sequences: %d\n", n;
    printf "  Total length: %d bp (%.2f Gb)\n", sum, sum/1e9;
}' "${GENOME_FASTA}.fai"

# Gene/mRNA/CDS counts from GFF
log_info "Annotation statistics:"
for feature in gene mRNA CDS exon; do
    count=$(awk -v f="${feature}" '$3==f' "${GENOME_GFF}" | wc -l)
    printf "  %s: %d\n" "${feature}" "${count}"
done

# Clean up compressed files to save space (optional - uncomment if needed)
# log_info "Cleaning up compressed files"
# rm -f "${SF_DIR}"/*.gz "${ATHA_DIR}"/*.gz "${OSAT_DIR}"/*.gz \
#        "${SLYC_DIR}"/*.gz "${LSAT_DIR}"/*.gz

log_done "All reference genomes downloaded and indexed"
log_info "Sunflower genome: ${GENOME_FASTA}"
log_info "HISAT2 index: ${HISAT2_INDEX}"
log_info "Comparative proteomes: ${REF_DIR}/{arabidopsis,rice,tomato,lettuce}/"
