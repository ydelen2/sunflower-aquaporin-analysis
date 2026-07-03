#!/bin/bash
#SBATCH --job-name=hisat2_align
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=1-174%20
#SBATCH --output=logs/05_hisat2_%A_%a.out
#SBATCH --error=logs/05_hisat2_%A_%a.err

# =============================================================================
# 05_hisat2_align.sh
# HISAT2 alignment of trimmed reads to HanXRQr2.0
# Array index = line number in all_samples.txt (1-based)
# Throttled to 20 concurrent to manage memory pressure
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

module purge
module load "$MOD_HISAT2"
module load "$MOD_SAMTOOLS"

init_dirs

# ---------------------------------------------------------------------------
# Resolve sample
# ---------------------------------------------------------------------------
SAMPLE=$(get_sample_by_index "$SLURM_ARRAY_TASK_ID")
[[ -z "$SAMPLE" ]] && die "No sample at index ${SLURM_ARRAY_TASK_ID}"

BP_ID=$(find_bioproject_for_srr "$SAMPLE")
[[ "$BP_ID" == "UNKNOWN" ]] && die "Cannot find BioProject for ${SAMPLE}"

log_msg "Aligning sample: ${SAMPLE} (BioProject: ${BP_ID})"

# ---------------------------------------------------------------------------
# Locate trimmed FASTQ files
# ---------------------------------------------------------------------------
TRIM_DIR="${TRIMMED_DIR}/${BP_ID}"

R1="${TRIM_DIR}/${SAMPLE}_1.trimmed.fastq.gz"
R2="${TRIM_DIR}/${SAMPLE}_2.trimmed.fastq.gz"
SE="${TRIM_DIR}/${SAMPLE}.trimmed.fastq.gz"

if [[ -f "$R1" && -f "$R2" ]]; then
    LAYOUT="PE"
elif [[ -f "$SE" ]]; then
    LAYOUT="SE"
else
    die "No trimmed FASTQ files found for ${SAMPLE} in ${TRIM_DIR}"
fi

# ---------------------------------------------------------------------------
# Verify HISAT2 index exists
# ---------------------------------------------------------------------------
[[ -f "${HISAT2_INDEX}.1.ht2" || -f "${HISAT2_INDEX}.1.ht2l" ]] \
    || die "HISAT2 index not found at ${HISAT2_INDEX}"

# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------
BAM_OUT="${BAM_DIR}/${BP_ID}"
mkdir -p "$BAM_OUT"

SORTED_BAM="${BAM_OUT}/${SAMPLE}.sorted.bam"
ALIGN_LOG="${BAM_OUT}/${SAMPLE}.hisat2.log"
FLAGSTAT="${BAM_OUT}/${SAMPLE}.flagstat.txt"

# Skip if sorted BAM already exists and is valid
if [[ -f "$SORTED_BAM" && -f "${SORTED_BAM}.bai" ]]; then
    log_msg "  SKIP (exists): ${SORTED_BAM}"
    exit 0
fi

# ---------------------------------------------------------------------------
# HISAT2 alignment
# ---------------------------------------------------------------------------
log_msg "  Running HISAT2..."

TMPDIR="${RNASEQ_DIR}/tmp_${SAMPLE}"
mkdir -p "$TMPDIR"

if [[ "$LAYOUT" == "PE" ]]; then
    hisat2 \
        -x "$HISAT2_INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        --dta \
        -p "$HISAT2_THREADS" \
        --new-summary \
        --summary-file "$ALIGN_LOG" \
        --rg-id "${SAMPLE}" \
        --rg "SM:${SAMPLE}" \
        --rg "PL:ILLUMINA" \
        --rg "LB:${BP_ID}" \
    2>> "${BAM_OUT}/${SAMPLE}.hisat2.stderr" \
    | samtools sort \
        -@ 4 \
        -m 2G \
        -T "${TMPDIR}/${SAMPLE}" \
        -o "$SORTED_BAM" \
        -
else
    hisat2 \
        -x "$HISAT2_INDEX" \
        -U "$SE" \
        --dta \
        -p "$HISAT2_THREADS" \
        --new-summary \
        --summary-file "$ALIGN_LOG" \
        --rg-id "${SAMPLE}" \
        --rg "SM:${SAMPLE}" \
        --rg "PL:ILLUMINA" \
        --rg "LB:${BP_ID}" \
    2>> "${BAM_OUT}/${SAMPLE}.hisat2.stderr" \
    | samtools sort \
        -@ 4 \
        -m 2G \
        -T "${TMPDIR}/${SAMPLE}" \
        -o "$SORTED_BAM" \
        -
fi

# ---------------------------------------------------------------------------
# Index BAM
# ---------------------------------------------------------------------------
log_msg "  Indexing BAM..."
samtools index -@ 4 "$SORTED_BAM"

# ---------------------------------------------------------------------------
# Flagstat
# ---------------------------------------------------------------------------
log_msg "  Running flagstat..."
samtools flagstat -@ 4 "$SORTED_BAM" > "$FLAGSTAT"

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
rm -rf "$TMPDIR"

# ---------------------------------------------------------------------------
# Quick alignment stats extraction
# ---------------------------------------------------------------------------
ALIGN_STATS="${STATS_DIR}/alignment_stats.tsv"

# Write header once
if [[ ! -s "$ALIGN_STATS" ]]; then
    echo -e "sample\tbioproject\tlayout\ttotal_reads\taligned_reads\talignment_rate\tuniquely_aligned" \
        > "$ALIGN_STATS" 2>/dev/null || true
fi

# Parse HISAT2 new-summary format
python3 -c "
import re, sys

with open('${ALIGN_LOG}') as f:
    text = f.read()

# Extract total reads
total_m = re.search(r'Total reads:\s+(\d+)', text) or re.search(r'(\d+) reads', text)
total = int(total_m.group(1)) if total_m else 0

# Extract alignment rate
rate_m = re.search(r'Overall alignment rate:\s+([\d.]+)%', text) or re.search(r'([\d.]+)% overall', text)
rate = float(rate_m.group(1)) if rate_m else 0.0

# Extract uniquely aligned
uniq_m = re.search(r'Aligned 1 time:\s+(\d+)', text) or re.search(r'(\d+).*aligned exactly 1 time', text)
uniq = int(uniq_m.group(1)) if uniq_m else 0

aligned = int(total * rate / 100)

print(f'${SAMPLE}\t${BP_ID}\t${LAYOUT}\t{total}\t{aligned}\t{rate}\t{uniq}')
" >> "$ALIGN_STATS"

log_msg "Done: ${SAMPLE} (alignment rate logged)"
