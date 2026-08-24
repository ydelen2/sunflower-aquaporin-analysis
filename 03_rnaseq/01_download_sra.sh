#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --array=0-5
#SBATCH --output=logs/01_download_%A_%a.out
#SBATCH --error=logs/01_download_%A_%a.err

# =============================================================================
# 01_download_sra.sh
# Download SRA data for each BioProject via SLURM array
# Array index 0-5 maps to the 6 BioProjects defined in config.sh
# Prerequisites: Run 02_get_srr_ids.sh first to populate accession_lists/
# =============================================================================

set -eo pipefail

SCRIPT_DIR="/work/dweikat/ydelen2/aquaporin_study"
source "/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/config.sh"

module purge

init_dirs

# ---------------------------------------------------------------------------
# Resolve BioProject for this array task
# ---------------------------------------------------------------------------
TASK_IDX=${SLURM_ARRAY_TASK_ID}
ENTRY="${BIOPROJECTS[$TASK_IDX]}"
parse_bioproject "$ENTRY"

ACC_FILE="${ACCESSION_DIR}/${BP_ID}_srr.txt"
OUT_DIR="${FASTQ_DIR}/${BP_ID}"
PREFETCH_DIR="${RNASEQ_DIR}/sra_cache/${BP_ID}"

mkdir -p "$OUT_DIR" "$PREFETCH_DIR"

log_msg "Task ${TASK_IDX}: Downloading ${BP_ID} (${BP_NAME})"

# Validate accession list exists
[[ -f "$ACC_FILE" ]] || die "Accession file not found: ${ACC_FILE}. Run 02_get_srr_ids.sh first."

TOTAL=$(wc -l < "$ACC_FILE")
log_msg "  ${TOTAL} accessions to download"

# ---------------------------------------------------------------------------
# Download loop: prefetch -> validate -> fasterq-dump
# ---------------------------------------------------------------------------
DOWNLOADED=0
FAILED=0
SKIPPED=0

while IFS= read -r SRR; do
    [[ -z "$SRR" ]] && continue

    # Skip if FASTQ files already exist
    if ls "${OUT_DIR}/${SRR}"_*.fastq.gz 1>/dev/null 2>&1; then
        log_msg "  SKIP (exists): ${SRR}"
        SKIPPED=$((SKIPPED+1))
        continue
    fi

    log_msg "  Downloading: ${SRR}"

    # Step 1: prefetch
    if ! prefetch "$SRR" \
            --output-directory "$PREFETCH_DIR" \
            --max-size 50G \
            --progress 2>&1; then
        log_msg "  FAILED prefetch: ${SRR}"
        FAILED=$((FAILED+1))
        continue
    fi

    # Step 2: vdb-validate
    SRA_FILE="${PREFETCH_DIR}/${SRR}/${SRR}.sra"
    if [[ -f "$SRA_FILE" ]]; then
        if ! vdb-validate "$SRA_FILE" 2>&1; then
            log_msg "  FAILED validation: ${SRR}"
            FAILED=$((FAILED+1))
            continue
        fi
    fi

    # Step 3: fasterq-dump (split into R1/R2 for paired-end)
    if ! fasterq-dump "$SRR" \
            --outdir "$OUT_DIR" \
            --split-3 \
            --threads 1 \
            --temp "${RNASEQ_DIR}/sra_cache/tmp_${SRR}" \
            --progress 2>&1; then
        log_msg "  FAILED fasterq-dump: ${SRR}"
        FAILED=$((FAILED+1))
        continue
    fi

    # Step 4: Compress FASTQ files
    log_msg "  Compressing: ${SRR}"
    for fq in "${OUT_DIR}/${SRR}"*.fastq; do
        if [[ -f "$fq" ]]; then
            pigz -p 4 "$fq"
        fi
    done

    # Cleanup prefetch cache for this sample
    rm -rf "${PREFETCH_DIR}/${SRR}" "${RNASEQ_DIR}/sra_cache/tmp_${SRR}"

    DOWNLOADED=$((DOWNLOADED+1))
    log_msg "  DONE: ${SRR} (${DOWNLOADED}/${TOTAL})"

done < "$ACC_FILE"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log_msg "=== Download Summary for ${BP_ID} ==="
log_msg "  Total accessions: ${TOTAL}"
log_msg "  Downloaded: ${DOWNLOADED}"
log_msg "  Skipped (existing): ${SKIPPED}"
log_msg "  Failed: ${FAILED}"

if [[ "$FAILED" -gt 0 ]]; then
    log_msg "WARNING: ${FAILED} downloads failed for ${BP_ID}"
    exit 1
fi

log_msg "Done: ${BP_ID}"
