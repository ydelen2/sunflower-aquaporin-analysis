#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/04_multiqc_%j.out
#SBATCH --error=logs/04_multiqc_%j.err

# =============================================================================
# 04_multiqc.sh
# Aggregate all QC reports (FastQC raw, FastQC trimmed, fastp) into
# a unified MultiQC report.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

module purge
module load "$MOD_MULTIQC"

init_dirs

log_msg "Running MultiQC..."

# ---------------------------------------------------------------------------
# MultiQC on raw FastQC + fastp outputs
# ---------------------------------------------------------------------------
multiqc \
    "$FASTQC_RAW_DIR" \
    "$FASTQC_TRIM_DIR" \
    "$FASTP_DIR" \
    --outdir "$MULTIQC_DIR" \
    --filename "multiqc_qc_report" \
    --title "Aquaporin Study - QC Summary" \
    --force \
    --no-data-dir

log_msg "MultiQC QC report: ${MULTIQC_DIR}/multiqc_qc_report.html"

# ---------------------------------------------------------------------------
# If alignment stats exist, run a second MultiQC for alignment
# ---------------------------------------------------------------------------
if [[ -d "$BAM_DIR" ]] && ls "${BAM_DIR}"/*/*.hisat2.log 1>/dev/null 2>&1; then
    log_msg "Running MultiQC on alignment logs..."

    multiqc \
        "$BAM_DIR" \
        "$COUNTS_DIR" \
        --outdir "$MULTIQC_DIR" \
        --filename "multiqc_alignment_report" \
        --title "Aquaporin Study - Alignment Summary" \
        --force \
        --no-data-dir

    log_msg "MultiQC alignment report: ${MULTIQC_DIR}/multiqc_alignment_report.html"
fi

log_msg "Done."
