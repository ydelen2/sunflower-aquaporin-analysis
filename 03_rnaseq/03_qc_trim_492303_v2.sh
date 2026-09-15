#!/bin/bash
#SBATCH --job-name=trim_492303
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-96%20
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/logs/03_trim_492303_%A_%a.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/logs/03_trim_492303_%A_%a.err

# =============================================================================
# 03_qc_trim_492303_v2.sh
# FastQC + fastp for every run of PRJNA492303 (flooding, 96 runs).
# The first analysis trimmed and aligned only 46 of the 96 runs, an arbitrary
# block of consecutive accessions; the full experiment is a balanced
# 2 genotypes x 2 treatments x 2 tissues x 3 ages x 4 replicates design.
# Parameters are identical to 03_qc_trim.sh. Runs already trimmed are skipped.
# Array index = line number in accession_lists/PRJNA492303_srr.txt.
# =============================================================================
set -eo pipefail
source /work/dweikat/ydelen2/aquaporin_study/03_rnaseq/config.sh

module purge
module load "$MOD_FASTQC"
module load "$MOD_FASTP"
init_dirs

BP_ID="PRJNA492303"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACCESSION_DIR}/${BP_ID}_srr.txt")
[[ -z "$SAMPLE" ]] && die "No run at index ${SLURM_ARRAY_TASK_ID}"
log_msg "Sample ${SAMPLE} (${BP_ID}), array index ${SLURM_ARRAY_TASK_ID}"

RAW_DIR="${FASTQ_DIR}/${BP_ID}"
R1="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
R2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"
[[ -s "$R1" && -s "$R2" ]] || die "Raw FASTQ missing for ${SAMPLE} in ${RAW_DIR}"

TRIM_OUT="${TRIMMED_DIR}/${BP_ID}"
QC_RAW_OUT="${FASTQC_RAW_DIR}/${BP_ID}"
QC_TRIM_OUT="${FASTQC_TRIM_DIR}/${BP_ID}"
FASTP_OUT="${FASTP_DIR}/${BP_ID}"
mkdir -p "$TRIM_OUT" "$QC_RAW_OUT" "$QC_TRIM_OUT" "$FASTP_OUT"

TRIM_R1="${TRIM_OUT}/${SAMPLE}_1.trimmed.fastq.gz"
TRIM_R2="${TRIM_OUT}/${SAMPLE}_2.trimmed.fastq.gz"
FASTP_HTML="${FASTP_OUT}/${SAMPLE}_fastp.html"
FASTP_JSON="${FASTP_OUT}/${SAMPLE}_fastp.json"

if [[ -s "$TRIM_R1" && -s "$TRIM_R2" && -s "$FASTP_JSON" ]]; then
    log_msg "  SKIP ${SAMPLE}: trimmed files already present"
    exit 0
fi

if [[ ! -s "${QC_RAW_OUT}/${SAMPLE}_1_fastqc.zip" ]]; then
    log_msg "  FastQC on raw reads"
    fastqc -t "$FASTP_THREADS" -o "$QC_RAW_OUT" "$R1" "$R2"
fi

log_msg "  fastp"
fastp \
    --in1 "$R1" \
    --in2 "$R2" \
    --out1 "$TRIM_R1" \
    --out2 "$TRIM_R2" \
    --qualified_quality_phred "$FASTP_QUAL" \
    --length_required "$FASTP_MINLEN" \
    --detect_adapter_for_pe \
    --thread "$FASTP_THREADS" \
    --html "$FASTP_HTML" \
    --json "$FASTP_JSON" \
    --report_title "${SAMPLE}"

log_msg "  FastQC on trimmed reads"
fastqc -t "$FASTP_THREADS" -o "$QC_TRIM_OUT" "$TRIM_R1" "$TRIM_R2"

log_msg "Done ${SAMPLE}"
