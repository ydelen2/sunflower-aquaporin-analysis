#!/bin/bash
#SBATCH --job-name=hisat2_492303
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=1-96%20
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/logs/05_hisat2_492303_%A_%a.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/logs/05_hisat2_492303_%A_%a.err

# =============================================================================
# 05_hisat2_492303_v2.sh
# HISAT2 alignment of every PRJNA492303 run (96), same parameters as
# 05_hisat2_align.sh. Runs that already have a sorted, indexed BAM are skipped.
# Array index = line number in accession_lists/PRJNA492303_srr.txt.
# =============================================================================
set -eo pipefail
source /work/dweikat/ydelen2/aquaporin_study/03_rnaseq/config.sh

module purge
module load "$MOD_HISAT2"
module load "$MOD_SAMTOOLS"
init_dirs

BP_ID="PRJNA492303"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACCESSION_DIR}/${BP_ID}_srr.txt")
[[ -z "$SAMPLE" ]] && die "No run at index ${SLURM_ARRAY_TASK_ID}"
log_msg "Sample ${SAMPLE} (${BP_ID}), array index ${SLURM_ARRAY_TASK_ID}"

TRIM_DIR="${TRIMMED_DIR}/${BP_ID}"
R1="${TRIM_DIR}/${SAMPLE}_1.trimmed.fastq.gz"
R2="${TRIM_DIR}/${SAMPLE}_2.trimmed.fastq.gz"
[[ -s "$R1" && -s "$R2" ]] || die "Trimmed FASTQ missing for ${SAMPLE} in ${TRIM_DIR}"

BAM_OUT="${BAM_DIR}/${BP_ID}"
mkdir -p "$BAM_OUT"
SORTED_BAM="${BAM_OUT}/${SAMPLE}.sorted.bam"
ALIGN_LOG="${BAM_OUT}/${SAMPLE}.hisat2.log"
FLAGSTAT="${BAM_OUT}/${SAMPLE}.flagstat.txt"

if [[ -s "$SORTED_BAM" && -s "${SORTED_BAM}.bai" && -s "$ALIGN_LOG" ]]; then
    log_msg "  SKIP ${SAMPLE}: BAM already present"
    exit 0
fi

TMPDIR="${RNASEQ_DIR}/tmp_${SAMPLE}"
mkdir -p "$TMPDIR"

log_msg "  hisat2"
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

log_msg "  index and flagstat"
samtools index -@ 4 "$SORTED_BAM"
samtools flagstat -@ 4 "$SORTED_BAM" > "$FLAGSTAT"
rm -rf "$TMPDIR"

log_msg "Done ${SAMPLE}"
