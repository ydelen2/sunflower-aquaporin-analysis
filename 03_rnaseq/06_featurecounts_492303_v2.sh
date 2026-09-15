#!/bin/bash
#SBATCH --job-name=fc_492303
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --array=1-96%20
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/logs/06_fc_492303_%A_%a.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/logs/06_fc_492303_%A_%a.err

# =============================================================================
# 06_featurecounts_492303_v2.sh
# Gene-level counts for PRJNA492303, one BAM per array task, with the settings of
# 06_featurecounts.sh (-p --countReadPairs -s 0 -g gene_name -t exon). Counting
# the 96 BAMs in one process took about 28 min per BAM because paired reads are
# re-sorted internally, so the work is split across tasks; 06b_merge_counts_492303_v2.sh
# joins the per-sample files into one matrix afterwards.
# Output: counts/PRJNA492303_fc/<run>.txt (+ .summary)
# =============================================================================
set -eo pipefail
source /work/dweikat/ydelen2/aquaporin_study/03_rnaseq/config.sh

module purge
module load "$MOD_SUBREAD"
init_dirs

BP_ID="PRJNA492303"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACCESSION_DIR}/${BP_ID}_srr.txt")
[[ -z "$SAMPLE" ]] && die "No run at index ${SLURM_ARRAY_TASK_ID}"
[[ -f "$GTF_FILE" ]] || die "GTF annotation not found: ${GTF_FILE}"
BAM="${BAM_DIR}/${BP_ID}/${SAMPLE}.sorted.bam"
[[ -s "$BAM" && -s "${BAM}.bai" ]] || die "BAM missing for ${SAMPLE}"

OUT_DIR="${COUNTS_DIR}/${BP_ID}_fc"
mkdir -p "$OUT_DIR"
OUT="${OUT_DIR}/${SAMPLE}.txt"
if [[ -s "$OUT" && -s "${OUT}.summary" ]]; then
    log_msg "  SKIP ${SAMPLE}: counts already present"
    exit 0
fi

TMP="${RNASEQ_DIR}/tmp_fc_${SAMPLE}"
mkdir -p "$TMP"
log_msg "featureCounts ${SAMPLE}"
featureCounts \
    -a "$GTF_FILE" \
    -o "$OUT" \
    -T "$FEATURECOUNTS_THREADS" \
    -p --countReadPairs \
    -s 0 \
    -g gene_name \
    -t exon \
    --tmpDir "$TMP" \
    "$BAM"
rm -rf "$TMP"
grep -E "Assigned|Unassigned_NoFeatures|Unassigned_MultiMapping|Unassigned_Unmapped" "${OUT}.summary" | tr '\n' ' '; echo
log_msg "Done ${SAMPLE}"
