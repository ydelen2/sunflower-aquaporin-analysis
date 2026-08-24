#!/bin/bash
#SBATCH --job-name=qc_trim
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-249%30
#SBATCH --output=logs/03_qc_trim_%A_%a.out
#SBATCH --error=logs/03_qc_trim_%A_%a.err

# =============================================================================
# 03_qc_trim.sh
# Per-sample QC and adapter/quality trimming
# Array index = line number in all_samples.txt (1-based)
# Throttled to 30 concurrent tasks to avoid I/O saturation
# =============================================================================

set -eo pipefail

SCRIPT_DIR="/work/dweikat/ydelen2/aquaporin_study"
source "/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/config.sh"

module purge
module load "$MOD_FASTQC"
module load "$MOD_FASTP"

init_dirs

# ---------------------------------------------------------------------------
# Resolve sample for this array task
# ---------------------------------------------------------------------------
SAMPLE=$(get_sample_by_index "$SLURM_ARRAY_TASK_ID")
[[ -z "$SAMPLE" ]] && die "No sample at index ${SLURM_ARRAY_TASK_ID}"

BP_ID=$(find_bioproject_for_srr "$SAMPLE")
[[ "$BP_ID" == "UNKNOWN" ]] && die "Cannot find BioProject for ${SAMPLE}"

log_msg "Processing sample: ${SAMPLE} (BioProject: ${BP_ID})"

# ---------------------------------------------------------------------------
# Locate input FASTQ files
# ---------------------------------------------------------------------------
RAW_DIR="${FASTQ_DIR}/${BP_ID}"

# Detect paired-end vs single-end
R1="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
R2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"
SE="${RAW_DIR}/${SAMPLE}.fastq.gz"

if [[ -f "$R1" && -f "$R2" ]]; then
    LAYOUT="PE"
    log_msg "  Layout: paired-end"
elif [[ -f "$SE" ]]; then
    LAYOUT="SE"
    log_msg "  Layout: single-end"
else
    die "No FASTQ files found for ${SAMPLE} in ${RAW_DIR}"
fi

# ---------------------------------------------------------------------------
# Output directories (per-BioProject subdirs for organization)
# ---------------------------------------------------------------------------
TRIM_OUT="${TRIMMED_DIR}/${BP_ID}"
QC_RAW_OUT="${FASTQC_RAW_DIR}/${BP_ID}"
QC_TRIM_OUT="${FASTQC_TRIM_DIR}/${BP_ID}"
FASTP_OUT="${FASTP_DIR}/${BP_ID}"

mkdir -p "$TRIM_OUT" "$QC_RAW_OUT" "$QC_TRIM_OUT" "$FASTP_OUT"

# ---------------------------------------------------------------------------
# Step 1: FastQC on raw reads
# ---------------------------------------------------------------------------
log_msg "  Running FastQC on raw reads..."

if [[ "$LAYOUT" == "PE" ]]; then
    fastqc -t "$FASTP_THREADS" -o "$QC_RAW_OUT" "$R1" "$R2"
else
    fastqc -t "$FASTP_THREADS" -o "$QC_RAW_OUT" "$SE"
fi

# ---------------------------------------------------------------------------
# Step 2: fastp trimming
# ---------------------------------------------------------------------------
log_msg "  Running fastp..."

TRIM_R1="${TRIM_OUT}/${SAMPLE}_1.trimmed.fastq.gz"
TRIM_R2="${TRIM_OUT}/${SAMPLE}_2.trimmed.fastq.gz"
TRIM_SE="${TRIM_OUT}/${SAMPLE}.trimmed.fastq.gz"
FASTP_HTML="${FASTP_OUT}/${SAMPLE}_fastp.html"
FASTP_JSON="${FASTP_OUT}/${SAMPLE}_fastp.json"

if [[ "$LAYOUT" == "PE" ]]; then
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
else
    fastp \
        --in1 "$SE" \
        --out1 "$TRIM_SE" \
        --qualified_quality_phred "$FASTP_QUAL" \
        --length_required "$FASTP_MINLEN" \
        --thread "$FASTP_THREADS" \
        --html "$FASTP_HTML" \
        --json "$FASTP_JSON" \
        --report_title "${SAMPLE}"
fi

# ---------------------------------------------------------------------------
# Step 3: FastQC on trimmed reads
# ---------------------------------------------------------------------------
log_msg "  Running FastQC on trimmed reads..."

if [[ "$LAYOUT" == "PE" ]]; then
    fastqc -t "$FASTP_THREADS" -o "$QC_TRIM_OUT" "$TRIM_R1" "$TRIM_R2"
else
    fastqc -t "$FASTP_THREADS" -o "$QC_TRIM_OUT" "$TRIM_SE"
fi

# ---------------------------------------------------------------------------
# Step 4: Extract quick stats from fastp JSON
# ---------------------------------------------------------------------------
STATS_FILE="${STATS_DIR}/qc_stats.tsv"

# Write header once (race-safe: check if file is empty)
if [[ ! -s "$STATS_FILE" ]]; then
    echo -e "sample\tbioproject\tlayout\traw_reads\ttrimmed_reads\tpct_passed\tgc_before\tgc_after\tq30_before\tq30_after" \
        > "$STATS_FILE" 2>/dev/null || true
fi

# Parse fastp JSON with python (available on most systems)
python3 -c "
import json, sys
with open('${FASTP_JSON}') as f:
    d = json.load(f)
s = d['summary']
raw = s['before_filtering']['total_reads']
trim = s['after_filtering']['total_reads']
pct = round(100.0 * trim / raw, 2) if raw > 0 else 0
gc_b = round(s['before_filtering']['gc_content'] * 100, 2)
gc_a = round(s['after_filtering']['gc_content'] * 100, 2)
q30_b = round(s['before_filtering']['q30_rate'] * 100, 2)
q30_a = round(s['after_filtering']['q30_rate'] * 100, 2)
print(f'${SAMPLE}\t${BP_ID}\t${LAYOUT}\t{raw}\t{trim}\t{pct}\t{gc_b}\t{gc_a}\t{q30_b}\t{q30_a}')
" >> "$STATS_FILE"

log_msg "Done: ${SAMPLE}"
