#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/06_featurecounts_%j.out
#SBATCH --error=logs/06_featurecounts_%j.err

# =============================================================================
# 06_featurecounts.sh
# Gene-level read counting with featureCounts (subread)
# Processes ALL BAM files in a single run for a unified count matrix
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

module purge
module load "$MOD_SUBREAD"

init_dirs

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
[[ -f "$GTF_FILE" ]] || die "GTF annotation not found: ${GTF_FILE}"

# Collect all sorted BAM files
BAM_FILES=()
while IFS= read -r SAMPLE; do
    BP_ID=$(find_bioproject_for_srr "$SAMPLE")
    BAM="${BAM_DIR}/${BP_ID}/${SAMPLE}.sorted.bam"
    if [[ -f "$BAM" ]]; then
        BAM_FILES+=("$BAM")
    else
        log_msg "WARNING: Missing BAM for ${SAMPLE}: ${BAM}"
    fi
done < "${RNASEQ_DIR}/all_samples.txt"

NBAMS=${#BAM_FILES[@]}
log_msg "Found ${NBAMS} BAM files for counting"

[[ "$NBAMS" -gt 0 ]] || die "No BAM files found"

# ---------------------------------------------------------------------------
# Determine library layout from first BAM
# Check if the majority of reads are paired
# ---------------------------------------------------------------------------
FIRST_BAM="${BAM_FILES[0]}"
module load "$MOD_SAMTOOLS"
PAIRED_COUNT=$(samtools view -c -f 1 "$FIRST_BAM" 2>/dev/null | head -1 || echo 0)

PAIR_FLAG=""
if [[ "$PAIRED_COUNT" -gt 0 ]]; then
    PAIR_FLAG="-p --countReadPairs"
    log_msg "Detected paired-end data -> using -p --countReadPairs"
else
    log_msg "Detected single-end data"
fi

# ---------------------------------------------------------------------------
# featureCounts
# ---------------------------------------------------------------------------
RAW_COUNTS="${COUNTS_DIR}/gene_counts_raw.txt"
COUNTS_SUMMARY="${COUNTS_DIR}/gene_counts_raw.txt.summary"

log_msg "Running featureCounts..."

featureCounts \
    -a "$GTF_FILE" \
    -o "$RAW_COUNTS" \
    -T "$FEATURECOUNTS_THREADS" \
    $PAIR_FLAG \
    -s 0 \
    -g gene_id \
    -t exon \
    --tmpDir "${RNASEQ_DIR}/tmp_fc" \
    "${BAM_FILES[@]}"

log_msg "Raw counts written to: ${RAW_COUNTS}"

# ---------------------------------------------------------------------------
# Clean up column headers (strip path, keep sample ID)
# ---------------------------------------------------------------------------
CLEAN_COUNTS="${COUNTS_DIR}/gene_counts_clean.txt"

head -1 "$RAW_COUNTS" > "$CLEAN_COUNTS"  # comment line

# Process header: replace full BAM paths with sample names
python3 -c "
import re, sys

with open('${RAW_COUNTS}') as f:
    lines = f.readlines()

# Skip comment line (line 0), process header (line 1)
header = lines[1].strip().split('\t')
clean_header = []
for col in header:
    # Extract sample name from path like .../SRR12345.sorted.bam
    m = re.search(r'([SED]RR\d+)\.sorted\.bam', col)
    if m:
        clean_header.append(m.group(1))
    else:
        clean_header.append(col)

print('\t'.join(clean_header))

# Data lines
for line in lines[2:]:
    print(line.rstrip())
" >> "$CLEAN_COUNTS"

log_msg "Clean count matrix: ${CLEAN_COUNTS}"

# ---------------------------------------------------------------------------
# Generate TPM/FPKM normalization script
# ---------------------------------------------------------------------------
NORM_SCRIPT="${SCRIPTS_DIR}/normalize_counts.R"

cat > "$NORM_SCRIPT" << 'RSCRIPT'
#!/usr/bin/env Rscript
# =============================================================================
# normalize_counts.R
# Compute TPM and FPKM from featureCounts output
# Usage: Rscript normalize_counts.R <featurecounts_output> <output_dir>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript normalize_counts.R <featurecounts_file> <output_dir>")
}

counts_file <- args[1]
output_dir  <- args[2]

# Read featureCounts output (skip first comment line)
raw <- read.delim(counts_file, comment.char = "#", header = TRUE,
                  check.names = FALSE)

# Columns: Geneid, Chr, Start, End, Strand, Length, <sample1>, <sample2>, ...
gene_info <- raw[, 1:6]
count_mat <- as.matrix(raw[, 7:ncol(raw)])
rownames(count_mat) <- raw$Geneid
gene_lengths <- raw$Length  # effective gene length in bp

# Clean sample names from BAM paths
colnames(count_mat) <- gsub(".*/", "", colnames(count_mat))
colnames(count_mat) <- gsub("\\.sorted\\.bam$", "", colnames(count_mat))

cat("Loaded", nrow(count_mat), "genes x", ncol(count_mat), "samples\n")

# --- TPM calculation ---
# TPM_i = (count_i / length_i) / sum(count_j / length_j) * 1e6
rpk <- sweep(count_mat, 1, gene_lengths / 1000, "/")  # reads per kilobase
scaling_factors <- colSums(rpk)
tpm <- sweep(rpk, 2, scaling_factors / 1e6, "/")

# --- FPKM calculation ---
# FPKM_i = count_i / (length_i_kb * total_reads_M)
lib_sizes <- colSums(count_mat)
rpkm <- sweep(count_mat, 1, gene_lengths / 1000, "/")
fpkm <- sweep(rpkm, 2, lib_sizes / 1e6, "/")

# --- Write outputs ---
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tpm_out <- data.frame(Geneid = rownames(tpm), tpm, check.names = FALSE)
fpkm_out <- data.frame(Geneid = rownames(fpkm), fpkm, check.names = FALSE)
raw_out <- data.frame(Geneid = rownames(count_mat), count_mat, check.names = FALSE)

write.table(raw_out, file.path(output_dir, "gene_counts_raw_matrix.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tpm_out, file.path(output_dir, "gene_tpm_matrix.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fpkm_out, file.path(output_dir, "gene_fpkm_matrix.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Output files written to:", output_dir, "\n")
cat("  gene_counts_raw_matrix.tsv\n")
cat("  gene_tpm_matrix.tsv\n")
cat("  gene_fpkm_matrix.tsv\n")
RSCRIPT

chmod +x "$NORM_SCRIPT"
log_msg "Normalization script: ${NORM_SCRIPT}"

# ---------------------------------------------------------------------------
# Run normalization
# ---------------------------------------------------------------------------
module load "$MOD_MINIFORGE"

log_msg "Running TPM/FPKM normalization..."
Rscript "$NORM_SCRIPT" "$RAW_COUNTS" "$COUNTS_DIR"

log_msg "Done."
