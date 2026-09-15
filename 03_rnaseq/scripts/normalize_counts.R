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
