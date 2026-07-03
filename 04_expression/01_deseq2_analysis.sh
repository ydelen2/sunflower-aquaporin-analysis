#!/bin/bash
#SBATCH --job-name=deseq2_aqp
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/deseq2_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/deseq2_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 01_deseq2_analysis.sh
# DESeq2 differential expression analysis for sunflower aquaporin study
# Handles multiple BioProjects with separate analysis per project
###############################################################################

set -euo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
EXPR_DIR="${PROJ_DIR}/04_expression"
COUNTS_DIR="${PROJ_DIR}/03_rnaseq/counts"
RESULTS_DIR="${EXPR_DIR}/results/deseq2"
LOG_DIR="${EXPR_DIR}/logs"
R_LIBS="${PROJ_DIR}/R_libs"

mkdir -p "${RESULTS_DIR}" "${LOG_DIR}"

echo "=== DESeq2 Analysis ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-local}"

module purge
module load R/4.1

# ---- Embedded R script ----
Rscript --no-save --no-restore - <<'RSCRIPT_EOF'

###############################################################################
# DESeq2 Differential Expression Analysis
###############################################################################

# -- Paths --
PROJ_DIR   <- "/work/dweikat/ydelen2/aquaporin_study"
COUNTS_DIR <- file.path(PROJ_DIR, "03_rnaseq/counts")
RESULTS_DIR <- file.path(PROJ_DIR, "04_expression/results/deseq2")
R_LIBS     <- file.path(PROJ_DIR, "R_libs")
AQP_LIST   <- file.path(PROJ_DIR, "02_gene_family/results/verified_aquaporins.txt")

.libPaths(c(R_LIBS, .libPaths()))

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(tidyverse)
    library(pheatmap)
})

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)

cat("R session info:\n")
cat("R version:", R.version.string, "\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")
cat("Threads:", parallel::detectCores(), "\n\n")

# Register parallel backend for DESeq2
BPPARAM <- BiocParallel::MulticoreParam(workers = 8)
register(BPPARAM)

###############################################################################
# 1. Load featureCounts output
###############################################################################

count_files <- list.files(COUNTS_DIR, pattern = "\\.counts$|\\.txt$",
                          full.names = TRUE, recursive = TRUE)

if (length(count_files) == 0) {
    stop("No count files found in ", COUNTS_DIR)
}

cat("Found", length(count_files), "count files\n")

# featureCounts summary format: first column = Geneid, then sample columns
# If single merged counts matrix exists, load it; otherwise merge individual files.
merged_file <- file.path(COUNTS_DIR, "all_counts.txt")

if (file.exists(merged_file)) {
    cat("Loading merged count matrix:", merged_file, "\n")
    raw <- read.delim(merged_file, comment.char = "#", row.names = 1)
    # featureCounts includes Chr, Start, End, Strand, Length before counts
    meta_cols <- c("Chr", "Start", "End", "Strand", "Length")
    gene_lengths <- raw$Length  # save for TPM later
    names(gene_lengths) <- rownames(raw)
    count_mat <- as.matrix(raw[, !(colnames(raw) %in% meta_cols)])
} else {
    cat("No merged file found. Loading individual count files and merging.\n")
    count_list <- list()
    gene_lengths <- NULL

    for (f in count_files) {
        tab <- read.delim(f, comment.char = "#", row.names = 1)
        meta_cols <- c("Chr", "Start", "End", "Strand", "Length")
        if (is.null(gene_lengths) && "Length" %in% colnames(tab)) {
            gene_lengths <- tab$Length
            names(gene_lengths) <- rownames(tab)
        }
        counts_only <- tab[, !(colnames(tab) %in% meta_cols), drop = FALSE]
        count_list[[f]] <- counts_only
    }
    count_mat <- do.call(cbind, count_list)
    count_mat <- as.matrix(count_mat)
}

# Clean sample names: strip path prefixes and .bam/.sorted suffixes
colnames(count_mat) <- gsub(".*\\/", "", colnames(count_mat))
colnames(count_mat) <- gsub("\\.sorted\\.bam$|\\.bam$", "", colnames(count_mat))

cat("Count matrix dimensions:", nrow(count_mat), "genes x",
    ncol(count_mat), "samples\n")
cat("Sample names:\n")
cat(paste(" ", colnames(count_mat)), sep = "\n")
cat("\n")

# Save gene lengths for TPM calculation downstream
if (!is.null(gene_lengths)) {
    saveRDS(gene_lengths, file.path(RESULTS_DIR, "gene_lengths.rds"))
}

###############################################################################
# 2. Build sample metadata from filenames
###############################################################################
# Expected naming convention: {BioProject}_{StressType}_{Tissue}_{Genotype}_{Rep}
# Adjust parsing logic below to match actual naming scheme.

sample_names <- colnames(count_mat)

# Attempt to parse structured sample names
parse_sample_name <- function(sname) {
    parts <- strsplit(sname, "_")[[1]]
    if (length(parts) >= 5) {
        data.frame(
            sample     = sname,
            bioproject = parts[1],
            stress     = parts[2],
            tissue     = parts[3],
            genotype   = parts[4],
            replicate  = parts[5],
            stringsAsFactors = FALSE
        )
    } else if (length(parts) >= 3) {
        # Minimal: project_stress_rep
        data.frame(
            sample     = sname,
            bioproject = parts[1],
            stress     = parts[2],
            tissue     = "leaf",
            genotype   = "HA",
            replicate  = parts[3],
            stringsAsFactors = FALSE
        )
    } else {
        data.frame(
            sample     = sname,
            bioproject = "unknown",
            stress     = "unknown",
            tissue     = "unknown",
            genotype   = "unknown",
            replicate  = "1",
            stringsAsFactors = FALSE
        )
    }
}

coldata <- do.call(rbind, lapply(sample_names, parse_sample_name))
rownames(coldata) <- coldata$sample

# If a separate sample_metadata.tsv exists, prefer it
meta_file <- file.path(PROJ_DIR, "03_rnaseq/sample_metadata.tsv")
if (file.exists(meta_file)) {
    cat("Loading external metadata from:", meta_file, "\n")
    coldata <- read.delim(meta_file, row.names = 1, stringsAsFactors = FALSE)
    # Ensure row order matches count matrix
    common <- intersect(colnames(count_mat), rownames(coldata))
    if (length(common) == 0) stop("No overlap between count matrix columns and metadata rows")
    count_mat <- count_mat[, common]
    coldata   <- coldata[common, , drop = FALSE]
}

cat("Metadata:\n")
print(str(coldata))

# Determine condition column (stress/condition)
cond_col <- if ("stress" %in% colnames(coldata)) "stress" else
            if ("condition" %in% colnames(coldata)) "condition" else
            stop("Metadata must have 'stress' or 'condition' column")

coldata[[cond_col]] <- factor(coldata[[cond_col]])

# Save metadata
write.table(coldata, file.path(RESULTS_DIR, "sample_metadata.tsv"),
            sep = "\t", quote = FALSE)
saveRDS(count_mat, file.path(RESULTS_DIR, "raw_counts.rds"))

###############################################################################
# 3. DESeq2 analysis per BioProject
###############################################################################

run_deseq2_project <- function(counts, metadata, project_id, cond_col) {

    out_dir <- file.path(RESULTS_DIR, project_id)
    plot_dir <- file.path(out_dir, "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    cat("\n========================================\n")
    cat("Processing BioProject:", project_id, "\n")
    cat("Samples:", nrow(metadata), "\n")
    cat("Conditions:", paste(levels(factor(metadata[[cond_col]])), collapse = ", "), "\n")
    cat("========================================\n")

    metadata[[cond_col]] <- factor(metadata[[cond_col]])

    # Identify control condition
    control_terms <- c("control", "ctrl", "ck", "mock", "normal", "unstressed")
    conditions <- levels(metadata[[cond_col]])
    control_cond <- conditions[tolower(conditions) %in% control_terms]

    if (length(control_cond) == 0) {
        cat("WARNING: No control condition auto-detected. Using first level as reference.\n")
        control_cond <- conditions[1]
    } else {
        control_cond <- control_cond[1]
    }
    cat("Reference (control) condition:", control_cond, "\n")

    metadata[[cond_col]] <- relevel(metadata[[cond_col]], ref = control_cond)

    # Build design formula
    # If multiple tissues/genotypes present, include them
    design_formula <- as.formula(paste0("~ ", cond_col))
    if (length(unique(metadata$tissue)) > 1 && "tissue" %in% colnames(metadata)) {
        design_formula <- as.formula(paste0("~ tissue + ", cond_col))
    }
    cat("Design formula:", deparse(design_formula), "\n")

    # Pre-filter: remove genes with < 10 total counts
    keep <- rowSums(counts) >= 10
    cat("Genes before filtering:", nrow(counts), "\n")
    counts <- counts[keep, ]
    cat("Genes after filtering:", nrow(counts), "\n")

    # Construct DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData   = metadata,
                                  design    = design_formula)

    # Run DESeq2
    dds <- DESeq(dds, parallel = TRUE)
    saveRDS(dds, file.path(out_dir, "dds.rds"))

    # Normalized counts (for downstream use)
    norm_counts <- counts(dds, normalized = TRUE)
    write.table(norm_counts, file.path(out_dir, "normalized_counts.tsv"),
                sep = "\t", quote = FALSE)
    saveRDS(norm_counts, file.path(out_dir, "normalized_counts.rds"))

    # VST for visualization
    vsd <- vst(dds, blind = FALSE)
    saveRDS(vsd, file.path(out_dir, "vsd.rds"))

    # ---- PCA Plot ----
    pca_data <- plotPCA(vsd, intgroup = cond_col, returnData = TRUE)
    pct_var  <- round(100 * attr(pca_data, "percentVar"))

    p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[cond_col]])) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        ggtitle(paste0("PCA - ", project_id)) +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5))

    ggsave(file.path(plot_dir, "pca_plot.pdf"), p_pca, width = 8, height = 6)
    ggsave(file.path(plot_dir, "pca_plot.png"), p_pca, width = 8, height = 6, dpi = 300)

    # ---- Per-contrast results ----
    stress_conditions <- setdiff(levels(metadata[[cond_col]]), control_cond)
    all_results <- list()

    for (stress in stress_conditions) {
        contrast_name <- paste0(stress, "_vs_", control_cond)
        cat("  Contrast:", contrast_name, "\n")

        res <- results(dds, contrast = c(cond_col, stress, control_cond),
                       alpha = 0.05, pAdjustMethod = "BH")
        res <- res[order(res$padj), ]

        # Summary
        cat("    Total tested:", sum(!is.na(res$padj)), "\n")
        cat("    DEGs (padj<0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
        cat("    Up:", sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE), "\n")
        cat("    Down:", sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE), "\n")

        res_df <- as.data.frame(res)
        res_df$gene_id <- rownames(res_df)
        res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

        write.table(res_df, file.path(out_dir, paste0("DEG_", contrast_name, ".tsv")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        saveRDS(res, file.path(out_dir, paste0("DEG_", contrast_name, ".rds")))

        # DEG list (padj < 0.05, |log2FC| > 1)
        sig <- res_df[!is.na(res_df$padj) &
                      res_df$padj < 0.05 &
                      abs(res_df$log2FoldChange) > 1, ]
        write.table(sig, file.path(out_dir, paste0("DEG_sig_", contrast_name, ".tsv")),
                    sep = "\t", row.names = FALSE, quote = FALSE)

        all_results[[contrast_name]] <- res

        # ---- MA Plot ----
        pdf(file.path(plot_dir, paste0("MA_", contrast_name, ".pdf")), width = 8, height = 6)
        plotMA(res, main = paste0("MA Plot: ", contrast_name), ylim = c(-5, 5))
        dev.off()

        # ---- Volcano Plot ----
        vol_df <- as.data.frame(res)
        vol_df$significant <- ifelse(!is.na(vol_df$padj) &
                                     vol_df$padj < 0.05 &
                                     abs(vol_df$log2FoldChange) > 1,
                                     ifelse(vol_df$log2FoldChange > 0, "Up", "Down"),
                                     "NS")
        vol_df$significant <- factor(vol_df$significant, levels = c("Up", "Down", "NS"))

        p_vol <- ggplot(vol_df, aes(x = log2FoldChange,
                                    y = -log10(pvalue),
                                    color = significant)) +
            geom_point(alpha = 0.5, size = 1) +
            scale_color_manual(values = c("Up" = "firebrick3",
                                          "Down" = "steelblue3",
                                          "NS" = "grey60")) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
            geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
            xlim(c(-10, 10)) +
            ggtitle(paste0("Volcano: ", contrast_name)) +
            xlab("log2 Fold Change") +
            ylab("-log10(p-value)") +
            theme_bw(base_size = 14) +
            theme(plot.title = element_text(hjust = 0.5),
                  legend.title = element_blank())

        ggsave(file.path(plot_dir, paste0("volcano_", contrast_name, ".pdf")),
               p_vol, width = 8, height = 6)
        ggsave(file.path(plot_dir, paste0("volcano_", contrast_name, ".png")),
               p_vol, width = 8, height = 6, dpi = 300)
    }

    saveRDS(all_results, file.path(out_dir, "all_contrast_results.rds"))

    # ---- DEG summary table ----
    summary_df <- data.frame(
        contrast = names(all_results),
        total_tested = sapply(all_results, function(r) sum(!is.na(r$padj))),
        deg_005 = sapply(all_results, function(r) sum(r$padj < 0.05, na.rm = TRUE)),
        deg_005_lfc1 = sapply(all_results, function(r)
            sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE)),
        up_005_lfc1 = sapply(all_results, function(r)
            sum(r$padj < 0.05 & r$log2FoldChange > 1, na.rm = TRUE)),
        down_005_lfc1 = sapply(all_results, function(r)
            sum(r$padj < 0.05 & r$log2FoldChange < -1, na.rm = TRUE))
    )
    write.table(summary_df, file.path(out_dir, "DEG_summary.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    cat("\nProject", project_id, "complete.\n")
    return(invisible(list(dds = dds, vsd = vsd, results = all_results)))
}

###############################################################################
# 4. Run per BioProject
###############################################################################

projects <- unique(coldata$bioproject)
cat("\nBioProjects detected:", paste(projects, collapse = ", "), "\n\n")

project_results <- list()

for (proj in projects) {
    idx <- coldata$bioproject == proj
    proj_counts <- count_mat[, idx]
    proj_meta   <- coldata[idx, , drop = FALSE]

    tryCatch({
        project_results[[proj]] <- run_deseq2_project(
            counts     = proj_counts,
            metadata   = proj_meta,
            project_id = proj,
            cond_col   = cond_col
        )
    }, error = function(e) {
        cat("ERROR in project", proj, ":", conditionMessage(e), "\n")
    })
}

###############################################################################
# 5. Cross-project normalized counts merge (for downstream scripts)
###############################################################################

# Merge all normalized counts (note: cross-project normalization is separate)
all_norm <- list()
for (proj in names(project_results)) {
    nf <- file.path(RESULTS_DIR, proj, "normalized_counts.rds")
    if (file.exists(nf)) all_norm[[proj]] <- readRDS(nf)
}

if (length(all_norm) > 0) {
    common_genes <- Reduce(intersect, lapply(all_norm, rownames))
    merged_norm <- do.call(cbind, lapply(all_norm, function(x) x[common_genes, ]))
    write.table(merged_norm, file.path(RESULTS_DIR, "all_normalized_counts.tsv"),
                sep = "\t", quote = FALSE)
    saveRDS(merged_norm, file.path(RESULTS_DIR, "all_normalized_counts.rds"))
    cat("\nMerged normalized counts:", nrow(merged_norm), "genes x",
        ncol(merged_norm), "samples\n")
}

cat("\n=== DESeq2 analysis complete ===\n")
cat("Results in:", RESULTS_DIR, "\n")
sessionInfo()

RSCRIPT_EOF

echo "=== DESeq2 analysis finished: $(date) ==="
