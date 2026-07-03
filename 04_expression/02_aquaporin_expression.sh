#!/bin/bash
#SBATCH --job-name=aqp_expr
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/aqp_expr_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/aqp_expr_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 02_aquaporin_expression.sh
# Aquaporin-focused expression profiling with ComplexHeatmap
###############################################################################

set -euo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
EXPR_DIR="${PROJ_DIR}/04_expression"
LOG_DIR="${EXPR_DIR}/logs"

mkdir -p "${LOG_DIR}"

echo "=== Aquaporin Expression Analysis ==="
echo "Start: $(date)"
echo "Node: $(hostname)"

module purge
module load R/4.1

Rscript --no-save --no-restore - <<'RSCRIPT_EOF'

###############################################################################
# Aquaporin Expression Profiling
###############################################################################

PROJ_DIR    <- "/work/dweikat/ydelen2/aquaporin_study"
DESEQ_DIR   <- file.path(PROJ_DIR, "04_expression/results/deseq2")
RESULTS_DIR <- file.path(PROJ_DIR, "04_expression/results/aquaporin_expression")
AQP_LIST    <- file.path(PROJ_DIR, "02_gene_family/results/verified_aquaporins.txt")
R_LIBS      <- file.path(PROJ_DIR, "R_libs")

.libPaths(c(R_LIBS, .libPaths()))

suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(tidyverse)
    library(pheatmap)
    library(RColorBrewer)
})

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(RESULTS_DIR, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("R version:", R.version.string, "\n")
cat("ComplexHeatmap version:", as.character(packageVersion("ComplexHeatmap")), "\n\n")

###############################################################################
# 1. Load aquaporin gene list
###############################################################################

if (!file.exists(AQP_LIST)) stop("Aquaporin gene list not found: ", AQP_LIST)

aqp_genes <- read.delim(AQP_LIST, header = TRUE, stringsAsFactors = FALSE)
cat("Loaded", nrow(aqp_genes), "aquaporin genes\n")
cat("Columns:", paste(colnames(aqp_genes), collapse = ", "), "\n")

# Expected columns: gene_id, gene_name (or symbol), subfamily (PIP/TIP/NIP/SIP/XIP)
# Detect subfamily column
subfamily_col <- intersect(c("subfamily", "Subfamily", "family", "subgroup", "type"),
                           colnames(aqp_genes))
if (length(subfamily_col) == 0) {
    cat("WARNING: No subfamily column found. Inferring from gene name.\n")
    # Infer subfamily from gene name pattern (e.g., HaPIP1;1 -> PIP)
    aqp_genes$subfamily <- gsub(".*?(PIP|TIP|NIP|SIP|XIP).*", "\\1",
                                aqp_genes[, grep("name|symbol|gene", colnames(aqp_genes),
                                                 ignore.case = TRUE, value = TRUE)[1]])
    subfamily_col <- "subfamily"
} else {
    subfamily_col <- subfamily_col[1]
}

# Detect gene_id column
id_col <- intersect(c("gene_id", "Gene_ID", "GeneID", "id"), colnames(aqp_genes))
if (length(id_col) == 0) {
    id_col <- colnames(aqp_genes)[1]
    cat("Using first column as gene ID:", id_col, "\n")
} else {
    id_col <- id_col[1]
}

cat("Subfamily distribution:\n")
print(table(aqp_genes[[subfamily_col]]))
cat("\n")

###############################################################################
# 2. Load normalized counts and gene lengths
###############################################################################

# Try merged normalized counts first
norm_file <- file.path(DESEQ_DIR, "all_normalized_counts.rds")
if (!file.exists(norm_file)) {
    # Fall back to first project
    proj_dirs <- list.dirs(DESEQ_DIR, recursive = FALSE)
    norm_candidates <- file.path(proj_dirs, "normalized_counts.rds")
    norm_file <- norm_candidates[file.exists(norm_candidates)][1]
}
if (is.na(norm_file) || !file.exists(norm_file)) {
    stop("No normalized count matrix found in DESeq2 results")
}

norm_counts <- readRDS(norm_file)
cat("Normalized count matrix:", nrow(norm_counts), "genes x",
    ncol(norm_counts), "samples\n")

# Load sample metadata
meta_file <- file.path(DESEQ_DIR, "sample_metadata.tsv")
if (file.exists(meta_file)) {
    sample_meta <- read.delim(meta_file, stringsAsFactors = FALSE)
} else {
    cat("WARNING: No sample metadata found. Creating minimal metadata from column names.\n")
    sample_meta <- data.frame(sample = colnames(norm_counts),
                              stringsAsFactors = FALSE)
    rownames(sample_meta) <- sample_meta$sample
}

# Gene lengths for TPM
length_file <- file.path(DESEQ_DIR, "gene_lengths.rds")
gene_lengths <- if (file.exists(length_file)) readRDS(length_file) else NULL

###############################################################################
# 3. Calculate TPM and FPKM for aquaporin genes
###############################################################################

# Match aquaporin genes to count matrix
aqp_ids <- aqp_genes[[id_col]]
found   <- aqp_ids[aqp_ids %in% rownames(norm_counts)]
missing <- aqp_ids[!aqp_ids %in% rownames(norm_counts)]

cat("Aquaporin genes found in count matrix:", length(found), "/", length(aqp_ids), "\n")
if (length(missing) > 0) {
    cat("Missing genes:", paste(head(missing, 10), collapse = ", "),
        if (length(missing) > 10) paste0("... (", length(missing) - 10, " more)") else "", "\n")
}

aqp_norm <- norm_counts[found, , drop = FALSE]

# TPM calculation
if (!is.null(gene_lengths)) {
    # Load raw counts for TPM (TPM is calculated from raw, not normalized)
    raw_file <- file.path(DESEQ_DIR, "raw_counts.rds")
    if (file.exists(raw_file)) {
        raw_counts <- readRDS(raw_file)
    } else {
        cat("Using normalized counts for RPK (raw counts not available)\n")
        raw_counts <- norm_counts
    }

    # RPK = reads / (gene length in kb)
    common_genes <- intersect(rownames(raw_counts), names(gene_lengths))
    rpk <- raw_counts[common_genes, ] / (gene_lengths[common_genes] / 1000)

    # TPM = RPK / sum(RPK) * 1e6
    tpm <- apply(rpk, 2, function(x) x / sum(x) * 1e6)

    # FPKM: (reads * 1e9) / (total_reads * gene_length)
    total_reads <- colSums(raw_counts[common_genes, ])
    fpkm <- sweep(raw_counts[common_genes, ], 2, total_reads, "/") * 1e9
    fpkm <- sweep(fpkm, 1, gene_lengths[common_genes], "/")

    # Extract AQP-specific
    aqp_found_len <- intersect(found, common_genes)
    aqp_tpm  <- tpm[aqp_found_len, , drop = FALSE]
    aqp_fpkm <- fpkm[aqp_found_len, , drop = FALSE]

    write.table(aqp_tpm, file.path(RESULTS_DIR, "aquaporin_TPM.tsv"),
                sep = "\t", quote = FALSE)
    write.table(aqp_fpkm, file.path(RESULTS_DIR, "aquaporin_FPKM.tsv"),
                sep = "\t", quote = FALSE)

    # Full genome TPM/FPKM
    write.table(tpm, file.path(RESULTS_DIR, "all_genes_TPM.tsv"),
                sep = "\t", quote = FALSE)
    cat("TPM and FPKM matrices saved\n")
} else {
    cat("Gene lengths not available. Skipping TPM/FPKM. Using normalized counts.\n")
    aqp_tpm <- aqp_norm
}

# Save normalized AQP counts
write.table(aqp_norm, file.path(RESULTS_DIR, "aquaporin_normalized_counts.tsv"),
            sep = "\t", quote = FALSE)
saveRDS(aqp_norm, file.path(RESULTS_DIR, "aquaporin_normalized_counts.rds"))

###############################################################################
# 4. ComplexHeatmap: AQP expression across conditions
###############################################################################

# Prepare expression matrix for heatmap (use TPM if available, else normalized)
heatmap_mat <- if (exists("aqp_tpm")) aqp_tpm else aqp_norm
heatmap_mat <- log2(heatmap_mat + 1)  # log2 transform

# Map gene IDs to subfamily info
aqp_info <- aqp_genes[aqp_genes[[id_col]] %in% rownames(heatmap_mat), ]
rownames(aqp_info) <- aqp_info[[id_col]]
aqp_info <- aqp_info[rownames(heatmap_mat), ]

# Order by subfamily
subfamily_order <- c("PIP", "TIP", "NIP", "SIP", "XIP")
aqp_info$subfamily_f <- factor(aqp_info[[subfamily_col]],
                                levels = intersect(subfamily_order,
                                                   unique(aqp_info[[subfamily_col]])))
ord <- order(aqp_info$subfamily_f)
heatmap_mat <- heatmap_mat[ord, ]
aqp_info    <- aqp_info[ord, ]

# Row annotation: subfamily
subfamily_colors <- c(PIP = "#E41A1C", TIP = "#377EB8", NIP = "#4DAF4A",
                      SIP = "#984EA3", XIP = "#FF7F00")
present_subfamilies <- unique(as.character(aqp_info$subfamily_f))
subfamily_colors <- subfamily_colors[names(subfamily_colors) %in% present_subfamilies]

row_ha <- rowAnnotation(
    Subfamily = aqp_info$subfamily_f,
    col = list(Subfamily = subfamily_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "top"
)

# Column annotation: stress type, tissue (if available in metadata)
if ("stress" %in% colnames(sample_meta) || "condition" %in% colnames(sample_meta)) {
    cond_col_name <- if ("stress" %in% colnames(sample_meta)) "stress" else "condition"

    # Ensure metadata rows match heatmap columns
    col_meta <- sample_meta[match(colnames(heatmap_mat), rownames(sample_meta)), , drop = FALSE]

    stress_vals <- unique(col_meta[[cond_col_name]])
    n_stress <- length(stress_vals)
    stress_colors <- setNames(brewer.pal(max(3, min(n_stress, 12)), "Set3")[1:n_stress],
                              stress_vals)

    col_ha_args <- list(
        Stress = col_meta[[cond_col_name]],
        col = list(Stress = stress_colors),
        show_annotation_name = TRUE,
        annotation_name_side = "left"
    )

    if ("tissue" %in% colnames(col_meta)) {
        tissue_vals <- unique(col_meta$tissue)
        n_tissue <- length(tissue_vals)
        tissue_colors <- setNames(brewer.pal(max(3, min(n_tissue, 8)), "Set2")[1:n_tissue],
                                  tissue_vals)
        col_ha_args$Tissue <- col_meta$tissue
        col_ha_args$col$Tissue <- tissue_colors
    }

    col_ha <- do.call(HeatmapAnnotation, col_ha_args)
} else {
    col_ha <- NULL
}

# Row labels: use gene name if available, else gene_id
name_col <- intersect(c("gene_name", "name", "symbol", "Name"), colnames(aqp_info))
row_labels <- if (length(name_col) > 0) aqp_info[[name_col[1]]] else rownames(heatmap_mat)

# Scale rows (z-score)
heatmap_scaled <- t(scale(t(heatmap_mat)))
heatmap_scaled[is.nan(heatmap_scaled)] <- 0

# Color scale
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Cluster within subfamilies: split by subfamily
ht <- Heatmap(
    heatmap_scaled,
    name = "Z-score\n(log2 TPM)",
    col = col_fun,
    row_labels = row_labels,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,
    row_split = aqp_info$subfamily_f,
    row_gap = unit(2, "mm"),
    cluster_row_slices = FALSE,
    row_title_rot = 0,
    left_annotation = row_ha,
    top_annotation = col_ha,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    border = TRUE,
    heatmap_legend_param = list(direction = "vertical")
)

# PDF output
pdf(file.path(plot_dir, "aquaporin_expression_heatmap.pdf"),
    width = max(12, ncol(heatmap_mat) * 0.3 + 4),
    height = max(10, nrow(heatmap_mat) * 0.25 + 3))
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(10, 10, 10, 20), "mm"))
dev.off()

# PNG output
png(file.path(plot_dir, "aquaporin_expression_heatmap.png"),
    width = max(12, ncol(heatmap_mat) * 0.3 + 4),
    height = max(10, nrow(heatmap_mat) * 0.25 + 3),
    units = "in", res = 300)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(10, 10, 10, 20), "mm"))
dev.off()

cat("Main heatmap saved\n")

###############################################################################
# 5. Mean expression heatmap (conditions as columns, not individual samples)
###############################################################################

if (exists("cond_col_name") && cond_col_name %in% colnames(sample_meta)) {
    # Calculate mean expression per condition
    conditions <- unique(col_meta[[cond_col_name]])
    mean_expr <- sapply(conditions, function(cond) {
        cols <- rownames(col_meta)[col_meta[[cond_col_name]] == cond]
        cols <- intersect(cols, colnames(heatmap_mat))
        if (length(cols) > 1) rowMeans(heatmap_mat[, cols]) else heatmap_mat[, cols]
    })
    colnames(mean_expr) <- conditions

    # Scale
    mean_scaled <- t(scale(t(mean_expr)))
    mean_scaled[is.nan(mean_scaled)] <- 0

    ht_mean <- Heatmap(
        mean_scaled,
        name = "Z-score",
        col = col_fun,
        row_labels = row_labels,
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 11),
        column_names_rot = 45,
        row_split = aqp_info$subfamily_f,
        row_gap = unit(2, "mm"),
        cluster_row_slices = FALSE,
        row_title_rot = 0,
        left_annotation = row_ha,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        border = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", mean_expr[i, j]),
                      x, y, gp = gpar(fontsize = 6, col = "black"))
        }
    )

    pdf(file.path(plot_dir, "aquaporin_mean_expression_heatmap.pdf"),
        width = max(8, ncol(mean_expr) * 1.2 + 4),
        height = max(10, nrow(mean_expr) * 0.25 + 3))
    draw(ht_mean, heatmap_legend_side = "right",
         padding = unit(c(10, 10, 10, 20), "mm"))
    dev.off()

    png(file.path(plot_dir, "aquaporin_mean_expression_heatmap.png"),
        width = max(8, ncol(mean_expr) * 1.2 + 4),
        height = max(10, nrow(mean_expr) * 0.25 + 3),
        units = "in", res = 300)
    draw(ht_mean, heatmap_legend_side = "right",
         padding = unit(c(10, 10, 10, 20), "mm"))
    dev.off()

    cat("Mean expression heatmap saved\n")
}

###############################################################################
# 6. Per-stress expression barplots
###############################################################################

if (exists("cond_col_name")) {
    # Long format for ggplot
    expr_long <- as.data.frame(heatmap_mat) %>%
        rownames_to_column("gene_id") %>%
        pivot_longer(-gene_id, names_to = "sample", values_to = "log2_expr") %>%
        left_join(aqp_info %>% rownames_to_column("gene_id") %>%
                  select(gene_id, !!sym(subfamily_col)),
                  by = "gene_id") %>%
        left_join(col_meta %>% rownames_to_column("sample") %>%
                  select(sample, !!sym(cond_col_name)),
                  by = "sample")

    if (length(name_col) > 0) {
        name_map <- setNames(aqp_info[[name_col[1]]], rownames(aqp_info))
        expr_long$gene_label <- name_map[expr_long$gene_id]
    } else {
        expr_long$gene_label <- expr_long$gene_id
    }

    # Per-stress barplots
    stresses <- unique(expr_long[[cond_col_name]])

    for (stress in stresses) {
        sub <- expr_long[expr_long[[cond_col_name]] == stress, ]
        if (nrow(sub) == 0) next

        mean_sub <- sub %>%
            group_by(gene_id, gene_label, !!sym(subfamily_col)) %>%
            summarise(mean_expr = mean(log2_expr, na.rm = TRUE),
                      se = sd(log2_expr, na.rm = TRUE) / sqrt(n()),
                      .groups = "drop") %>%
            arrange(!!sym(subfamily_col), desc(mean_expr))

        mean_sub$gene_label <- factor(mean_sub$gene_label,
                                       levels = mean_sub$gene_label)

        p <- ggplot(mean_sub, aes(x = gene_label, y = mean_expr,
                                   fill = !!sym(subfamily_col))) +
            geom_col(width = 0.7) +
            geom_errorbar(aes(ymin = mean_expr - se, ymax = mean_expr + se),
                          width = 0.3) +
            scale_fill_manual(values = subfamily_colors) +
            labs(title = paste0("Aquaporin Expression: ", stress),
                 x = NULL, y = "log2(TPM + 1)", fill = "Subfamily") +
            theme_bw(base_size = 12) +
            theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
                  plot.title = element_text(hjust = 0.5))

        ggsave(file.path(plot_dir, paste0("barplot_", stress, ".pdf")),
               p, width = max(8, nrow(mean_sub) * 0.3), height = 6)
        ggsave(file.path(plot_dir, paste0("barplot_", stress, ".png")),
               p, width = max(8, nrow(mean_sub) * 0.3), height = 6, dpi = 300)
    }

    cat("Per-stress barplots saved\n")

    # Subfamily-level summary across stresses
    subfamily_summary <- expr_long %>%
        group_by(!!sym(subfamily_col), !!sym(cond_col_name)) %>%
        summarise(mean_expr = mean(log2_expr, na.rm = TRUE),
                  se = sd(log2_expr, na.rm = TRUE) / sqrt(n()),
                  .groups = "drop")

    p_sub <- ggplot(subfamily_summary,
                    aes(x = !!sym(cond_col_name), y = mean_expr,
                        fill = !!sym(subfamily_col))) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        geom_errorbar(aes(ymin = mean_expr - se, ymax = mean_expr + se),
                      position = position_dodge(width = 0.8), width = 0.3) +
        scale_fill_manual(values = subfamily_colors) +
        labs(title = "Aquaporin Subfamily Expression Across Stresses",
             x = "Condition", y = "Mean log2(TPM + 1)", fill = "Subfamily") +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5))

    ggsave(file.path(plot_dir, "subfamily_expression_summary.pdf"), p_sub,
           width = 10, height = 7)
    ggsave(file.path(plot_dir, "subfamily_expression_summary.png"), p_sub,
           width = 10, height = 7, dpi = 300)
}

cat("\n=== Aquaporin expression analysis complete ===\n")
cat("Results in:", RESULTS_DIR, "\n")
sessionInfo()

RSCRIPT_EOF

echo "=== Aquaporin expression analysis finished: $(date) ==="
