#!/bin/bash
#SBATCH --job-name=wgcna_v3
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/wgcna_v3_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/wgcna_v3_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 03_wgcna_analysis_v3.sh
# WGCNA co-expression network analysis of PRJNA869183.
#
#   1. All 57 libraries are used. The sample tree is only plotted; samples
#      are removed only if goodSamplesGenes flags them.
#   2. Traits are coded from the treatment each library received, with the
#      pre-treatment and untreated libraries (cold_0h, heat_0h, drought_normal,
#      control) coded 0 for every trait.
#        cold        cold_8h, cold_16h, cold_32h
#        heat        heat_4h, heat_8h, heat_16h, heat_32h
#        drought     drought_moderate, drought_severe
#        rewatering  drought_water_deficit (water regain after drought, per NCBI BioSample)
#        salt        salt_NaCl
#        SA          salt_SA1 to salt_SA4 (SA applied after salt stress, per NCBI BioSample)
#   3. Results go to 04_expression/results/wgcna_v3/.
# Network settings: variance-stabilized counts, genes above the
# first quartile of MAD (aquaporin genes always kept), signed hybrid network,
# soft-thresholding power from pickSoftThreshold (scale-free fit R^2 > 0.85),
# signed TOM, minModuleSize 30, mergeCutHeight 0.25.
###############################################################################

set -eo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
EXPR_DIR="${PROJ_DIR}/04_expression"
LOG_DIR="${EXPR_DIR}/logs"

mkdir -p "${LOG_DIR}"

echo "=== WGCNA Analysis ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "Memory: $(free -h | head -2)"

module purge
module load miniforge/24.5
eval "$(conda shell.bash hook)"
conda activate /work/dweikat/ydelen2/aquaporin_study/conda_envs/aqp_env

Rscript --no-save --no-restore - <<'RSCRIPT_EOF'

###############################################################################
# WGCNA Co-expression Network Analysis
###############################################################################

PROJ_DIR    <- "/work/dweikat/ydelen2/aquaporin_study"
DESEQ_DIR   <- file.path(PROJ_DIR, "04_expression/results/deseq2")
RESULTS_DIR <- file.path(PROJ_DIR, "04_expression/results/wgcna_v3")
AQP_LIST    <- file.path(PROJ_DIR, "02_gene_family/results/verified_aquaporins.txt")
R_LIBS      <- file.path(PROJ_DIR, "R_libs")

.libPaths(c(R_LIBS, .libPaths()))

suppressPackageStartupMessages({
    library(WGCNA)
    library(DESeq2)
    library(ggplot2)
    library(tidyverse)
    library(clusterProfiler)
    library(pheatmap)
})

# WGCNA settings
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 16)
# Override cor with WGCNA version to avoid conflicts
cor <- WGCNA::cor
enableWGCNAThreads(nThreads = 16)

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(RESULTS_DIR, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("WGCNA version:", as.character(packageVersion("WGCNA")), "\n")
cat("Threads:", 16, "\n\n")

###############################################################################
# 1. Load data - PRJNA869183 project
###############################################################################

# Prefer PRJNA869183-specific DESeq2 output
target_project <- "PRJNA869183"

vsd_file <- file.path(DESEQ_DIR, target_project, "vsd.rds")
dds_file <- file.path(DESEQ_DIR, target_project, "dds.rds")
norm_file <- file.path(DESEQ_DIR, target_project, "normalized_counts.rds")
meta_file <- file.path(DESEQ_DIR, "sample_metadata.tsv")

# Load VST data
if (file.exists(vsd_file)) {
    cat("Loading VST data for", target_project, "\n")
    vsd <- readRDS(vsd_file)
    expr_mat <- assay(vsd)
} else if (file.exists(norm_file)) {
    cat("VST not found. Using normalized counts with log2 transform.\n")
    norm_counts <- readRDS(norm_file)
    expr_mat <- log2(norm_counts + 1)
} else {
    stop("No expression data found for ", target_project)
}

cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

# Load metadata
if (file.exists(meta_file)) {
    full_meta <- read.delim(meta_file, stringsAsFactors = FALSE)
    # make sample names the row names so that metadata can be matched to count columns
    id_candidates <- intersect(c("sample", "sample_id", "run", "Run", "SRR", "sample_name"), colnames(full_meta))
    if (length(id_candidates) > 0) {
        rownames(full_meta) <- make.unique(as.character(full_meta[[id_candidates[1]]]))
    } else if (!any(colnames(expr_mat) %in% rownames(full_meta))) {
        first_match <- which(sapply(full_meta, function(x) any(as.character(x) %in% colnames(expr_mat))))[1]
        if (!is.na(first_match)) rownames(full_meta) <- make.unique(as.character(full_meta[[first_match]]))
    }
    cat("Metadata samples matching expression columns:", sum(colnames(expr_mat) %in% rownames(full_meta)), "/", ncol(expr_mat), "\n")
    # Filter to PRJNA869183
    if ("bioproject" %in% colnames(full_meta)) {
        sample_meta <- full_meta[full_meta$bioproject == target_project, ]
    } else {
        sample_meta <- full_meta[rownames(full_meta) %in% colnames(expr_mat), ]
    }
} else {
    stop("Sample metadata not found: ", meta_file)
}

# Ensure alignment
common_samples <- intersect(colnames(expr_mat), rownames(sample_meta))
if (length(common_samples) < 6) stop("Too few samples for WGCNA: ", length(common_samples))
expr_mat    <- expr_mat[, common_samples]
sample_meta <- sample_meta[common_samples, ]

cat("Using", ncol(expr_mat), "samples for WGCNA\n")
cat("Conditions:", paste(unique(sample_meta$stress), collapse = ", "), "\n\n")

# Load aquaporin gene list
aqp_genes <- read.delim(AQP_LIST, header = TRUE, stringsAsFactors = FALSE)
id_col <- intersect(c("gene_id", "Gene_ID", "GeneID", "id"), colnames(aqp_genes))
if (length(id_col) == 0) id_col <- colnames(aqp_genes)[1] else id_col <- id_col[1]
aqp_ids <- aqp_genes[[id_col]]

###############################################################################
# 2. Gene filtering for WGCNA
###############################################################################

# Keep genes with sufficient variation (top 75% by MAD)
gene_mad <- apply(expr_mat, 1, mad)
mad_threshold <- quantile(gene_mad, 0.25)
keep <- gene_mad > mad_threshold & gene_mad > 0

# Always keep aquaporin genes regardless of MAD
aqp_in_mat <- aqp_ids[aqp_ids %in% rownames(expr_mat)]
keep[aqp_in_mat] <- TRUE

expr_filtered <- expr_mat[keep, ]
cat("Genes after MAD filtering:", nrow(expr_filtered), "\n")
cat("Aquaporin genes retained:", sum(aqp_in_mat %in% rownames(expr_filtered)),
    "/", length(aqp_in_mat), "\n\n")

# WGCNA expects samples as rows, genes as columns
datExpr <- t(expr_filtered)

###############################################################################
# 3. Sample clustering and outlier detection
###############################################################################

sampleTree <- hclust(dist(datExpr), method = "average")

pdf(file.path(plot_dir, "sample_clustering.pdf"), width = 12, height = 8)
par(cex = 0.8)
plot(sampleTree, main = "Sample Clustering (Average Linkage)",
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1, cex.main = 1.5)
dev.off()

# Samples are not removed on the basis of the cluster tree. The two top-level
# branches are reported so that the grouping can be inspected.
top2 <- cutree(sampleTree, k = 2)
cat("Top-level branches of the sample tree (conditions per branch):\n")
print(table(sample_meta$condition, top2))
cat("Samples kept:", nrow(datExpr), "\n\n")

# Check for genes and samples with too many missing values or zero variance
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    sample_meta <- sample_meta[gsg$goodSamples, ]
    cat("Removed bad samples/genes. New dims:", nrow(datExpr), "x", ncol(datExpr), "\n")
}

###############################################################################
# 4. Soft threshold selection
###############################################################################

powers <- c(1:20, seq(22, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                          networkType = "signed hybrid",
                          verbose = 5)

# Plot scale-free topology fit
pdf(file.path(plot_dir, "soft_threshold.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

# Scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     main = "Scale Independence",
     type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red", lty = 2)

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean Connectivity",
     type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = 0.9, col = "red")
dev.off()

# Select power: use estimated or default to 12
selected_power <- sft$powerEstimate
if (is.na(selected_power) || selected_power < 1) {
    selected_power <- 12
    cat("WARNING: No power estimate from pickSoftThreshold. Using default:", selected_power, "\n")
} else {
    cat("Selected soft threshold power:", selected_power, "\n")
}

saveRDS(sft, file.path(RESULTS_DIR, "soft_threshold.rds"))
write.table(sft$fitIndices, file.path(RESULTS_DIR, "soft_threshold_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

###############################################################################
# 5. Network construction (blockwiseModules)
###############################################################################

cat("\nConstructing network with blockwiseModules...\n")
cat("This may take a while. Power:", selected_power, "\n")

net <- blockwiseModules(
    datExpr,
    power                = selected_power,
    networkType          = "signed hybrid",
    TOMType              = "signed",
    minModuleSize        = 30,
    reassignThreshold    = 0,
    mergeCutHeight       = 0.25,
    numericLabels        = TRUE,
    pamRespectsDendro     = FALSE,
    saveTOMs             = FALSE,
    saveTOMFileBase       = file.path(RESULTS_DIR, "TOM"),
    maxBlockSize         = 20000,
    verbose              = 5,
    nThreads             = 16
)

cat("Modules detected:", max(net$colors), "\n")
cat("Module sizes:\n")
print(table(net$colors))

saveRDS(net, file.path(RESULTS_DIR, "wgcna_network.rds"))

# Convert numeric labels to colors
moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)

# Module dendrogram
pdf(file.path(plot_dir, "module_dendrogram.pdf"), width = 16, height = 8)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()

###############################################################################
# 6. Module-trait correlation
###############################################################################

# Binary trait matrix coded from the treatment of each library
trait_def <- list(
    cold       = c("cold_8h", "cold_16h", "cold_32h"),
    heat       = c("heat_4h", "heat_8h", "heat_16h", "heat_32h"),
    drought    = c("drought_moderate", "drought_severe"),
    rewatering = c("drought_water_deficit"),
    salt       = c("salt_NaCl"),
    SA         = c("salt_SA1", "salt_SA2", "salt_SA3", "salt_SA4")
)
untreated <- c("cold_0h", "heat_0h", "drought_normal", "control")
stopifnot(all(sample_meta$condition %in% c(unlist(trait_def), untreated)))
traitData <- sapply(trait_def, function(conds) as.numeric(sample_meta$condition %in% conds))
rownames(traitData) <- rownames(datExpr)
conditions <- colnames(traitData)
cat("\nLibraries per trait:\n"); print(colSums(traitData))
trait_of_sample <- apply(traitData, 1, function(x) if (any(x == 1)) names(trait_def)[x == 1] else "untreated")
cat("Trait coding by condition:\n"); print(table(sample_meta$condition, trait_of_sample))

cat("\nTrait matrix:\n")
print(head(traitData))

# Module eigengenes
MEs <- net$MEs
# Reorder to match moduleColors
MEs <- orderMEs(MEs)

# Correlate eigengenes with traits
moduleTraitCor  <- cor(MEs, traitData, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Heatmap of module-trait correlations
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPval, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

pdf(file.path(plot_dir, "module_trait_correlation.pdf"), width = 10, height = 12)
par(mar = c(6, 10, 3, 2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traitData),
               yLabels = colnames(MEs),
               ySymbols = gsub("ME", "", colnames(MEs)),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()

saveRDS(list(cor = moduleTraitCor, pval = moduleTraitPval),
        file.path(RESULTS_DIR, "module_trait_correlation.rds"))

# Save module-trait table
mt_df <- as.data.frame(moduleTraitCor) %>%
    rownames_to_column("module") %>%
    pivot_longer(-module, names_to = "trait", values_to = "correlation") %>%
    left_join(
        as.data.frame(moduleTraitPval) %>%
            rownames_to_column("module") %>%
            pivot_longer(-module, names_to = "trait", values_to = "pvalue"),
        by = c("module", "trait")
    )
write.table(mt_df, file.path(RESULTS_DIR, "module_trait_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

###############################################################################
# 7. Aquaporin module membership
###############################################################################

cat("\n=== Aquaporin Gene Module Analysis ===\n")

# Which modules contain aquaporin genes?
aqp_in_network <- aqp_ids[aqp_ids %in% colnames(datExpr)]
cat("Aquaporin genes in network:", length(aqp_in_network), "\n")

aqp_modules <- data.frame(
    gene_id = aqp_in_network,
    module_number = net$colors[aqp_in_network],
    module_color  = moduleColors[aqp_in_network],
    stringsAsFactors = FALSE
)

# Merge with aquaporin info
aqp_info <- aqp_genes[match(aqp_in_network, aqp_genes[[id_col]]), ]
aqp_modules <- cbind(aqp_modules, aqp_info[, setdiff(colnames(aqp_info), id_col)])

# Module membership (kME) and gene significance (GS)
# kME = correlation of gene expression with module eigengene
kME_all <- as.data.frame(cor(datExpr, MEs, use = "p"))
colnames(kME_all) <- gsub("ME", "kME_", colnames(kME_all))

# GS for each trait
GS_all <- list()
for (trait in colnames(traitData)) {
    gs <- as.data.frame(cor(datExpr, traitData[, trait], use = "p"))
    colnames(gs) <- paste0("GS_", trait)
    gs_pval <- as.data.frame(corPvalueStudent(gs[[1]], nrow(datExpr)))
    colnames(gs_pval) <- paste0("GS_pval_", trait)
    GS_all[[trait]] <- cbind(gs, gs_pval)
}
GS_combined <- do.call(cbind, GS_all)

# AQP-specific kME and GS
aqp_kME <- kME_all[aqp_in_network, , drop = FALSE]
aqp_GS  <- GS_combined[aqp_in_network, , drop = FALSE]

aqp_full <- cbind(aqp_modules, aqp_kME, aqp_GS)
write.table(aqp_full, file.path(RESULTS_DIR, "aquaporin_module_membership.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nAquaporin module distribution:\n")
print(table(aqp_modules$module_color))

# Modules containing AQP genes
aqp_module_colors <- unique(aqp_modules$module_color)
cat("Modules containing aquaporin genes:", paste(aqp_module_colors, collapse = ", "), "\n")

###############################################################################
# 8. Hub gene identification in AQP-containing modules
###############################################################################

cat("\n=== Hub Gene Identification ===\n")

hub_results <- list()

for (mod_color in aqp_module_colors) {
    mod_genes <- names(moduleColors)[moduleColors == mod_color]
    mod_me_col <- paste0("ME", which(labels2colors(0:max(net$colors)) == mod_color) - 1)

    if (!(mod_me_col %in% colnames(MEs))) {
        # Try alternative naming
        mod_me_col <- paste0("ME", mod_color)
        if (!(mod_me_col %in% colnames(MEs))) {
            cat("  Skipping module", mod_color, "- eigengene not found\n")
            next
        }
    }

    # kME for this module
    kme_col <- paste0("kME_", gsub("ME", "", mod_me_col))
    if (kme_col %in% colnames(kME_all)) {
        mod_kme <- kME_all[mod_genes, kme_col, drop = FALSE]
    } else {
        mod_kme <- data.frame(kME = cor(datExpr[, mod_genes], MEs[, mod_me_col], use = "p"))
        colnames(mod_kme) <- "kME"
    }

    colnames(mod_kme)[1] <- "kME"
    mod_kme$gene_id <- rownames(mod_kme)
    mod_kme$is_aquaporin <- mod_kme$gene_id %in% aqp_ids
    mod_kme <- mod_kme[order(-abs(mod_kme[[1]])), ]

    # Top 20 hub genes
    hubs <- head(mod_kme, 20)
    hubs$module <- mod_color
    hub_results[[mod_color]] <- hubs

    cat("\nModule:", mod_color, "(", length(mod_genes), "genes,",
        sum(mod_genes %in% aqp_ids), "AQP)\n")
    cat("  Top 5 hubs:\n")
    print(head(hubs, 5))
}

if (length(hub_results) > 0) {
    hub_df <- do.call(rbind, hub_results)
    write.table(hub_df, file.path(RESULTS_DIR, "hub_genes_aqp_modules.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
}

###############################################################################
# 9. Run summary
###############################################################################

param_lines <- c(
    paste("samples", nrow(datExpr), sep = "\t"),
    paste("genes", ncol(datExpr), sep = "\t"),
    paste("soft_threshold_power", selected_power, sep = "\t"),
    paste("network_type", "signed hybrid", sep = "\t"),
    paste("TOM_type", "signed", sep = "\t"),
    paste("min_module_size", 30, sep = "\t"),
    paste("merge_cut_height", 0.25, sep = "\t"),
    paste("modules_excluding_unassigned", max(net$colors), sep = "\t"),
    paste("genes_unassigned", sum(net$colors == 0), sep = "\t"),
    paste("aquaporin_genes_in_network", length(aqp_in_network), sep = "\t"),
    paste("aquaporin_genes_unassigned", sum(net$colors[aqp_in_network] == 0), sep = "\t"),
    paste("modules_with_aquaporins_excluding_unassigned", length(setdiff(unique(net$colors[aqp_in_network]), 0)), sep = "\t")
)
writeLines(param_lines, file.path(RESULTS_DIR, "wgcna_parameters.tsv"))
cat("\n=== Run summary ===\n"); cat(param_lines, sep = "\n")
cat("\nModule sizes:\n"); print(table(net$colors))
cat("\nEnd:", format(Sys.time()), "\n")
sessionInfo()
RSCRIPT_EOF

echo "=== WGCNA v3 finished: $(date) ==="
