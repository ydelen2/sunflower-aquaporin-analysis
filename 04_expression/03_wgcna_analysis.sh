#!/bin/bash
#SBATCH --job-name=wgcna_aqp
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/wgcna_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/wgcna_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 03_wgcna_analysis.sh
# WGCNA co-expression network analysis - PRJNA869183 (4 stresses)
###############################################################################

set -euo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
EXPR_DIR="${PROJ_DIR}/04_expression"
LOG_DIR="${EXPR_DIR}/logs"

mkdir -p "${LOG_DIR}"

echo "=== WGCNA Analysis ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "Memory: $(free -h | head -2)"

module purge
module load R/4.1

Rscript --no-save --no-restore - <<'RSCRIPT_EOF'

###############################################################################
# WGCNA Co-expression Network Analysis
###############################################################################

PROJ_DIR    <- "/work/dweikat/ydelen2/aquaporin_study"
DESEQ_DIR   <- file.path(PROJ_DIR, "04_expression/results/deseq2")
RESULTS_DIR <- file.path(PROJ_DIR, "04_expression/results/wgcna")
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
# Draw cut line at reasonable height for outlier detection
cut_height <- max(sampleTree$height) * 0.85
abline(h = cut_height, col = "red", lty = 2)
dev.off()

# Remove outliers if any
clust <- cutreeStatic(sampleTree, cutHeight = cut_height, minSize = 3)
keep_samples <- clust == 1
if (sum(!keep_samples) > 0) {
    cat("Removing", sum(!keep_samples), "outlier samples:\n")
    cat(paste("  ", rownames(datExpr)[!keep_samples]), sep = "\n")
    datExpr <- datExpr[keep_samples, ]
    sample_meta <- sample_meta[keep_samples, ]
}
cat("Samples after outlier removal:", nrow(datExpr), "\n\n")

# Check for genes with zero variance after outlier removal
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
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
     ylab = "Scale Free Topology Model Fit (signed R²)",
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
    saveTOMs             = TRUE,
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

# Create binary trait matrix
cond_col_name <- if ("stress" %in% colnames(sample_meta)) "stress" else "condition"
conditions <- unique(sample_meta[[cond_col_name]])

traitData <- sapply(conditions, function(cond) {
    as.numeric(sample_meta[[cond_col_name]] == cond)
})
rownames(traitData) <- rownames(datExpr)
colnames(traitData) <- conditions

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
# 9. GO/KEGG enrichment of AQP-containing modules
###############################################################################

cat("\n=== Module Enrichment Analysis ===\n")

# Attempt organism database for sunflower
# Since H. annuus is not in OrgDb, use Arabidopsis orthologs or eggNOG annotations
# Check if custom annotation exists
annot_file <- file.path(PROJ_DIR, "02_gene_family/results/eggnog_annotations.tsv")
ath_ortho_file <- file.path(PROJ_DIR, "02_gene_family/results/arabidopsis_orthologs.tsv")

use_custom_annot <- FALSE
gene2go <- NULL
gene2kegg <- NULL

if (file.exists(annot_file)) {
    cat("Loading eggNOG annotations for enrichment\n")
    eggnog <- read.delim(annot_file, stringsAsFactors = FALSE, comment.char = "#")

    # Parse GO terms from eggNOG (column names vary)
    go_col <- grep("GOs|GO_terms", colnames(eggnog), value = TRUE, ignore.case = TRUE)
    id_eggnog <- grep("query|gene", colnames(eggnog), value = TRUE, ignore.case = TRUE)[1]

    if (length(go_col) > 0 && !is.na(id_eggnog)) {
        eggnog_clean <- eggnog[eggnog[[go_col[1]]] != "" & eggnog[[go_col[1]]] != "-", ]
        gene2go <- eggnog_clean %>%
            select(gene = !!sym(id_eggnog), GO = !!sym(go_col[1])) %>%
            separate_rows(GO, sep = ",") %>%
            filter(grepl("^GO:", GO))
        use_custom_annot <- TRUE
    }

    kegg_col <- grep("KEGG|kegg_ko", colnames(eggnog), value = TRUE, ignore.case = TRUE)
    if (length(kegg_col) > 0 && !is.na(id_eggnog)) {
        eggnog_kegg <- eggnog[eggnog[[kegg_col[1]]] != "" & eggnog[[kegg_col[1]]] != "-", ]
        gene2kegg <- eggnog_kegg %>%
            select(gene = !!sym(id_eggnog), KEGG = !!sym(kegg_col[1])) %>%
            separate_rows(KEGG, sep = ",")
    }

} else if (file.exists(ath_ortho_file)) {
    cat("Loading Arabidopsis ortholog mapping for enrichment\n")
    ath_map <- read.delim(ath_ortho_file, stringsAsFactors = FALSE)
    # Will use enrichGO with org.At.tair.db
    tryCatch({
        library(org.At.tair.db)
        cat("Using org.At.tair.db for Arabidopsis-based enrichment\n")
    }, error = function(e) {
        cat("org.At.tair.db not available. Skipping GO enrichment.\n")
    })
}

# Run enrichment per AQP-containing module
enrich_results <- list()

for (mod_color in aqp_module_colors) {
    mod_genes <- names(moduleColors)[moduleColors == mod_color]
    cat("\nEnrichment for module:", mod_color, "(", length(mod_genes), "genes)\n")

    if (use_custom_annot && !is.null(gene2go)) {
        # Custom enricher with TERM2GENE
        tryCatch({
            ego <- enricher(
                gene     = mod_genes,
                universe = colnames(datExpr),
                TERM2GENE = gene2go[, c("GO", "gene")],
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1
            )
            if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
                enrich_results[[paste0(mod_color, "_GO")]] <- ego
                write.table(as.data.frame(ego),
                            file.path(RESULTS_DIR, paste0("enrichGO_module_", mod_color, ".tsv")),
                            sep = "\t", row.names = FALSE, quote = FALSE)

                pdf(file.path(plot_dir, paste0("enrichGO_dotplot_", mod_color, ".pdf")),
                    width = 10, height = 8)
                print(dotplot(ego, showCategory = 20,
                              title = paste0("GO Enrichment: Module ", mod_color)))
                dev.off()
            }
        }, error = function(e) {
            cat("  GO enrichment error:", conditionMessage(e), "\n")
        })

        # KEGG enrichment
        if (!is.null(gene2kegg)) {
            tryCatch({
                ekegg <- enricher(
                    gene     = mod_genes,
                    universe = colnames(datExpr),
                    TERM2GENE = gene2kegg[, c("KEGG", "gene")],
                    pvalueCutoff = 0.05
                )
                if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
                    enrich_results[[paste0(mod_color, "_KEGG")]] <- ekegg
                    write.table(as.data.frame(ekegg),
                                file.path(RESULTS_DIR, paste0("enrichKEGG_module_", mod_color, ".tsv")),
                                sep = "\t", row.names = FALSE, quote = FALSE)
                }
            }, error = function(e) {
                cat("  KEGG enrichment error:", conditionMessage(e), "\n")
            })
        }

    } else if (exists("ath_map")) {
        # Arabidopsis ortholog-based enrichment
        mod_ath <- ath_map[ath_map[[1]] %in% mod_genes, ]
        ath_col <- grep("ath|arabidopsis|tair", colnames(ath_map), ignore.case = TRUE, value = TRUE)
        if (length(ath_col) > 0 && exists("org.At.tair.db")) {
            ath_ids <- unique(mod_ath[[ath_col[1]]])
            bg_ath  <- unique(ath_map[[ath_col[1]]])

            tryCatch({
                ego <- enrichGO(
                    gene     = ath_ids,
                    universe = bg_ath,
                    OrgDb    = org.At.tair.db,
                    keyType  = "TAIR",
                    ont      = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1
                )
                if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
                    enrich_results[[paste0(mod_color, "_GO_BP")]] <- ego
                    write.table(as.data.frame(ego),
                                file.path(RESULTS_DIR, paste0("enrichGO_module_", mod_color, ".tsv")),
                                sep = "\t", row.names = FALSE, quote = FALSE)
                }
            }, error = function(e) {
                cat("  Arabidopsis GO enrichment error:", conditionMessage(e), "\n")
            })
        }
    } else {
        cat("  No annotation source available. Skipping enrichment.\n")
    }
}

if (length(enrich_results) > 0) {
    saveRDS(enrich_results, file.path(RESULTS_DIR, "module_enrichment_results.rds"))
}

###############################################################################
# 10. Export network edges for Cytoscape
###############################################################################

cat("\n=== Exporting Network for Cytoscape ===\n")

# Export edges for AQP-containing modules only (full network is too large)
for (mod_color in aqp_module_colors) {
    mod_genes <- names(moduleColors)[moduleColors == mod_color]

    if (length(mod_genes) > 2000) {
        cat("Module", mod_color, "too large (", length(mod_genes),
            "genes). Exporting top 500 by kME.\n")
        kme_col_name <- paste0("kME_", which(labels2colors(0:max(net$colors)) == mod_color) - 1)
        if (kme_col_name %in% colnames(kME_all)) {
            top_genes <- rownames(kME_all[mod_genes, ])[
                order(-abs(kME_all[mod_genes, kme_col_name]))[1:500]
            ]
            mod_genes <- top_genes
        } else {
            mod_genes <- mod_genes[1:500]
        }
    }

    # Get TOM for this module
    # TOM files saved by blockwiseModules
    tom_file <- list.files(RESULTS_DIR, pattern = "TOM.*\\.RData$", full.names = TRUE)

    if (length(tom_file) > 0) {
        cat("Loading TOM from:", tom_file[1], "\n")
        load(tom_file[1])  # loads TOM object

        if (exists("TOM")) {
            # Gene indices in TOM
            block_genes <- net$blockGenes[[1]]
            gene_names_block <- colnames(datExpr)[block_genes]

            mod_idx <- which(gene_names_block %in% mod_genes)

            if (length(mod_idx) > 1) {
                mod_tom <- TOM[mod_idx, mod_idx]
                rownames(mod_tom) <- gene_names_block[mod_idx]
                colnames(mod_tom) <- gene_names_block[mod_idx]

                # Export edges above threshold
                threshold <- quantile(mod_tom[upper.tri(mod_tom)], 0.95)

                edges <- which(mod_tom > threshold & upper.tri(mod_tom), arr.ind = TRUE)
                if (nrow(edges) > 0) {
                    edge_df <- data.frame(
                        source = rownames(mod_tom)[edges[, 1]],
                        target = colnames(mod_tom)[edges[, 2]],
                        weight = mod_tom[edges],
                        stringsAsFactors = FALSE
                    )

                    # Mark aquaporin genes
                    edge_df$source_aqp <- edge_df$source %in% aqp_ids
                    edge_df$target_aqp <- edge_df$target %in% aqp_ids

                    write.table(edge_df,
                                file.path(RESULTS_DIR,
                                          paste0("cytoscape_edges_", mod_color, ".tsv")),
                                sep = "\t", row.names = FALSE, quote = FALSE)
                    cat("Module", mod_color, ":", nrow(edge_df), "edges exported\n")
                }
            }
        }
    } else {
        cat("TOM file not found. Computing adjacency for Cytoscape export.\n")

        # Compute adjacency directly (smaller subset)
        mod_expr <- datExpr[, mod_genes]
        adj <- adjacency(mod_expr, power = selected_power, type = "signed hybrid")

        # Top edges
        threshold <- quantile(adj[upper.tri(adj)], 0.95)
        edges <- which(adj > threshold & upper.tri(adj), arr.ind = TRUE)

        if (nrow(edges) > 0) {
            edge_df <- data.frame(
                source = rownames(adj)[edges[, 1]],
                target = colnames(adj)[edges[, 2]],
                weight = adj[edges],
                source_aqp = rownames(adj)[edges[, 1]] %in% aqp_ids,
                target_aqp = colnames(adj)[edges[, 2]] %in% aqp_ids,
                stringsAsFactors = FALSE
            )
            write.table(edge_df,
                        file.path(RESULTS_DIR,
                                  paste0("cytoscape_edges_", mod_color, ".tsv")),
                        sep = "\t", row.names = FALSE, quote = FALSE)
            cat("Module", mod_color, ":", nrow(edge_df), "edges exported\n")
        }
    }
}

# Node attribute table for Cytoscape
node_attrs <- data.frame(
    gene_id = colnames(datExpr),
    module_color = moduleColors[colnames(datExpr)],
    is_aquaporin = colnames(datExpr) %in% aqp_ids,
    stringsAsFactors = FALSE
)
node_attrs <- cbind(node_attrs, kME_all[colnames(datExpr), ])

write.table(node_attrs, file.path(RESULTS_DIR, "cytoscape_node_attributes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

###############################################################################
# 11. Summary
###############################################################################

cat("\n=== WGCNA Analysis Summary ===\n")
cat("Total genes in network:", ncol(datExpr), "\n")
cat("Total samples:", nrow(datExpr), "\n")
cat("Soft threshold power:", selected_power, "\n")
cat("Number of modules:", max(net$colors), "\n")
cat("Aquaporin genes in network:", length(aqp_in_network), "\n")
cat("AQP-containing modules:", paste(aqp_module_colors, collapse = ", "), "\n")
cat("Results in:", RESULTS_DIR, "\n\n")

sessionInfo()

RSCRIPT_EOF

echo "=== WGCNA analysis finished: $(date) ==="
