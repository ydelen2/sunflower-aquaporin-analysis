#!/bin/bash
#SBATCH --job-name=enrichment_aqp
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/enrichment_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/enrichment_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 04_go_kegg_enrichment.sh
# GO/KEGG enrichment analysis for stress-specific DEGs
###############################################################################

set -euo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
EXPR_DIR="${PROJ_DIR}/04_expression"
LOG_DIR="${EXPR_DIR}/logs"

mkdir -p "${LOG_DIR}"

echo "=== GO/KEGG Enrichment Analysis ==="
echo "Start: $(date)"
echo "Node: $(hostname)"

module purge
module load R/4.1

Rscript --no-save --no-restore - <<'RSCRIPT_EOF'

###############################################################################
# GO/KEGG Enrichment Analysis for Stress DEGs
###############################################################################

PROJ_DIR    <- "/work/dweikat/ydelen2/aquaporin_study"
DESEQ_DIR   <- file.path(PROJ_DIR, "04_expression/results/deseq2")
RESULTS_DIR <- file.path(PROJ_DIR, "04_expression/results/enrichment")
AQP_LIST    <- file.path(PROJ_DIR, "02_gene_family/results/verified_aquaporins.txt")
R_LIBS      <- file.path(PROJ_DIR, "R_libs")

.libPaths(c(R_LIBS, .libPaths()))

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
    library(tidyverse)
})

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(RESULTS_DIR, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("clusterProfiler version:", as.character(packageVersion("clusterProfiler")), "\n\n")

###############################################################################
# 1. Load annotation source
###############################################################################

# Priority: eggNOG annotations > Arabidopsis orthologs > KEGG API
annot_file     <- file.path(PROJ_DIR, "02_gene_family/results/eggnog_annotations.tsv")
ath_ortho_file <- file.path(PROJ_DIR, "02_gene_family/results/arabidopsis_orthologs.tsv")

# Annotation containers
TERM2GENE_GO  <- NULL
TERM2NAME_GO  <- NULL
TERM2GENE_KEGG <- NULL
TERM2NAME_KEGG <- NULL
use_orgdb      <- FALSE
ath_map        <- NULL

if (file.exists(annot_file)) {
    cat("=== Loading eggNOG-mapper annotations ===\n")
    eggnog <- read.delim(annot_file, stringsAsFactors = FALSE, comment.char = "#")
    cat("eggNOG entries:", nrow(eggnog), "\n")
    cat("Columns:", paste(colnames(eggnog), collapse = ", "), "\n")

    # Identify column names (eggNOG v2 output has varying headers)
    id_col_egg <- grep("^query$|^#query$|^gene$", colnames(eggnog),
                       ignore.case = TRUE, value = TRUE)
    if (length(id_col_egg) == 0) id_col_egg <- colnames(eggnog)[1]
    else id_col_egg <- id_col_egg[1]

    go_col <- grep("GOs$|GO_terms$", colnames(eggnog), ignore.case = TRUE, value = TRUE)
    kegg_col <- grep("KEGG_ko$|KEGG_Pathway$", colnames(eggnog), ignore.case = TRUE, value = TRUE)
    kegg_path_col <- grep("KEGG_Pathway$", colnames(eggnog), ignore.case = TRUE, value = TRUE)

    # Build GO TERM2GENE
    if (length(go_col) > 0) {
        go_data <- eggnog[, c(id_col_egg, go_col[1])]
        colnames(go_data) <- c("gene", "GO")
        go_data <- go_data[go_data$GO != "" & go_data$GO != "-", ]
        TERM2GENE_GO <- go_data %>%
            separate_rows(GO, sep = ",") %>%
            filter(grepl("^GO:", GO)) %>%
            mutate(GO = trimws(GO)) %>%
            select(term = GO, gene = gene)
        cat("GO annotations:", nrow(TERM2GENE_GO), "gene-term pairs\n")
        cat("Unique GO terms:", length(unique(TERM2GENE_GO$term)), "\n")
    }

    # Build KEGG TERM2GENE (from KO or pathway)
    if (length(kegg_path_col) > 0) {
        kegg_data <- eggnog[, c(id_col_egg, kegg_path_col[1])]
        colnames(kegg_data) <- c("gene", "KEGG")
        kegg_data <- kegg_data[kegg_data$KEGG != "" & kegg_data$KEGG != "-", ]
        TERM2GENE_KEGG <- kegg_data %>%
            separate_rows(KEGG, sep = ",") %>%
            mutate(KEGG = trimws(KEGG)) %>%
            filter(KEGG != "") %>%
            select(term = KEGG, gene = gene)
        cat("KEGG annotations:", nrow(TERM2GENE_KEGG), "gene-term pairs\n")
    } else if (length(kegg_col) > 0) {
        kegg_data <- eggnog[, c(id_col_egg, kegg_col[1])]
        colnames(kegg_data) <- c("gene", "KO")
        kegg_data <- kegg_data[kegg_data$KO != "" & kegg_data$KO != "-", ]
        TERM2GENE_KEGG <- kegg_data %>%
            separate_rows(KO, sep = ",") %>%
            mutate(KO = trimws(KO)) %>%
            filter(KO != "") %>%
            select(term = KO, gene = gene)
        cat("KEGG KO annotations:", nrow(TERM2GENE_KEGG), "gene-term pairs\n")
    }

} else if (file.exists(ath_ortho_file)) {
    cat("=== Loading Arabidopsis ortholog mapping ===\n")
    ath_map <- read.delim(ath_ortho_file, stringsAsFactors = FALSE)
    cat("Ortholog pairs:", nrow(ath_map), "\n")
    cat("Columns:", paste(colnames(ath_map), collapse = ", "), "\n")

    # Detect columns
    han_col <- grep("sunflower|helianthus|han|query|gene_id", colnames(ath_map),
                    ignore.case = TRUE, value = TRUE)
    ath_col <- grep("ath|arabidopsis|tair|ortholog|target", colnames(ath_map),
                    ignore.case = TRUE, value = TRUE)

    if (length(han_col) == 0) han_col <- colnames(ath_map)[1]
    else han_col <- han_col[1]
    if (length(ath_col) == 0) ath_col <- colnames(ath_map)[2]
    else ath_col <- ath_col[1]

    cat("Mapping:", han_col, "->", ath_col, "\n")

    # Try loading org.At.tair.db
    tryCatch({
        library(org.At.tair.db)
        use_orgdb <- TRUE
        cat("org.At.tair.db loaded successfully\n")
    }, error = function(e) {
        cat("org.At.tair.db not available.\n")
        cat("Attempting to use org.Hs.eg.db as fallback (not ideal).\n")
    })

} else {
    cat("WARNING: No annotation source found.\n")
    cat("Expected one of:\n")
    cat("  ", annot_file, "\n")
    cat("  ", ath_ortho_file, "\n")
    cat("Will attempt KEGG API enrichment only.\n")
}

###############################################################################
# 2. Collect DEG lists from all contrasts
###############################################################################

cat("\n=== Loading DEG lists ===\n")

# Search all project subdirectories for DEG files
project_dirs <- list.dirs(DESEQ_DIR, recursive = FALSE)
deg_files <- list.files(DESEQ_DIR, pattern = "^DEG_sig_.*\\.tsv$",
                        full.names = TRUE, recursive = TRUE)

cat("Found", length(deg_files), "significant DEG files\n")

if (length(deg_files) == 0) {
    # Try non-sig files and apply cutoffs
    deg_files <- list.files(DESEQ_DIR, pattern = "^DEG_.*_vs_.*\\.tsv$",
                            full.names = TRUE, recursive = TRUE)
    deg_files <- deg_files[!grepl("DEG_sig_", deg_files)]
    cat("Using full DEG result files:", length(deg_files), "\n")
}

# Load all DEG lists
deg_lists       <- list()  # All significant DEGs
deg_lists_up    <- list()  # Upregulated
deg_lists_down  <- list()  # Downregulated
deg_full        <- list()  # Full results (for GSEA)

universe_genes <- NULL

for (f in deg_files) {
    contrast_name <- gsub("^DEG_sig_|^DEG_|\\.tsv$", "", basename(f))
    # Extract project name from path
    proj_name <- basename(dirname(f))
    if (proj_name != "deseq2") {
        contrast_label <- paste0(proj_name, ":", contrast_name)
    } else {
        contrast_label <- contrast_name
    }

    df <- read.delim(f, stringsAsFactors = FALSE)
    if (nrow(df) == 0) {
        cat("  ", contrast_label, ": 0 DEGs, skipping\n")
        next
    }

    # Detect gene_id column
    gene_col <- intersect(c("gene_id", "Gene_ID", "GeneID"), colnames(df))
    if (length(gene_col) == 0) {
        if ("Row.names" %in% colnames(df)) gene_col <- "Row.names"
        else gene_col <- colnames(df)[1]
    } else {
        gene_col <- gene_col[1]
    }

    # If full results (not pre-filtered), apply cutoffs
    if (grepl("^DEG_[^s]", basename(f)) && "padj" %in% colnames(df)) {
        if (is.null(universe_genes)) universe_genes <- df[[gene_col]]
        sig <- df[!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1, ]
    } else {
        sig <- df
    }

    if (nrow(sig) == 0) {
        cat("  ", contrast_label, ": 0 significant DEGs after filtering\n")
        next
    }

    deg_lists[[contrast_label]] <- sig[[gene_col]]

    if ("log2FoldChange" %in% colnames(sig)) {
        deg_lists_up[[contrast_label]]   <- sig[[gene_col]][sig$log2FoldChange > 0]
        deg_lists_down[[contrast_label]] <- sig[[gene_col]][sig$log2FoldChange < 0]
    }

    # For GSEA: ranked gene list
    if ("log2FoldChange" %in% colnames(df) && "padj" %in% colnames(df)) {
        ranked <- df[!is.na(df$padj), ]
        gene_rank <- setNames(-log10(ranked$padj) * sign(ranked$log2FoldChange),
                              ranked[[gene_col]])
        gene_rank <- sort(gene_rank, decreasing = TRUE)
        deg_full[[contrast_label]] <- gene_rank
    }

    cat("  ", contrast_label, ":", length(deg_lists[[contrast_label]]),
        "DEGs (", length(deg_lists_up[[contrast_label]]), "up,",
        length(deg_lists_down[[contrast_label]]), "down)\n")
}

if (length(deg_lists) == 0) {
    stop("No DEG lists loaded. Check DESeq2 output in ", DESEQ_DIR)
}

# Universe: all tested genes
if (is.null(universe_genes)) {
    # Load from raw counts
    raw_file <- file.path(DESEQ_DIR, "raw_counts.rds")
    if (file.exists(raw_file)) {
        universe_genes <- rownames(readRDS(raw_file))
    }
}
cat("\nUniverse size:", length(universe_genes), "genes\n\n")

###############################################################################
# 3. GO Enrichment (ORA) per stress contrast
###############################################################################

cat("=== GO Over-Representation Analysis ===\n")

go_results     <- list()
go_results_up  <- list()
go_results_down <- list()

run_go_enrichment <- function(gene_list, label, universe = NULL) {
    if (length(gene_list) < 5) {
        cat("  ", label, ": too few genes (", length(gene_list), "), skipping\n")
        return(NULL)
    }

    result <- NULL

    if (!is.null(TERM2GENE_GO)) {
        # eggNOG-based enrichment
        tryCatch({
            result <- enricher(
                gene         = gene_list,
                universe     = universe,
                TERM2GENE    = TERM2GENE_GO,
                TERM2NAME    = TERM2NAME_GO,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1,
                pAdjustMethod = "BH",
                minGSSize    = 10,
                maxGSSize    = 500
            )
        }, error = function(e) {
            cat("  enricher error for", label, ":", conditionMessage(e), "\n")
        })

    } else if (use_orgdb && !is.null(ath_map)) {
        # Arabidopsis ortholog-based enrichment
        han_col_name <- grep("sunflower|helianthus|han|query|gene_id",
                             colnames(ath_map), ignore.case = TRUE, value = TRUE)[1]
        ath_col_name <- grep("ath|arabidopsis|tair|ortholog|target",
                             colnames(ath_map), ignore.case = TRUE, value = TRUE)[1]

        ath_genes <- ath_map[[ath_col_name]][ath_map[[han_col_name]] %in% gene_list]
        ath_genes <- unique(ath_genes[!is.na(ath_genes) & ath_genes != ""])

        ath_universe <- NULL
        if (!is.null(universe)) {
            ath_universe <- ath_map[[ath_col_name]][ath_map[[han_col_name]] %in% universe]
            ath_universe <- unique(ath_universe[!is.na(ath_universe)])
        }

        if (length(ath_genes) < 3) {
            cat("  ", label, ": too few Arabidopsis orthologs mapped\n")
            return(NULL)
        }

        # Clean TAIR IDs (remove transcript suffix)
        ath_genes <- gsub("\\.\\d+$", "", ath_genes)
        if (!is.null(ath_universe)) ath_universe <- gsub("\\.\\d+$", "", ath_universe)

        for (ont in c("BP", "MF", "CC")) {
            tryCatch({
                ego <- enrichGO(
                    gene          = ath_genes,
                    universe      = ath_universe,
                    OrgDb         = org.At.tair.db,
                    keyType       = "TAIR",
                    ont           = ont,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1
                )
                if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
                    if (is.null(result)) {
                        result <- ego
                    }
                    # Save per-ontology
                    write.table(as.data.frame(ego),
                                file.path(RESULTS_DIR,
                                          paste0("GO_", ont, "_", gsub(":", "_", label), ".tsv")),
                                sep = "\t", row.names = FALSE, quote = FALSE)
                }
            }, error = function(e) {
                cat("  enrichGO", ont, "error:", conditionMessage(e), "\n")
            })
        }
    }

    return(result)
}

# Run for each contrast
for (contrast in names(deg_lists)) {
    cat("\nContrast:", contrast, "\n")

    go_results[[contrast]] <- run_go_enrichment(
        deg_lists[[contrast]], paste0(contrast, "_all"), universe_genes)

    if (contrast %in% names(deg_lists_up)) {
        go_results_up[[contrast]] <- run_go_enrichment(
            deg_lists_up[[contrast]], paste0(contrast, "_up"), universe_genes)
    }
    if (contrast %in% names(deg_lists_down)) {
        go_results_down[[contrast]] <- run_go_enrichment(
            deg_lists_down[[contrast]], paste0(contrast, "_down"), universe_genes)
    }
}

# Save all GO results
go_notnull <- go_results[!sapply(go_results, is.null)]
saveRDS(go_results, file.path(RESULTS_DIR, "go_enrichment_all.rds"))

# Per-contrast plots
for (contrast in names(go_notnull)) {
    ego <- go_notnull[[contrast]]
    if (nrow(as.data.frame(ego)) == 0) next

    safe_name <- gsub(":", "_", contrast)

    write.table(as.data.frame(ego),
                file.path(RESULTS_DIR, paste0("GO_", safe_name, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # Dotplot
    tryCatch({
        p <- dotplot(ego, showCategory = 20,
                     title = paste0("GO Enrichment: ", contrast))
        ggsave(file.path(plot_dir, paste0("GO_dotplot_", safe_name, ".pdf")),
               p, width = 10, height = 8)
        ggsave(file.path(plot_dir, paste0("GO_dotplot_", safe_name, ".png")),
               p, width = 10, height = 8, dpi = 300)
    }, error = function(e) {
        cat("  dotplot error:", conditionMessage(e), "\n")
    })

    # Barplot
    tryCatch({
        p <- barplot(ego, showCategory = 15,
                     title = paste0("GO Enrichment: ", contrast))
        ggsave(file.path(plot_dir, paste0("GO_barplot_", safe_name, ".pdf")),
               p, width = 10, height = 7)
    }, error = function(e) {
        cat("  barplot error:", conditionMessage(e), "\n")
    })
}

###############################################################################
# 4. KEGG Enrichment per contrast
###############################################################################

cat("\n=== KEGG Enrichment ===\n")

kegg_results <- list()

for (contrast in names(deg_lists)) {
    gene_list <- deg_lists[[contrast]]
    if (length(gene_list) < 5) next

    safe_name <- gsub(":", "_", contrast)

    if (!is.null(TERM2GENE_KEGG)) {
        tryCatch({
            ekegg <- enricher(
                gene         = gene_list,
                universe     = universe_genes,
                TERM2GENE    = TERM2GENE_KEGG,
                TERM2NAME    = TERM2NAME_KEGG,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1
            )
            if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
                kegg_results[[contrast]] <- ekegg
                write.table(as.data.frame(ekegg),
                            file.path(RESULTS_DIR, paste0("KEGG_", safe_name, ".tsv")),
                            sep = "\t", row.names = FALSE, quote = FALSE)

                p <- dotplot(ekegg, showCategory = 20,
                             title = paste0("KEGG: ", contrast))
                ggsave(file.path(plot_dir, paste0("KEGG_dotplot_", safe_name, ".pdf")),
                       p, width = 10, height = 8)
                ggsave(file.path(plot_dir, paste0("KEGG_dotplot_", safe_name, ".png")),
                       p, width = 10, height = 8, dpi = 300)
            }
        }, error = function(e) {
            cat("  KEGG enricher error for", contrast, ":", conditionMessage(e), "\n")
        })

    } else if (use_orgdb && !is.null(ath_map)) {
        # Map to Arabidopsis and use KEGG API
        han_col_name <- grep("sunflower|helianthus|han|query|gene_id",
                             colnames(ath_map), ignore.case = TRUE, value = TRUE)[1]
        ath_col_name <- grep("ath|arabidopsis|tair|ortholog|target",
                             colnames(ath_map), ignore.case = TRUE, value = TRUE)[1]

        ath_genes <- ath_map[[ath_col_name]][ath_map[[han_col_name]] %in% gene_list]
        ath_genes <- unique(gsub("\\.\\d+$", "", ath_genes))
        ath_genes <- ath_genes[!is.na(ath_genes) & ath_genes != ""]

        # Convert TAIR to ENTREZ for KEGG
        tryCatch({
            library(org.At.tair.db)
            id_map <- bitr(ath_genes, fromType = "TAIR", toType = "ENTREZID",
                           OrgDb = org.At.tair.db)

            ekegg <- enrichKEGG(
                gene         = id_map$ENTREZID,
                organism     = "ath",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1
            )
            if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
                kegg_results[[contrast]] <- ekegg
                write.table(as.data.frame(ekegg),
                            file.path(RESULTS_DIR, paste0("KEGG_", safe_name, ".tsv")),
                            sep = "\t", row.names = FALSE, quote = FALSE)

                p <- dotplot(ekegg, showCategory = 20,
                             title = paste0("KEGG: ", contrast))
                ggsave(file.path(plot_dir, paste0("KEGG_dotplot_", safe_name, ".pdf")),
                       p, width = 10, height = 8)
            }
        }, error = function(e) {
            cat("  KEGG API error for", contrast, ":", conditionMessage(e), "\n")
        })
    }
}

saveRDS(kegg_results, file.path(RESULTS_DIR, "kegg_enrichment_all.rds"))

###############################################################################
# 5. compareCluster: cross-stress comparison
###############################################################################

cat("\n=== Cross-Stress Comparison ===\n")

if (length(deg_lists) >= 2) {
    # compareCluster for GO
    if (!is.null(TERM2GENE_GO)) {
        tryCatch({
            cc_go <- compareCluster(
                geneClusters  = deg_lists,
                fun           = "enricher",
                TERM2GENE     = TERM2GENE_GO,
                TERM2NAME     = TERM2NAME_GO,
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1
            )

            if (!is.null(cc_go) && nrow(as.data.frame(cc_go)) > 0) {
                saveRDS(cc_go, file.path(RESULTS_DIR, "compareCluster_GO.rds"))
                write.table(as.data.frame(cc_go),
                            file.path(RESULTS_DIR, "compareCluster_GO.tsv"),
                            sep = "\t", row.names = FALSE, quote = FALSE)

                p <- dotplot(cc_go, showCategory = 10,
                             title = "GO Enrichment Comparison Across Stresses") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                ggsave(file.path(plot_dir, "compareCluster_GO_dotplot.pdf"),
                       p, width = 14, height = 10)
                ggsave(file.path(plot_dir, "compareCluster_GO_dotplot.png"),
                       p, width = 14, height = 10, dpi = 300)
                cat("compareCluster GO plot saved\n")
            }
        }, error = function(e) {
            cat("compareCluster GO error:", conditionMessage(e), "\n")
        })
    } else if (use_orgdb) {
        # Using Arabidopsis orthologs for compareCluster
        tryCatch({
            han_col_name <- grep("sunflower|helianthus|han|query|gene_id",
                                 colnames(ath_map), ignore.case = TRUE, value = TRUE)[1]
            ath_col_name <- grep("ath|arabidopsis|tair|ortholog|target",
                                 colnames(ath_map), ignore.case = TRUE, value = TRUE)[1]

            ath_deg_lists <- lapply(deg_lists, function(genes) {
                ath <- ath_map[[ath_col_name]][ath_map[[han_col_name]] %in% genes]
                unique(gsub("\\.\\d+$", "", ath[!is.na(ath) & ath != ""]))
            })
            ath_deg_lists <- ath_deg_lists[sapply(ath_deg_lists, length) >= 5]

            if (length(ath_deg_lists) >= 2) {
                cc_go <- compareCluster(
                    geneClusters  = ath_deg_lists,
                    fun           = "enrichGO",
                    OrgDb         = org.At.tair.db,
                    keyType       = "TAIR",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1
                )

                if (!is.null(cc_go) && nrow(as.data.frame(cc_go)) > 0) {
                    saveRDS(cc_go, file.path(RESULTS_DIR, "compareCluster_GO.rds"))
                    write.table(as.data.frame(cc_go),
                                file.path(RESULTS_DIR, "compareCluster_GO.tsv"),
                                sep = "\t", row.names = FALSE, quote = FALSE)

                    p <- dotplot(cc_go, showCategory = 10,
                                 title = "GO BP Comparison Across Stresses") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                    ggsave(file.path(plot_dir, "compareCluster_GO_dotplot.pdf"),
                           p, width = 14, height = 10)
                    ggsave(file.path(plot_dir, "compareCluster_GO_dotplot.png"),
                           p, width = 14, height = 10, dpi = 300)
                }
            }
        }, error = function(e) {
            cat("compareCluster (Arabidopsis) error:", conditionMessage(e), "\n")
        })
    }

    # compareCluster for KEGG
    if (!is.null(TERM2GENE_KEGG)) {
        tryCatch({
            cc_kegg <- compareCluster(
                geneClusters  = deg_lists,
                fun           = "enricher",
                TERM2GENE     = TERM2GENE_KEGG,
                TERM2NAME     = TERM2NAME_KEGG,
                pvalueCutoff  = 0.05
            )
            if (!is.null(cc_kegg) && nrow(as.data.frame(cc_kegg)) > 0) {
                saveRDS(cc_kegg, file.path(RESULTS_DIR, "compareCluster_KEGG.rds"))
                write.table(as.data.frame(cc_kegg),
                            file.path(RESULTS_DIR, "compareCluster_KEGG.tsv"),
                            sep = "\t", row.names = FALSE, quote = FALSE)

                p <- dotplot(cc_kegg, showCategory = 10,
                             title = "KEGG Pathway Comparison Across Stresses") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                ggsave(file.path(plot_dir, "compareCluster_KEGG_dotplot.pdf"),
                       p, width = 14, height = 10)
                ggsave(file.path(plot_dir, "compareCluster_KEGG_dotplot.png"),
                       p, width = 14, height = 10, dpi = 300)
            }
        }, error = function(e) {
            cat("compareCluster KEGG error:", conditionMessage(e), "\n")
        })
    }

    # Up/Down comparison
    if (length(deg_lists_up) >= 2 && !is.null(TERM2GENE_GO)) {
        tryCatch({
            # Combine up and down with labels
            combined_lists <- c(
                setNames(deg_lists_up, paste0(names(deg_lists_up), "_UP")),
                setNames(deg_lists_down, paste0(names(deg_lists_down), "_DOWN"))
            )
            combined_lists <- combined_lists[sapply(combined_lists, length) >= 5]

            if (length(combined_lists) >= 2) {
                cc_updown <- compareCluster(
                    geneClusters = combined_lists,
                    fun          = "enricher",
                    TERM2GENE    = TERM2GENE_GO,
                    TERM2NAME    = TERM2NAME_GO,
                    pvalueCutoff = 0.05
                )
                if (!is.null(cc_updown) && nrow(as.data.frame(cc_updown)) > 0) {
                    saveRDS(cc_updown, file.path(RESULTS_DIR, "compareCluster_GO_updown.rds"))

                    p <- dotplot(cc_updown, showCategory = 8,
                                 title = "GO: Up vs Down Regulated Across Stresses") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
                    ggsave(file.path(plot_dir, "compareCluster_GO_updown.pdf"),
                           p, width = 16, height = 12)
                }
            }
        }, error = function(e) {
            cat("Up/Down compareCluster error:", conditionMessage(e), "\n")
        })
    }
}

###############################################################################
# 6. GSEA (Gene Set Enrichment Analysis)
###############################################################################

cat("\n=== GSEA ===\n")

if (!is.null(TERM2GENE_GO) && length(deg_full) > 0) {
    gsea_results <- list()

    for (contrast in names(deg_full)) {
        gene_rank <- deg_full[[contrast]]
        gene_rank <- gene_rank[!is.na(gene_rank) & is.finite(gene_rank)]
        gene_rank <- gene_rank[!duplicated(names(gene_rank))]

        if (length(gene_rank) < 100) next

        safe_name <- gsub(":", "_", contrast)
        cat("GSEA:", contrast, "(", length(gene_rank), "ranked genes)\n")

        tryCatch({
            gsea_res <- GSEA(
                geneList     = gene_rank,
                TERM2GENE    = TERM2GENE_GO,
                TERM2NAME    = TERM2NAME_GO,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                minGSSize    = 15,
                maxGSSize    = 500,
                verbose      = FALSE
            )
            if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
                gsea_results[[contrast]] <- gsea_res
                write.table(as.data.frame(gsea_res),
                            file.path(RESULTS_DIR, paste0("GSEA_", safe_name, ".tsv")),
                            sep = "\t", row.names = FALSE, quote = FALSE)

                # Ridge plot
                tryCatch({
                    p <- ridgeplot(gsea_res, showCategory = 15) +
                        ggtitle(paste0("GSEA Ridge: ", contrast))
                    ggsave(file.path(plot_dir, paste0("GSEA_ridge_", safe_name, ".pdf")),
                           p, width = 12, height = 10)
                }, error = function(e) {
                    cat("  ridgeplot error:", conditionMessage(e), "\n")
                })

                # Enrichment score plot for top terms
                tryCatch({
                    top_terms <- head(as.data.frame(gsea_res)$ID, 4)
                    for (i in seq_along(top_terms)) {
                        p <- gseaplot2(gsea_res, geneSetID = top_terms[i],
                                       title = top_terms[i])
                        ggsave(file.path(plot_dir,
                                         paste0("GSEA_plot_", safe_name, "_", i, ".pdf")),
                               p, width = 10, height = 7)
                    }
                }, error = function(e) {
                    cat("  gseaplot error:", conditionMessage(e), "\n")
                })
            }
        }, error = function(e) {
            cat("  GSEA error for", contrast, ":", conditionMessage(e), "\n")
        })
    }

    if (length(gsea_results) > 0) {
        saveRDS(gsea_results, file.path(RESULTS_DIR, "gsea_results_all.rds"))
    }
}

###############################################################################
# 7. Aquaporin-specific enrichment context
###############################################################################

cat("\n=== Aquaporin Gene Enrichment Context ===\n")

aqp_genes <- read.delim(AQP_LIST, header = TRUE, stringsAsFactors = FALSE)
aqp_id_col <- intersect(c("gene_id", "Gene_ID", "GeneID", "id"), colnames(aqp_genes))
if (length(aqp_id_col) == 0) aqp_id_col <- colnames(aqp_genes)[1] else aqp_id_col <- aqp_id_col[1]
aqp_ids <- aqp_genes[[aqp_id_col]]

# Check which AQP genes are DEGs per stress
aqp_deg_summary <- data.frame(
    contrast = names(deg_lists),
    total_DEGs = sapply(deg_lists, length),
    aqp_DEGs = sapply(deg_lists, function(x) sum(x %in% aqp_ids)),
    aqp_up = sapply(names(deg_lists), function(n) {
        if (n %in% names(deg_lists_up)) sum(deg_lists_up[[n]] %in% aqp_ids) else 0
    }),
    aqp_down = sapply(names(deg_lists), function(n) {
        if (n %in% names(deg_lists_down)) sum(deg_lists_down[[n]] %in% aqp_ids) else 0
    }),
    aqp_genes_list = sapply(deg_lists, function(x) {
        ids <- x[x %in% aqp_ids]
        paste(ids, collapse = ";")
    })
)

write.table(aqp_deg_summary, file.path(RESULTS_DIR, "aquaporin_DEG_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Aquaporin DEG summary:\n")
print(aqp_deg_summary[, 1:5])

###############################################################################
# Summary
###############################################################################

cat("\n=== Enrichment Analysis Complete ===\n")
cat("Results in:", RESULTS_DIR, "\n")
cat("GO results:", sum(!sapply(go_results, is.null)), "/", length(go_results), "contrasts\n")
cat("KEGG results:", sum(!sapply(kegg_results, is.null)), "/", length(kegg_results), "contrasts\n")

sessionInfo()

RSCRIPT_EOF

echo "=== GO/KEGG enrichment analysis finished: $(date) ==="
