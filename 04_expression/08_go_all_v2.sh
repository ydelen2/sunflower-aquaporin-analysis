#!/bin/bash
#SBATCH --job-name=go_all_v2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/go_all_v2_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/go_all_v2_%j.err

# =============================================================================
# 08_go_all_v2.sh
# GO over-representation analysis for the 16 stress contrasts used in the
# manuscript, run once with one explicit gene universe so that every contrast is
# tested against the same background. The first run had taken its universe from
# whichever results table was listed first; here the universe is the union of all
# genes with a DESeq2 result in any of the 16 contrasts, mapped to Arabidopsis
# identifiers. Settings as in 04_go_kegg_enrichment.sh (org.At.tair.db, BP/MF/CC,
# BH, p < 0.05, q < 0.1).
# Input : 04_expression/results/aquaporin_deg_v2/contrast_map.tsv (file -> label)
# Output: 04_expression/results/enrichment/GO_all_contrasts_v2.tsv
#         04_expression/results/enrichment/GO_all_contrasts_v2_universe.txt
# =============================================================================
set -eo pipefail
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
mkdir -p "${PROJ_DIR}/04_expression/logs" "${PROJ_DIR}/04_expression/results/enrichment"

module purge
module load miniforge/24.5
eval "$(conda shell.bash hook)"
conda activate /work/dweikat/ydelen2/aquaporin_study/conda_envs/aqp_env

Rscript --no-save --no-restore - <<'RSCRIPT_EOF'
suppressPackageStartupMessages({ library(clusterProfiler); library(org.At.tair.db) })
PROJ_DIR  <- "/work/dweikat/ydelen2/aquaporin_study"
DESEQ_DIR <- file.path(PROJ_DIR, "04_expression/results/deseq2")
MAP       <- file.path(PROJ_DIR, "04_expression/results/aquaporin_deg_v2/contrast_map.tsv")
OUT       <- file.path(PROJ_DIR, "04_expression/results/enrichment/GO_all_contrasts_v2.tsv")
OUT_UNI   <- file.path(PROJ_DIR, "04_expression/results/enrichment/GO_all_contrasts_v2_universe.txt")

ortho_file <- file.path(PROJ_DIR, "02_gene_family/results/arabidopsis_orthologs.tsv")
if (!file.exists(ortho_file)) ortho_file <- file.path(PROJ_DIR, "02_gene_family/results/arabidopsis_orthologs_v2.tsv")
ath_map <- read.delim(ortho_file, stringsAsFactors = FALSE)
han_col <- grep("sunflower|helianthus|han|query|gene_id", colnames(ath_map), ignore.case = TRUE, value = TRUE)
ath_col <- grep("ath|arabidopsis|tair|ortholog|target", colnames(ath_map), ignore.case = TRUE, value = TRUE)
han_col <- if (length(han_col)) han_col[1] else colnames(ath_map)[1]
ath_col <- if (length(ath_col)) ath_col[1] else colnames(ath_map)[2]
to_tair <- function(genes) {
    g <- ath_map[[ath_col]][ath_map[[han_col]] %in% genes]
    g <- unique(g[!is.na(g) & g != ""])
    gsub("\\.\\d+$", "", g)
}

cmap <- read.delim(MAP, stringsAsFactors = FALSE)
cat("contrasts:", nrow(cmap), "\n")
full <- lapply(cmap$file, function(f) read.delim(file.path(DESEQ_DIR, f), stringsAsFactors = FALSE))
universe <- unique(unlist(lapply(full, function(d) d$gene_id[!is.na(d$padj)])))
cat("universe (genes with a DESeq2 result in any contrast):", length(universe), "\n")
ath_universe <- to_tair(universe)
cat("universe mapped to TAIR:", length(ath_universe), "\n")
writeLines(universe, OUT_UNI)

# manuscript labels for the contrast column (project prefix : short contrast)
short <- c(Cold = "869:cold", Heat = "869:heat", Drought_869 = "869:drought", Salt_NaCl = "869:salt", Rehydration = "869:rehyd",
           Drought_7d = "797:d7", Drought_14d = "797:d14", Drought_21d = "797:d21", PEG_72h = "1041:PEG", Flooding = "492:flood",
           Sclerotinia = "908:scler", Orobanche_A = "850:pA", Orobanche_B = "850:pB", Orobanche_C = "850:pC", Orobanche_D = "850:pD", Orobanche_E = "850:pE")

res <- list()
for (i in seq_len(nrow(cmap))) {
    d <- full[[i]]
    sig <- d$gene_id[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) > 1]
    lab <- short[[cmap$label[i]]]
    ath_genes <- to_tair(sig)
    cat(sprintf("%-12s DEGs %5d  mapped %5d\n", lab, length(sig), length(ath_genes)))
    if (length(ath_genes) < 3) next
    for (ont in c("BP", "MF", "CC")) {
        ego <- tryCatch(enrichGO(gene = ath_genes, universe = ath_universe, OrgDb = org.At.tair.db,
                                 keyType = "TAIR", ont = ont, pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05, qvalueCutoff = 0.1),
                        error = function(e) { cat("enrichGO error", lab, ont, ":", conditionMessage(e), "\n"); NULL })
        if (!is.null(ego)) {
            df <- as.data.frame(ego)
            if (nrow(df)) { df$contrast <- lab; df$ontology <- ont; res[[paste(lab, ont)]] <- df }
        }
    }
}
out <- do.call(rbind, res)
out <- out[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "contrast", "ontology")]
write.table(out, OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat("written", OUT, "with", nrow(out), "rows\n")
print(table(out$contrast))
RSCRIPT_EOF
