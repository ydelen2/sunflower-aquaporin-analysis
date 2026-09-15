#!/bin/bash
#SBATCH --job-name=deseq2_492303
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/deseq2_492303_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/deseq2_492303_%j.err

# =============================================================================
# 06_deseq2_492303_full_v2.sh
# DESeq2 for PRJNA492303 (flooding) on all 96 runs with the design that the
# experiment actually has: ~ genotype + tissue + age + condition
# (HA351/RHA428, leaf/root, 16/24/30 days, control/flooding, 4 replicates).
# Replaces the 46-run analysis whose design (~ tissue + condition) was
# confounded by an unbalanced subset of runs.
# Inputs : 03_rnaseq/counts/PRJNA492303_full_counts_clean.txt
#          03_rnaseq/accession_lists/PRJNA492303_metadata_full.tsv
# Output : 04_expression/results/deseq2/PRJNA492303_fixed/  (same file layout as
#          the other *_fixed directories, so 04_aquaporin_deg_table_v2.sh picks it up)
# =============================================================================
set -eo pipefail
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
mkdir -p "${PROJ_DIR}/04_expression/logs"

module purge
module load miniforge/24.5
eval "$(conda shell.bash hook)"
conda activate /work/dweikat/ydelen2/aquaporin_study/conda_envs/aqp_env

Rscript --no-save --no-restore - <<'RSCRIPT_EOF'
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
PROJ_DIR <- "/work/dweikat/ydelen2/aquaporin_study"
counts_file <- file.path(PROJ_DIR, "03_rnaseq/counts/PRJNA492303_full_counts_clean.txt")
meta_file   <- file.path(PROJ_DIR, "03_rnaseq/accession_lists/PRJNA492303_metadata_full.tsv")
out_dir     <- file.path(PROJ_DIR, "04_expression/results/deseq2/PRJNA492303_fixed")
plot_dir    <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
meta   <- read.delim(meta_file, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
common <- intersect(colnames(counts), rownames(meta))
cat("count columns:", ncol(counts), "| metadata rows:", nrow(meta), "| in common:", length(common), "\n")
if (length(common) != 96) stop("expected 96 runs in common, found ", length(common))
counts <- counts[, common]
meta   <- meta[common, ]

meta$genotype  <- factor(meta$genotype)
meta$tissue    <- factor(meta$tissue, levels = c("leaf", "root"))
meta$age       <- factor(meta$age, levels = c("d16", "d24", "d30"))
meta$condition <- relevel(factor(meta$condition), ref = "control")
print(table(meta$genotype, meta$condition, meta$tissue))

keep <- rowSums(counts) >= 10
cat("genes before filtering:", nrow(counts), "| after:", sum(keep), "\n")
counts <- counts[keep, ]

BiocParallel::register(BiocParallel::MulticoreParam(workers = 8))
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = meta,
                              design = ~ genotype + tissue + age + condition)
dds <- DESeq(dds, parallel = TRUE)
saveRDS(dds, file.path(out_dir, "dds.rds"))

norm_counts <- counts(dds, normalized = TRUE)
write.table(norm_counts, file.path(out_dir, "normalized_counts.tsv"), sep = "\t", quote = FALSE)
saveRDS(norm_counts, file.path(out_dir, "normalized_counts.rds"))
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(out_dir, "vsd.rds"))
write.table(meta, file.path(out_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pca <- plotPCA(vsd, intgroup = c("condition", "tissue", "genotype"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color = condition, shape = tissue)) + geom_point(size = 3) +
     facet_wrap(~ genotype) +
     xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) +
     theme_bw(base_size = 12)
ggsave(file.path(plot_dir, "pca_plot.pdf"), p, width = 9, height = 5)
ggsave(file.path(plot_dir, "pca_plot.png"), p, width = 9, height = 5, dpi = 300)

contrast_name <- "flooding_vs_control"
res <- results(dds, contrast = c("condition", "flooding", "control"), alpha = 0.05, pAdjustMethod = "BH")
res <- res[order(res$padj), ]
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
write.table(res_df, file.path(out_dir, paste0("DEG_", contrast_name, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(res, file.path(out_dir, paste0("DEG_", contrast_name, ".rds")))
sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
write.table(sig, file.path(out_dir, paste0("DEG_sig_", contrast_name, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

summary_df <- data.frame(contrast = contrast_name,
                         total_tested = sum(!is.na(res$padj)),
                         deg_005 = sum(res$padj < 0.05, na.rm = TRUE),
                         deg_005_lfc1 = nrow(sig),
                         up_005_lfc1 = sum(sig$log2FoldChange > 0),
                         down_005_lfc1 = sum(sig$log2FoldChange < 0))
write.table(summary_df, file.path(out_dir, "DEG_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(summary_df)

pdf(file.path(plot_dir, paste0("MA_", contrast_name, ".pdf")), width = 8, height = 6)
plotMA(res, main = paste0("MA plot: ", contrast_name), ylim = c(-5, 5))
dev.off()

# design note for the records
writeLines(c("PRJNA492303 full design: 96 runs, ~ genotype + tissue + age + condition",
             "contrast: condition flooding vs control",
             paste0("run on ", format(Sys.time(), "%Y-%m-%d"))),
           file.path(out_dir, "DESIGN_NOTE.txt"))
cat("done\n")
RSCRIPT_EOF
