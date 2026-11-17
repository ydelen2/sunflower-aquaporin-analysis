#!/bin/bash
#SBATCH --job-name=deseq2_v3
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/deseq2_v3_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/deseq2_v3_%j.err

# =============================================================================
# 10_deseq2_designs_v3.sh
# Documents and tests the DESeq2 designs behind the 16 stress contrasts.
#
# Part A  Re-fits every design from the raw count matrix and compares the
#         result with the stored per-contrast tables
#         (04_expression/results/deseq2/<project>[_fixed]/DEG_*.tsv). For
#         PRJNA869183, whose grouping was not kept in a script, several
#         candidate groupings are fitted and the one that reproduces the
#         stored table is reported.
# Part B  Tests whether the treatment response depends on genotype, tissue or
#         plant age in the datasets that have these factors: likelihood-ratio
#         tests of the interaction terms and treatment contrasts within each
#         genotype, tissue or age (Wald tests).
# Part C  Applies one Benjamini-Hochberg correction to all tests of the 16
#         contrasts together and compares it with the per-contrast correction.
#
# Inputs
#   04_expression/results/deseq2/raw_counts.rds, sample_metadata.tsv
#   03_rnaseq/counts/PRJNA492303_full_counts_clean.txt
#   03_rnaseq/accession_lists/PRJNA492303_metadata_full.tsv
#   04_expression/results/aquaporin_deg_v2/contrast_map.tsv
#   02_gene_family/results/verified_aquaporins.txt
# Outputs (04_expression/results/deseq2_designs_v3/)
#   reproduction_report.tsv       candidate design x stored contrast: agreement
#   design_cells.tsv              libraries per factor combination in each dataset
#   interaction_tests_summary.tsv genome-wide and aquaporin counts per test
#   interaction_tests_aquaporin.tsv   LRT results for the aquaporin genes
#   stratified_contrasts_summary.tsv  genome-wide and aquaporin counts per stratum
#   stratified_contrasts_aquaporin.tsv  Wald results for the aquaporin genes
#   joint_fdr_aquaporin.tsv       per-contrast and joint adjusted P values
#   joint_fdr_summary.tsv
# =============================================================================
set -eo pipefail
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
mkdir -p "${PROJ_DIR}/04_expression/logs" "${PROJ_DIR}/04_expression/results/deseq2_designs_v3"

module purge
module load miniforge/24.5
eval "$(conda shell.bash hook)"
conda activate /work/dweikat/ydelen2/aquaporin_study/conda_envs/aqp_env

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RS="${PROJ_DIR}/04_expression/results/deseq2_designs_v3/deseq2_designs_v3.R"
cat > "${RS}" <<'RSCRIPT_EOF'
suppressPackageStartupMessages({ library(DESeq2); library(BiocParallel) })
PROJ  <- "/work/dweikat/ydelen2/aquaporin_study"
DESEQ <- file.path(PROJ, "04_expression/results/deseq2")
OUT   <- file.path(PROJ, "04_expression/results/deseq2_designs_v3")
.libPaths(c(file.path(PROJ, "R_libs"), .libPaths()))
BiocParallel::register(SerialParam())
cat("R:", R.version.string, "| DESeq2:", as.character(packageVersion("DESeq2")), "\n")

raw  <- readRDS(file.path(DESEQ, "raw_counts.rds"))
meta <- read.delim(file.path(DESEQ, "sample_metadata.tsv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
cat("raw counts:", nrow(raw), "genes x", ncol(raw), "libraries | metadata rows:", nrow(meta), "\n")
aqp <- read.delim(file.path(PROJ, "02_gene_family/results/verified_aquaporins.txt"), stringsAsFactors = FALSE)$gene_id
cat("aquaporin genes:", length(aqp), "\n")

# ------------------------------------------------------------------ helpers
fit <- function(cts, coldata, design, test = "Wald", reduced = NULL) {
    cts <- round(as.matrix(cts)); storage.mode(cts) <- "integer"
    stopifnot(identical(colnames(cts), rownames(coldata)))
    keep <- rowSums(cts) >= 10
    dds <- DESeqDataSetFromMatrix(countData = cts[keep, , drop = FALSE], colData = coldata, design = design)
    cat(format(Sys.time(), "%H:%M:%S"), "fit:", ncol(cts), "libraries,", sum(keep), "genes,", test, deparse(design), "\n")
    if (test == "LRT") DESeq(dds, test = "LRT", reduced = reduced, quiet = TRUE)
    else DESeq(dds, quiet = TRUE)
}
wald <- function(dds, factor_name, num, den) {
    as.data.frame(results(dds, contrast = c(factor_name, num, den), alpha = 0.05, pAdjustMethod = "BH"))
}
n_deg <- function(df) {
    s <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1
    c(total = sum(s), up = sum(s & df$log2FoldChange > 0), down = sum(s & df$log2FoldChange < 0))
}
read_pub <- function(path) {
    d <- read.delim(path, stringsAsFactors = FALSE)
    if (!("gene_id" %in% colnames(d))) { d$gene_id <- rownames(d) }
    rownames(d) <- d$gene_id
    d
}
compare <- function(new, pub_path, dataset, contrast, candidate) {
    pub <- read_pub(pub_path)
    g <- intersect(rownames(new), rownames(pub))
    ok <- g[!is.na(new[g, "log2FoldChange"]) & !is.na(pub[g, "log2FoldChange"])]
    dl <- abs(new[ok, "log2FoldChange"] - pub[ok, "log2FoldChange"])
    okp <- g[!is.na(new[g, "pvalue"]) & !is.na(pub[g, "pvalue"]) & new[g, "pvalue"] > 0 & pub[g, "pvalue"] > 0]
    a <- n_deg(new); b <- n_deg(pub)
    data.frame(dataset = dataset, contrast = contrast, candidate = candidate,
               genes_new = nrow(new), genes_published = nrow(pub), genes_common = length(g),
               tested_new = sum(!is.na(new$padj)), tested_published = sum(!is.na(pub$padj)),
               deg_new = a[["total"]], up_new = a[["up"]], down_new = a[["down"]],
               deg_published = b[["total"]], up_published = b[["up"]], down_published = b[["down"]],
               lfc_correlation = round(cor(new[ok, "log2FoldChange"], pub[ok, "log2FoldChange"]), 6),
               lfc_max_abs_diff = signif(max(dl), 4), lfc_median_abs_diff = signif(median(dl), 4),
               log10p_correlation = round(cor(log10(new[okp, "pvalue"]), log10(pub[okp, "pvalue"])), 6),
               exact = max(dl) < 1e-6 && a[["total"]] == b[["total"]] && a[["up"]] == b[["up"]],
               stringsAsFactors = FALSE)
}
report <- list(); add_report <- function(x) report[[length(report) + 1]] <<- x
cells <- list()
add_cells <- function(dataset, coldata, cols) {
    tb <- as.data.frame(table(coldata[, cols, drop = FALSE]), stringsAsFactors = FALSE)
    tb <- tb[tb$Freq > 0, , drop = FALSE]
    for (cc in c("genotype", "tissue", "age", "condition")) if (!(cc %in% colnames(tb))) tb[[cc]] <- ""
    cells[[length(cells) + 1]] <<- data.frame(dataset = dataset, genotype = tb$genotype, tissue = tb$tissue, age = tb$age,
                                              condition = tb$condition, libraries = tb$Freq, stringsAsFactors = FALSE)
}
int_sum <- list(); int_aqp <- list(); str_sum <- list(); str_aqp <- list()
run_lrt <- function(dataset, cts, coldata, full, reduced, term) {
    r <- tryCatch({
        d <- fit(cts, coldata, full, test = "LRT", reduced = reduced)
        as.data.frame(results(d, alpha = 0.05, pAdjustMethod = "BH"))
    }, error = function(e) { cat("LRT failed:", dataset, term, conditionMessage(e), "\n"); NULL })
    if (is.null(r)) return(invisible(NULL))
    sig <- !is.na(r$padj) & r$padj < 0.05
    a <- intersect(aqp, rownames(r))
    int_sum[[length(int_sum) + 1]] <<- data.frame(dataset = dataset, term = term, full_design = paste(deparse(full), collapse = ""),
        reduced_design = paste(deparse(reduced), collapse = ""), genes_tested = sum(!is.na(r$padj)), genes_significant = sum(sig),
        aquaporins_tested = sum(!is.na(r[a, "padj"])), aquaporins_significant = sum(sig[match(a, rownames(r))]), stringsAsFactors = FALSE)
    int_aqp[[length(int_aqp) + 1]] <<- data.frame(dataset = dataset, term = term, gene_id = a, baseMean = r[a, "baseMean"],
        stat = r[a, "stat"], pvalue = r[a, "pvalue"], padj = r[a, "padj"], stringsAsFactors = FALSE)
    cat(sprintf("  LRT %-14s %-22s tested %6d  significant %6d  aquaporins %d of %d\n", dataset, term, sum(!is.na(r$padj)), sum(sig),
                sum(sig[match(a, rownames(r))]), sum(!is.na(r[a, "padj"]))))
}
run_strata <- function(dataset, cts, coldata, design, factor_name, pairs, stratum_type) {
    d <- tryCatch(fit(cts, coldata, design), error = function(e) { cat("stratified fit failed:", dataset, conditionMessage(e), "\n"); NULL })
    if (is.null(d)) return(invisible(NULL))
    for (p in pairs) {
        r <- wald(d, factor_name, p[["num"]], p[["den"]])
        k <- n_deg(r); a <- intersect(aqp, rownames(r)); ra <- r[a, ]
        sa <- !is.na(ra$padj) & ra$padj < 0.05 & abs(ra$log2FoldChange) > 1
        str_sum[[length(str_sum) + 1]] <<- data.frame(dataset = dataset, stratum_type = stratum_type, stratum = p[["stratum"]],
            contrast = p[["label"]], design = paste(deparse(design), collapse = ""), genes_tested = sum(!is.na(r$padj)),
            deg = k[["total"]], up = k[["up"]], down = k[["down"]], aquaporin_deg = sum(sa),
            aquaporin_up = sum(sa & ra$log2FoldChange > 0), aquaporin_down = sum(sa & ra$log2FoldChange < 0), stringsAsFactors = FALSE)
        str_aqp[[length(str_aqp) + 1]] <<- data.frame(dataset = dataset, stratum_type = stratum_type, stratum = p[["stratum"]],
            contrast = p[["label"]], gene_id = a, baseMean = ra$baseMean, log2FoldChange = ra$log2FoldChange, lfcSE = ra$lfcSE,
            pvalue = ra$pvalue, padj = ra$padj, stringsAsFactors = FALSE)
        cat(sprintf("  strata %-14s %-10s %-34s DEG %6d  aquaporin DEG %d (%d up, %d down)\n", dataset, p[["stratum"]], p[["label"]],
                    k[["total"]], sum(sa), sum(sa & ra$log2FoldChange > 0), sum(sa & ra$log2FoldChange < 0)))
    }
}

# =============================================================== PRJNA869183
cat("\n=== PRJNA869183: candidate groupings ===\n")
m <- meta[meta$project_id == "PRJNA869183", ]; s <- intersect(colnames(raw), m$sample_id); m <- m[s, ]
cat("libraries:", length(s), "\n"); print(table(m$condition))
add_cells("PRJNA869183", data.frame(condition = m$condition), "condition")
cold_t <- c("cold_8h", "cold_16h", "cold_32h"); heat_t <- c("heat_4h", "heat_8h", "heat_16h", "heat_32h")
sa <- c("salt_SA1", "salt_SA2", "salt_SA3", "salt_SA4")
grouping <- function(cond, drought_stress, pool_sa, pool_temp = TRUE) {
    g <- cond
    if (pool_temp) { g[cond %in% cold_t] <- "cold_stress"; g[cond %in% heat_t] <- "heat_stress" }
    g[cond %in% drought_stress] <- "drought_stress"
    g[cond == "drought_normal"] <- "drought_control"
    if (pool_sa) g[cond %in% sa] <- "salt_SA"
    g
}
G <- list(
    G1_drought_MD_SD_WD__SA_separate = grouping(m$condition, c("drought_moderate", "drought_severe", "drought_water_deficit"), FALSE),
    G2_drought_MD_SD__SA_separate    = grouping(m$condition, c("drought_moderate", "drought_severe"), FALSE),
    G3_drought_MD_SD__SA_pooled      = grouping(m$condition, c("drought_moderate", "drought_severe"), TRUE),
    G4_drought_MD_SD_WD__SA_pooled   = grouping(m$condition, c("drought_moderate", "drought_severe", "drought_water_deficit"), TRUE),
    G5_drought_levels_separate__SA_separate = grouping(m$condition, character(0), FALSE),
    G6_drought_levels_separate__SA_pooled   = grouping(m$condition, character(0), TRUE),
    G7_all_19_conditions = { g <- m$condition; g[g == "drought_normal"] <- "drought_control"; g }
)
pub869 <- list(
    cold        = list(file = "PRJNA869183_fixed/DEG_cold_stress_vs_cold_0h.tsv", num = "cold_stress", den = "cold_0h"),
    heat        = list(file = "PRJNA869183_fixed/DEG_heat_stress_vs_heat_0h.tsv", num = "heat_stress", den = "heat_0h"),
    drought     = list(file = "PRJNA869183_fixed/DEG_drought_stress_vs_drought_control.tsv", num = "drought_stress", den = "drought_control"),
    rehydration = list(file = "PRJNA869183_fixed/DEG_rehydrated_vs_severe.tsv", num = "drought_water_deficit", den = "drought_severe"),
    salt        = list(file = "PRJNA869183_fixed/DEG_salt_NaCl_vs_control.tsv", num = "salt_NaCl", den = "control")
)
for (gn in names(G)) {
    cd <- data.frame(group = factor(G[[gn]]), row.names = s)
    cat("fitting", gn, "levels:", nlevels(cd$group), "\n")
    d <- tryCatch(fit(raw[, s], cd, ~ group), error = function(e) { cat("  failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(d)) next
    for (cn in names(pub869)) {
        p <- pub869[[cn]]
        if (!all(c(p$num, p$den) %in% levels(cd$group))) next
        add_report(compare(wald(d, "group", p$num, p$den), file.path(DESEQ, p$file), "PRJNA869183", cn, paste0("all57__", gn)))
    }
}
# two-group and stress-family subsets
subsets <- list(
    cold__two_groups        = list(con = "cold", keep = c(cold_t, "cold_0h"), map = list(cold_stress = cold_t, cold_0h = "cold_0h")),
    heat__two_groups        = list(con = "heat", keep = c(heat_t, "heat_0h"), map = list(heat_stress = heat_t, heat_0h = "heat_0h")),
    drought__MD_SD_vs_ND    = list(con = "drought", keep = c("drought_moderate", "drought_severe", "drought_normal"),
                                   map = list(drought_stress = c("drought_moderate", "drought_severe"), drought_control = "drought_normal")),
    drought__MD_SD_WD_vs_ND = list(con = "drought", keep = c("drought_moderate", "drought_severe", "drought_water_deficit", "drought_normal"),
                                   map = list(drought_stress = c("drought_moderate", "drought_severe", "drought_water_deficit"), drought_control = "drought_normal")),
    drought__family_MD_SD_pooled = list(con = "drought", keep = c("drought_moderate", "drought_severe", "drought_water_deficit", "drought_normal"),
                                   map = list(drought_stress = c("drought_moderate", "drought_severe"), drought_control = "drought_normal", drought_water_deficit = "drought_water_deficit")),
    rehydration__two_groups = list(con = "rehydration", keep = c("drought_water_deficit", "drought_severe"),
                                   map = list(drought_water_deficit = "drought_water_deficit", drought_severe = "drought_severe")),
    rehydration__drought_family = list(con = "rehydration", keep = c("drought_moderate", "drought_severe", "drought_water_deficit", "drought_normal"),
                                   map = list(drought_water_deficit = "drought_water_deficit", drought_severe = "drought_severe",
                                              drought_moderate = "drought_moderate", drought_normal = "drought_normal")),
    salt__two_groups        = list(con = "salt", keep = c("salt_NaCl", "control"), map = list(salt_NaCl = "salt_NaCl", control = "control")),
    salt__family_SA_separate = list(con = "salt", keep = c("salt_NaCl", "control", sa),
                                   map = list(salt_NaCl = "salt_NaCl", control = "control", salt_SA1 = "salt_SA1", salt_SA2 = "salt_SA2", salt_SA3 = "salt_SA3", salt_SA4 = "salt_SA4"))
)
for (sn in names(subsets)) {
    x <- subsets[[sn]]; ss <- s[m$condition %in% x$keep]
    g <- m[ss, "condition"]; for (lv in names(x$map)) g[m[ss, "condition"] %in% x$map[[lv]]] <- lv
    cd <- data.frame(group = factor(g), row.names = ss)
    d <- tryCatch(fit(raw[, ss], cd, ~ group), error = function(e) { cat("  failed:", sn, conditionMessage(e), "\n"); NULL })
    if (is.null(d)) next
    p <- pub869[[x$con]]
    add_report(compare(wald(d, "group", p$num, p$den), file.path(DESEQ, p$file), "PRJNA869183", x$con, paste0("subset__", sn)))
}

# =============================================================== PRJNA797473
cat("\n=== PRJNA797473 ===\n")
m <- meta[meta$project_id == "PRJNA797473", ]; s <- intersect(colnames(raw), m$sample_id); m <- m[s, ]
cd <- data.frame(genotype = factor(m$genotype), condition = relevel(factor(m$condition), ref = "control"), row.names = s)
print(table(cd$genotype, cd$condition)); add_cells("PRJNA797473", cd, c("genotype", "condition"))
d <- fit(raw[, s], cd, ~ genotype + condition)
for (lv in setdiff(levels(cd$condition), "control"))
    add_report(compare(wald(d, "condition", lv, "control"), file.path(DESEQ, "PRJNA797473_fixed", paste0("DEG_", lv, "_vs_control.tsv")),
                       "PRJNA797473", lv, "~genotype+condition"))
run_lrt("PRJNA797473", raw[, s], cd, ~ genotype + condition + genotype:condition, ~ genotype + condition, "genotype:condition")
cd$group <- factor(paste(cd$genotype, cd$condition, sep = "_"))
prs <- list(); for (gt in levels(cd$genotype)) for (lv in setdiff(levels(cd$condition), "control"))
    prs[[length(prs) + 1]] <- list(stratum = gt, label = paste0(lv, "_vs_control"), num = paste(gt, lv, sep = "_"), den = paste(gt, "control", sep = "_"))
run_strata("PRJNA797473", raw[, s], cd, ~ group, "group", prs, "genotype")

# ============================================================== PRJNA1041959
cat("\n=== PRJNA1041959 ===\n")
m <- meta[meta$project_id == "PRJNA1041959", ]; s <- intersect(colnames(raw), m$sample_id); m <- m[s, ]
cd <- data.frame(genotype = factor(m$genotype), tissue = factor(m$tissue), condition = relevel(factor(m$condition), ref = "control"), row.names = s)
print(table(cd$genotype, cd$tissue, cd$condition)); add_cells("PRJNA1041959", cd, c("genotype", "tissue", "condition"))
d <- fit(raw[, s], cd, ~ genotype + tissue + condition)
trt <- setdiff(levels(cd$condition), "control")
add_report(compare(wald(d, "condition", trt, "control"), file.path(DESEQ, "PRJNA1041959_fixed", paste0("DEG_", trt, "_vs_control.tsv")),
                   "PRJNA1041959", trt, "~genotype+tissue+condition"))
full <- ~ genotype + tissue + condition + tissue:condition + genotype:condition
run_lrt("PRJNA1041959", raw[, s], cd, full, ~ genotype + tissue + condition + genotype:condition, "tissue:condition")
run_lrt("PRJNA1041959", raw[, s], cd, full, ~ genotype + tissue + condition + tissue:condition, "genotype:condition")
cd$tc <- factor(paste(cd$tissue, cd$condition, sep = "_")); cd$gc <- factor(paste(cd$genotype, cd$condition, sep = "_"))
run_strata("PRJNA1041959", raw[, s], cd, ~ genotype + tc, "tc",
           lapply(levels(cd$tissue), function(x) list(stratum = x, label = paste0(trt, "_vs_control"), num = paste(x, trt, sep = "_"), den = paste(x, "control", sep = "_"))), "tissue")
run_strata("PRJNA1041959", raw[, s], cd, ~ tissue + gc, "gc",
           lapply(levels(cd$genotype), function(x) list(stratum = x, label = paste0(trt, "_vs_control"), num = paste(x, trt, sep = "_"), den = paste(x, "control", sep = "_"))), "genotype")

# =============================================================== PRJNA908908
cat("\n=== PRJNA908908 ===\n")
m <- meta[meta$project_id == "PRJNA908908", ]; s <- intersect(colnames(raw), m$sample_id); m <- m[s, ]
cd <- data.frame(genotype = factor(m$genotype), condition = relevel(factor(m$condition), ref = "control"), row.names = s)
print(table(cd$genotype, cd$condition)); add_cells("PRJNA908908", cd, c("genotype", "condition"))
d <- fit(raw[, s], cd, ~ genotype + condition)
trt <- setdiff(levels(cd$condition), "control")
add_report(compare(wald(d, "condition", trt, "control"), file.path(DESEQ, "PRJNA908908", "DEG_sclerotinia_vs_control.tsv"),
                   "PRJNA908908", trt, "~genotype+condition"))
run_lrt("PRJNA908908", raw[, s], cd, ~ genotype + condition + genotype:condition, ~ genotype + condition, "genotype:condition")
cd$group <- factor(paste(cd$genotype, cd$condition, sep = "_"))
run_strata("PRJNA908908", raw[, s], cd, ~ group, "group",
           lapply(levels(cd$genotype), function(x) list(stratum = x, label = paste0(trt, "_vs_control"), num = paste(x, trt, sep = "_"), den = paste(x, "control", sep = "_"))), "genotype")

# =============================================================== PRJNA850121
cat("\n=== PRJNA850121 ===\n")
m <- meta[meta$project_id == "PRJNA850121", ]; s <- intersect(colnames(raw), m$sample_id); m <- m[s, ]
cd <- data.frame(genotype = factor(ifelse(m$genotype == "K", "HZ2399", "SQ25")), condition = relevel(factor(m$condition), ref = "control"), row.names = s)
print(table(m$genotype, m$condition)); print(table(cd$genotype, cd$condition)); add_cells("PRJNA850121", cd, c("genotype", "condition"))
d <- fit(raw[, s], cd, ~ genotype + condition)
stages <- setdiff(levels(cd$condition), "control")
for (lv in stages)
    add_report(compare(wald(d, "condition", lv, "control"), file.path(DESEQ, "PRJNA850121_fixed", paste0("DEG_", lv, "_vs_control.tsv")),
                       "PRJNA850121", lv, "~genotype_fixed+condition"))
run_lrt("PRJNA850121", raw[, s], cd, ~ genotype + condition + genotype:condition, ~ genotype + condition, "genotype:condition")
cd$group <- factor(paste(cd$genotype, cd$condition, sep = "_"))
prs <- list(); for (gt in levels(cd$genotype)) for (lv in stages)
    prs[[length(prs) + 1]] <- list(stratum = gt, label = paste0(lv, "_vs_control"), num = paste(gt, lv, sep = "_"), den = paste(gt, "control", sep = "_"))
run_strata("PRJNA850121", raw[, s], cd, ~ group, "group", prs, "genotype")

# =============================================================== PRJNA492303
cat("\n=== PRJNA492303 (96 libraries) ===\n")
fc <- read.delim(file.path(PROJ, "03_rnaseq/counts/PRJNA492303_full_counts_clean.txt"), row.names = 1, check.names = FALSE)
fm <- read.delim(file.path(PROJ, "03_rnaseq/accession_lists/PRJNA492303_metadata_full.tsv"), stringsAsFactors = FALSE)
rownames(fm) <- fm$sample_id; s <- intersect(colnames(fc), fm$sample_id); fm <- fm[s, ]
cd <- data.frame(genotype = factor(fm$genotype), tissue = factor(fm$tissue), age = factor(fm$age),
                 condition = relevel(factor(fm$condition), ref = "control"), row.names = s)
print(table(cd$genotype, cd$tissue, cd$condition)); add_cells("PRJNA492303", cd, c("genotype", "tissue", "age", "condition"))
d <- fit(fc[, s], cd, ~ genotype + tissue + age + condition)
add_report(compare(wald(d, "condition", "flooding", "control"), file.path(DESEQ, "PRJNA492303_fixed", "DEG_flooding_vs_control.tsv"),
                   "PRJNA492303", "flooding", "~genotype+tissue+age+condition"))
full <- ~ genotype + tissue + age + condition + tissue:condition + genotype:condition + age:condition
run_lrt("PRJNA492303", fc[, s], cd, full, ~ genotype + tissue + age + condition + genotype:condition + age:condition, "tissue:condition")
run_lrt("PRJNA492303", fc[, s], cd, full, ~ genotype + tissue + age + condition + tissue:condition + age:condition, "genotype:condition")
run_lrt("PRJNA492303", fc[, s], cd, full, ~ genotype + tissue + age + condition + tissue:condition + genotype:condition, "age:condition")
cd$tc <- factor(paste(cd$tissue, cd$condition, sep = "_")); cd$gc <- factor(paste(cd$genotype, cd$condition, sep = "_"))
cd$ac <- factor(paste(cd$age, cd$condition, sep = "_"))
mk <- function(lv) lapply(lv, function(x) list(stratum = x, label = "flooding_vs_control", num = paste(x, "flooding", sep = "_"), den = paste(x, "control", sep = "_")))
run_strata("PRJNA492303", fc[, s], cd, ~ genotype + age + tc, "tc", mk(levels(cd$tissue)), "tissue")
run_strata("PRJNA492303", fc[, s], cd, ~ tissue + age + gc, "gc", mk(levels(cd$genotype)), "genotype")
run_strata("PRJNA492303", fc[, s], cd, ~ genotype + tissue + ac, "ac", mk(levels(cd$age)), "age")

# ------------------------------------------------------------------ write A, B
wr <- function(x, name) write.table(x, file.path(OUT, name), sep = "\t", quote = FALSE, row.names = FALSE)
rep_df <- do.call(rbind, report); wr(rep_df, "reproduction_report.tsv")
wr(do.call(rbind, cells), "design_cells.tsv")
wr(do.call(rbind, int_sum), "interaction_tests_summary.tsv"); wr(do.call(rbind, int_aqp), "interaction_tests_aquaporin.tsv")
wr(do.call(rbind, str_sum), "stratified_contrasts_summary.tsv"); wr(do.call(rbind, str_aqp), "stratified_contrasts_aquaporin.tsv")
cat("\n=== Reproduction report ===\n")
print(rep_df[, c("dataset", "contrast", "candidate", "deg_new", "deg_published", "lfc_correlation", "lfc_max_abs_diff", "log10p_correlation", "exact")], row.names = FALSE)

# =============================================================== Part C
cat("\n=== Joint correction across the 16 contrasts ===\n")
cmap <- read.delim(file.path(PROJ, "04_expression/results/aquaporin_deg_v2/contrast_map.tsv"), stringsAsFactors = FALSE)
allp <- list()
for (i in seq_len(nrow(cmap))) {
    x <- read_pub(file.path(DESEQ, cmap$file[i]))
    x <- x[!is.na(x$padj), c("gene_id", "log2FoldChange", "pvalue", "padj")]
    x$contrast <- cmap$label[i]; rownames(x) <- NULL; allp[[i]] <- x
}
allp <- do.call(rbind, allp)
allp$padj_joint <- p.adjust(allp$pvalue, method = "BH")
cat("tests entering the joint correction:", nrow(allp), "| contrasts:", length(unique(allp$contrast)), "\n")
A <- allp[allp$gene_id %in% aqp, ]
A$de_within <- A$padj < 0.05 & abs(A$log2FoldChange) > 1
A$de_joint  <- A$padj_joint < 0.05 & abs(A$log2FoldChange) > 1
A$de_bonferroni_contrasts <- A$padj < 0.05 / length(unique(allp$contrast)) & abs(A$log2FoldChange) > 1
wr(A[, c("gene_id", "contrast", "log2FoldChange", "pvalue", "padj", "padj_joint", "de_within", "de_joint", "de_bonferroni_contrasts")], "joint_fdr_aquaporin.tsv")
cnt <- function(flag) c(calls = sum(flag), genes_ge1 = length(unique(A$gene_id[flag])), genes_ge2 = sum(table(A$gene_id[flag]) >= 2))
js <- rbind(per_contrast_BH = cnt(A$de_within), joint_BH_all_tests = cnt(A$de_joint), bonferroni_over_contrasts = cnt(A$de_bonferroni_contrasts))
js <- data.frame(correction = rownames(js), js, aquaporin_genes_total = length(aqp), row.names = NULL)
gw <- allp$padj < 0.05 & abs(allp$log2FoldChange) > 1; gj <- allp$padj_joint < 0.05 & abs(allp$log2FoldChange) > 1
js$genome_wide_calls <- c(sum(gw), sum(gj), sum(allp$padj < 0.05 / length(unique(allp$contrast)) & abs(allp$log2FoldChange) > 1))
wr(js, "joint_fdr_summary.tsv"); print(js, row.names = FALSE)
cat("\nDone:", format(Sys.time()), "\n")
sessionInfo()
RSCRIPT_EOF
Rscript --no-save --no-restore "${RS}"
echo "=== 10_deseq2_designs_v3 finished: $(date) ==="
