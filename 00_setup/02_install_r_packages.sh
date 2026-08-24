#!/usr/bin/env bash
#SBATCH --job-name=install_r_pkgs
#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/install_r_packages_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/install_r_packages_%j.err
# =============================================================================
# 02_install_r_packages.sh - Install R/Bioconductor packages to local library
# Submit: sbatch 02_install_r_packages.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

log_step "Installing R/Bioconductor packages"

# Load R module
load_r

# Create local R library
mkdir -p "${R_LIBS_USER}"
export R_LIBS_USER

log_info "R version: $(R --version | head -1)"
log_info "Library path: ${R_LIBS_USER}"

# ---------------------------------------------------------------------------
# Install packages via Rscript
# ---------------------------------------------------------------------------
Rscript --no-save --no-restore - <<'REOF'

# Set library path
lib_path <- Sys.getenv("R_LIBS_USER")
.libPaths(c(lib_path, .libPaths()))
cat("Library paths:\n")
print(.libPaths())

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = lib_path)
}
library(BiocManager)

# ---------------------------------------------------------------------------
# Bioconductor packages
# ---------------------------------------------------------------------------
bioc_packages <- c(
    "DESeq2",
    "WGCNA",
    "clusterProfiler",
    "ComplexHeatmap",
    "enrichplot",
    "org.At.tair.db",
    "GenomicFeatures",
    "rtracklayer",
    "tximport",
    "AnnotationDbi"
)

# ---------------------------------------------------------------------------
# CRAN packages
# ---------------------------------------------------------------------------
cran_packages <- c(
    "pheatmap",
    "ggplot2",
    "tidyverse",
    "RColorBrewer",
    "circlize",
    "ape",
    "ggtree",
    "ggrepel",
    "scales",
    "reshape2",
    "gridExtra",
    "VennDiagram",
    "UpSetR",
    "dendextend"
)

# ---------------------------------------------------------------------------
# Install function with error handling
# ---------------------------------------------------------------------------
install_pkg <- function(pkg, lib) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("[SKIP] %s already installed\n", pkg))
        return(invisible(NULL))
    }
    cat(sprintf("[INSTALL] %s ...\n", pkg))
    tryCatch({
        BiocManager::install(pkg, lib = lib, update = FALSE, ask = FALSE,
                             Ncpus = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
        cat(sprintf("[OK] %s installed successfully\n", pkg))
    }, error = function(e) {
        cat(sprintf("[FAIL] %s - %s\n", pkg, conditionMessage(e)))
    })
}

# Install all
for (pkg in c(bioc_packages, cran_packages)) {
    install_pkg(pkg, lib_path)
}

# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
cat("\n=== Verification ===\n")
all_pkgs <- c(bioc_packages, cran_packages)
installed <- sapply(all_pkgs, requireNamespace, quietly = TRUE)
status <- ifelse(installed, "OK", "MISSING")
result <- data.frame(Package = all_pkgs, Status = status, stringsAsFactors = FALSE)
print(result)

n_ok <- sum(installed)
n_total <- length(all_pkgs)
cat(sprintf("\nInstalled %d / %d packages\n", n_ok, n_total))

if (n_ok < n_total) {
    cat("WARNING: Some packages failed to install. Check errors above.\n")
    quit(status = 1)
}
cat("All packages installed successfully.\n")

REOF

log_done "R package installation complete"
