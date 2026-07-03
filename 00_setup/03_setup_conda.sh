#!/usr/bin/env bash
#SBATCH --job-name=setup_conda
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/setup_conda_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/setup_conda_%j.err
# =============================================================================
# 03_setup_conda.sh — Create conda environment with Python tools
# Submit: sbatch 03_setup_conda.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

log_step "Setting up conda environment: ${CONDA_ENV_NAME}"

# ---------------------------------------------------------------------------
# Load miniforge module
# ---------------------------------------------------------------------------
load_miniforge

# Initialize conda for this shell session
eval "$(conda shell.bash hook)"

# ---------------------------------------------------------------------------
# Remove existing environment if present (clean rebuild)
# ---------------------------------------------------------------------------
if conda env list | grep -q "${CONDA_ENV_PREFIX}"; then
    log_warn "Existing environment found at ${CONDA_ENV_PREFIX} — removing"
    conda env remove --prefix "${CONDA_ENV_PREFIX}" --yes
fi

mkdir -p "$(dirname "${CONDA_ENV_PREFIX}")"

# ---------------------------------------------------------------------------
# Create environment with core Python stack
# ---------------------------------------------------------------------------
log_step "Creating conda environment at ${CONDA_ENV_PREFIX}"

conda create --prefix "${CONDA_ENV_PREFIX}" --yes \
    -c conda-forge -c bioconda \
    python=3.10 \
    biopython \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scipy \
    scikit-learn \
    requests

# ---------------------------------------------------------------------------
# Activate and install additional bioinformatics tools
# ---------------------------------------------------------------------------
conda activate "${CONDA_ENV_PREFIX}"

log_step "Installing bioinformatics tools"

# MCScanX — collinearity / synteny detection
conda install --yes -c bioconda mcscanx

# JCVI (python library for MCScan-based synteny visualization)
conda install --yes -c bioconda jcvi

# AGAT — GFF/GTF manipulation toolkit
conda install --yes -c bioconda agat

# KaKs_Calculator 2.0 — synonymous/nonsynonymous substitution rates
conda install --yes -c bioconda kaks-calculator

# ParaAT — parallel alignment tool for Ka/Ks
conda install --yes -c bioconda paraat

# SRA toolkit — for downloading FASTQ from NCBI
conda install --yes -c bioconda sra-tools

# InterProScan helper (optional — large install, skip if HPC has module)
# conda install --yes -c bioconda interproscan

# ---------------------------------------------------------------------------
# pip-only packages (if not on conda)
# ---------------------------------------------------------------------------
pip install --no-cache-dir \
    pygenometracks \
    gffutils

# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
log_step "Verifying installed packages"

python -c "
import sys
packages = [
    'Bio',         # biopython
    'pandas',
    'numpy',
    'matplotlib',
    'seaborn',
    'scipy',
    'sklearn',
    'gffutils',
]
missing = []
for pkg in packages:
    try:
        __import__(pkg)
        print(f'  [OK] {pkg}')
    except ImportError:
        print(f'  [FAIL] {pkg}')
        missing.append(pkg)

if missing:
    print(f'\\nMissing packages: {missing}')
    sys.exit(1)
print('\\nAll Python packages verified.')
"

# Check command-line tools
for tool in MCScanX agat_convert_sp_gff2gtf.pl fasterq-dump; do
    if command -v "${tool}" &>/dev/null; then
        log_info "[OK] ${tool} found"
    else
        log_warn "[MISSING] ${tool} not found in PATH"
    fi
done

conda deactivate

log_done "Conda environment ready at: ${CONDA_ENV_PREFIX}"
log_info "Activate with: conda activate ${CONDA_ENV_PREFIX}"
