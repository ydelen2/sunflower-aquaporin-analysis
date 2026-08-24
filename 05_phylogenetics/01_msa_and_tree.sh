#!/bin/bash
#SBATCH --job-name=aqp_phylo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/05_msa_tree_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/05_msa_tree_%j.err
#SBATCH --mail-type=END,FAIL

# ============================================================================
# 01_msa_and_tree.sh
# Multiple sequence alignment and phylogenetic tree for aquaporin gene family
# Combines sunflower AQPs with reference AQPs from Arabidopsis, rice, tomato,
# lettuce. Aligns with MAFFT L-INS-i, trims with trimAl, builds ML tree
# with IQ-TREE (ModelFinder + ultrafast bootstrap + SH-aLRT).
# ============================================================================

set -eo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
source "${PROJ_DIR}/config.sh"

WORK_DIR="${PROJ_DIR}/05_phylogenetics"
RESULTS_DIR="${WORK_DIR}/results"
REF_DIR="${PROJ_DIR}/01_references"
GENE_FAM_DIR="${PROJ_DIR}/02_gene_family/results"
LOG_DIR="${PROJ_DIR}/logs"

mkdir -p "${RESULTS_DIR}" "${LOG_DIR}"

# Sunflower aquaporin protein sequences from gene family identification step
SUNFLOWER_AQP="${GENE_FAM_DIR}/aquaporin_proteins.fa"

# Reference aquaporin sequences (one FASTA per species, pre-curated)
# Expected: Arabidopsis (~35), rice (~33), tomato (~30), lettuce (~35)
AT_REF="${REF_DIR}/arabidopsis_aquaporins.fa"
OS_REF="${REF_DIR}/rice_aquaporins.fa"
SL_REF="${REF_DIR}/tomato_aquaporins.fa"
LS_REF="${REF_DIR}/lettuce_aquaporins.fa"

# ---------------------------------------------------------------------------
# Modules
# ---------------------------------------------------------------------------
module purge
module load mafft/7.526
module load iqtree/2.2
module load miniforge/24.5

echo "================================================================"
echo "Pipeline: MSA + Phylogenetic Tree"
echo "Started: $(date)"
echo "Job ID:  ${SLURM_JOB_ID}"
echo "Nodes:   ${SLURM_JOB_NODELIST}"
echo "CPUs:    ${SLURM_NTASKS_PER_NODE}"
echo "================================================================"

# ---------------------------------------------------------------------------
# Step 1: Combine all aquaporin sequences into a single FASTA
# ---------------------------------------------------------------------------
COMBINED="${RESULTS_DIR}/all_aquaporins_combined.fa"

echo "[$(date '+%H:%M:%S')] Combining aquaporin sequences..."

# Validate that all input files exist
for f in "${SUNFLOWER_AQP}" "${AT_REF}" "${OS_REF}" "${SL_REF}" "${LS_REF}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Required file not found: ${f}" >&2
        exit 1
    fi
    echo "  Found: ${f} ($(grep -c '^>' "${f}") sequences)"
done

cat "${SUNFLOWER_AQP}" "${AT_REF}" "${OS_REF}" "${SL_REF}" "${LS_REF}" > "${COMBINED}"

TOTAL_SEQ=$(grep -c '^>' "${COMBINED}")
echo "  Total combined sequences: ${TOTAL_SEQ}"

# Sanity check: remove duplicate headers if any
awk '/^>/{h=$0; if(seen[h]++){next}} {print}' "${COMBINED}" > "${COMBINED}.dedup"
DEDUP_SEQ=$(grep -c '^>' "${COMBINED}.dedup")
if [[ "${DEDUP_SEQ}" -ne "${TOTAL_SEQ}" ]]; then
    echo "  WARNING: Removed $((TOTAL_SEQ - DEDUP_SEQ)) duplicate sequences"
    mv "${COMBINED}.dedup" "${COMBINED}"
else
    rm "${COMBINED}.dedup"
fi

# ---------------------------------------------------------------------------
# Step 2: Multiple Sequence Alignment with MAFFT
# ---------------------------------------------------------------------------
MSA_RAW="${RESULTS_DIR}/all_aquaporins_mafft.fa"

echo "[$(date '+%H:%M:%S')] Running MAFFT L-INS-i alignment..."

# L-INS-i: most accurate for <200 sequences with one alignable domain.
# Falls back to FFT-NS-2 automatically if sequence count is too high.
mafft \
    --localpair \
    --maxiterate 1000 \
    --thread "${SLURM_NTASKS_PER_NODE}" \
    --reorder \
    "${COMBINED}" > "${MSA_RAW}" 2> "${LOG_DIR}/mafft_${SLURM_JOB_ID}.log"

echo "  Alignment completed: ${MSA_RAW}"
echo "  Alignment length: $(head -2 "${MSA_RAW}" | tail -1 | wc -c) characters"

# ---------------------------------------------------------------------------
# Step 3: Trim alignment with trimAl
# ---------------------------------------------------------------------------
MSA_TRIMMED="${RESULTS_DIR}/all_aquaporins_trimmed.fa"

echo "[$(date '+%H:%M:%S')] Trimming alignment with trimAl..."

# Install trimAl via conda if not available
TRIMAL_BIN=$(which trimal 2>/dev/null || echo "")
if [[ -z "${TRIMAL_BIN}" ]]; then
    CONDA_ENV="${PROJ_DIR}/envs/trimal_env"
    if [[ ! -d "${CONDA_ENV}" ]]; then
        echo "  Installing trimAl via conda..."
        conda create -y -p "${CONDA_ENV}" -c bioconda trimal 2>&1 | tail -3
    fi
    conda activate "${CONDA_ENV}"
fi

trimal \
    -in "${MSA_RAW}" \
    -out "${MSA_TRIMMED}" \
    -htmlout "${RESULTS_DIR}/trimming_report.html" \
    -automated1

TRIMMED_LEN=$(head -2 "${MSA_TRIMMED}" | tail -1 | tr -d '\n' | wc -c)
RAW_LEN=$(head -2 "${MSA_RAW}" | tail -1 | tr -d '\n' | wc -c)
echo "  Alignment trimmed: ${RAW_LEN} -> ${TRIMMED_LEN} positions"

# Deactivate conda env if it was activated
conda deactivate 2>/dev/null || true

# ---------------------------------------------------------------------------
# Step 4: Phylogenetic tree with IQ-TREE
# ---------------------------------------------------------------------------
TREE_PREFIX="${RESULTS_DIR}/aquaporin_tree"

echo "[$(date '+%H:%M:%S')] Running IQ-TREE phylogenetic analysis..."

# ModelFinder Plus: automatic model selection
# UFBoot2: ultrafast bootstrap 1000 replicates
# SH-aLRT: SH-like approximate likelihood ratio test 1000 replicates
# -bnni: reduce UFBoot overestimation via NNI optimization
iqtree2 \
    -s "${MSA_TRIMMED}" \
    --prefix "${TREE_PREFIX}" \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -bnni \
    -nt "${SLURM_NTASKS_PER_NODE}" \
    --seed 42 \
    2>&1 | tee "${LOG_DIR}/iqtree_${SLURM_JOB_ID}.log"

echo "[$(date '+%H:%M:%S')] IQ-TREE completed."
echo "  Best model: $(grep 'Best-fit model' "${TREE_PREFIX}.iqtree" | head -1)"
echo "  Tree file:  ${TREE_PREFIX}.treefile"
echo "  Consensus:  ${TREE_PREFIX}.contree"

# ---------------------------------------------------------------------------
# Step 5: Generate alignment statistics summary
# ---------------------------------------------------------------------------
STATS_FILE="${RESULTS_DIR}/alignment_stats.txt"

echo "[$(date '+%H:%M:%S')] Generating alignment statistics..."

cat > "${STATS_FILE}" <<STATS
Aquaporin Phylogenetic Analysis Summary
========================================
Date: $(date)
Job ID: ${SLURM_JOB_ID}

Input Sequences
---------------
Sunflower AQPs: $(grep -c '^>' "${SUNFLOWER_AQP}")
Arabidopsis:    $(grep -c '^>' "${AT_REF}")
Rice:           $(grep -c '^>' "${OS_REF}")
Tomato:         $(grep -c '^>' "${SL_REF}")
Lettuce:        $(grep -c '^>' "${LS_REF}")
Total:          ${TOTAL_SEQ}

Alignment
---------
Raw alignment length:     ${RAW_LEN} aa
Trimmed alignment length: ${TRIMMED_LEN} aa

IQ-TREE Results
---------------
$(grep 'Best-fit model' "${TREE_PREFIX}.iqtree" | head -1)
$(grep 'Log-likelihood of the tree' "${TREE_PREFIX}.iqtree" | head -1)
$(grep 'Total tree length' "${TREE_PREFIX}.iqtree" | head -1)

Output Files
------------
Combined FASTA:    ${COMBINED}
Raw alignment:     ${MSA_RAW}
Trimmed alignment: ${MSA_TRIMMED}
Tree file:         ${TREE_PREFIX}.treefile
Consensus tree:    ${TREE_PREFIX}.contree
IQ-TREE log:       ${TREE_PREFIX}.iqtree
STATS

echo "================================================================"
echo "Pipeline completed: $(date)"
echo "================================================================"
