#!/bin/bash
#SBATCH --job-name=aqp_phylo_v2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/05_msa_tree_v2_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/05_msa_tree_v2_%j.err
#SBATCH --mail-type=END,FAIL

# ============================================================================
# 01_msa_and_tree_v2.sh
# Multiple sequence alignment and phylogenetic tree for aquaporin gene family
#
# Difference from the first version: reference aquaporins are now identified
# by scanning each reference proteome with the PF00230 (MIP) HMM profile,
# instead of downloading a fixed list of UniProt accessions. Many of those
# accessions are no longer served by UniProt, which left the rice set with
# only 14 of 33 sequences and the Arabidopsis set with 32 of 35.
#
# Alignment (MAFFT L-INS-i), trimming (trimAl -automated1) and tree inference
# (IQ-TREE, ModelFinder + UFBoot + SH-aLRT) are unchanged.
# ============================================================================

set -eo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
source "${PROJ_DIR}/config.sh"

WORK_DIR="${PROJ_DIR}/05_phylogenetics"
RESULTS_DIR="${WORK_DIR}/results_v2"
REF_DIR="${PROJ_DIR}/01_references"
GENE_FAM_DIR="${PROJ_DIR}/02_gene_family/results"
LOG_DIR="${PROJ_DIR}/logs"

mkdir -p "${RESULTS_DIR}" "${LOG_DIR}"

# Sunflower aquaporin protein sequences from gene family identification step
SUNFLOWER_AQP="${GENE_FAM_DIR}/aquaporin_proteins.fa"

# PF00230 profile downloaded by 02_gene_family/01_hmm_search.sh
HMM_FILE="${PROJ_DIR}/02_gene_family/01_hmm_search/PF00230.hmm"

# HMM E-value threshold, same as the sunflower search in 02_gene_family
HMM_EVALUE="1e-5"

# ---------------------------------------------------------------------------
# Modules
# The pinned versions are the ones used for the original analysis. If the
# cluster has retired a version, fall back to the unversioned module so the
# job does not die at submission time.
# ---------------------------------------------------------------------------
module purge

# Try each candidate module name in turn and stop at the first that loads.
try_modules() {
    local m
    for m in "$@"; do
        if module load "${m}" 2>/dev/null; then
            echo "  loaded module: ${m}"
            return 0
        fi
    done
    echo "  note: no module loaded from: $*"
    return 0
}

try_modules "hmmer/3.4"      "hmmer"      "HMMER/3.4"   "HMMER"
try_modules "mafft/7.526"    "mafft"      "MAFFT/7.526" "MAFFT"
try_modules "iqtree/2.2"     "iqtree"     "iqtree/2"    "IQ-TREE/2.2" "IQTREE"
try_modules "miniforge/24.5" "miniforge"  "anaconda"    "miniconda"

# If a tool is not provided as a module on this cluster, install it into a
# dedicated conda environment under the project directory. Same approach the
# script already uses for trimAl.
ensure_tool() {
    local binary="$1" pkg="$2" envname="$3"

    command -v "${binary}" &>/dev/null && return 0

    local env_path="${PROJ_DIR}/envs/${envname}"
    if [[ -x "${env_path}/bin/${binary}" ]]; then
        export PATH="${env_path}/bin:${PATH}"
        echo "  ${binary}: using existing conda env ${env_path}"
        return 0
    fi

    if ! command -v conda &>/dev/null; then
        echo "ERROR: ${binary} is not available and conda was not found either." >&2
        echo "       Load a conda module first, or install ${pkg} manually." >&2
        return 1
    fi

    echo "  ${binary} not on PATH; installing ${pkg} into ${env_path}..."
    mkdir -p "${PROJ_DIR}/envs"
    conda create -y -p "${env_path}" -c conda-forge -c bioconda "${pkg}" 2>&1 | tail -5
    export PATH="${env_path}/bin:${PATH}"

    command -v "${binary}" &>/dev/null
}

ensure_tool hmmsearch hmmer  hmmer_env  || exit 1
ensure_tool mafft     mafft  mafft_env  || exit 1

echo "================================================================"
echo "Pipeline: MSA + Phylogenetic Tree (v2, HMM-based reference sets)"
echo "Started: $(date)"
echo "Job ID:  ${SLURM_JOB_ID}"
echo "CPUs:    ${SLURM_NTASKS_PER_NODE}"
echo "================================================================"

# IQ-TREE ships as iqtree2 or iqtree depending on the build
if ! command -v iqtree2 &>/dev/null && ! command -v iqtree &>/dev/null; then
    ensure_tool iqtree2 iqtree iqtree_env || ensure_tool iqtree iqtree iqtree_env || exit 1
fi

# Resolve the executables we need
HMMSEARCH_BIN=$(command -v hmmsearch || true)
MAFFT_BIN=$(command -v mafft || true)
IQTREE_BIN=$(command -v iqtree2 || command -v iqtree || true)

for pair in "hmmsearch:${HMMSEARCH_BIN}" "mafft:${MAFFT_BIN}" "iqtree:${IQTREE_BIN}"; do
    name="${pair%%:*}"
    path="${pair#*:}"
    if [[ -z "${path}" ]]; then
        echo "ERROR: ${name} could not be provided by a module or by conda." >&2
        echo "       Check what this cluster offers with: module spider ${name}" >&2
        exit 1
    fi
    echo "  ${name}: ${path}"
done

# ---------------------------------------------------------------------------
# Step 0: Resolve reference proteomes and check inputs
# The pipeline stores the NCBI downloads as <species>/protein.faa; the names
# in config.sh are kept as a fallback.
# ---------------------------------------------------------------------------
resolve_proteome() {
    local label="$1"
    shift
    local candidate
    for candidate in "$@"; do
        if [[ -s "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done
    echo "ERROR: no proteome found for ${label}. Tried: $*" >&2
    return 1
}

ATHA_FA=$(resolve_proteome "Arabidopsis" "${REF_DIR}/arabidopsis/protein.faa" "${ATHA_PROT:-}") || exit 1
OSAT_FA=$(resolve_proteome "rice"        "${REF_DIR}/rice/protein.faa"        "${OSAT_PROT:-}") || exit 1
SLYC_FA=$(resolve_proteome "tomato"      "${REF_DIR}/tomato/protein.faa"      "${SLYC_PROT:-}") || exit 1
LSAT_FA=$(resolve_proteome "lettuce"     "${REF_DIR}/lettuce/protein.faa"     "${LSAT_PROT:-}") || exit 1

# PF00230 profile: download it if the gene family step did not leave one behind
if [[ ! -s "${HMM_FILE}" ]]; then
    echo "  PF00230 profile not found at ${HMM_FILE}, downloading..."
    mkdir -p "$(dirname "${HMM_FILE}")"
    wget -q -O "${HMM_FILE}" \
        "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00230?annotation=hmm" || true
    if file "${HMM_FILE}" 2>/dev/null | grep -q gzip; then
        mv "${HMM_FILE}" "${HMM_FILE}.gz"
        gunzip -f "${HMM_FILE}.gz"
    fi
fi

for f in "${SUNFLOWER_AQP}" "${HMM_FILE}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: required input missing: ${f}" >&2
        exit 1
    fi
done

echo "[$(date '+%H:%M:%S')] Inputs present."
echo "  Sunflower aquaporins: $(grep -c '^>' "${SUNFLOWER_AQP}")"
echo "  Arabidopsis proteome: ${ATHA_FA}"
echo "  Rice proteome       : ${OSAT_FA}"
echo "  Tomato proteome     : ${SLYC_FA}"
echo "  Lettuce proteome    : ${LSAT_FA}"

# ---------------------------------------------------------------------------
# Step 1: Identify aquaporins in each reference proteome with PF00230
# ---------------------------------------------------------------------------
COMBINED="${RESULTS_DIR}/all_aquaporins_combined.fa"
COUNTS_FILE="${RESULTS_DIR}/sequence_counts.tsv"

echo "[$(date '+%H:%M:%S')] Scanning reference proteomes with PF00230..."

# Sunflower sequences go in first, tagged Ha
awk '/^>/ {print ">Ha_" substr($1,2); next} {print}' "${SUNFLOWER_AQP}" > "${COMBINED}"

printf "species\ttag\tproteome\tsequences\n" > "${COUNTS_FILE}"
printf "Helianthus annuus\tHa\t%s\t%s\n" \
    "$(basename "${SUNFLOWER_AQP}")" "$(grep -c '^>' "${SUNFLOWER_AQP}")" >> "${COUNTS_FILE}"

scan_proteome() {
    local tag="$1"
    local proteome="$2"
    local label="$3"

    local domtbl="${RESULTS_DIR}/${tag}_PF00230.domtblout"
    local ids="${RESULTS_DIR}/${tag}_hits.txt"
    local faa="${RESULTS_DIR}/${tag}_aquaporins.fa"

    echo "  ${label}: scanning $(grep -c '^>' "${proteome}") proteins..."

    "${HMMSEARCH_BIN}" \
        --domtblout "${domtbl}" \
        -E "${HMM_EVALUE}" \
        --cpu "${SLURM_NTASKS_PER_NODE}" \
        "${HMM_FILE}" \
        "${proteome}" \
        > /dev/null

    grep -v '^#' "${domtbl}" | awk '{print $1}' | sort -u > "${ids}" || true

    # Pull the matching sequences out of the proteome
    python3 - "${proteome}" "${ids}" "${faa}" <<'PYCODE'
import sys

proteome, id_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3]

with open(id_file) as fh:
    wanted = {line.strip() for line in fh if line.strip()}

kept = {}
name = None
buf = []

def flush():
    if name is not None:
        acc = name.split()[0]
        if acc in wanted:
            kept[acc] = ''.join(buf)

with open(proteome) as fh:
    for line in fh:
        if line.startswith('>'):
            flush()
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    flush()

with open(out_file, 'w') as out:
    for acc, seq in kept.items():
        out.write('>%s\n' % acc)
        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60] + '\n')

print('    kept %d sequences' % len(kept))
PYCODE

    awk -v t="${tag}" '/^>/ {print ">" t "_" substr($1,2); next} {print}' "${faa}" >> "${COMBINED}"

    local n
    n=$(grep -c '^>' "${faa}" || true)
    printf "%s\t%s\t%s\t%s\n" "${label}" "${tag}" "$(basename "${proteome}")" "${n}" >> "${COUNTS_FILE}"
    echo "  ${label}: ${n} aquaporin sequences"
}

scan_proteome "At" "${ATHA_FA}" "Arabidopsis thaliana"
scan_proteome "Os" "${OSAT_FA}" "Oryza sativa"
scan_proteome "Sl" "${SLYC_FA}" "Solanum lycopersicum"
scan_proteome "Ls" "${LSAT_FA}" "Lactuca sativa"

TOTAL_SEQ=$(grep -c '^>' "${COMBINED}" || true)
echo "  Total combined sequences: ${TOTAL_SEQ}"

# Drop duplicate headers if any slipped through
awk '/^>/{h=$0; if(seen[h]++){skip=1; next} skip=0} skip!=1 {print}' "${COMBINED}" > "${COMBINED}.dedup"
DEDUP_SEQ=$(grep -c '^>' "${COMBINED}.dedup" || true)
if [[ "${DEDUP_SEQ}" -ne "${TOTAL_SEQ}" ]]; then
    echo "  WARNING: removed $((TOTAL_SEQ - DEDUP_SEQ)) duplicate headers"
    mv "${COMBINED}.dedup" "${COMBINED}"
    TOTAL_SEQ="${DEDUP_SEQ}"
else
    rm "${COMBINED}.dedup"
fi

echo
echo "  Sequence counts per species:"
column -t -s $'\t' "${COUNTS_FILE}" | sed 's/^/    /'
echo

# ---------------------------------------------------------------------------
# Step 2: Multiple Sequence Alignment with MAFFT
# ---------------------------------------------------------------------------
MSA_RAW="${RESULTS_DIR}/all_aquaporins_mafft.fa"

echo "[$(date '+%H:%M:%S')] Running MAFFT L-INS-i alignment..."

"${MAFFT_BIN}" \
    --localpair \
    --maxiterate 1000 \
    --thread "${SLURM_NTASKS_PER_NODE}" \
    --reorder \
    "${COMBINED}" > "${MSA_RAW}" 2> "${LOG_DIR}/mafft_v2_${SLURM_JOB_ID}.log"

echo "  Alignment completed: ${MSA_RAW}"

# ---------------------------------------------------------------------------
# Step 3: Trim alignment with trimAl
# ---------------------------------------------------------------------------
MSA_TRIMMED="${RESULTS_DIR}/all_aquaporins_trimmed.fa"

echo "[$(date '+%H:%M:%S')] Trimming alignment with trimAl..."

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

conda deactivate 2>/dev/null || true

# ---------------------------------------------------------------------------
# Step 4: Phylogenetic tree with IQ-TREE
# ---------------------------------------------------------------------------
TREE_PREFIX="${RESULTS_DIR}/aquaporin_tree"

echo "[$(date '+%H:%M:%S')] Running IQ-TREE phylogenetic analysis..."

"${IQTREE_BIN}" \
    -s "${MSA_TRIMMED}" \
    --prefix "${TREE_PREFIX}" \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -bnni \
    -nt "${SLURM_NTASKS_PER_NODE}" \
    --seed 42 \
    2>&1 | tee "${LOG_DIR}/iqtree_v2_${SLURM_JOB_ID}.log"

echo "[$(date '+%H:%M:%S')] IQ-TREE completed."

# ---------------------------------------------------------------------------
# Step 5: Summary for the manuscript
# ---------------------------------------------------------------------------
STATS_FILE="${RESULTS_DIR}/alignment_stats.txt"

ALN_SEQS=$(grep -m1 'Input data:' "${TREE_PREFIX}.iqtree" || true)
BEST_MODEL=$(grep -m1 'Best-fit model' "${TREE_PREFIX}.iqtree" || true)

cat > "${STATS_FILE}" <<STATS
Aquaporin Phylogenetic Analysis Summary (v2, HMM-based reference sets)
======================================================================
Date: $(date)
Job ID: ${SLURM_JOB_ID}

Input Sequences
---------------
$(column -t -s $'\t' "${COUNTS_FILE}")
Total: ${TOTAL_SEQ}

Alignment
---------
${ALN_SEQS}

IQ-TREE Results
---------------
${BEST_MODEL}
$(grep -m1 'Log-likelihood of the tree' "${TREE_PREFIX}.iqtree" || true)
$(grep -m1 'Total tree length' "${TREE_PREFIX}.iqtree" || true)

Output Files
------------
Combined FASTA:    ${COMBINED}
Raw alignment:     ${MSA_RAW}
Trimmed alignment: ${MSA_TRIMMED}
Tree file:         ${TREE_PREFIX}.treefile
Consensus tree:    ${TREE_PREFIX}.contree
IQ-TREE report:    ${TREE_PREFIX}.iqtree
Per-species counts: ${COUNTS_FILE}
STATS

echo
echo "================================================================"
cat "${STATS_FILE}"
echo "================================================================"
echo "Pipeline completed: $(date)"
