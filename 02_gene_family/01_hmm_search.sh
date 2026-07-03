#!/bin/bash
#SBATCH --job-name=aqp_hmmsearch
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/01_hmm_search_%j.out
#SBATCH --error=logs/01_hmm_search_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 01_hmm_search.sh
# HMM-based identification of aquaporin (MIP family) genes in sunflower
# Uses PF00230 HMM profile from Pfam/InterPro
# Reference: HanXRQr2.0 (GCF_002127325.2)
###############################################################################

set -euo pipefail

# ── Source project configuration ─────────────────────────────────────────────
CONFIG="/work/dweikat/ydelen2/aquaporin_study/config.sh"
if [[ ! -f "${CONFIG}" ]]; then
    echo "FATAL: config file not found: ${CONFIG}" >&2
    exit 1
fi
source "${CONFIG}"

# ── Module loads ─────────────────────────────────────────────────────────────
module purge
module load hmmer/3.4
module load samtools/1.23

# ── Directory setup ──────────────────────────────────────────────────────────
WORK_DIR="${PROJ_DIR}/02_gene_family/01_hmm_search"
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${WORK_DIR}" "${LOG_DIR}"

LOGFILE="${LOG_DIR}/hmm_search_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "================================================================="
echo "  HMM Search for Aquaporins (PF00230 / MIP)"
echo "  Started: $(date)"
echo "  Job ID:  ${SLURM_JOB_ID:-local}"
echo "================================================================="

# ── Input validation ─────────────────────────────────────────────────────────
PROTEOME="${REF_DIR}/GCF_002127325.2_HanXRQr2.0_protein.faa"
if [[ ! -f "${PROTEOME}" ]]; then
    echo "FATAL: Proteome file not found: ${PROTEOME}" >&2
    echo "  Download from NCBI: datasets download genome accession GCF_002127325.2 --include protein" >&2
    exit 1
fi

echo "Proteome: ${PROTEOME}"
echo "Protein count: $(grep -c '^>' "${PROTEOME}")"

# ── Step 1: Obtain PF00230 HMM profile ──────────────────────────────────────
HMM_FILE="${WORK_DIR}/PF00230.hmm"

if [[ ! -f "${HMM_FILE}" ]]; then
    echo ""
    echo "[Step 1] Downloading PF00230 HMM profile from InterPro/Pfam..."

    # Try InterPro endpoint first (Pfam now hosted at InterPro)
    wget -q -O "${HMM_FILE}" \
        "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00230?annotation=hmm" \
        2>/dev/null || true

    # Fallback: Pfam legacy
    if [[ ! -s "${HMM_FILE}" ]]; then
        echo "  InterPro download failed, trying Pfam legacy endpoint..."
        wget -q -O "${HMM_FILE}.gz" \
            "https://pfam.xfam.org/family/PF00230/hmm" \
            2>/dev/null || true
        if [[ -f "${HMM_FILE}.gz" ]]; then
            gunzip -f "${HMM_FILE}.gz" 2>/dev/null || mv "${HMM_FILE}.gz" "${HMM_FILE}"
        fi
    fi

    # Validate HMM file
    if [[ ! -s "${HMM_FILE}" ]]; then
        echo "FATAL: Failed to download PF00230 HMM profile." >&2
        echo "  Manual download:" >&2
        echo "    Go to https://www.ebi.ac.uk/interpro/entry/pfam/PF00230/" >&2
        echo "    Download the HMM file and place at: ${HMM_FILE}" >&2
        exit 1
    fi

    # Press the HMM (prepare binary format for hmmsearch)
    hmmpress "${HMM_FILE}"
    echo "  HMM profile ready."
else
    echo "[Step 1] PF00230 HMM profile already exists: ${HMM_FILE}"
fi

# ── Step 2: Run hmmsearch ────────────────────────────────────────────────────
echo ""
echo "[Step 2] Running hmmsearch against sunflower proteome..."
echo "  E-value threshold: 1e-5"
echo "  CPUs: ${SLURM_NTASKS_PER_NODE:-8}"

HMMSEARCH_OUT="${WORK_DIR}/hmmsearch_results.out"
HMMSEARCH_TBL="${WORK_DIR}/hmmsearch_results.tbl"
HMMSEARCH_DOM="${WORK_DIR}/hmmsearch_results.domtbl"

hmmsearch \
    --cpu "${SLURM_NTASKS_PER_NODE:-8}" \
    -E 1e-5 \
    --domE 1e-5 \
    --tblout "${HMMSEARCH_TBL}" \
    --domtblout "${HMMSEARCH_DOM}" \
    -o "${HMMSEARCH_OUT}" \
    "${HMM_FILE}" \
    "${PROTEOME}"

echo "  hmmsearch complete."

# ── Step 3: Parse results ────────────────────────────────────────────────────
echo ""
echo "[Step 3] Parsing hmmsearch results..."

# Extract unique protein IDs from table output (skip comment lines)
HIT_IDS="${WORK_DIR}/hmm_hit_ids.txt"
grep -v '^#' "${HMMSEARCH_TBL}" \
    | awk '{print $1}' \
    | sort -u \
    > "${HIT_IDS}"

N_HITS=$(wc -l < "${HIT_IDS}")
echo "  Unique protein hits: ${N_HITS}"

if [[ "${N_HITS}" -eq 0 ]]; then
    echo "WARNING: No hits found. Check proteome file format and HMM profile." >&2
    exit 1
fi

# Create detailed hit table: protein_id, E-value, score, bias, best_domain_E
HIT_TABLE="${WORK_DIR}/hmm_hit_table.tsv"
echo -e "protein_id\tE_value\tscore\tbias\tbest_dom_evalue\tbest_dom_score" > "${HIT_TABLE}"
grep -v '^#' "${HMMSEARCH_TBL}" \
    | awk 'BEGIN{OFS="\t"} {print $1, $5, $6, $7, $8, $9}' \
    | sort -k3,3nr \
    >> "${HIT_TABLE}"

echo "  Hit table: ${HIT_TABLE}"

# ── Step 4: Extract hit sequences ───────────────────────────────────────────
echo ""
echo "[Step 4] Extracting hit protein sequences..."

HIT_FASTA="${WORK_DIR}/hmm_hit_sequences.faa"

# Index proteome if needed
if [[ ! -f "${PROTEOME}.fai" ]]; then
    samtools faidx "${PROTEOME}"
fi

# Extract sequences using samtools faidx
# Handle NCBI-style headers: samtools faidx uses the first word before space
samtools faidx "${PROTEOME}" $(cat "${HIT_IDS}") > "${HIT_FASTA}" 2>/dev/null || {
    # Fallback: use awk-based extraction for complex headers
    echo "  samtools faidx extraction had issues, using awk fallback..."
    awk 'BEGIN {
        while ((getline line < "'"${HIT_IDS}"'") > 0) ids[line] = 1
    }
    /^>/ {
        id = substr($1, 2)
        found = (id in ids)
    }
    found' "${PROTEOME}" > "${HIT_FASTA}"
}

N_EXTRACTED=$(grep -c '^>' "${HIT_FASTA}" || echo 0)
echo "  Extracted sequences: ${N_EXTRACTED}"

# ── Step 5: Parse domain architecture from domtblout ─────────────────────────
echo ""
echo "[Step 5] Parsing domain architecture..."

DOMAIN_ARCH="${WORK_DIR}/hmm_domain_architecture.tsv"
echo -e "protein_id\tdomain_num\thmm_from\thmm_to\tali_from\tali_to\tenv_from\tenv_to\tdom_evalue\tdom_score" \
    > "${DOMAIN_ARCH}"
grep -v '^#' "${HMMSEARCH_DOM}" \
    | awk 'BEGIN{OFS="\t"} {print $1, $10, $16, $17, $18, $19, $20, $21, $13, $14}' \
    >> "${DOMAIN_ARCH}"

echo "  Domain architecture: ${DOMAIN_ARCH}"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "================================================================="
echo "  HMM Search Summary"
echo "================================================================="
echo "  HMM profile:           PF00230 (MIP)"
echo "  Proteome:              $(basename "${PROTEOME}")"
echo "  E-value threshold:     1e-5"
echo "  Total protein hits:    ${N_HITS}"
echo "  Sequences extracted:   ${N_EXTRACTED}"
echo ""
echo "  Output files:"
echo "    ${HMMSEARCH_TBL}"
echo "    ${HMMSEARCH_DOM}"
echo "    ${HIT_IDS}"
echo "    ${HIT_TABLE}"
echo "    ${HIT_FASTA}"
echo "    ${DOMAIN_ARCH}"
echo ""
echo "  Completed: $(date)"
echo "================================================================="
