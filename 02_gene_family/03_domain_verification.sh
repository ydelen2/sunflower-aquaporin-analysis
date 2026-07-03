#!/bin/bash
#SBATCH --job-name=aqp_domain_ver
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/03_domain_verification_%j.out
#SBATCH --error=logs/03_domain_verification_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 03_domain_verification.sh
# Domain verification of aquaporin candidates:
#   1. hmmscan against Pfam-A for MIP domain (PF00230) confirmation
#   2. Transmembrane helix prediction (aquaporins should have 6 TM helices)
#   3. Filtering: MIP domain required AND ≥4 TM domains
###############################################################################

set -euo pipefail

CONFIG="/work/dweikat/ydelen2/aquaporin_study/config.sh"
if [[ ! -f "${CONFIG}" ]]; then
    echo "FATAL: config file not found: ${CONFIG}" >&2
    exit 1
fi
source "${CONFIG}"

module purge
module load hmmer/3.4
module load miniforge/24.5

# ── Directories ──────────────────────────────────────────────────────────────
WORK_DIR="${PROJ_DIR}/02_gene_family/03_domain_verification"
PARENT_DIR="${PROJ_DIR}/02_gene_family"
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${WORK_DIR}" "${LOG_DIR}"

LOGFILE="${LOG_DIR}/domain_verification_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "================================================================="
echo "  Domain Verification of Aquaporin Candidates"
echo "  Started: $(date)"
echo "  Job ID:  ${SLURM_JOB_ID:-local}"
echo "================================================================="

# ── Input validation ─────────────────────────────────────────────────────────
CANDIDATE_FASTA="${PARENT_DIR}/candidate_aquaporins.faa"
CANDIDATE_IDS="${PARENT_DIR}/candidate_aquaporin_ids.txt"

for f in "${CANDIDATE_FASTA}" "${CANDIDATE_IDS}"; do
    if [[ ! -f "${f}" ]]; then
        echo "FATAL: Required input not found: ${f}" >&2
        echo "  Run 02_blast_validation.sh first." >&2
        exit 1
    fi
done

N_CANDIDATES=$(grep -c '^>' "${CANDIDATE_FASTA}")
echo "Input candidates: ${N_CANDIDATES}"

# ── Step 1: Pfam-A database check ───────────────────────────────────────────
echo ""
echo "[Step 1] Checking Pfam-A database..."

# Common Pfam-A database locations on HPC clusters
PFAM_DB=""
PFAM_SEARCH_PATHS=(
    "${REF_DIR}/pfam/Pfam-A.hmm"
    "/common/databases/pfam/Pfam-A.hmm"
    "/work/common/pfam/Pfam-A.hmm"
    "${PROJ_DIR}/databases/pfam/Pfam-A.hmm"
)

for dbpath in "${PFAM_SEARCH_PATHS[@]}"; do
    if [[ -f "${dbpath}" ]]; then
        PFAM_DB="${dbpath}"
        break
    fi
done

if [[ -z "${PFAM_DB}" ]]; then
    echo "  Pfam-A database not found in standard locations."
    echo "  Downloading Pfam-A.hmm (this is ~1.5 GB)..."

    PFAM_DB_DIR="${PROJ_DIR}/databases/pfam"
    mkdir -p "${PFAM_DB_DIR}"
    PFAM_DB="${PFAM_DB_DIR}/Pfam-A.hmm"

    if [[ ! -f "${PFAM_DB}" ]]; then
        wget -q -O "${PFAM_DB}.gz" \
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
        gunzip "${PFAM_DB}.gz"
        hmmpress "${PFAM_DB}"
    fi
fi

echo "  Pfam-A database: ${PFAM_DB}"

# Verify Pfam-A is pressed
if [[ ! -f "${PFAM_DB}.h3i" ]]; then
    echo "  Pressing Pfam-A database..."
    hmmpress "${PFAM_DB}"
fi

# ── Step 2: hmmscan against Pfam-A ──────────────────────────────────────────
echo ""
echo "[Step 2] Running hmmscan against Pfam-A database..."

HMMSCAN_OUT="${WORK_DIR}/hmmscan_pfam.out"
HMMSCAN_TBL="${WORK_DIR}/hmmscan_pfam.tbl"
HMMSCAN_DOM="${WORK_DIR}/hmmscan_pfam.domtbl"

hmmscan \
    --cpu "${SLURM_NTASKS_PER_NODE:-8}" \
    -E 1e-5 \
    --domE 1e-5 \
    --tblout "${HMMSCAN_TBL}" \
    --domtblout "${HMMSCAN_DOM}" \
    -o "${HMMSCAN_OUT}" \
    "${PFAM_DB}" \
    "${CANDIDATE_FASTA}"

echo "  hmmscan complete."

# ── Step 3: Parse MIP domain results ────────────────────────────────────────
echo ""
echo "[Step 3] Parsing MIP domain (PF00230) results..."

# Extract proteins with MIP domain hit
MIP_HITS="${WORK_DIR}/mip_domain_hits.tsv"
echo -e "protein_id\tdomain\taccession\tdom_evalue\tdom_score\tali_from\tali_to" \
    > "${MIP_HITS}"

grep -v '^#' "${HMMSCAN_DOM}" \
    | awk '$2 == "PF00230" || $2 ~ /^PF00230\./ {
        OFS="\t"
        print $4, $1, $2, $13, $14, $18, $19
    }' \
    >> "${MIP_HITS}"

MIP_PROTEIN_IDS="${WORK_DIR}/mip_positive_ids.txt"
grep -v '^#' "${HMMSCAN_DOM}" \
    | awk '$2 == "PF00230" || $2 ~ /^PF00230\./ {print $4}' \
    | sort -u \
    > "${MIP_PROTEIN_IDS}"

N_MIP=$(wc -l < "${MIP_PROTEIN_IDS}")
echo "  Proteins with MIP domain: ${N_MIP}"

# Also extract all domain annotations for candidates
ALL_DOMAINS="${WORK_DIR}/all_domain_annotations.tsv"
echo -e "protein_id\tdomain_name\taccession\tdom_evalue\tdom_score\tali_from\tali_to\tdescription" \
    > "${ALL_DOMAINS}"

grep -v '^#' "${HMMSCAN_DOM}" \
    | awk 'BEGIN{OFS="\t"} {
        # Reconstruct description from remaining fields
        desc = ""
        for (i=23; i<=NF; i++) desc = (desc == "" ? $i : desc " " $i)
        print $4, $1, $2, $13, $14, $18, $19, desc
    }' \
    >> "${ALL_DOMAINS}"

echo "  Full domain annotations: ${ALL_DOMAINS}"

# ── Step 4: Transmembrane helix prediction ───────────────────────────────────
echo ""
echo "[Step 4] Predicting transmembrane helices..."

# Use a Python-based TM prediction since TMHMM may not be available
# This uses a hydrophobicity-based sliding window approach as a proxy
# For publication, DeepTMHMM or TMHMM should be run separately

TM_SCRIPT="${WORK_DIR}/predict_tm.py"
cat > "${TM_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""
Transmembrane helix prediction using Kyte-Doolittle hydrophobicity.
This is a simple screening method. For publication, use DeepTMHMM or TMHMM.

Aquaporins are expected to have 6 TM helices.
"""
import sys
import os

# Kyte-Doolittle hydrophobicity scale
KD_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    'X': 0.0, 'U': 0.0, 'B': -3.5, 'Z': -3.5, 'J': 4.15
}

def parse_fasta(fasta_file):
    """Parse FASTA file, yield (header, sequence) tuples."""
    header = None
    seq_parts = []
    with open(fasta_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if header is not None:
        yield header, ''.join(seq_parts)

def predict_tm_helices(sequence, window=19, threshold=1.6, min_helix_len=15, merge_dist=5):
    """
    Predict TM helices using sliding window hydrophobicity.
    Returns list of (start, end) tuples (1-indexed).
    """
    if len(sequence) < window:
        return []

    # Calculate hydrophobicity profile
    scores = []
    for i in range(len(sequence) - window + 1):
        segment = sequence[i:i+window]
        score = sum(KD_SCALE.get(aa, 0) for aa in segment) / window
        scores.append(score)

    # Find regions above threshold
    in_tm = False
    tm_regions = []
    start = 0
    for i, score in enumerate(scores):
        if score >= threshold and not in_tm:
            start = i
            in_tm = True
        elif score < threshold and in_tm:
            # TM region ended
            end = i + window - 1
            region_len = end - start
            if region_len >= min_helix_len:
                tm_regions.append((start + 1, end))  # 1-indexed
            in_tm = False
    # Handle case where TM extends to end
    if in_tm:
        end = len(sequence)
        region_len = end - start
        if region_len >= min_helix_len:
            tm_regions.append((start + 1, end))

    # Merge close TM helices
    merged = []
    for region in tm_regions:
        if merged and region[0] - merged[-1][1] <= merge_dist:
            merged[-1] = (merged[-1][0], region[1])
        else:
            merged.append(region)

    return merged

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.faa> <output.tsv>", file=sys.stderr)
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_file = sys.argv[2]

    with open(output_file, 'w') as out:
        out.write("protein_id\tnum_tm_helices\ttm_positions\tsequence_length\n")

        for prot_id, sequence in parse_fasta(input_fasta):
            tm_helices = predict_tm_helices(sequence)
            n_tm = len(tm_helices)
            tm_pos = ';'.join(f"{s}-{e}" for s, e in tm_helices) if tm_helices else 'none'
            out.write(f"{prot_id}\t{n_tm}\t{tm_pos}\t{len(sequence)}\n")

    # Print summary to stderr
    counts = {}
    with open(output_file) as fh:
        next(fh)  # skip header
        for line in fh:
            n = int(line.strip().split('\t')[1])
            counts[n] = counts.get(n, 0) + 1

    print("TM helix count distribution:", file=sys.stderr)
    for n_tm in sorted(counts.keys()):
        print(f"  {n_tm} TM: {counts[n_tm]} proteins", file=sys.stderr)

if __name__ == '__main__':
    main()
PYEOF

TM_RESULTS="${WORK_DIR}/tm_helix_predictions.tsv"
python3 "${TM_SCRIPT}" "${CANDIDATE_FASTA}" "${TM_RESULTS}"

echo "  TM predictions: ${TM_RESULTS}"

# Extract IDs with ≥4 TM helices
TM_PASS_IDS="${WORK_DIR}/tm_pass_ids.txt"
awk -F'\t' 'NR > 1 && $2 >= 4 {print $1}' "${TM_RESULTS}" \
    | sort -u > "${TM_PASS_IDS}"
echo "  Proteins with ≥4 TM helices: $(wc -l < "${TM_PASS_IDS}")"

# ── Step 5: Apply combined filter ────────────────────────────────────────────
echo ""
echo "[Step 5] Applying combined filter: MIP domain AND ≥4 TM helices..."

VERIFIED_IDS="${WORK_DIR}/verified_aquaporin_ids.txt"
comm -12 <(sort "${MIP_PROTEIN_IDS}") <(sort "${TM_PASS_IDS}") \
    > "${VERIFIED_IDS}"

N_VERIFIED=$(wc -l < "${VERIFIED_IDS}")
echo "  Verified aquaporins: ${N_VERIFIED}"

# ── Step 6: Create verified output files ─────────────────────────────────────
echo ""
echo "[Step 6] Creating verified aquaporin output files..."

VERIFIED_FASTA="${WORK_DIR}/verified_aquaporins.faa"
awk 'BEGIN {
    while ((getline line < "'"${VERIFIED_IDS}"'") > 0) ids[line] = 1
}
/^>/ {
    id = substr($1, 2)
    found = (id in ids)
}
found' "${CANDIDATE_FASTA}" > "${VERIFIED_FASTA}"

# Create comprehensive verification table
VERIFIED_TABLE="${WORK_DIR}/verified_aquaporins.tsv"
echo -e "protein_id\tmip_domain\tnum_tm_helices\ttm_positions\tseq_length\tmip_evalue\tmip_score" \
    > "${VERIFIED_TABLE}"

while IFS= read -r pid; do
    # Get TM info
    tm_line=$(awk -F'\t' -v id="${pid}" '$1 == id {print $2"\t"$3"\t"$4}' "${TM_RESULTS}")
    n_tm=$(echo "${tm_line}" | cut -f1)
    tm_pos=$(echo "${tm_line}" | cut -f2)
    seq_len=$(echo "${tm_line}" | cut -f3)

    # Get MIP domain info (best hit)
    mip_info=$(awk -F'\t' -v id="${pid}" 'NR > 1 && $1 == id {print $4"\t"$5; exit}' "${MIP_HITS}")
    mip_eval=$(echo "${mip_info}" | cut -f1)
    mip_score=$(echo "${mip_info}" | cut -f2)

    echo -e "${pid}\tyes\t${n_tm}\t${tm_pos}\t${seq_len}\t${mip_eval:-NA}\t${mip_score:-NA}"
done < "${VERIFIED_IDS}" >> "${VERIFIED_TABLE}"

# Identify proteins that failed verification
FAILED_IDS="${WORK_DIR}/failed_verification.tsv"
echo -e "protein_id\thas_mip\tnum_tm\treason" > "${FAILED_IDS}"

while IFS= read -r pid; do
    if ! grep -qxF "${pid}" "${VERIFIED_IDS}"; then
        has_mip="no"
        grep -qxF "${pid}" "${MIP_PROTEIN_IDS}" 2>/dev/null && has_mip="yes"

        n_tm=0
        tm_val=$(awk -F'\t' -v id="${pid}" '$1 == id {print $2}' "${TM_RESULTS}")
        [[ -n "${tm_val}" ]] && n_tm="${tm_val}"

        reason=""
        [[ "${has_mip}" == "no" ]] && reason="no_MIP_domain"
        [[ "${n_tm}" -lt 4 ]] && reason="${reason:+${reason};}insufficient_TM_helices(${n_tm})"

        echo -e "${pid}\t${has_mip}\t${n_tm}\t${reason}"
    fi
done < "${CANDIDATE_IDS}" >> "${FAILED_IDS}"

N_FAILED=$(( $(wc -l < "${FAILED_IDS}") - 1 ))
echo "  Failed verification: ${N_FAILED}"

# ── Step 7: Copy final outputs to parent directory ───────────────────────────
echo ""
echo "[Step 7] Copying final outputs..."

cp "${VERIFIED_IDS}" "${PARENT_DIR}/verified_aquaporin_ids.txt"
cp "${VERIFIED_FASTA}" "${PARENT_DIR}/verified_aquaporins.faa"
cp "${VERIFIED_TABLE}" "${PARENT_DIR}/verified_aquaporins.tsv"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "================================================================="
echo "  Domain Verification Summary"
echo "================================================================="
echo "  Input candidates:         ${N_CANDIDATES}"
echo "  MIP domain positive:      ${N_MIP}"
echo "  ≥4 TM helices:            $(wc -l < "${TM_PASS_IDS}")"
echo "  Verified (MIP + TM):      ${N_VERIFIED}"
echo "  Failed verification:      ${N_FAILED}"
echo ""
echo "  NOTE: TM prediction here uses Kyte-Doolittle hydrophobicity."
echo "  For publication, run DeepTMHMM (https://dtu.biolib.com/DeepTMHMM)"
echo "  or TMHMM-2.0 on verified_aquaporins.faa for accurate TM topology."
echo ""
echo "  Output files:"
echo "    ${VERIFIED_TABLE}"
echo "    ${VERIFIED_FASTA}"
echo "    ${VERIFIED_IDS}"
echo "    ${FAILED_IDS}"
echo ""
echo "  Completed: $(date)"
echo "================================================================="
