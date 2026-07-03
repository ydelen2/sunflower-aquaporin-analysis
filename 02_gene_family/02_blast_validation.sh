#!/bin/bash
#SBATCH --job-name=aqp_blast_val
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/02_blast_validation_%j.out
#SBATCH --error=logs/02_blast_validation_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 02_blast_validation.sh
# BLASTP validation of aquaporin candidates using known Arabidopsis (35) and
# rice (33) aquaporins as queries against sunflower proteome.
# Merges results with HMM hits, removes redundancy.
###############################################################################

set -euo pipefail

CONFIG="/work/dweikat/ydelen2/aquaporin_study/config.sh"
if [[ ! -f "${CONFIG}" ]]; then
    echo "FATAL: config file not found: ${CONFIG}" >&2
    exit 1
fi
source "${CONFIG}"

module purge
module load blast/2.17

# ── Directories ──────────────────────────────────────────────────────────────
WORK_DIR="${PROJ_DIR}/02_gene_family/02_blast_validation"
HMM_DIR="${PROJ_DIR}/02_gene_family/01_hmm_search"
LOG_DIR="${WORK_DIR}/logs"
QUERY_DIR="${WORK_DIR}/query_sequences"
mkdir -p "${WORK_DIR}" "${LOG_DIR}" "${QUERY_DIR}"

LOGFILE="${LOG_DIR}/blast_validation_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "================================================================="
echo "  BLASTP Validation of Aquaporin Candidates"
echo "  Started: $(date)"
echo "  Job ID:  ${SLURM_JOB_ID:-local}"
echo "================================================================="

# ── Input validation ─────────────────────────────────────────────────────────
PROTEOME="${REF_DIR}/GCF_002127325.2_HanXRQr2.0_protein.faa"
HMM_HITS="${HMM_DIR}/hmm_hit_ids.txt"

for f in "${PROTEOME}" "${HMM_HITS}"; do
    if [[ ! -f "${f}" ]]; then
        echo "FATAL: Required input not found: ${f}" >&2
        exit 1
    fi
done

echo "Proteome: ${PROTEOME}"
echo "HMM hits: $(wc -l < "${HMM_HITS}") proteins"

# ── Step 1: Prepare known aquaporin query sequences ──────────────────────────
echo ""
echo "[Step 1] Preparing Arabidopsis and rice aquaporin query sequences..."

# Create Arabidopsis aquaporin accessions file
# 35 genes: 13 PIP (PIP1;1-PIP1;5, PIP2;1-PIP2;8), 10 TIP (TIP1;1-TIP1;3,
# TIP2;1-TIP2;3, TIP3;1-TIP3;2, TIP4;1, TIP5;1), 9 NIP (NIP1;1-NIP1;2,
# NIP2;1, NIP3;1, NIP4;1-NIP4;2, NIP5;1, NIP6;1, NIP7;1), 3 SIP (SIP1;1-SIP1;2, SIP2;1)
AT_AQP_FILE="${QUERY_DIR}/arabidopsis_aquaporins.txt"
cat > "${AT_AQP_FILE}" << 'ATEOF'
# Arabidopsis thaliana aquaporins - TAIR/UniProt accessions
# Subfamily	Gene_name	UniProt_ID
PIP	AtPIP1;1	Q06611
PIP	AtPIP1;2	P61837
PIP	AtPIP1;3	Q8LFP7
PIP	AtPIP1;4	Q39196
PIP	AtPIP1;5	Q8VZ49
PIP	AtPIP2;1	P43286
PIP	AtPIP2;2	P30302
PIP	AtPIP2;3	P43287
PIP	AtPIP2;4	Q9FF53
PIP	AtPIP2;5	Q9SV31
PIP	AtPIP2;6	Q9ZV07
PIP	AtPIP2;7	O22768
PIP	AtPIP2;8	Q9FIZ4
TIP	AtTIP1;1	P25818
TIP	AtTIP1;2	Q41963
TIP	AtTIP1;3	Q9C9Z1
TIP	AtTIP2;1	Q41951
TIP	AtTIP2;2	Q38857
TIP	AtTIP2;3	Q9FGL2
TIP	AtTIP3;1	Q08733
TIP	AtTIP3;2	Q9ZV07
TIP	AtTIP4;1	Q9LKJ5
TIP	AtTIP5;1	Q9C5A0
NIP	AtNIP1;1	Q8VWS5
NIP	AtNIP1;2	Q9FMN2
NIP	AtNIP2;1	Q9FPR0
NIP	AtNIP3;1	Q9SNS7
NIP	AtNIP4;1	Q9FKS6
NIP	AtNIP4;2	Q9FKS5
NIP	AtNIP5;1	Q9LNV3
NIP	AtNIP6;1	Q9LIK0
NIP	AtNIP7;1	Q8GZQ4
SIP	AtSIP1;1	Q9FN05
SIP	AtSIP1;2	Q9M8T1
SIP	AtSIP2;1	Q9LPV5
ATEOF

# Create rice aquaporin accessions file
# 33 genes in Oryza sativa
OS_AQP_FILE="${QUERY_DIR}/rice_aquaporins.txt"
cat > "${OS_AQP_FILE}" << 'OSEOF'
# Oryza sativa aquaporins
# Subfamily	Gene_name	UniProt_ID
PIP	OsPIP1;1	Q7XUA0
PIP	OsPIP1;2	Q6YZU4
PIP	OsPIP1;3	Q7XUA2
PIP	OsPIP2;1	Q65XA0
PIP	OsPIP2;2	Q6Z2T3
PIP	OsPIP2;3	Q7Y1W8
PIP	OsPIP2;4	Q6H4L2
PIP	OsPIP2;5	Q8H3Q1
PIP	OsPIP2;6	Q6YZR7
PIP	OsPIP2;7	Q6YZR6
PIP	OsPIP2;8	Q69NK3
TIP	OsTIP1;1	Q93YL5
TIP	OsTIP1;2	Q10MU8
TIP	OsTIP2;1	Q40735
TIP	OsTIP2;2	Q7XPJ3
TIP	OsTIP3;1	Q8GWC5
TIP	OsTIP3;2	Q10LA2
TIP	OsTIP4;1	Q6ERH5
TIP	OsTIP4;2	Q75GK6
TIP	OsTIP4;3	Q7F2Y3
TIP	OsTIP5;1	Q6L523
NIP	OsNIP1;1	Q7F730
NIP	OsNIP1;2	Q0DFC1
NIP	OsNIP1;3	Q75HM1
NIP	OsNIP1;4	Q7EZZ0
NIP	OsNIP2;1	Q6Z5J3
NIP	OsNIP2;2	Q6ZFJ9
NIP	OsNIP3;1	Q653R6
NIP	OsNIP3;2	Q7F1V5
NIP	OsNIP3;3	Q651J3
NIP	OsNIP4;1	Q7F2S9
SIP	OsSIP1;1	Q6K7N6
SIP	OsSIP2;1	Q6ZBL5
OSEOF

# Download protein sequences from UniProt
AT_QUERY="${QUERY_DIR}/arabidopsis_aqp.faa"
OS_QUERY="${QUERY_DIR}/rice_aqp.faa"

download_uniprot_sequences() {
    local accession_file="$1"
    local output_fasta="$2"
    local species_label="$3"

    if [[ -f "${output_fasta}" ]] && [[ $(grep -c '^>' "${output_fasta}") -gt 5 ]]; then
        echo "  ${species_label} sequences already downloaded: $(grep -c '^>' "${output_fasta}") proteins"
        return 0
    fi

    > "${output_fasta}"

    local count=0
    local fail_count=0
    while IFS=$'\t' read -r subfamily gene_name uniprot_id; do
        [[ "${subfamily}" =~ ^#.* ]] && continue
        [[ -z "${uniprot_id}" ]] && continue

        local tmp_seq="${QUERY_DIR}/tmp_${uniprot_id}.faa"
        if wget -q -O "${tmp_seq}" \
            "https://rest.uniprot.org/uniprotkb/${uniprot_id}.fasta" 2>/dev/null; then
            if [[ -s "${tmp_seq}" ]]; then
                # Rename header to include subfamily and gene name
                sed "1s/^>.*/>sp|${uniprot_id}|${gene_name} [${subfamily}]/" "${tmp_seq}" >> "${output_fasta}"
                count=$((count + 1))
            else
                fail_count=$((fail_count + 1))
            fi
        else
            fail_count=$((fail_count + 1))
        fi
        rm -f "${tmp_seq}"
        sleep 0.5  # Rate limiting
    done < "${accession_file}"

    echo "  ${species_label}: ${count} sequences downloaded, ${fail_count} failed"

    if [[ "${count}" -lt 10 ]]; then
        echo "WARNING: Low sequence count for ${species_label}. Check UniProt accessions." >&2
    fi
}

download_uniprot_sequences "${AT_AQP_FILE}" "${AT_QUERY}" "Arabidopsis"
download_uniprot_sequences "${OS_AQP_FILE}" "${OS_QUERY}" "Rice"

# Combine query sequences
COMBINED_QUERY="${QUERY_DIR}/combined_aqp_query.faa"
cat "${AT_QUERY}" "${OS_QUERY}" > "${COMBINED_QUERY}"
echo "  Combined query: $(grep -c '^>' "${COMBINED_QUERY}") sequences"

# ── Step 2: Create BLAST database ────────────────────────────────────────────
echo ""
echo "[Step 2] Creating BLAST database from sunflower proteome..."

BLASTDB_DIR="${WORK_DIR}/blastdb"
mkdir -p "${BLASTDB_DIR}"
BLASTDB="${BLASTDB_DIR}/HanXRQr2_proteome"

if [[ ! -f "${BLASTDB}.pdb" ]]; then
    makeblastdb \
        -in "${PROTEOME}" \
        -dbtype prot \
        -out "${BLASTDB}" \
        -title "HanXRQr2.0_proteome" \
        -parse_seqids
    echo "  BLAST database created."
else
    echo "  BLAST database already exists."
fi

# ── Step 3: Run BLASTP ───────────────────────────────────────────────────────
echo ""
echo "[Step 3] Running BLASTP searches..."

# BLASTP with Arabidopsis queries
AT_BLAST_OUT="${WORK_DIR}/blastp_arabidopsis.outfmt6"
echo "  Arabidopsis queries..."
blastp \
    -query "${AT_QUERY}" \
    -db "${BLASTDB}" \
    -out "${AT_BLAST_OUT}" \
    -evalue 1e-10 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -num_threads "${SLURM_NTASKS_PER_NODE:-8}" \
    -max_target_seqs 50 \
    -seg yes
echo "    Hits: $(wc -l < "${AT_BLAST_OUT}")"

# BLASTP with rice queries
OS_BLAST_OUT="${WORK_DIR}/blastp_rice.outfmt6"
echo "  Rice queries..."
blastp \
    -query "${OS_QUERY}" \
    -db "${BLASTDB}" \
    -out "${OS_BLAST_OUT}" \
    -evalue 1e-10 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -num_threads "${SLURM_NTASKS_PER_NODE:-8}" \
    -max_target_seqs 50 \
    -seg yes
echo "    Hits: $(wc -l < "${OS_BLAST_OUT}")"

# ── Step 4: Extract BLAST hit IDs ────────────────────────────────────────────
echo ""
echo "[Step 4] Extracting and merging hit IDs..."

BLAST_HIT_IDS="${WORK_DIR}/blast_hit_ids.txt"
cat "${AT_BLAST_OUT}" "${OS_BLAST_OUT}" \
    | awk -F'\t' '{print $2}' \
    | sort -u \
    > "${BLAST_HIT_IDS}"

echo "  Unique BLAST hits: $(wc -l < "${BLAST_HIT_IDS}")"

# ── Step 5: Merge with HMM results, deduplicate ─────────────────────────────
echo ""
echo "[Step 5] Merging HMM and BLAST results..."

MERGED_IDS="${WORK_DIR}/merged_candidate_ids.txt"
cat "${HMM_HITS}" "${BLAST_HIT_IDS}" \
    | sort -u \
    > "${MERGED_IDS}"

echo "  HMM-only hits:     $(comm -23 <(sort "${HMM_HITS}") <(sort "${BLAST_HIT_IDS}") | wc -l)"
echo "  BLAST-only hits:   $(comm -13 <(sort "${HMM_HITS}") <(sort "${BLAST_HIT_IDS}") | wc -l)"
echo "  Shared hits:       $(comm -12 <(sort "${HMM_HITS}") <(sort "${BLAST_HIT_IDS}") | wc -l)"
echo "  Total merged:      $(wc -l < "${MERGED_IDS}")"

# Create Venn-style classification
CLASSIFICATION="${WORK_DIR}/candidate_classification.tsv"
echo -e "protein_id\thmm_hit\tblast_hit\tconfidence" > "${CLASSIFICATION}"

while IFS= read -r pid; do
    hmm_flag=0
    blast_flag=0
    grep -qxF "${pid}" "${HMM_HITS}" 2>/dev/null && hmm_flag=1
    grep -qxF "${pid}" "${BLAST_HIT_IDS}" 2>/dev/null && blast_flag=1

    if [[ "${hmm_flag}" -eq 1 && "${blast_flag}" -eq 1 ]]; then
        confidence="high"
    else
        confidence="medium"
    fi

    echo -e "${pid}\t${hmm_flag}\t${blast_flag}\t${confidence}"
done < "${MERGED_IDS}" >> "${CLASSIFICATION}"

# ── Step 6: Extract candidate sequences ──────────────────────────────────────
echo ""
echo "[Step 6] Extracting merged candidate sequences..."

CANDIDATE_FASTA="${WORK_DIR}/candidate_aquaporins.faa"

# Build extraction command
awk 'BEGIN {
    while ((getline line < "'"${MERGED_IDS}"'") > 0) ids[line] = 1
}
/^>/ {
    id = substr($1, 2)
    found = (id in ids)
}
found' "${PROTEOME}" > "${CANDIDATE_FASTA}"

N_EXTRACTED=$(grep -c '^>' "${CANDIDATE_FASTA}" || echo 0)
echo "  Extracted: ${N_EXTRACTED} sequences"

# ── Step 7: Create BLAST best-hit annotation ─────────────────────────────────
echo ""
echo "[Step 7] Creating best-hit annotation (subfamily hints from queries)..."

BEST_HITS="${WORK_DIR}/best_blast_hits.tsv"
echo -e "sunflower_id\tbest_hit_query\tsubfamily_hint\tpident\tevalue\tbitscore\tqcovs" \
    > "${BEST_HITS}"

# For each sunflower candidate, find best BLAST hit from known aquaporins
cat "${AT_BLAST_OUT}" "${OS_BLAST_OUT}" \
    | sort -t$'\t' -k2,2 -k12,12nr \
    | awk -F'\t' '!seen[$2]++ {
        # Parse subfamily from query header (format: sp|ID|GeneName [SUBFAMILY])
        split($1, a, "|")
        gene = a[3]
        sub("\\[", "", gene)
        sub("\\]", "", gene)
        # Extract subfamily from bracketed text
        match($1, /\[([A-Z]+)\]/, m)
        subfamily = m[1]
        if (subfamily == "") subfamily = "unknown"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $2, $1, subfamily, $3, $11, $12, $15
    }' >> "${BEST_HITS}"

echo "  Best-hit annotations: $(( $(wc -l < "${BEST_HITS}") - 1 )) proteins"

# ── Step 8: Create final gene list ───────────────────────────────────────────
echo ""
echo "[Step 8] Creating final candidate aquaporin gene list..."

FINAL_LIST="${WORK_DIR}/final_aquaporin_candidates.txt"
cp "${MERGED_IDS}" "${FINAL_LIST}"

# Also create a symlink/copy in parent directory for downstream scripts
PARENT_OUTPUT="${PROJ_DIR}/02_gene_family/candidate_aquaporin_ids.txt"
cp "${FINAL_LIST}" "${PARENT_OUTPUT}"
cp "${CANDIDATE_FASTA}" "${PROJ_DIR}/02_gene_family/candidate_aquaporins.faa"
cp "${CLASSIFICATION}" "${PROJ_DIR}/02_gene_family/candidate_classification.tsv"

echo "  Final candidate list: ${FINAL_LIST}"
echo "  Copied to: ${PARENT_OUTPUT}"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "================================================================="
echo "  BLASTP Validation Summary"
echo "================================================================="
echo "  Arabidopsis query proteins:  $(grep -c '^>' "${AT_QUERY}" || echo 0)"
echo "  Rice query proteins:         $(grep -c '^>' "${OS_QUERY}" || echo 0)"
echo "  BLAST E-value threshold:     1e-10"
echo ""
echo "  Arabidopsis BLAST hits:      $(wc -l < "${AT_BLAST_OUT}")"
echo "  Rice BLAST hits:             $(wc -l < "${OS_BLAST_OUT}")"
echo "  Unique BLAST hit proteins:   $(wc -l < "${BLAST_HIT_IDS}")"
echo ""
echo "  HMM hits (from step 01):     $(wc -l < "${HMM_HITS}")"
echo "  Merged candidates (total):   $(wc -l < "${MERGED_IDS}")"
echo "  High-confidence (both):      $(awk -F'\t' '$4=="high"' "${CLASSIFICATION}" | wc -l)"
echo ""
echo "  Completed: $(date)"
echo "================================================================="
