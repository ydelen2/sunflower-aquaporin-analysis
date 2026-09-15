#!/bin/bash
#SBATCH --job-name=aqp_domain_v2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/03_domain_verification_v2_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/03_domain_verification_v2_%j.err
#SBATCH --mail-type=END,FAIL

# ============================================================================
# 03_domain_verification_v2.sh
# Domain verification of aquaporin candidates, second version.
#
# Difference from the first version: transmembrane helices are predicted with
# DeepTMHMM (Hallgren et al. 2022) instead of the Kyte-Doolittle sliding
# window that the first version used as a screening proxy. The proxy
# under-counted helices in many full-length aquaporins, so the "at least 4
# TM helices" rule removed about 30 annotated aquaporin genes that carry a
# complete MIP domain.
#
# Criteria kept from the first version:
#   1. MIP domain (Pfam PF00230) present, domain E-value <= 1e-5
#   2. at least 4 predicted TM helices
#
# Outputs go to 02_gene_family/03_domain_verification_v2/. The three files
# that downstream scripts read from 02_gene_family/ (verified_aquaporin_ids.txt,
# verified_aquaporins.faa, verified_aquaporins.tsv) are replaced, after the
# first-version copies are saved with a .v1_backup suffix.
#
# DeepTMHMM runs through the BioLib client (pip package pybiolib); the
# computation itself happens on the BioLib servers, so the node needs
# outbound internet access. If the job fails at that step, run the TM step
# on the login node with:
#   cd 02_gene_family/03_domain_verification_v2 && biolib run DTU/DeepTMHMM --fasta candidates_clean.faa
# and resubmit; the script picks up an existing biolib_results/TMRs.gff3.
# ============================================================================

set -eo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
source "${PROJ_DIR}/config.sh"

PARENT_DIR="${PROJ_DIR}/02_gene_family"
WORK_DIR="${PARENT_DIR}/03_domain_verification_v2"
LOG_DIR="${PROJ_DIR}/logs"
mkdir -p "${WORK_DIR}" "${LOG_DIR}"

CANDIDATE_FASTA="${PARENT_DIR}/candidate_aquaporins.faa"
CANDIDATE_IDS="${PARENT_DIR}/candidate_aquaporin_ids.txt"
HMM_FILE="${PARENT_DIR}/01_hmm_search/PF00230.hmm"

MIP_EVALUE="1e-5"
MIN_TM=4

# ---------------------------------------------------------------------------
# Modules and tools
# ---------------------------------------------------------------------------
module purge

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

try_modules "hmmer/3.4"      "hmmer"     "HMMER/3.4" "HMMER"
try_modules "miniforge/24.5" "miniforge" "anaconda"  "miniconda"

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
        return 1
    fi
    echo "  ${binary} not on PATH; installing ${pkg} into ${env_path}..."
    mkdir -p "${PROJ_DIR}/envs"
    conda create -y -p "${env_path}" -c conda-forge -c bioconda "${pkg}" 2>&1 | tail -5
    export PATH="${env_path}/bin:${PATH}"
    command -v "${binary}" &>/dev/null
}

ensure_tool hmmsearch hmmer hmmer_env || exit 1

# BioLib client for DeepTMHMM: a small conda env with python and pip
ensure_biolib() {
    command -v biolib &>/dev/null && return 0
    local env_path="${PROJ_DIR}/envs/tm_env"
    if [[ -x "${env_path}/bin/biolib" ]]; then
        export PATH="${env_path}/bin:${PATH}"
        echo "  biolib: using existing conda env ${env_path}"
        return 0
    fi
    if ! command -v conda &>/dev/null; then
        echo "ERROR: conda not found; cannot install the BioLib client." >&2
        return 1
    fi
    echo "  installing pybiolib into ${env_path}..."
    mkdir -p "${PROJ_DIR}/envs"
    conda create -y -p "${env_path}" -c conda-forge "python=3.10" pip 2>&1 | tail -3
    "${env_path}/bin/pip" install --quiet pybiolib
    export PATH="${env_path}/bin:${PATH}"
    command -v biolib &>/dev/null
}

echo "================================================================"
echo "Domain verification v2 (PF00230 + DeepTMHMM)"
echo "Started: $(date)"
echo "Job ID:  ${SLURM_JOB_ID:-local}"
echo "================================================================"

for f in "${CANDIDATE_FASTA}" "${CANDIDATE_IDS}" "${HMM_FILE}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required input not found: ${f}" >&2
        exit 1
    fi
done
N_CANDIDATES=$(grep -c '^>' "${CANDIDATE_FASTA}" || true)
echo "Input candidates: ${N_CANDIDATES}"

# ---------------------------------------------------------------------------
# Step 1: MIP domain (PF00230) with hmmsearch
# ---------------------------------------------------------------------------
echo ""
echo "[Step 1] hmmsearch with the PF00230 profile..."

DOMTBL="${WORK_DIR}/pf00230_candidates.domtbl"
hmmsearch \
    --cpu "${SLURM_NTASKS_PER_NODE:-8}" \
    -E "${MIP_EVALUE}" \
    --domE "${MIP_EVALUE}" \
    --domtblout "${DOMTBL}" \
    -o "${WORK_DIR}/pf00230_candidates.out" \
    "${HMM_FILE}" "${CANDIDATE_FASTA}"

MIP_HITS="${WORK_DIR}/mip_domain_hits.tsv"
printf 'protein_id\tdom_evalue\tdom_score\tali_from\tali_to\tdomain_len\n' > "${MIP_HITS}"
grep -v '^#' "${DOMTBL}" \
    | awk 'BEGIN{OFS="\t"} {print $1, $13, $14, $18, $19, ($19-$18+1)}' \
    | sort -k1,1 -k3,3gr \
    | awk -F'\t' '!seen[$1]++' \
    >> "${MIP_HITS}"

MIP_IDS="${WORK_DIR}/mip_positive_ids.txt"
awk -F'\t' 'NR > 1 {print $1}' "${MIP_HITS}" | sort -u > "${MIP_IDS}"
N_MIP=$(wc -l < "${MIP_IDS}")
echo "  Proteins with a PF00230 domain (E <= ${MIP_EVALUE}): ${N_MIP}"

# ---------------------------------------------------------------------------
# Step 2: transmembrane helices with DeepTMHMM
# ---------------------------------------------------------------------------
echo ""
echo "[Step 2] Transmembrane helix prediction with DeepTMHMM..."

# DeepTMHMM wants plain headers; keep only the accession
CLEAN_FASTA="${WORK_DIR}/candidates_clean.faa"
awk '/^>/ {print ">" substr($1, 2); next} {print}' "${CANDIDATE_FASTA}" > "${CLEAN_FASTA}"

TMR_GFF="${WORK_DIR}/biolib_results/TMRs.gff3"
if [[ -s "${TMR_GFF}" ]]; then
    echo "  existing DeepTMHMM output found, reusing: ${TMR_GFF}"
else
    ensure_biolib || exit 1
    cd "${WORK_DIR}"
    biolib run DTU/DeepTMHMM --fasta "${CLEAN_FASTA}"
    cd "${PROJ_DIR}"
    if [[ ! -s "${TMR_GFF}" ]]; then
        echo "ERROR: DeepTMHMM produced no TMRs.gff3 in ${WORK_DIR}/biolib_results/." >&2
        echo "       Check network access from the node, or run the biolib command on the login node." >&2
        exit 1
    fi
fi

# Count TMhelix segments per protein from the GFF3
TM_RESULTS="${WORK_DIR}/tm_helix_predictions.tsv"
python3 - "${TMR_GFF}" "${CLEAN_FASTA}" "${TM_RESULTS}" << 'PYEOF'
import sys, re
gff, fasta, out = sys.argv[1:4]
lengths, order = {}, []
seq_id = None
for line in open(fasta):
    line = line.strip()
    if line.startswith('>'):
        seq_id = line[1:].split()[0]
        order.append(seq_id); lengths[seq_id] = 0
    elif seq_id:
        lengths[seq_id] += len(line)
helices = {k: [] for k in order}
for line in open(gff):
    if line.startswith('#') or not line.strip():
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 4:
        continue
    pid, feat, start, end = f[0], f[1], f[2], f[3]
    if feat == 'TMhelix' and pid in helices:
        helices[pid].append((int(start), int(end)))
with open(out, 'w') as fh:
    fh.write('protein_id\tnum_tm_helices\ttm_positions\tsequence_length\n')
    for pid in order:
        h = sorted(helices[pid])
        pos = ';'.join('%d-%d' % (s, e) for s, e in h) if h else 'none'
        fh.write('%s\t%d\t%s\t%d\n' % (pid, len(h), pos, lengths[pid]))
counts = {}
for pid in order:
    n = len(helices[pid]); counts[n] = counts.get(n, 0) + 1
print('  TM helix count distribution (DeepTMHMM):')
for n in sorted(counts):
    print('    %d TM: %d proteins' % (n, counts[n]))
if sum(len(v) for v in helices.values()) == 0:
    sys.stderr.write('ERROR: no TMhelix features were parsed from %s; check the GFF3 format.\n' % gff)
    sys.exit(1)
PYEOF

TM_PASS_IDS="${WORK_DIR}/tm_pass_ids.txt"
awk -F'\t' -v m="${MIN_TM}" 'NR > 1 && $2 >= m {print $1}' "${TM_RESULTS}" | sort -u > "${TM_PASS_IDS}"
echo "  Proteins with >= ${MIN_TM} TM helices: $(wc -l < "${TM_PASS_IDS}")"

# ---------------------------------------------------------------------------
# Step 3: combined filter
# ---------------------------------------------------------------------------
echo ""
echo "[Step 3] Combined filter: PF00230 domain AND >= ${MIN_TM} TM helices..."

VERIFIED_IDS="${WORK_DIR}/verified_aquaporin_ids.txt"
comm -12 <(sort "${MIP_IDS}") <(sort "${TM_PASS_IDS}") > "${VERIFIED_IDS}"
N_VERIFIED=$(wc -l < "${VERIFIED_IDS}")
echo "  Verified aquaporin proteins: ${N_VERIFIED}"

VERIFIED_FASTA="${WORK_DIR}/verified_aquaporins.faa"
awk 'BEGIN { while ((getline line < "'"${VERIFIED_IDS}"'") > 0) ids[line] = 1 }
     /^>/ { id = substr($1, 2); found = (id in ids) }
     found' "${CANDIDATE_FASTA}" > "${VERIFIED_FASTA}"

VERIFIED_TABLE="${WORK_DIR}/verified_aquaporins.tsv"
python3 - "${VERIFIED_IDS}" "${TM_RESULTS}" "${MIP_HITS}" "${VERIFIED_TABLE}" << 'PYEOF'
import sys
ids_f, tm_f, mip_f, out_f = sys.argv[1:5]
tm = {}
for i, line in enumerate(open(tm_f)):
    if i == 0: continue
    f = line.rstrip('\n').split('\t'); tm[f[0]] = f[1:4]
mip = {}
for i, line in enumerate(open(mip_f)):
    if i == 0: continue
    f = line.rstrip('\n').split('\t'); mip.setdefault(f[0], f[1:6])
with open(out_f, 'w') as out:
    out.write('protein_id\tmip_domain\tnum_tm_helices\ttm_positions\tseq_length\tmip_evalue\tmip_score\tmip_ali_from\tmip_ali_to\n')
    for line in open(ids_f):
        pid = line.strip()
        if not pid: continue
        t = tm.get(pid, ['NA', 'NA', 'NA']); m = mip.get(pid, ['NA'] * 5)
        out.write('\t'.join([pid, 'yes', t[0], t[1], t[2], m[0], m[1], m[2], m[3]]) + '\n')
PYEOF

FAILED="${WORK_DIR}/failed_verification.tsv"
python3 - "${CANDIDATE_IDS}" "${VERIFIED_IDS}" "${MIP_IDS}" "${TM_RESULTS}" "${FAILED}" << 'PYEOF'
import sys
cand_f, ver_f, mip_f, tm_f, out_f = sys.argv[1:6]
ver = {l.strip() for l in open(ver_f) if l.strip()}
mip = {l.strip() for l in open(mip_f) if l.strip()}
tm = {}
for i, line in enumerate(open(tm_f)):
    if i == 0: continue
    f = line.rstrip('\n').split('\t'); tm[f[0]] = (int(f[1]), f[3])
n = 0
with open(out_f, 'w') as out:
    out.write('protein_id\thas_mip\tnum_tm\tseq_length\treason\n')
    for line in open(cand_f):
        pid = line.strip()
        if not pid or pid in ver: continue
        has = 'yes' if pid in mip else 'no'
        ntm, ln = tm.get(pid, (0, 'NA'))
        reasons = []
        if has == 'no': reasons.append('no_MIP_domain')
        if ntm < 4: reasons.append('insufficient_TM_helices(%d)' % ntm)
        out.write('%s\t%s\t%d\t%s\t%s\n' % (pid, has, ntm, ln, ';'.join(reasons)))
        n += 1
print('  Failed verification: %d' % n)
PYEOF

# ---------------------------------------------------------------------------
# Step 4: comparison with the first version, then replace the shared outputs
# ---------------------------------------------------------------------------
echo ""
echo "[Step 4] Comparing with the first-version gene set and updating shared outputs..."

for f in verified_aquaporin_ids.txt verified_aquaporins.faa verified_aquaporins.tsv; do
    if [[ -f "${PARENT_DIR}/${f}" && ! -f "${PARENT_DIR}/${f}.v1_backup" ]]; then
        cp "${PARENT_DIR}/${f}" "${PARENT_DIR}/${f}.v1_backup"
        echo "  saved first-version copy: ${f}.v1_backup"
    fi
done

OLD_IDS="${PARENT_DIR}/verified_aquaporin_ids.txt.v1_backup"
COMPARE="${WORK_DIR}/comparison_v1_vs_v2.tsv"
if [[ -f "${OLD_IDS}" ]]; then
    {
        printf 'protein_id\tstatus\n'
        comm -13 <(sort "${OLD_IDS}") <(sort "${VERIFIED_IDS}") | sed 's/$/\tnew_in_v2/'
        comm -23 <(sort "${OLD_IDS}") <(sort "${VERIFIED_IDS}") | sed 's/$/\tdropped_in_v2/'
        comm -12 <(sort "${OLD_IDS}") <(sort "${VERIFIED_IDS}") | sed 's/$/\tin_both/'
    } > "${COMPARE}"
    echo "  v1 proteins: $(wc -l < "${OLD_IDS}") | v2 proteins: ${N_VERIFIED}"
    echo "  new in v2:     $(grep -c 'new_in_v2' "${COMPARE}" || true)"
    echo "  dropped in v2: $(grep -c 'dropped_in_v2' "${COMPARE}" || true)"
    echo "  comparison table: ${COMPARE}"
fi

cp "${VERIFIED_IDS}"   "${PARENT_DIR}/verified_aquaporin_ids.txt"
cp "${VERIFIED_FASTA}" "${PARENT_DIR}/verified_aquaporins.faa"
cp "${VERIFIED_TABLE}" "${PARENT_DIR}/verified_aquaporins.tsv"

# 02_gene_family/results/ holds symlinks into the first-version folder
# (aquaporin_proteins.fa -> ../03_domain_verification/verified_aquaporins.faa,
#  verified_aquaporins.txt -> ../03_domain_verification/verified_aquaporin_ids.txt)
# and a copy of the verified table. Downstream scripts read these, so point
# them at the v2 files. verified_aquaporins.txt is rebuilt later as a real
# gene table by 07_build_gene_tables_v2.sh.
RES_DIR="${PARENT_DIR}/results"
mkdir -p "${RES_DIR}"
for f in aquaporin_proteins.fa verified_aquaporins.txt verified_aquaporins.tsv; do
    if [[ -e "${RES_DIR}/${f}" || -L "${RES_DIR}/${f}" ]] && [[ ! -e "${RES_DIR}/${f}.v1_backup" ]]; then
        cp -L "${RES_DIR}/${f}" "${RES_DIR}/${f}.v1_backup" 2>/dev/null || true
    fi
done
ln -sfn "../03_domain_verification_v2/verified_aquaporins.faa"    "${RES_DIR}/aquaporin_proteins.fa"
ln -sfn "../03_domain_verification_v2/verified_aquaporin_ids.txt" "${RES_DIR}/verified_aquaporins.txt"
cp "${VERIFIED_TABLE}" "${RES_DIR}/verified_aquaporins.tsv"
echo "  results/ links refreshed: aquaporin_proteins.fa, verified_aquaporins.txt, verified_aquaporins.tsv"

echo ""
echo "================================================================"
echo "Summary"
echo "  candidates:                 ${N_CANDIDATES}"
echo "  PF00230 positive:           ${N_MIP}"
echo "  >= ${MIN_TM} TM (DeepTMHMM):       $(wc -l < "${TM_PASS_IDS}")"
echo "  verified (both criteria):   ${N_VERIFIED}"
echo "  outputs:                    ${WORK_DIR}/"
echo "  shared outputs updated in:  ${PARENT_DIR}/"
echo "Finished: $(date)"
echo "================================================================"
