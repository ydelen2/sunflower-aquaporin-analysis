#!/usr/bin/env bash
#SBATCH --job-name=synteny_v3
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/synteny_v3_%j.log
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/synteny_v3_%j.err

# =============================================================================
# 05_synteny_one_protein_per_gene_v3.sh
# MCScanX collinearity with one protein per gene in the BLASTP searches.
#
# The BLASTP searches use only the proteins of the MCScanX position files (one
# protein per gene), with a separate search for every direction, so that hits
# to other isoforms of a gene or to proteins of the same species do not take
# up the five target sequences per query:
#
#   Ha_self_v3   sunflower against sunflower
#   Ha_Ls_v3     sunflower against lettuce and lettuce against sunflower
#   Ha_At_v3     sunflower against Arabidopsis and Arabidopsis against sunflower
#
# BLASTP: E <= 1e-10, five target sequences per query. MCScanX: default
# parameters, minimum block size 5. The position files of
# 07_synteny/data are used unchanged.
#
# Outputs
#   07_synteny/blast_v3/*.blast
#   07_synteny/mcscanx/Ha_self_v3, Ha_Ls_v3, Ha_At_v3 (collinearity files)
# =============================================================================
set -euo pipefail

PROJ_DIR="/work/dweikat/ydelen2/aquaporin_study"
SYN="${PROJ_DIR}/07_synteny"
REF="${PROJ_DIR}/01_references"
HA_PROT="${REF}/sunflower/protein.faa"
AT_PROT="${REF}/arabidopsis/protein_ncbi.faa"
LS_PROT="${REF}/lettuce/protein_ncbi.faa"
[[ -f "${AT_PROT}" ]] || AT_PROT="${REF}/arabidopsis/protein.faa"
[[ -f "${LS_PROT}" ]] || LS_PROT="${REF}/lettuce/protein.faa"
B="${SYN}/blast_v3"
THREADS="${SLURM_NTASKS_PER_NODE:-16}"
mkdir -p "${B}" "${PROJ_DIR}/logs"

module load blast/2.17
module load miniforge/24.5
set +u
conda activate "${PROJ_DIR}/conda_envs/aquaporin_env"
set -u
command -v MCScanX > /dev/null || { echo "MCScanX not found" >&2; exit 1; }
blastp -version | head -1
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

for f in "${HA_PROT}" "${AT_PROT}" "${LS_PROT}" "${SYN}/data/Ha.gff" "${SYN}/data/At.gff" "${SYN}/data/Ls.gff"; do
    [[ -s "${f}" ]] || { echo "missing input: ${f}" >&2; exit 1; }
done

# one protein per gene: the proteins named in the MCScanX position files
python3 - "${SYN}/data" "${B}" "${HA_PROT}" "${AT_PROT}" "${LS_PROT}" << 'PYEOF'
import os, sys
data_dir, out_dir = sys.argv[1], sys.argv[2]
for sp, faa in zip(("Ha", "At", "Ls"), sys.argv[3:6]):
    keep = {line.split()[1] for line in open(os.path.join(data_dir, sp + ".gff")) if line.strip()}
    n = 0; write = False
    with open(os.path.join(out_dir, sp + ".one_per_gene.faa"), "w") as out:
        for line in open(faa):
            if line.startswith(">"):
                write = line[1:].split()[0] in keep
                n += write
            if write:
                out.write(line)
    print("%s: proteins in the position file %d, sequences written %d" % (sp, len(keep), n))
    if n != len(keep):
        sys.exit("not every protein of %s.gff was found in %s" % (sp, faa))
PYEOF

for sp in Ha At Ls; do
    makeblastdb -in "${B}/${sp}.one_per_gene.faa" -dbtype prot -out "${B}/${sp}_db" > /dev/null
done

run_blast() {   # query species, database species
    local out="${B}/$1_vs_$2.blast"
    if [[ -s "${out}.done" ]]; then log "  $1 vs $2 already finished"; return 0; fi
    log "  BLASTP $1 vs $2"
    blastp -query "${B}/$1.one_per_gene.faa" -db "${B}/$2_db" -out "${out}" \
           -evalue 1e-10 -num_threads "${THREADS}" -outfmt 6 -max_target_seqs 5
    wc -l < "${out}" > "${out}.done"
    log "  $1 vs $2: $(cat "${out}.done") hits"
}
run_blast Ha Ha
run_blast Ha Ls
run_blast Ls Ha
run_blast Ha At
run_blast At Ha

run_mcscanx() {   # run name, position files, BLAST files
    local name="$1" gffs="$2" blasts="$3" d="${SYN}/mcscanx/$1"
    mkdir -p "${d}"
    cat ${gffs} > "${d}/${name}.gff"
    cat ${blasts} > "${d}/${name}.blast"
    ( cd "${d}" && MCScanX "${name}" -s 5 > "${name}.mcscanx.log" 2>&1 )
    [[ -s "${d}/${name}.collinearity" ]] || { echo "MCScanX gave no collinearity file for ${name}" >&2; exit 1; }
    log "  ${name}: $(grep 'Number of collinear genes' "${d}/${name}.collinearity") | blocks $(grep -c '^## Alignment' "${d}/${name}.collinearity")"
}
log "MCScanX"
run_mcscanx Ha_self_v3 "${SYN}/data/Ha.gff" "${B}/Ha_vs_Ha.blast"
run_mcscanx Ha_Ls_v3 "${SYN}/data/Ha.gff ${SYN}/data/Ls.gff" "${B}/Ha_vs_Ls.blast ${B}/Ls_vs_Ha.blast"
run_mcscanx Ha_At_v3 "${SYN}/data/Ha.gff ${SYN}/data/At.gff" "${B}/Ha_vs_At.blast ${B}/At_vs_Ha.blast"
log "finished"
