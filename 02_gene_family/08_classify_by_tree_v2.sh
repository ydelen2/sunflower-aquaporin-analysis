#!/bin/bash
#SBATCH --job-name=aqp_classify
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/08_classify_v2_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/08_classify_v2_%j.err
#SBATCH --mail-type=END,FAIL

# ============================================================================
# 08_classify_by_tree_v2.sh
# Subfamily assignment of the sunflower aquaporins by phylogenetic placement.
# Runs classify_by_tree.py on the IQ-TREE tree from 05_phylogenetics/results_v2
# (five species, PF00230-derived reference sets). Each sunflower protein is
# assigned the subfamily carried by the majority of its 5 nearest reference
# tips (patristic distance). Reference subfamilies come from the RefSeq
# product names of the reference proteomes.
#
# Run after 01_msa_and_tree_v2.sh. Needs classify_by_tree.py in
# 02_gene_family/ (copied together with this script).
# ============================================================================

set -eo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
source "${PROJ_DIR}/config.sh"

REF_DIR="${PROJ_DIR}/01_references"
TREE="${PROJ_DIR}/05_phylogenetics/results_v2/aquaporin_tree.treefile"
OUT_DIR="${PROJ_DIR}/02_gene_family/08_classification_v2"
SCRIPT="${PROJ_DIR}/02_gene_family/classify_by_tree.py"
SUNFLOWER_FAA="${PROJ_DIR}/02_gene_family/candidate_aquaporins.faa"
K=5

mkdir -p "${OUT_DIR}" "${PROJ_DIR}/logs"

module purge
module load miniforge/24.5 2>/dev/null || module load miniforge 2>/dev/null || true

ENV_PATH="${PROJ_DIR}/envs/biopython_env"
if [[ ! -x "${ENV_PATH}/bin/python" ]]; then
    echo "  creating conda env with biopython at ${ENV_PATH}..."
    mkdir -p "${PROJ_DIR}/envs"
    conda create -y -p "${ENV_PATH}" -c conda-forge "python=3.10" biopython 2>&1 | tail -3
fi
export PATH="${ENV_PATH}/bin:${PATH}"

# Reference proteomes: same resolution order as 01_msa_and_tree_v2.sh
resolve_proteome() {
    local first="$1" second="$2"
    if [[ -s "${first}" ]]; then echo "${first}"; return 0; fi
    if [[ -n "${second}" && -s "${second}" ]]; then echo "${second}"; return 0; fi
    return 1
}
ATHA=$(resolve_proteome "${REF_DIR}/arabidopsis/protein.faa" "${ATHA_PROT:-}") || { echo "ERROR: Arabidopsis proteome not found" >&2; exit 1; }
OSAT=$(resolve_proteome "${REF_DIR}/rice/protein.faa"        "${OSAT_PROT:-}") || { echo "ERROR: rice proteome not found" >&2; exit 1; }
SLYC=$(resolve_proteome "${REF_DIR}/tomato/protein.faa"      "${SLYC_PROT:-}") || { echo "ERROR: tomato proteome not found" >&2; exit 1; }
LSAT=$(resolve_proteome "${REF_DIR}/lettuce/protein.faa"     "${LSAT_PROT:-}") || { echo "ERROR: lettuce proteome not found" >&2; exit 1; }

for f in "${TREE}" "${SCRIPT}" "${SUNFLOWER_FAA}"; do
    [[ -s "${f}" ]] || { echo "ERROR: missing input ${f}" >&2; exit 1; }
done

echo "================================================================"
echo "Subfamily classification by tree placement (k = ${K})"
echo "tree: ${TREE}"
echo "================================================================"

python "${SCRIPT}" \
    --tree "${TREE}" \
    --ref "${ATHA}" --ref "${OSAT}" --ref "${SLYC}" --ref "${LSAT}" \
    --sunflower "${SUNFLOWER_FAA}" \
    --k "${K}" \
    --out "${OUT_DIR}/subfamily_by_tree.tsv"

echo ""
echo "Agreement between tree-based call and NCBI product name:"
awk -F'\t' 'NR > 1 { total++; if ($2 == $7) agree++; else print "  disagreement: " $1 "\t" $2 " (tree)\t" $7 " (name)\t" $8 }
            END { printf "  %d of %d agree\n", agree, total }' "${OUT_DIR}/subfamily_by_tree.tsv"

echo ""
echo "Output: ${OUT_DIR}/subfamily_by_tree.tsv"
echo "Finished: $(date)"
