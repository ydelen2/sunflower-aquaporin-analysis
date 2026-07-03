#!/usr/bin/env bash
# =============================================================================
# 01_create_dirs.sh — Create the full project directory structure
# Run locally on SWAN login node (not a SLURM job)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

log_step "Creating project directory structure under ${PROJ_DIR}"

# ---------------------------------------------------------------------------
# Core directory tree
# ---------------------------------------------------------------------------
directories=(
    "${SETUP_DIR}"
    "${REF_DIR}/sunflower/hisat2_index"
    "${REF_DIR}/arabidopsis"
    "${REF_DIR}/rice"
    "${REF_DIR}/tomato"
    "${REF_DIR}/lettuce"
    "${REF_DIR}/aquaporin_hmm"
    "${GENE_FAM_DIR}/blast_results"
    "${GENE_FAM_DIR}/hmmer_results"
    "${GENE_FAM_DIR}/identified_genes"
    "${GENE_FAM_DIR}/sequences"
    "${GENE_FAM_DIR}/domain_analysis"
    "${FASTQ_DIR}"
    "${QC_DIR}/raw"
    "${QC_DIR}/trimmed"
    "${QC_DIR}/multiqc"
    "${TRIMMED_DIR}"
    "${ALIGNED_DIR}"
    "${COUNTS_DIR}"
    "${EXPR_DIR}/deseq2"
    "${EXPR_DIR}/wgcna"
    "${EXPR_DIR}/heatmaps"
    "${PHYLO_DIR}/alignments"
    "${PHYLO_DIR}/trees"
    "${CIS_DIR}/promoters"
    "${CIS_DIR}/meme_results"
    "${SYNTENY_DIR}/mcscanx"
    "${SYNTENY_DIR}/kaks"
    "${STRUCT_DIR}/alphafold"
    "${STRUCT_DIR}/models"
    "${FIG_DIR}/publication"
    "${FIG_DIR}/supplementary"
    "${LOG_DIR}"
    "${PROJ_DIR}/R_libs"
    "${PROJ_DIR}/conda_envs"
    "${PROJ_DIR}/metadata"
    "${PROJ_DIR}/scripts"
)

# Per-BioProject FASTQ subdirectories
for bp in "${BIOPROJECTS[@]}"; do
    directories+=("${FASTQ_DIR}/${bp}")
done

for dir in "${directories[@]}"; do
    mkdir -p "${dir}"
    log_info "Created: ${dir}"
done

# ---------------------------------------------------------------------------
# Create placeholder metadata files
# ---------------------------------------------------------------------------
log_step "Creating sample metadata templates"

for bp in "${BIOPROJECTS[@]}"; do
    meta_file="${PROJ_DIR}/metadata/${bp}_samples.tsv"
    if [[ ! -f "${meta_file}" ]]; then
        printf "sample_id\tsra_accession\tbioproject\tcondition\ttissue\tgenotype\ttime_point\treplicate\n" \
            > "${meta_file}"
        log_info "Created metadata template: ${meta_file}"
    fi
done

# PRJNA869183 metadata header with specific conditions
cat > "${PROJ_DIR}/metadata/PRJNA869183_samples.tsv" <<'EOF'
sample_id	sra_accession	bioproject	condition	tissue	genotype	time_point	replicate
# cold_8h_rep1	SRR21038001	PRJNA869183	cold_4C	leaf	seedling	8h	1
# Fill in actual accessions after SRA query
EOF

# PRJNA492303 metadata
cat > "${PROJ_DIR}/metadata/PRJNA492303_samples.tsv" <<'EOF'
sample_id	sra_accession	bioproject	condition	tissue	genotype	time_point	replicate
# flooding_HA351_leaf_rep1	SRR7837001	PRJNA492303	flooding	leaf	HA351_resistant	NA	1
# flooding_RHA428_root_rep1	SRR7837005	PRJNA492303	flooding	root	RHA428_susceptible	NA	1
# Fill in actual accessions after SRA query
EOF

# ---------------------------------------------------------------------------
# Create a master sample sheet
# ---------------------------------------------------------------------------
cat > "${PROJ_DIR}/metadata/bioproject_summary.tsv" <<'EOF'
bioproject	stress_type	stress_category	tissue	genotypes	n_samples	description
PRJNA869183	cold_heat_drought_salt	abiotic	leaf	seedling	~45	Cold 4°C (8/16/32h), Heat 39°C (4/8/16/32h), Drought PEG, Salt 150mM NaCl
PRJNA492303	flooding	abiotic	leaf,root	HA351,RHA428	48	Flooding stress; resistant vs susceptible; 4 reps
PRJNA1041959	drought_PEG	abiotic	leaf,root	K55,K58	24	Drought PEG; sensitive vs tolerant genotypes
PRJNA797473	drought	abiotic	leaf	unknown	~9	Drought 0/7/14 days
PRJNA908908	sclerotinia	biotic	leaf	unknown	~12	Sclerotinia sclerotiorum infection
PRJNA850121	orobanche	biotic	root	unknown	36	Orobanche cumana parasitism
EOF

log_done "Directory structure created successfully"
log_info "Project root: ${PROJ_DIR}"
log_info "Run 'tree ${PROJ_DIR} -d -L 3' to visualize"
