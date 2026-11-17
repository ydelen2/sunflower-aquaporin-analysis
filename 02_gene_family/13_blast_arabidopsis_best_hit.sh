#!/bin/bash
#SBATCH --job-name=blast_ath
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=02_gene_family/logs/blast_ath2_%j.out
#SBATCH --error=02_gene_family/logs/blast_ath2_%j.err

# =============================================================================
# 13_blast_arabidopsis_best_hit.sh
# BLASTP of all sunflower proteins against the Arabidopsis RefSeq proteome, one
# best subject per query (E <= 1e-5). The output is the source of the
# LOC-to-TAIR table used for the GO enrichment (Table S4,
# results/02_gene_family/results/arabidopsis_orthologs_v2.tsv).
#
# Gene-level table from this output: a sunflower gene is kept when one of its
# proteins has a hit with identity above 30 % and query coverage above 50 %;
# where the proteins of a gene have different subjects, the hit with the
# highest bit score is taken, and the subject protein is replaced by its TAIR
# locus identifier. Hits to organelle-encoded proteins and to proteins without
# a nuclear TAIR locus identifier, and sunflower organelle genes, are not
# mapped. This gene-level step is not part of the script.
#
# Output
#   02_gene_family/results/sunflower_to_arabidopsis_blast.tsv
#   columns: qseqid sseqid pident length evalue bitscore qcovs
# =============================================================================

module load blast/2.17

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
HA_PROT="${PROJ_DIR}/01_references/GCF_002127325.2_HanXRQr2.0_protein.faa"
ATH_PROT="${PROJ_DIR}/01_references/arabidopsis/protein.faa"
BLASTDB="${PROJ_DIR}/01_references/arabidopsis/ath_blastdb"
OUT="${PROJ_DIR}/02_gene_family/results/sunflower_to_arabidopsis_blast.tsv"

echo "=== BLASTp: Sunflower vs Arabidopsis ==="
echo "Start: $(date)"
echo "Query: $HA_PROT"
echo "DB source: $ATH_PROT"

# Build Arabidopsis BLAST DB
makeblastdb -in "$ATH_PROT" -dbtype prot -out "$BLASTDB" 2>&1
echo "BLAST DB built"

# Run BLASTp
blastp -query "$HA_PROT" -db "$BLASTDB" \
    -out "$OUT" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs" \
    -evalue 1e-5 -max_target_seqs 1 -num_threads 8 2>&1

echo "Hits: $(wc -l < $OUT)"
echo "=== BLAST finished: $(date) ==="
