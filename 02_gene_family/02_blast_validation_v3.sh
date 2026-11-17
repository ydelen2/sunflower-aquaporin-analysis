#!/bin/bash
#SBATCH --job-name=aqp_blast_v3
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/02_blast_validation_v3_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/02_blast_validation_v3_%j.err

# =============================================================================
# 02_blast_validation_v3.sh
# BLASTP of the Arabidopsis and rice aquaporins against the sunflower proteome,
# as a cross-check of the HMM search for aquaporin candidates.
#
# The queries are the Arabidopsis and rice proteins of the reference set of
# the phylogenetic analysis (RefSeq proteins with a PF00230 domain,
# 05_phylogenetics/results_v2); the number of subjects is not limited.
#
# Inputs
#   05_phylogenetics/results_v2/all_aquaporins_combined.fa  (records At_*, Os_*)
#   01_references/sunflower/protein.faa                     sunflower proteome
#   02_gene_family/candidate_aquaporins.faa                 proteins found by the HMM search
#   02_gene_family/verified_aquaporin_ids.txt               the 99 verified proteins
# Outputs (02_gene_family/02_blast_validation_v3/)
#   query_aquaporins.faa      the query sequences
#   blastp_hits.tsv           all hits at E <= 1e-10
#   blastp_subjects.tsv       one row per sunflower protein hit, with its best
#                             query and whether the HMM search found it
#   blastp_summary.txt        counts
# =============================================================================
set -eo pipefail
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
OUT="${PROJ_DIR}/02_gene_family/02_blast_validation_v3"
REFSET="${PROJ_DIR}/05_phylogenetics/results_v2/all_aquaporins_combined.fa"
PROTEOME="${PROJ_DIR}/01_references/sunflower/protein.faa"
HMM_HITS="${PROJ_DIR}/02_gene_family/candidate_aquaporins.faa"
VERIFIED="${PROJ_DIR}/02_gene_family/verified_aquaporin_ids.txt"
mkdir -p "${OUT}/db" "${PROJ_DIR}/logs"

module purge
module load blast/2.17

for f in "${REFSET}" "${PROTEOME}" "${HMM_HITS}" "${VERIFIED}"; do
    [[ -s "${f}" ]] || { echo "missing input: ${f}" >&2; exit 1; }
done

# query set: the Arabidopsis and rice records of the reference set
awk '/^>/{keep = ($0 ~ /^>(At|Os)_/)} keep' "${REFSET}" > "${OUT}/query_aquaporins.faa"
echo "query sequences: $(grep -c '^>' "${OUT}/query_aquaporins.faa") (Arabidopsis $(grep -c '^>At_' "${OUT}/query_aquaporins.faa"), rice $(grep -c '^>Os_' "${OUT}/query_aquaporins.faa"))"
echo "sunflower proteins: $(grep -c '^>' "${PROTEOME}")"

makeblastdb -in "${PROTEOME}" -dbtype prot -out "${OUT}/db/sunflower_proteome" > /dev/null
blastp -version | head -1
blastp -query "${OUT}/query_aquaporins.faa" -db "${OUT}/db/sunflower_proteome" -evalue 1e-10 \
       -num_threads "${SLURM_NTASKS_PER_NODE:-8}" -max_target_seqs 100000 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
       > "${OUT}/blastp_hits.tsv"
echo "hits: $(wc -l < "${OUT}/blastp_hits.tsv")"

python3 - "${OUT}" "${PROTEOME}" "${HMM_HITS}" "${VERIFIED}" << 'PYEOF'
import collections, os, sys
out, proteome, hmm_f, ver_f = sys.argv[1:5]
desc = {}
for line in open(proteome):
    if line.startswith(">"):
        p = line[1:].rstrip("\n").split(None, 1)
        desc[p[0]] = p[1] if len(p) > 1 else ""
hmm = {line[1:].split()[0] for line in open(hmm_f) if line.startswith(">")}
verified = {l.strip() for l in open(ver_f) if l.strip()}
best = {}
nq = collections.defaultdict(set)
for line in open(os.path.join(out, "blastp_hits.tsv")):
    f = line.rstrip("\n").split("\t")
    s = f[1]
    nq[s].add(f[0])
    if s not in best or float(f[11]) > float(best[s][11]):
        best[s] = f
with open(os.path.join(out, "blastp_subjects.tsv"), "w") as fh:
    fh.write("protein_id\tlength_aa\tqueries_with_a_hit\tbest_query\tidentity\tevalue\tbitscore\tquery_coverage\tfound_by_HMM_search\tverified_aquaporin\tproduct\n")
    for s in sorted(best, key=lambda k: -float(best[k][11])):
        f = best[s]
        fh.write("\t".join([s, f[13], str(len(nq[s])), f[0], f[2], f[10], f[11], f[14], "yes" if s in hmm else "no", "yes" if s in verified else "no", desc.get(s, "")]) + "\n")
subj = set(best)
lines = ["query sequences with a hit: %d" % len({f[0] for f in best.values()} | {q for v in nq.values() for q in v}),
         "sunflower proteins hit at E <= 1e-10: %d" % len(subj),
         "  found by the HMM search: %d of %d HMM hits" % (len(subj & hmm), len(hmm)),
         "  not found by the HMM search: %d" % len(subj - hmm),
         "HMM hits without a BLASTP hit: %d %s" % (len(hmm - subj), " ".join(sorted(hmm - subj))),
         "verified aquaporins with a BLASTP hit: %d of %d" % (len(verified & subj), len(verified)),
         "verified aquaporins without a BLASTP hit: %s" % (" ".join(sorted(verified - subj)) or "none")]
for s in sorted(subj - hmm, key=lambda k: -float(best[k][11])):
    f = best[s]
    lines.append("  BLASTP only: %s %s aa, best query %s, identity %s, E %s, query coverage %s, %s" % (s, f[13], f[0], f[2], f[10], f[14], desc.get(s, "")))
open(os.path.join(out, "blastp_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
PYEOF
echo "Finished: $(date)"
