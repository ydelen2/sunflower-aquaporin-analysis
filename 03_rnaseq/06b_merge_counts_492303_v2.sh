#!/bin/bash
#SBATCH --job-name=fcmerge_492303
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/logs/06b_fcmerge_492303_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/03_rnaseq/logs/06b_fcmerge_492303_%j.err

# =============================================================================
# 06b_merge_counts_492303_v2.sh
# Joins the per-run featureCounts files of PRJNA492303 into one matrix
# (gene_id + one column per run, run order as in the accession list) and writes
# an assignment summary per run.
# Output: counts/PRJNA492303_full_counts_clean.txt
#         counts/PRJNA492303_full_counts_summary.tsv
# =============================================================================
set -eo pipefail
source /work/dweikat/ydelen2/aquaporin_study/03_rnaseq/config.sh
BP_ID="PRJNA492303"
FC_DIR="${COUNTS_DIR}/${BP_ID}_fc"
LIST="${ACCESSION_DIR}/${BP_ID}_srr.txt"
CLEAN="${COUNTS_DIR}/${BP_ID}_full_counts_clean.txt"
SUMM="${COUNTS_DIR}/${BP_ID}_full_counts_summary.tsv"

python3 - "$FC_DIR" "$LIST" "$CLEAN" "$SUMM" <<'PY'
import sys, os
fc_dir, lst, clean, summ = sys.argv[1:5]
runs = [l.strip() for l in open(lst) if l.strip()]
genes = None; cols = {}; summary = []
for r in runs:
    f = os.path.join(fc_dir, r + ".txt")
    if not os.path.exists(f):
        sys.exit("missing counts for " + r)
    ids, vals = [], []
    with open(f) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            p = line.rstrip("\n").split("\t")
            ids.append(p[0]); vals.append(p[-1])
    if genes is None:
        genes = ids
    elif ids != genes:
        sys.exit("gene order differs in " + r)
    cols[r] = vals
    st = {}
    with open(f + ".summary") as fh:
        next(fh)
        for line in fh:
            k, v = line.rstrip("\n").split("\t")[:2]; st[k] = int(v)
    tot = sum(st.values()); summary.append((r, st.get("Assigned", 0), tot, st.get("Unassigned_NoFeatures", 0), st.get("Unassigned_MultiMapping", 0), st.get("Unassigned_Unmapped", 0)))
with open(clean, "w") as out:
    out.write("gene_id\t" + "\t".join(runs) + "\n")
    for i, g in enumerate(genes):
        out.write(g + "\t" + "\t".join(cols[r][i] for r in runs) + "\n")
with open(summ, "w") as out:
    out.write("run\tassigned\ttotal\tpct_assigned\tno_features\tmulti_mapping\tunmapped\n")
    for r, a, t, nf, mm, um in summary:
        out.write("%s\t%d\t%d\t%.1f\t%d\t%d\t%d\n" % (r, a, t, 100.0 * a / t if t else 0, nf, mm, um))
print("matrix:", clean, "| genes:", len(genes), "| runs:", len(runs))
pcts = sorted(100.0 * a / t for _, a, t, *_ in summary)
print("assigned percent: min %.1f, median %.1f, max %.1f" % (pcts[0], pcts[len(pcts) // 2], pcts[-1]))
PY
log_msg "Done: ${CLEAN}"
