#!/bin/bash
#SBATCH --job-name=aqp_deg_table
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/aqp_deg_table_v2_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/aqp_deg_table_v2_%j.err
#SBATCH --mail-type=END,FAIL

# ============================================================================
# 04_aquaporin_deg_table_v2.sh
# Aquaporin differential-expression matrix across all DESeq2 contrasts.
# In the first analysis this table (Table S4) was assembled from the
# per-contrast DESeq2 result files in an unscripted step; this script makes it
# reproducible. DESeq2 itself is not re-run: the genome-wide results in
# 04_expression/results/deseq2/<project>/DEG_<contrast>.tsv are unchanged.
#
# Inputs
#   02_gene_family/results/verified_aquaporins.txt   (07_build_gene_tables_v2.sh)
#   04_expression/results/deseq2/*/DEG_*_vs_*.tsv     (01_deseq2_analysis.sh)
# Outputs (04_expression/results/aquaporin_deg_v2/ and 09_figures/supplementary/)
#   Table_S4_aqp_deg_all_contrasts.tsv   gene x contrast log2FC and padj
#   aquaporin_deg_counts.tsv             up / down counts per contrast
#   aquaporin_deg_long.tsv               one row per gene x contrast (significant only)
#   broadly_responsive_genes.tsv         genes ranked by number of contrasts
#   contrast_map.tsv                     file name -> label used in the manuscript
# ============================================================================

set -eo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
DESEQ="${PROJ_DIR}/04_expression/results/deseq2"
AQP="${PROJ_DIR}/02_gene_family/results/verified_aquaporins.txt"
OUT="${PROJ_DIR}/04_expression/results/aquaporin_deg_v2"
SUPP="${PROJ_DIR}/09_figures/supplementary"
mkdir -p "${OUT}" "${SUPP}" "${PROJ_DIR}/04_expression/logs"

[[ -s "${AQP}" ]] || { echo "ERROR: ${AQP} missing (run 07_build_gene_tables_v2.sh first)" >&2; exit 1; }
[[ -d "${DESEQ}" ]] || { echo "ERROR: ${DESEQ} missing" >&2; exit 1; }
if [[ -f "${SUPP}/Table_S4_aqp_deg_all_contrasts.tsv" && ! -f "${SUPP}/Table_S4_aqp_deg_all_contrasts.tsv.v1_backup" ]]; then
    cp "${SUPP}/Table_S4_aqp_deg_all_contrasts.tsv" "${SUPP}/Table_S4_aqp_deg_all_contrasts.tsv.v1_backup"
fi

python3 - "${DESEQ}" "${AQP}" "${OUT}" "${SUPP}" << 'PYEOF'
import sys, os, glob, re, collections, math
deseq, aqp_f, out, supp = sys.argv[1:5]
LFC_CUT, PADJ_CUT = 1.0, 0.05

# labels used in the manuscript for the 16 contrasts; anything not listed keeps project__contrast
LABELS = {
    ("PRJNA869183", "cold"): "Cold", ("PRJNA869183", "heat"): "Heat", ("PRJNA869183", "drought"): "Drought_869",
    ("PRJNA869183", "salt"): "Salt_NaCl", ("PRJNA869183", "rehydrat"): "Rehydration",
    ("PRJNA797473", "7d"): "Drought_7d", ("PRJNA797473", "14d"): "Drought_14d", ("PRJNA797473", "21d"): "Drought_21d",
    ("PRJNA1041959", "peg"): "PEG_72h", ("PRJNA492303", "flood"): "Flooding", ("PRJNA908908", "sclero"): "Sclerotinia",
}
ORDER = ["Cold", "Heat", "Drought_869", "Salt_NaCl", "Rehydration", "Drought_7d", "Drought_14d", "Drought_21d",
         "PEG_72h", "Flooding", "Sclerotinia", "Orobanche_A", "Orobanche_B", "Orobanche_C", "Orobanche_D", "Orobanche_E"]

def label_for(project, contrast):
    c = contrast.lower()
    stress = c.split("_vs_")[0]
    if project == "PRJNA850121":
        m = re.search(r"(?:^|_)([a-e])$", stress)          # parasitized_A ... parasitized_E
        if m:
            return "Orobanche_" + m.group(1).upper()
    for (p, key), lab in LABELS.items():
        if p == project and key in stress:
            return lab
    # Orobanche stages are single letters, drought time points end in d
    m = re.match(r"^(?:stage_?)?([a-e])(?:_|$)", stress)
    if project == "PRJNA850121" and m:
        return "Orobanche_" + m.group(1).upper()
    m = re.match(r"^(?:drought_?|d)?(\d+)d", stress)
    if project == "PRJNA797473" and m:
        return "Drought_%sd" % m.group(1)
    return project + "__" + contrast

genes = []
for i, line in enumerate(open(aqp_f)):
    f = line.rstrip("\n").split("\t")
    if i == 0:
        h = f; continue
    genes.append(dict(zip(h, f)))
gid = [g["gene_id"] for g in genes]
info = {g["gene_id"]: g for g in genes}
print("aquaporin genes:", len(gid))

files = sorted(glob.glob(os.path.join(deseq, "*", "DEG_*_vs_*.tsv")))
files = [f for f in files if "DEG_sig_" not in os.path.basename(f)]
# the corrected designs live in <project>_fixed directories; when one exists the
# uncorrected <project> directory is ignored (Additional file 1: Table S10)
dirs = {os.path.basename(os.path.dirname(f)) for f in files}
fixed = {d[:-6] for d in dirs if d.endswith("_fixed")}
files = [f for f in files if os.path.basename(os.path.dirname(f)) not in fixed]
print("project directories used:", ", ".join(sorted({os.path.basename(os.path.dirname(f)) for f in files})))
print("contrast result files:", len(files))
table = collections.defaultdict(dict)   # gene -> label -> (lfc, padj)
labels = []
cmap = []
for f in files:
    project = os.path.basename(os.path.dirname(f)).replace("_fixed", "")
    contrast = os.path.basename(f)[4:-4]
    lab = label_for(project, contrast)
    cmap.append((f, project, contrast, lab))
    if lab in labels:
        print("WARNING: duplicate label", lab, "from", f); lab = project + "__" + contrast
    labels.append(lab)
    with open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {c: k for k, c in enumerate(hdr)}
        gcol = ci.get("gene_id", 0); lcol = ci["log2FoldChange"]; pcol = ci["padj"]
        for line in fh:
            r = line.rstrip("\n").split("\t")
            g = r[gcol]
            if g in info:
                try:
                    l = float(r[lcol]) if r[lcol] not in ("NA", "") else float("nan")
                    p = float(r[pcol]) if r[pcol] not in ("NA", "") else float("nan")
                except ValueError:
                    l = p = float("nan")
                table[g][lab] = (l, p)
labels = [l for l in ORDER if l in labels] + [l for l in labels if l not in ORDER]
with open(os.path.join(out, "contrast_map.tsv"), "w") as fh:
    fh.write("file\tproject\tcontrast\tlabel\n")
    for f, p, c, l in cmap:
        fh.write("%s\t%s\t%s\t%s\n" % (os.path.relpath(f, deseq), p, c, l))

def fmt(x):
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else ("%.6g" % x)

# Table S4
s4 = os.path.join(supp, "Table_S4_aqp_deg_all_contrasts.tsv")
with open(s4, "w") as fh:
    fh.write("gene_id\tproduct\tsubfamily\t" + "\t".join("%s_LFC\t%s_padj" % (l, l) for l in labels) + "\n")
    for g in gid:
        row = [g, info[g].get("gene_name", ""), info[g].get("subfamily", "")]
        for l in labels:
            lfc, p = table[g].get(l, (float("nan"), float("nan")))
            row += [fmt(lfc), fmt(p)]
        fh.write("\t".join(row) + "\n")

# counts, long table, broadly responsive
counts = collections.OrderedDict((l, [0, 0]) for l in labels)
long_rows = []
per_gene = collections.Counter(); per_gene_dirs = collections.defaultdict(list)
for g in gid:
    for l in labels:
        lfc, p = table[g].get(l, (float("nan"), float("nan")))
        if not math.isnan(lfc) and not math.isnan(p) and abs(lfc) > LFC_CUT and p < PADJ_CUT:
            d = "up" if lfc > 0 else "down"
            counts[l][0 if d == "up" else 1] += 1
            long_rows.append((g, info[g].get("gene_name", ""), info[g].get("subfamily", ""), l, "%.3f" % lfc, "%.3g" % p, d))
            per_gene[g] += 1; per_gene_dirs[g].append("%s:%s" % (l, d))
with open(os.path.join(out, "aquaporin_deg_counts.tsv"), "w") as fh:
    fh.write("contrast\tup\tdown\ttotal\n")
    for l, (u, d) in counts.items():
        fh.write("%s\t%d\t%d\t%d\n" % (l, u, d, u + d))
with open(os.path.join(out, "aquaporin_deg_long.tsv"), "w") as fh:
    fh.write("gene_id\tproduct\tsubfamily\tcontrast\tlog2FC\tpadj\tdirection\n")
    for r in long_rows:
        fh.write("\t".join(r) + "\n")
with open(os.path.join(out, "broadly_responsive_genes.tsv"), "w") as fh:
    fh.write("gene_id\tproduct\tsubfamily\tchromosome\tn_contrasts\tcontrasts\n")
    for g, n in sorted(per_gene.items(), key=lambda kv: (-kv[1], kv[0])):
        fh.write("%s\t%s\t%s\t%s\t%d\t%s\n" % (g, info[g].get("gene_name", ""), info[g].get("subfamily", ""), info[g].get("chromosome", ""), n, ";".join(per_gene_dirs[g])))

de_any = sum(1 for g in gid if per_gene[g] > 0)
print("contrasts:", len(labels), "->", ", ".join(labels))
print("genes DE in at least one contrast: %d / %d (%.0f%%)" % (de_any, len(gid), 100.0 * de_any / len(gid)))
for l, (u, d) in counts.items():
    print("  %-14s up %2d  down %2d  total %2d" % (l, u, d, u + d))
print("genes DE in >= 5 contrasts:", sum(1 for g in gid if per_gene[g] >= 5))
print("top genes:", ", ".join("%s(%d)" % (g, n) for g, n in sorted(per_gene.items(), key=lambda kv: -kv[1])[:8]))
missing = [g for g in gid if not table[g]]
print("aquaporin genes absent from every DESeq2 result file:", len(missing), missing[:10])
PYEOF

echo "Outputs in ${OUT} and ${SUPP}/Table_S4_aqp_deg_all_contrasts.tsv"
echo "Finished: $(date)"
