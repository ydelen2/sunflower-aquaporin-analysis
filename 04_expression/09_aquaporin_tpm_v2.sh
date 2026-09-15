#!/bin/bash
#SBATCH --job-name=aqp_tpm_v2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/aqp_tpm_v2_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/04_expression/logs/aqp_tpm_v2_%j.err

# =============================================================================
# 09_aquaporin_tpm_v2.sh
# TPM of the aquaporin genes in every library of the six BioProjects (249
# libraries), computed from the featureCounts matrices: the main matrix
# (gene_counts_clean.txt) for five BioProjects and the 96-run matrix of
# PRJNA492303 (PRJNA492303_full_counts_clean.txt). TPM uses the gene lengths
# reported by featureCounts (sum of exon lengths) and is normalised per library
# over all annotated genes. Control libraries are then summarised per
# BioProject and tissue to give the baseline expression of each gene.
#
# Inputs
#   03_rnaseq/counts/gene_counts_clean.txt (counts); gene lengths from one per-run featureCounts
#   table of the same GTF (03_rnaseq/counts/PRJNA492303_fc/SRR7882942.txt, Length column)
#   03_rnaseq/counts/PRJNA492303_full_counts_clean.txt
#   04_expression/results/deseq2/sample_metadata.tsv
#   03_rnaseq/accession_lists/PRJNA492303_metadata_full.tsv
#   02_gene_family/results/verified_aquaporins.txt
# Outputs (04_expression/results/aquaporin_expression_v2/)
#   aquaporin_TPM_all_libraries.tsv   gene x library TPM (87 x 249)
#   library_metadata.tsv              BioProject, condition, tissue, genotype, age, control flag
#   aquaporin_TPM_control_means.tsv   mean TPM of control libraries per BioProject and tissue
#   tissue_summary.tsv                leaf/root comparison in the two BioProjects with both tissues
# =============================================================================
set -eo pipefail
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
mkdir -p "${PROJ_DIR}/04_expression/logs" "${PROJ_DIR}/04_expression/results/aquaporin_expression_v2"

python3 - "${PROJ_DIR}" <<'PY'
import sys, os, math, collections
P = sys.argv[1]
main_file = os.path.join(P, "03_rnaseq/counts/gene_counts_clean.txt")
raw_file  = os.path.join(P, "03_rnaseq/counts/PRJNA492303_fc/SRR7882942.txt")
fl_file   = os.path.join(P, "03_rnaseq/counts/PRJNA492303_full_counts_clean.txt")
meta_file = os.path.join(P, "04_expression/results/deseq2/sample_metadata.tsv")
fl_meta   = os.path.join(P, "03_rnaseq/accession_lists/PRJNA492303_metadata_full.tsv")
aqp_file  = os.path.join(P, "02_gene_family/results/verified_aquaporins.txt")
out_dir   = os.path.join(P, "04_expression/results/aquaporin_expression_v2")

def read_tsv(path, comment="#"):
    with open(path) as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip() and not l.startswith(comment)]
    return rows[0], rows[1:]

# gene lengths from a featureCounts table (Geneid Chr Start End Strand Length <sample>)
hdr, rows = read_tsv(raw_file)
li = hdr.index("Length")
length = {r[0]: int(r[li]) for r in rows}
# main count matrix with cleaned sample names (Geneid <samples>)
hdr, rows = read_tsv(main_file)
samples_main = hdr[1:]
counts = {s: {} for s in samples_main}
for r in rows:
    for j, s in enumerate(samples_main):
        counts[s][r[0]] = float(r[1 + j])
print("main matrix:", len(rows), "genes x", len(samples_main), "libraries; genes with a length:", len(length))

# metadata of the main matrix
h, m = read_tsv(meta_file)
col = {c: i for i, c in enumerate(h)}
meta = {}
for r in m:
    meta[r[col["sample_id"]]] = dict(project=r[col["project_id"]], condition=r[col["condition"]],
                                     tissue=r[col["tissue"]], genotype=r[col["genotype"]], age="")

# replace the PRJNA492303 subset by the 96-run matrix
old_flood = [s for s in samples_main if meta.get(s, {}).get("project") == "PRJNA492303"]
for s in old_flood:
    del counts[s]; del meta[s]
hdr2, rows2 = read_tsv(fl_file)
samples_fl = hdr2[1:]
for s in samples_fl:
    counts[s] = {}
for r in rows2:
    for j, s in enumerate(samples_fl):
        counts[s][r[0]] = float(r[1 + j])
h, m = read_tsv(fl_meta)
col = {c: i for i, c in enumerate(h)}
for r in m:
    meta[r[col["sample_id"]]] = dict(project="PRJNA492303", condition=r[col["condition"]], tissue=r[col["tissue"]],
                                     genotype=r[col["genotype"]], age=r[col["age"]])
genes_fl = set(r[0] for r in rows2)
missing = [g for g in length if g not in genes_fl]
print("PRJNA492303 matrix:", len(rows2), "genes x", len(samples_fl), "libraries; genes absent from it:", len(missing))
print("old PRJNA492303 libraries removed:", len(old_flood), "| libraries now:", len(counts))

# TPM per library over all genes with a length
old_set = set(old_flood)
samples = [s for s in samples_main if s not in old_set] + samples_fl
assert len(samples) == len(set(samples)) == len(counts), "duplicate library ids"
print("libraries in the TPM table:", len(samples))
tpm = {}
for s in samples:
    c = counts[s]
    rpk = {g: c.get(g, 0.0) / (length[g] / 1000.0) for g in length}
    tot = sum(rpk.values())
    tpm[s] = {g: v / tot * 1e6 for g, v in rpk.items()}

# aquaporin genes
h, m = read_tsv(aqp_file)
gi = h.index("gene_id") if "gene_id" in h else 0
aqp = [r[gi] for r in m]
absent = [g for g in aqp if g not in length]
print("aquaporin genes:", len(aqp), "| not in the count matrix:", absent)
aqp = [g for g in aqp if g in length]

# control flag: pre-treatment and untreated libraries
CONTROL = {"control", "cold_0h", "heat_0h", "drought_normal"}
for s in samples:
    meta[s]["control"] = "yes" if meta[s]["condition"] in CONTROL else "no"

with open(os.path.join(out_dir, "aquaporin_TPM_all_libraries.tsv"), "w") as fh:
    fh.write("gene_id\t" + "\t".join(samples) + "\n")
    for g in aqp:
        fh.write(g + "\t" + "\t".join("%.3f" % tpm[s][g] for s in samples) + "\n")
with open(os.path.join(out_dir, "library_metadata.tsv"), "w") as fh:
    fh.write("sample_id\tproject\tcondition\ttissue\tgenotype\tage\tcontrol\n")
    for s in samples:
        d = meta[s]
        fh.write("\t".join([s, d["project"], d["condition"], d["tissue"], d["genotype"], d["age"], d["control"]]) + "\n")

# control means per BioProject and tissue
groups = collections.OrderedDict()
for s in samples:
    d = meta[s]
    if d["control"] == "yes":
        groups.setdefault((d["project"], d["tissue"]), []).append(s)
keys = list(groups)
with open(os.path.join(out_dir, "aquaporin_TPM_control_means.tsv"), "w") as fh:
    fh.write("gene_id\t" + "\t".join("%s_%s" % k for k in keys) + "\n")
    fh.write("n_libraries\t" + "\t".join(str(len(groups[k])) for k in keys) + "\n")
    for g in aqp:
        fh.write(g + "\t" + "\t".join("%.3f" % (sum(tpm[s][g] for s in groups[k]) / len(groups[k])) for k in keys) + "\n")
print("control groups:", ["%s_%s n=%d" % (k[0], k[1], len(v)) for k, v in groups.items()])

# leaf versus root in the two BioProjects that sampled both tissues (control libraries)
with open(os.path.join(out_dir, "tissue_summary.tsv"), "w") as fh:
    fh.write("gene_id\tproject\tleaf_mean_TPM\troot_mean_TPM\tlog2_root_over_leaf\tn_leaf\tn_root\n")
    for proj in ("PRJNA1041959", "PRJNA492303"):
        lf = groups.get((proj, "leaf"), []); rt = groups.get((proj, "root"), [])
        for g in aqp:
            a = sum(tpm[s][g] for s in lf) / len(lf); b = sum(tpm[s][g] for s in rt) / len(rt)
            fh.write("%s\t%s\t%.3f\t%.3f\t%.3f\t%d\t%d\n" % (g, proj, a, b, math.log2((b + 0.1) / (a + 0.1)), len(lf), len(rt)))
print("written to", out_dir)
PY
