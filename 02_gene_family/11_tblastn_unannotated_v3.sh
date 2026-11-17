#!/bin/bash
#SBATCH --job-name=aqp_tblastn
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/11_tblastn_v3_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/11_tblastn_v3_%j.err

# =============================================================================
# 11_tblastn_unannotated_v3.sh
# Checks whether the genome sequence holds aquaporin-like loci that the RefSeq
# annotation does not model. The 99 aquaporin proteins of the inventory are
# searched against the genome with tblastn; the hits are merged into loci and
# compared with the annotated gene models.
#
# A locus is a set of tblastn HSPs on one strand of one sequence that lie
# within 10 kb of each other. For each locus the best query coverage (union of
# the HSPs of one query) is reported, and the locus is classified as
#   inventory     overlaps one of the 87 aquaporin genes
#   annotated     overlaps another annotated gene or pseudogene
#   unannotated   overlaps no annotated gene or pseudogene
#
# Inputs
#   01_references/sunflower/genomic.fna, genomic.gff
#   02_gene_family/verified_aquaporins.faa, results/verified_aquaporins.txt
# Outputs (02_gene_family/11_tblastn_v3/)
#   tblastn_hits.tsv          raw HSPs (outfmt 6 with qlen)
#   tblastn_loci.tsv          merged loci with class and overlapping genes
#   tblastn_loci_summary.tsv  counts per class and coverage bin
# =============================================================================
set -eo pipefail
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
OUT="${PROJ_DIR}/02_gene_family/11_tblastn_v3"
GENOME="${PROJ_DIR}/01_references/sunflower/genomic.fna"
GFF="${PROJ_DIR}/01_references/sunflower/genomic.gff"
QUERY="${PROJ_DIR}/02_gene_family/verified_aquaporins.faa"
GENES="${PROJ_DIR}/02_gene_family/results/verified_aquaporins.txt"
DB="${OUT}/genome_db/HanXRQr2"
mkdir -p "${OUT}/genome_db" "${PROJ_DIR}/logs"

module purge
module load blast/2.17

for f in "${GENOME}" "${GFF}" "${QUERY}" "${GENES}"; do
    [[ -s "${f}" ]] || { echo "missing input: ${f}" >&2; exit 1; }
done
echo "query proteins: $(grep -c '^>' "${QUERY}")"

if [[ ! -s "${DB}.nsq" && ! -s "${DB}.00.nsq" ]]; then
    echo "building nucleotide database: $(date)"
    makeblastdb -in "${GENOME}" -dbtype nucl -out "${DB}" -parse_seqids
fi

echo "tblastn: $(date)"
tblastn -query "${QUERY}" -db "${DB}" -evalue 1e-10 -num_threads 8 -max_target_seqs 5000 \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen" \
    > "${OUT}/tblastn_hits.tsv"
echo "HSPs: $(wc -l < "${OUT}/tblastn_hits.tsv")  $(date)"

python3 - "${OUT}" "${GFF}" "${GENES}" <<'PY'
import sys, os, re, collections
out, gff, genes_file = sys.argv[1:4]
GAP = 10000

# inventory genes
inv = set()
with open(genes_file) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    gi = hdr.index("gene_id") if "gene_id" in hdr else 0
    for line in fh:
        if line.strip():
            inv.add(line.rstrip("\n").split("\t")[gi])
print("inventory genes:", len(inv))

# annotated genes and pseudogenes
ann = collections.defaultdict(list)
n_gene = 0
with open(gff) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
            continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        gid = a.get("Name") or a.get("gene") or a.get("ID", "")
        biotype = a.get("gene_biotype", f[2])
        desc = a.get("description", "").replace("%2C", ",").replace("%3B", ";")
        ann[f[0]].append((int(f[3]), int(f[4]), gid, biotype, desc))
        n_gene += 1
for c in ann:
    ann[c].sort()
print("annotated gene and pseudogene features:", n_gene)

# HSPs
hsps = collections.defaultdict(list)
n = 0
with open(os.path.join(out, "tblastn_hits.tsv")) as fh:
    for line in fh:
        q, s, pid, ln, qs, qe, ss, se, ev, bs, ql = line.rstrip("\n").split("\t")
        if "|" in s:
            s = [x for x in s.split("|") if x][-1]
        ss, se = int(ss), int(se)
        strand = "+" if ss <= se else "-"
        lo, hi = min(ss, se), max(ss, se)
        hsps[(s, strand)].append((lo, hi, q, int(qs), int(qe), int(ql), float(pid), float(bs)))
        n += 1
print("HSPs read:", n)

def union_len(iv):
    iv = sorted(iv); tot = 0; cs, ce = iv[0]
    for a, b in iv[1:]:
        if a <= ce + 1:
            ce = max(ce, b)
        else:
            tot += ce - cs + 1; cs, ce = a, b
    return tot + ce - cs + 1

loci = []
for (chrom, strand), hs in hsps.items():
    hs.sort()
    cur = [hs[0]]; cur_end = hs[0][1]
    groups = []
    for h in hs[1:]:
        if h[0] - cur_end <= GAP:
            cur.append(h); cur_end = max(cur_end, h[1])
        else:
            groups.append(cur); cur = [h]; cur_end = h[1]
    groups.append(cur)
    for g in groups:
        start = min(h[0] for h in g); end = max(h[1] for h in g)
        byq = collections.defaultdict(list); ql = {}; ids = collections.defaultdict(list)
        for lo, hi, q, qs, qe, qlen, pid, bs in g:
            byq[q].append((qs, qe)); ql[q] = qlen; ids[q].append((pid, qe - qs + 1))
        best_q, best_cov = None, -1.0
        for q, iv in byq.items():
            cov = union_len(iv) / float(ql[q])
            if cov > best_cov:
                best_q, best_cov = q, cov
        wid = sum(p * l for p, l in ids[best_q]) / float(sum(l for p, l in ids[best_q]))
        ov = [(a, b, gid, bt, d) for a, b, gid, bt, d in ann.get(chrom, []) if a <= end and b >= start]
        inv_hit = [x for x in ov if x[2] in inv]
        if inv_hit:
            cls = "inventory"
        elif ov:
            cls = "annotated"
        else:
            cls = "unannotated"
        loci.append((chrom, start, end, strand, len(g), best_q, round(best_cov, 3), round(wid, 1), max(h[7] for h in g), cls,
                     ";".join("%s|%s|%s" % (x[2], x[3], x[4]) for x in ov)))

loci.sort(key=lambda x: (x[0], x[1]))
with open(os.path.join(out, "tblastn_loci.tsv"), "w") as fh:
    fh.write("sequence\tstart\tend\tstrand\tn_hsps\tbest_query\tbest_query_coverage\tidentity_of_best_query_pct\tmax_bitscore\tclass\toverlapping_genes\n")
    for l in loci:
        fh.write("\t".join(str(x) for x in l) + "\n")

bins = [("coverage_ge_0.8", 0.8, 2.0), ("coverage_0.5_to_0.8", 0.5, 0.8), ("coverage_0.2_to_0.5", 0.2, 0.5), ("coverage_lt_0.2", -1.0, 0.2)]
summ = collections.OrderedDict()
for cls in ("inventory", "annotated", "unannotated"):
    for name, lo, hi in bins:
        summ[(cls, name)] = sum(1 for l in loci if l[9] == cls and lo <= l[6] < hi)
with open(os.path.join(out, "tblastn_loci_summary.tsv"), "w") as fh:
    fh.write("class\tcoverage_bin\tloci\n")
    for (cls, name), v in summ.items():
        fh.write("%s\t%s\t%d\n" % (cls, name, v))
        print("%-12s %-22s %d" % (cls, name, v))
inv_found = set()
for l in loci:
    for part in l[10].split(";"):
        g = part.split("|")[0]
        if g in inv:
            inv_found.add(g)
print("inventory genes overlapped by at least one locus:", len(inv_found), "of", len(inv))
print("loci:", len(loci))
PY
echo "=== 11_tblastn_unannotated_v3 finished: $(date) ==="
