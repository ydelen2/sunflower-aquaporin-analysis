# -*- coding: utf-8 -*-
"""
Synteny support for the paralogous aquaporin gene pairs.

The duplication table (09_duplication_kaks_v3/duplication_kaks.tsv) classes a
pair by genomic distance only: tandem (same chromosome, gene starts <= 200 kb
apart), proximal (more than 200 kb and up to 1 Mb) and "segmental" (different chromosomes, or > 1 Mb apart).
Distance alone does not show that a pair arose by segmental or whole-genome
duplication, so the pairs of the third class are split here by the MCScanX
intra-genomic collinearity result:

  segmental_anchor    the two genes are an anchor pair of a collinear block
  segmental_in_block  the two genes are not an anchor pair, but lie at
                      corresponding positions of the two regions of one
                      collinear block (between matching flanking anchors,
                      with a tolerance of TOL anchors on each side)
  dispersed           no collinear block connects the two genes

Tandem and proximal pairs keep their class.

Inputs
  duplication_kaks.tsv            pair table with Ks
  Ha_self_v3.collinearity         MCScanX output (sunflower against itself)
  07_synteny/data/Ha.gff          MCScanX gene positions (one protein per gene)
  Ha_prot_to_gene.tsv             protein accession -> gene
Output
  duplication_synteny_class.tsv   the pair table with the columns
                                  synteny_class, block_id, block_anchors

usage: python 10_duplication_synteny_class_v3.py [dup.tsv] [collinearity] [gff] [prot_to_gene] [out.tsv]
"""
import bisect, collections, csv, os, re, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.environ.get("AQP_RESULTS", os.path.join(ROOT, "results"))
dup_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RES, "02_gene_family", "09_duplication_kaks_v3", "duplication_kaks.tsv")
col_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(RES, "07_synteny", "mcscanx", "Ha_self_v3", "Ha_self_v3.collinearity")
gff_path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(RES, "07_synteny", "data", "Ha.gff")
p2g_path = sys.argv[4] if len(sys.argv) > 4 else os.path.join(RES, "07_synteny", "data", "Ha_prot_to_gene.tsv")
out_path = sys.argv[5] if len(sys.argv) > 5 else os.path.join(RES, "02_gene_family", "10_duplication_synteny_v3", "duplication_synteny_class.tsv")
TOL = 2   # anchors of tolerance on each side when matching positions inside a block

# protein -> gene, gene -> position (the MCScanX gff holds one protein per gene)
p2g = {}
for line in open(p2g_path, encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 2:
        p2g[f[0]] = re.sub(r"^gene-", "", f[1])
gene_pos = {}
prot_pos = {}
for line in open(gff_path):
    c, p, s, e = line.split()[:4]
    prot_pos[p] = (c, (int(s) + int(e)) // 2)
    gene_pos[p2g.get(p, p)] = (c, (int(s) + int(e)) // 2)

# collinear blocks
blocks = []
cur = None
for line in open(col_path):
    if line.startswith("## Alignment"):
        m = re.search(r"Alignment (\d+): score=([\d.]+) e_value=(\S+) N=(\d+) (\S+)&(\S+) (\w+)", line)
        cur = {"id": int(m.group(1)), "n": int(m.group(4)), "c1": m.group(5), "c2": m.group(6), "orient": m.group(7), "pairs": []}
        blocks.append(cur)
    elif not line.startswith("#") and cur is not None:
        f = line.split("\t")
        cur["pairs"].append((f[1].strip(), f[2].strip()))
anchor = {}
for b in blocks:
    b["pos1"] = [prot_pos[a][1] for a, c in b["pairs"]]
    b["pos2"] = [prot_pos[c][1] for a, c in b["pairs"]]
    order = sorted(range(len(b["pairs"])), key=lambda i: b["pos1"][i])
    b["s1"] = [b["pos1"][i] for i in order]
    b["s2"] = [b["pos2"][i] for i in order]
    for a, c in b["pairs"]:
        anchor[frozenset((p2g.get(a, a), p2g.get(c, c)))] = b


def in_block(b, x_chr, x_pos, y_chr, y_pos):
    """x in region 1 and y in region 2 of block b, at corresponding positions."""
    if x_chr != b["c1"] or y_chr != b["c2"]:
        return False
    if not (b["s1"][0] <= x_pos <= b["s1"][-1]):
        return False
    i = bisect.bisect_left(b["s1"], x_pos)
    lo, hi = max(0, i - 1 - TOL), min(len(b["s1"]) - 1, i + TOL)
    seg = b["s2"][lo:hi + 1]
    return min(seg) <= y_pos <= max(seg)


def classify(g1, g2):
    key = frozenset((g1, g2))
    if key in anchor:
        return "segmental_anchor", anchor[key]
    if g1 not in gene_pos or g2 not in gene_pos:
        return "dispersed", None
    (c1, x), (c2, y) = gene_pos[g1], gene_pos[g2]
    hits = [b for b in blocks if in_block(b, c1, x, c2, y) or in_block(b, c2, y, c1, x)]
    if hits:
        return "segmental_in_block", max(hits, key=lambda b: b["n"])
    return "dispersed", None


rows = list(csv.DictReader(open(dup_path, encoding="utf-8"), delimiter="\t"))
hdr = list(rows[0].keys())
count = collections.Counter()
ks = collections.defaultdict(list)
sat = collections.Counter()
genes_by_class = collections.defaultdict(set)
for r in rows:
    if r["type"] == "segmental":
        cls, b = classify(r["gene1"], r["gene2"])
    else:
        cls, b = r["type"], None
    r["synteny_class"] = cls
    r["block_id"] = b["id"] if b else "NA"
    r["block_anchors"] = b["n"] if b else "NA"
    count[cls] += 1
    genes_by_class[cls].update((r["gene1"], r["gene2"]))
    try:
        ks[cls].append(float(r["Ks"]))
    except ValueError:
        if "saturat" in r.get("note", "").lower():
            sat[cls] += 1

os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
    cols = hdr + ["synteny_class", "block_id", "block_anchors"]
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

n = len(rows)
all_genes = {r["gene1"] for r in rows} | {r["gene2"] for r in rows}
print("pairs:", n, "| genes in pairs:", len(all_genes), "| collinear blocks:", len(blocks), "| anchor pairs:", sum(len(b["pairs"]) for b in blocks))
print("%-20s %5s %6s %10s %8s %8s %9s %6s" % ("class", "pairs", "pct", "Ks_median", "Ks_q1", "Ks_q3", "saturated", "genes"))
for cls in ("tandem", "proximal", "segmental_anchor", "segmental_in_block", "dispersed"):
    v = sorted(ks[cls])
    q = statistics.quantiles(v, n=4) if len(v) >= 4 else [float("nan")] * 3
    print("%-20s %5d %5.1f%% %10.2f %8.2f %8.2f %9d %6d" % (cls, count[cls], 100.0 * count[cls] / n, statistics.median(v) if v else float("nan"),
                                                           q[0], q[2], sat[cls], len(genes_by_class[cls])))
syn = ks["segmental_anchor"] + ks["segmental_in_block"]
syn_n = count["segmental_anchor"] + count["segmental_in_block"]
syn_genes = genes_by_class["segmental_anchor"] | genes_by_class["segmental_in_block"]
print("synteny-supported segmental pairs: %d (%.1f%%), median Ks %.2f, genes %d" % (syn_n, 100.0 * syn_n / n, statistics.median(syn), len(syn_genes)))
print("written", out_path)
