# -*- coding: utf-8 -*-
"""
Aquaporin gene pairs in the MCScanX collinearity files of the runs with one protein
per gene (05_synteny_one_protein_per_gene_v3.sh).

  intra-genomic (Ha_self_v3)   a pair is kept when both genes are in the aquaporin
                               inventory
  interspecies (Ha_Ls_v3,      a pair is kept when the sunflower gene is in the
  Ha_At_v3)                    inventory; the column partner_in_reference_set says
                               whether the partner protein is one of the reference
                               aquaporins of that species (PF00230 search,
                               05_phylogenetics/results_v2)

Each gene pair is written once, with the first collinear block in which it occurs.

Inputs
  07_synteny/mcscanx/<run>/<run>.collinearity
  07_synteny/data/<species>_prot_to_gene.tsv
  09_figures/supplementary/Table_S1_aquaporin_genes.tsv
  05_phylogenetics/results_v2/At_hits.txt, Ls_hits.txt
Outputs (07_synteny/results_v3/)
  Ha_self_intra_aqp_synteny.tsv, Ha_Ls_inter_aqp_synteny.tsv, Ha_At_inter_aqp_synteny.tsv
  collinearity_summary.tsv      blocks, anchor pairs and aquaporin numbers of each run

usage: python 06_reparse_synteny_v3.py
"""
import csv, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.environ.get("AQP_RESULTS", os.path.join(ROOT, "results"))
SYN = os.path.join(RES, "07_synteny")
OUT = os.path.join(SYN, "results_v3")
os.makedirs(OUT, exist_ok=True)

with open(os.path.join(RES, "09_figures", "supplementary", "Table_S1_aquaporin_genes.tsv"), encoding="utf-8") as fh:
    inventory = {r["loc_id"] for r in csv.DictReader(fh, delimiter="\t")}
p2g = {}
for sp in ("Ha", "At", "Ls"):
    for line in open(os.path.join(SYN, "data", sp + "_prot_to_gene.tsv"), encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) == 2:
            p2g[f[0]] = (sp, re.sub(r"^gene-", "", f[1]))
ref = {}
for sp in ("At", "Ls"):
    ref[sp] = {l.strip() for l in open(os.path.join(RES, "05_phylogenetics", "results_v2", sp + "_hits.txt")) if l.strip()}


def blocks(run):
    out, cur, stats = [], None, {}
    for line in open(os.path.join(SYN, "mcscanx", run, run + ".collinearity")):
        m = re.match(r"# Number of collinear genes: (\d+), Percentage: ([\d.]+)", line)
        if m:
            stats["collinear_genes"], stats["percent"] = m.group(1), m.group(2)
        m = re.match(r"# Number of all genes: (\d+)", line)
        if m:
            stats["all_genes"] = m.group(1)
        if line.startswith("## Alignment"):
            m = re.search(r"Alignment (\d+): .* N=(\d+) (\S+)&(\S+) ", line)
            cur = {"id": m.group(1), "n": int(m.group(2)), "c1": m.group(3), "c2": m.group(4), "pairs": []}
            out.append(cur)
        elif not line.startswith("#") and cur is not None:
            f = line.split("\t")
            cur["pairs"].append((f[1].strip(), f[2].strip()))
    return out, stats


HEADER = ["chr1", "chr2", "species1", "gene1", "prot1", "species2", "gene2", "prot2", "partner_in_reference_set", "block_id", "block_anchors"]
summary = [["run", "all_genes", "collinear_genes", "percent_collinear", "blocks", "blocks_between_genomes", "anchor_pairs", "aquaporin_pairs", "sunflower_aquaporin_genes", "partner_genes", "largest_block_anchors"]]
for run, other, name in (("Ha_self_v3", None, "Ha_self_intra_aqp_synteny.tsv"), ("Ha_Ls_v3", "Ls", "Ha_Ls_inter_aqp_synteny.tsv"), ("Ha_At_v3", "At", "Ha_At_inter_aqp_synteny.tsv")):
    bl, stats = blocks(run)
    rows, seen = [], set()
    for b in bl:
        between = b["c1"][:2] != b["c2"][:2]
        if (other is None) == between:
            continue
        for x, y in b["pairs"]:
            (sx, gx), (sy, gy) = p2g.get(x, ("?", x)), p2g.get(y, ("?", y))
            if other is None:
                if not (gx in inventory and gy in inventory):
                    continue
                flag = "NA"
            else:
                ha_gene, partner = (gx, y) if sx == "Ha" else (gy, x)
                if ha_gene not in inventory:
                    continue
                flag = "yes" if partner in ref[other] else "no"
            key = frozenset((gx, gy))
            if key in seen:
                continue
            seen.add(key)
            rows.append([b["c1"], b["c2"], sx, gx, x, sy, gy, y, flag, b["id"], str(b["n"])])
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    ha = {r[3] if r[2] == "Ha" else r[6] for r in rows} | ({r[6] for r in rows} if other is None else set())
    partners = set() if other is None else {r[6] if r[2] == "Ha" else r[3] for r in rows}
    n_between = sum(b["c1"][:2] != b["c2"][:2] for b in bl)
    summary.append([run, stats.get("all_genes", ""), stats.get("collinear_genes", ""), stats.get("percent", ""), str(len(bl)), str(n_between), str(sum(len(b["pairs"]) for b in bl)),
                    str(len(rows)), str(len(ha)), str(len(partners)), str(max(b["n"] for b in bl))])
    print("%-11s blocks %d (between genomes %d) | aquaporin pairs %d | sunflower aquaporin genes %d | partner genes %d | partners in the reference set %s" % (
        run, len(bl), n_between, len(rows), len(ha), len(partners), "NA" if other is None else "%d of %d" % (sum(r[8] == "yes" for r in rows), len(rows))))
with open(os.path.join(OUT, "collinearity_summary.tsv"), "w", encoding="utf-8", newline="\n") as fh:
    for r in summary:
        fh.write("\t".join(r) + "\n")
print("written", OUT)
