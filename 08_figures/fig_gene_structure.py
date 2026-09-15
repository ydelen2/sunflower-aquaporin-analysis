# -*- coding: utf-8 -*-
"""
Gene structure and conserved protein motifs of the sunflower aquaporins
(Additional file figure). Left: exon/intron structure of the representative
mRNA of each gene (CDS filled, UTR open, introns as lines), 5' to 3'.
Right: MEME motifs along the representative protein.

Inputs: 02_gene_family/05_gene_structure/aqp_genes.gff3, Table S1
(representative protein per gene, subfamily), MEME meme.xml
(02_gene_family/results/meme_protein or 05_phylogenetics/results/meme_motifs/meme_output).

usage: python fig_gene_structure.py [out_prefix]
"""
import os, re, sys, collections
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

import figlib as L

out_prefix = sys.argv[1] if len(sys.argv) > 1 else L.fig_path("FigS1_gene_structure_motifs")
s1 = L.load_table_s1()
p2g = L.load_protein_gene()
products = L.load_products()

# ---- GFF3: gene -> mRNAs -> exons/CDS; mRNA -> protein via CDS ID
gff = os.path.join(L.RESULTS, "02_gene_family", "05_gene_structure", "aqp_genes.gff3")
genes, mrna_gene, mrna_prot, exons, cdss, strand_of = {}, {}, {}, collections.defaultdict(list), collections.defaultdict(list), {}
for line in open(gff):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    typ, start, end, strand, attr = f[2], int(f[3]), int(f[4]), f[6], f[8]
    a = dict(kv.split("=", 1) for kv in attr.split(";") if "=" in kv)
    if typ == "gene":
        g = re.sub(r"^gene-", "", a["ID"]); genes[g] = (start, end); strand_of[g] = strand
    elif typ == "mRNA":
        mrna_gene[a["ID"]] = re.sub(r"^gene-", "", a["Parent"])
    elif typ == "exon":
        exons[a["Parent"]].append((start, end))
    elif typ == "CDS":
        cdss[a["Parent"]].append((start, end))
        mrna_prot[a["Parent"]] = re.sub(r"^cds-", "", a["ID"])

# representative protein per gene: Table S1 column when present, else longest CDS
rep = {}
for g, r in s1.items():
    if r.get("representative_protein"):
        rep[g] = r["representative_protein"]
prot_mrna = {p: m for m, p in mrna_prot.items()}
for g in s1:
    if g not in rep:
        cands = [(sum(e - s + 1 for s, e in cdss[m]), m) for m in mrna_gene if mrna_gene[m] == g]
        if cands:
            rep[g] = mrna_prot[max(cands)[1]]

# ---- MEME motifs
meme_paths = [os.path.join(L.RESULTS, "05_phylogenetics", "results", "meme_motifs", "meme_output", "meme.xml"),
              os.path.join(L.RESULTS, "02_gene_family", "results", "meme_protein", "meme.xml")]
meme = next((p for p in meme_paths if os.path.exists(p)), None)
sites = collections.defaultdict(list)   # protein -> [(motif_no, start, width)]
motif_w, motif_cons = {}, {}
if meme:
    root = ET.parse(meme).getroot()
    seq_name = {s.get("id"): s.get("name") for s in root.iter("sequence")}
    for m in root.iter("motif"):
        no = int(m.get("id").split("_")[1])
        if no <= 10:
            motif_w[no] = int(m.get("width")); motif_cons[no] = m.get("name")
    for ss in root.iter("scanned_sites"):
        name = seq_name.get(ss.get("sequence_id"))
        for s in ss.iter("scanned_site"):
            if float(s.get("pvalue")) < 1e-5 and int(s.get("motif_id").split("_")[1]) <= 10:
                sites[name].append((int(s.get("motif_id").split("_")[1]), int(s.get("position")), motif_w[int(s.get("motif_id").split("_")[1])]))
prot_len = {}
for line in open(os.path.join(L.RESULTS, "02_gene_family", "03_domain_verification", "tm_helix_predictions.tsv")):
    f = line.rstrip("\n").split("\t")
    if f[0] != "protein_id" and len(f) > 3:
        prot_len[f[0]] = int(f[3])
for pth in [os.path.join(L.RESULTS, "02_gene_family", "03_domain_verification_v2", "tm_helix_predictions.tsv")]:
    if os.path.exists(pth):
        for line in open(pth):
            f = line.rstrip("\n").split("\t")
            if f[0] != "protein_id" and len(f) > 3:
                prot_len[f[0]] = int(f[3])

# ---- rows ordered by subfamily then name
order = {"PIP": 0, "TIP": 1, "NIP": 2, "SIP": 3, "XIP": 4}
rows = sorted(s1.values(), key=lambda r: (order.get(r["subfamily"], 9), L.short_name(r.get("product", "")), r["loc_id"]))
n = len(rows)
motif_colors = plt.get_cmap("tab10")
fig_h = max(6, 0.17 * n + 1.5)
fig, (axg, axp) = plt.subplots(1, 2, figsize=(14, fig_h), gridspec_kw={"width_ratios": [1.25, 1.0]})
glens = sorted(genes[g][1] - genes[g][0] + 1 for g in s1 if g in genes)
max_glen = min(glens[-1], int(glens[int(0.95 * (len(glens) - 1))] * 1.15))   # clip the longest few genes
clipped = []
max_plen = max(prot_len.get(rep.get(g, ""), 300) for g in s1)

for i, r in enumerate(rows):
    g = r["loc_id"]; y = n - i
    col = L.SUBFAM_COLORS.get(r["subfamily"], "#7f7f7f")
    label = "%s %s" % (L.short_name(r.get("product", "")), g)
    axg.text(-0.01 * max_glen, y, label, ha="right", va="center", fontsize=4.6, color=col)
    m = prot_mrna.get(rep.get(g, ""))
    if g in genes and m:
        gs, ge = genes[g]; strand = strand_of[g]
        def rel(x):   # 5' -> 3' coordinate relative to gene start
            return (x - gs) if strand == "+" else (ge - x)
        ex = sorted((min(rel(s), rel(e)), max(rel(s), rel(e))) for s, e in exons[m])
        cd = sorted((min(rel(s), rel(e)), max(rel(s), rel(e))) for s, e in cdss[m])
        if ex[-1][1] > max_glen:
            clipped.append((g, ex[-1][1]))
            axg.text(max_glen * 1.005, y, "%.0f kb" % (ex[-1][1] / 1000.0), fontsize=4.2, va="center", color="#555555")
        axg.plot([ex[0][0], min(ex[-1][1], max_glen)], [y, y], color="#888888", lw=0.6, zorder=1)
        for s, e in ex:
            if s < max_glen:
                axg.add_patch(Rectangle((s, y - 0.28), min(e, max_glen) - s + 1, 0.56, facecolor="white", edgecolor=col, lw=0.5, zorder=2))
        for s, e in cd:
            if s < max_glen:
                axg.add_patch(Rectangle((s, y - 0.28), min(e, max_glen) - s + 1, 0.56, facecolor=col, edgecolor=col, lw=0.5, zorder=3))
    p = rep.get(g, "")
    pl = prot_len.get(p, 0)
    if pl:
        axp.plot([0, pl], [y, y], color="#888888", lw=0.8, zorder=1)
        for no, start, w in sorted(sites.get(p, [])):
            axp.add_patch(Rectangle((start, y - 0.3), w, 0.6, facecolor=motif_colors((no - 1) % 10), edgecolor="none", zorder=2))

for ax in (axg, axp):
    ax.set_ylim(0, n + 1); ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
axg.set_xlim(-0.62 * max_glen, max_glen * 1.12)
axg.set_xlabel("Gene length (kb), 5' to 3'", fontsize=8)
step = 2000 if max_glen <= 14000 else 5000
axg.set_xticks([k for k in range(0, int(max_glen) + 1, step)]); axg.set_xticklabels(["%d" % (k / 1000) for k in range(0, int(max_glen) + 1, step)], fontsize=7)
axg.set_title("A", fontsize=10, loc="left", fontweight="bold")
axp.set_xlim(-5, max_plen * 1.02); axp.set_xlabel("Protein length (aa)", fontsize=8); axp.tick_params(labelsize=7)
axp.set_title("B", fontsize=10, loc="left", fontweight="bold")
h1 = [Patch(facecolor=L.SUBFAM_COLORS[s], label=s) for s in ["PIP", "TIP", "NIP", "SIP", "XIP"] if any(r["subfamily"] == s for r in rows)]
h1 += [Patch(facecolor="white", edgecolor="black", label="UTR"), Patch(facecolor="black", label="CDS")]
axg.legend(handles=h1, loc="upper center", bbox_to_anchor=(0.5, -0.04), fontsize=6.5, frameon=False, ncol=7)
h2 = [Patch(facecolor=motif_colors((no - 1) % 10), label="Motif %d (%d aa)" % (no, motif_w[no])) for no in sorted(motif_w)]
axp.legend(handles=h2, loc="upper center", bbox_to_anchor=(0.5, -0.04), fontsize=6, frameon=False, ncol=5)
fig.tight_layout()
fig.savefig(out_prefix + ".pdf"); fig.savefig(out_prefix + ".png", dpi=300)
print("clipped genes:", clipped)
print("genes drawn:", n, "| proteins with motif sites:", sum(1 for r in rows if sites.get(rep.get(r["loc_id"], ""))), "| meme:", meme)
print("wrote", out_prefix + ".pdf/.png")
