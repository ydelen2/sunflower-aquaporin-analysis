# -*- coding: utf-8 -*-
"""
Aquaporin synteny figure.
  A  intra-genomic collinear aquaporin pairs in sunflower (circular layout,
     17 chromosomes, links coloured by subfamily)
  B  sunflower versus lettuce collinear aquaporin pairs (dual layout)
  C  sunflower versus Arabidopsis collinear aquaporin pairs (dual layout)

Inputs: 07_synteny/results_v2/*_aqp_synteny.tsv, 07_synteny/data/*.gff
(MCScanX gene positions), chromosome lengths, Table S1 (inventory and
subfamily). Intra-genomic pairs are drawn only when both genes are in the
inventory.

usage: python fig_synteny.py [out_prefix]
"""
import math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path
from matplotlib.lines import Line2D

import figlib as L

out_prefix = sys.argv[1] if len(sys.argv) > 1 else L.fig_path("Fig2_synteny")
SYN = os.path.join(L.RESULTS, "07_synteny")

s1 = L.load_table_s1()
p2g = L.load_protein_gene()
inv = set(s1)
sub_of_gene = {g: r["subfamily"] for g, r in s1.items()}
chrlen = L.load_chrom_lengths()

# gene positions (protein accession -> (chr label, midpoint))
pos = {}
for sp in ("Ha", "At", "Ls"):
    with open(os.path.join(SYN, "data", sp + ".gff")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            lab = L.chr_label(f[0])
            if lab:
                pos[f[1]] = (lab, (int(f[2]) + int(f[3])) / 2.0)

def load_pairs(name):
    rows = L.read_tsv(os.path.join(SYN, "results_v2", name))
    out = []
    for r in rows:
        out.append((r["prot1"], r["gene1"], r["species1"], r["prot2"], r["gene2"], r["species2"]))
    return out

intra = [p for p in load_pairs("Ha_self_intra_aqp_synteny.tsv") if p[1] in inv and p[4] in inv]
inter = {}
for sp, fn in (("Ls", "Ha_Ls_inter_aqp_synteny.tsv"), ("At", "Ha_At_inter_aqp_synteny.tsv")):
    rows = []
    for p1, g1, s1_, p2, g2, s2_ in load_pairs(fn):
        if s1_ == s2_:          # intra-genomic block inside a two-species run
            continue
        ha, other = ((p1, g1), (p2, g2)) if s1_ == "Ha" else ((p2, g2), (p1, g1))
        if ha[1] in inv and ha[0] in pos and other[0] in pos:
            rows.append((ha, other))
    inter[sp] = rows

print("intra AQP-AQP pairs:", len(intra), "| Ha-Ls:", len(inter["Ls"]), "| Ha-At:", len(inter["At"]))

fig = plt.figure(figsize=(12, 13))
gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1.0], hspace=0.12)
axA = fig.add_subplot(gs[0])
gsB = gs[1].subgridspec(2, 1, hspace=0.55)
axB = fig.add_subplot(gsB[0]); axC = fig.add_subplot(gsB[1])

# ------------------------------------------------------------------ panel A
ha_chr = ["Ha%d" % i for i in range(1, 18)]
total = sum(chrlen[c] for c in ha_chr)
gap = 0.012 * 2 * math.pi
avail = 2 * math.pi - gap * len(ha_chr)
start = {}
a = math.pi / 2
for c in ha_chr:
    span = avail * chrlen[c] / total
    start[c] = (a, a - span)   # clockwise
    a -= span + gap

def ang(c, bp):
    a0, a1 = start[c]
    return a0 + (a1 - a0) * bp / chrlen[c]

R = 1.0
axA.set_aspect("equal"); axA.axis("off"); axA.set_xlim(-1.32, 1.32); axA.set_ylim(-1.28, 1.28)
for c in ha_chr:
    a0, a1 = start[c]
    axA.add_patch(Wedge((0, 0), R + 0.06, math.degrees(a1), math.degrees(a0), width=0.06, facecolor="#d9d9d9", edgecolor="#555555", lw=0.6))
    am = (a0 + a1) / 2
    axA.text(1.15 * math.cos(am), 1.15 * math.sin(am), c, ha="center", va="center", fontsize=8, fontweight="bold")

# gene ticks
for xp, g in p2g.items():
    if g["gene_id"] in inv and xp in pos:
        c, mid = pos[xp]
        t = ang(c, mid)
        col = L.SUBFAM_COLORS.get(sub_of_gene.get(g["gene_id"], "other"), "#7f7f7f")
        axA.plot([(R + 0.06) * math.cos(t), (R + 0.09) * math.sin(t) * 0 + (R + 0.09) * math.cos(t)],
                 [(R + 0.06) * math.sin(t), (R + 0.09) * math.sin(t)], color=col, lw=1.2)

def chord(ax, t1, t2, col, r=R):
    x1, y1 = r * math.cos(t1), r * math.sin(t1)
    x2, y2 = r * math.cos(t2), r * math.sin(t2)
    path = Path([(x1, y1), (0, 0), (x2, y2)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
    ax.add_patch(PathPatch(path, facecolor="none", edgecolor=col, lw=1.0, alpha=0.75))

for p1, g1, _, p2, g2, _ in intra:
    if p1 in pos and p2 in pos:
        c1, m1 = pos[p1]; c2, m2 = pos[p2]
        chord(axA, ang(c1, m1), ang(c2, m2), L.SUBFAM_COLORS.get(sub_of_gene.get(g1, "other"), "#7f7f7f"))

axA.text(-1.30, 1.22, "A", fontsize=14, fontweight="bold")
h = [Line2D([0], [0], color=L.SUBFAM_COLORS[s], lw=3, label=s) for s in ["PIP", "TIP", "NIP", "SIP", "XIP"] if any(sub_of_gene.get(g) == s for g in inv)]
axA.legend(handles=h, loc="lower left", fontsize=8, frameon=False, title="Subfamily", title_fontsize=8)

# ------------------------------------------------------------ panels B, C
def dual(ax, other, letter, title):
    ax.axis("off")
    top = ha_chr
    bot = [c for c in sorted(chrlen) if c.startswith(other)]
    def layout(chrs, y, width=100.0, gapw=0.8):
        tot = sum(chrlen[c] for c in chrs); scale_ = (width - gapw * (len(chrs) - 1)) / tot
        x = 0.0; out = {}
        for c in chrs:
            w = chrlen[c] * scale_
            out[c] = (x, w)
            ax.add_patch(plt.Rectangle((x, y - 0.6), w, 1.2, facecolor="#d9d9d9", edgecolor="#555555", lw=0.6))
            ax.text(x + w / 2, y + (1.6 if y > 0 else -1.6), c, ha="center", va="center" if y > 0 else "top", fontsize=6.5)
            x += w + gapw
        return out
    ytop, ybot = 8.0, 0.0
    lt = layout(top, ytop); lb = layout(bot, ybot)
    def xof(lay, c, bp): x, w = lay[c]; return x + w * bp / chrlen[c]
    n = 0
    for (hp, hg), (op, og) in inter[other]:
        c1, m1 = pos[hp]; c2, m2 = pos[op]
        if c1 not in lt or c2 not in lb:
            continue
        x1 = xof(lt, c1, m1); x2 = xof(lb, c2, m2)
        col = L.SUBFAM_COLORS.get(sub_of_gene.get(hg, "other"), "#7f7f7f")
        path = Path([(x1, ytop - 0.6), (x1, (ytop + ybot) / 2), (x2, (ytop + ybot) / 2), (x2, ybot + 0.6)],
                    [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=col, lw=0.9, alpha=0.8))
        ax.plot([x1], [ytop - 0.6], "|", color=col, ms=6, mew=1.2)
        ax.plot([x2], [ybot + 0.6], "|", color=col, ms=6, mew=1.2)
        n += 1
    ax.set_xlim(-2, 102); ax.set_ylim(-3.5, 11)
    ax.text(-2, 10.6, letter, fontsize=14, fontweight="bold")
    ax.text(-1.5, ytop, "H. annuus", ha="right", va="center", fontsize=7, style="italic")
    ax.text(-1.5, ybot, L.SPECIES[other], ha="right", va="center", fontsize=7, style="italic")

dual(axB, "Ls", "B", "Sunflower versus lettuce collinear aquaporin pairs")
dual(axC, "At", "C", "Sunflower versus Arabidopsis collinear aquaporin pairs")

fig.savefig(out_prefix + ".pdf", bbox_inches="tight")
fig.savefig(out_prefix + ".png", dpi=300, bbox_inches="tight")
print("wrote", out_prefix + ".pdf/.png")
