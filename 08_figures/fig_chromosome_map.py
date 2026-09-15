# -*- coding: utf-8 -*-
"""
Chromosomal distribution of the sunflower aquaporin genes (Additional file
figure). Chromosomes are drawn to scale; each gene is a tick coloured by
subfamily and labelled with its short product name.

usage: python fig_chromosome_map.py [out_prefix]
"""
import os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch

import figlib as L

out_prefix = sys.argv[1] if len(sys.argv) > 1 else L.fig_path("FigS2_chromosome_map")
s1 = L.load_table_s1()
chrlen = L.load_chrom_lengths()
chrs = ["Ha%d" % i for i in range(1, 18)]

genes = []
for g, r in s1.items():
    chrom = L.HA_CHR.get(r.get("chromosome", ""), r.get("chr_name", ""))
    if chrom not in chrs:
        continue
    start = int(float(r["start"]))
    name = L.short_name(r.get("product", ""))
    genes.append((chrom, start, g, r["subfamily"], name))
print("genes placed:", len(genes))

fig, ax = plt.subplots(figsize=(19, 8.5))
maxlen = max(chrlen[c] for c in chrs) / 1e6
xgap = 1.3
for i, c in enumerate(chrs):
    x = i * xgap
    h = chrlen[c] / 1e6
    ax.add_patch(FancyBboxPatch((x - 0.12, 0), 0.24, h, boxstyle="round,pad=0,rounding_size=0.1",
                                facecolor="#e6e6e6", edgecolor="#555555", lw=0.8))
    ax.text(x, -4, c, ha="center", va="top", fontsize=9, fontweight="bold")
    gs = sorted([g for g in genes if g[0] == c], key=lambda g: g[1])
    n = len(gs)
    ax.text(x, h + 7, "n = %d" % n, ha="center", va="top", fontsize=6.5, color="#555555")
    # labels: spread to avoid overlap (minimum spacing in Mb)
    ys = [g[1] / 1e6 for g in gs]
    lab_y = []
    min_gap = maxlen * 0.022
    for y in ys:
        if lab_y and y - lab_y[-1] < min_gap:
            y = lab_y[-1] + min_gap
        lab_y.append(y)
    # shift down if the last label overshoots the chromosome
    if lab_y and lab_y[-1] > h + min_gap:
        shift = lab_y[-1] - (h + min_gap)
        lab_y = [y - shift for y in lab_y]
    for (chrom, start, g, sub, name), y, ly in zip(gs, ys, lab_y):
        col = L.SUBFAM_COLORS.get(sub, "#7f7f7f")
        ax.plot([x - 0.12, x + 0.12], [y, y], color=col, lw=1.6)
        ax.plot([x + 0.12, x + 0.30, x + 0.38], [y, ly, ly], color=col, lw=0.5)
        ax.text(x + 0.40, ly, "%s %s" % (name, g.replace("LOC", "")), va="center", fontsize=4.3, color=col)

ax.set_xlim(-0.6, len(chrs) * xgap)
ax.set_ylim(maxlen * 1.08, -12)
ax.set_ylabel("Position (Mb)", fontsize=9)
ax.set_xticks([])
for s in ("top", "right", "bottom"):
    ax.spines[s].set_visible(False)
h = [Patch(facecolor=L.SUBFAM_COLORS[s], label=s) for s in ["PIP", "TIP", "NIP", "SIP", "XIP"] if any(g[3] == s for g in genes)]
ax.legend(handles=h, loc="lower right", fontsize=8, frameon=False, title="Subfamily", title_fontsize=8)
fig.savefig(out_prefix + ".pdf", bbox_inches="tight")
fig.savefig(out_prefix + ".png", dpi=300, bbox_inches="tight")
print("wrote", out_prefix + ".pdf/.png")
