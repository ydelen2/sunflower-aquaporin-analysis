# -*- coding: utf-8 -*-
"""
Fig. 4: volcano plots of four stress contrasts with the aquaporin genes
highlighted. Needs the genome-wide DESeq2 result files
(04_expression/results/deseq2/<project>/DEG_*.tsv under results/).

usage: python fig_volcano.py
"""
import os, glob, math, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import figlib as L

s1 = L.load_table_s1()
D = os.path.join(L.RESULTS, "04_expression", "results", "deseq2")
PANELS = [("Cold", "PRJNA869183_fixed", "DEG_cold_stress_vs_cold_0h.tsv"),
          ("Drought 14 d", "PRJNA797473_fixed", "DEG_drought_14d_vs_control.tsv"),
          ("Salt (NaCl)", "PRJNA869183_fixed", "DEG_salt_NaCl_vs_control.tsv"),
          ("Sclerotinia", "PRJNA908908", None)]
fig, axes = plt.subplots(2, 2, figsize=(10, 8.5))
for ax, (title, proj, fn) in zip(axes.flat, PANELS):
    path = os.path.join(D, proj, fn) if fn else (glob.glob(os.path.join(D, proj, "DEG_*_vs_*.tsv")) or [None])[0]
    if not path or not os.path.exists(path):
        ax.text(0.5, 0.5, "file missing:\n%s/%s" % (proj, fn), ha="center", va="center", transform=ax.transAxes, fontsize=8); ax.set_title(title); continue
    rows = L.read_tsv(path)
    x, y, ax_x, ax_y, labels = [], [], [], [], []
    for r in rows:
        try:
            l, p = float(r["log2FoldChange"]), float(r["padj"])
        except (ValueError, KeyError):
            continue
        yy = -math.log10(max(p, 1e-300))
        if r["gene_id"] in s1:
            ax_x.append(l); ax_y.append(yy); labels.append((l, yy, r["gene_id"], abs(l) > 1 and p < 0.05))
        else:
            x.append(l); y.append(yy)
    x, y = np.array(x), np.array(y)
    sig_up = (x > 1) & (y > -math.log10(0.05)); sig_dn = (x < -1) & (y > -math.log10(0.05))
    ax.scatter(x[~(sig_up | sig_dn)], y[~(sig_up | sig_dn)], s=3, c="#cccccc", alpha=0.5, linewidths=0, rasterized=True)
    ax.scatter(x[sig_up], y[sig_up], s=3, c="#f5b7b1", alpha=0.7, linewidths=0, rasterized=True)
    ax.scatter(x[sig_dn], y[sig_dn], s=3, c="#aed6f1", alpha=0.7, linewidths=0, rasterized=True)
    ax.scatter(ax_x, ax_y, s=22, c="#d4ac0d", edgecolors="black", linewidths=0.5, zorder=5)
    # the three most extreme significant aquaporins are labelled; labels sit in a fixed column on
    # their own side of the panel (axes fraction) so they never overlap each other or leave the axes
    k = 0; side = {"left": 0, "right": 0}
    for l, yy, g, s in sorted(labels, key=lambda t: -abs(t[0])):
        if s and k < 3:
            sd = "left" if l < 0 else "right"
            ax.annotate("%s %s" % (L.short_name(s1[g]["product"]), g.replace("LOC", "")), (l, yy), fontsize=6,
                        xytext=(0.03 if sd == "left" else 0.97, 0.66 - 0.09 * side[sd]), textcoords="axes fraction", ha=sd, va="center",
                        arrowprops=dict(arrowstyle="-", lw=0.4, color="#555555", shrinkB=3))
            side[sd] += 1; k += 1
    ax.axvline(-1, ls="--", lw=0.6, color="black"); ax.axvline(1, ls="--", lw=0.6, color="black"); ax.axhline(-math.log10(0.05), ls="--", lw=0.6, color="black")
    n_up = sum(1 for l, yy, g, s in labels if s and l > 0); n_dn = sum(1 for l, yy, g, s in labels if s and l < 0)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("log2 fold change", fontsize=8); ax.set_ylabel("-log10 adjusted P", fontsize=8); ax.tick_params(labelsize=7)
    # x limits: the bulk of the genome-wide points (99.9th percentile) but never clipping an aquaporin
    lim = max(6, np.percentile(np.abs(x), 99.9), (max(abs(v) for v in ax_x) + 0.5) if ax_x else 0); ax.set_xlim(-lim, lim)
fig.tight_layout()
fig.savefig(L.fig_path("Fig4_volcano_plots.pdf")); fig.savefig(L.fig_path("Fig4_volcano_plots.png"), dpi=300)
print("wrote Fig4_volcano_plots.pdf/.png")
