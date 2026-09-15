# -*- coding: utf-8 -*-
"""
Maximum-likelihood tree of aquaporins from five species, fan layout.
Tip marker = species, tip label color = subfamily. Black dots mark nodes
with ultrafast bootstrap support >= 95.

usage: python fig_tree.py [treefile] [out_prefix]
"""
import math, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

import figlib as L

tree_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(L.RESULTS, "05_phylogenetics", "results_v2", "aquaporin_tree.treefile")
out_prefix = sys.argv[2] if len(sys.argv) > 2 else L.fig_path("Fig1_phylogenetic_tree")

tree = Phylo.read(tree_path, "newick")
tree.root_at_midpoint()
tree.ladderize()

products = L.load_products()
sub_ha = L.load_subfamily_by_protein()
ref_names = L.load_ref_names()

tips = tree.get_terminals()
n = len(tips)

# ---- layout: angle per tip, radius = cumulative branch length ---------------
angle = {}
for i, t in enumerate(tips):
    angle[t] = 2 * math.pi * i / n
for cl in tree.get_nonterminals(order="postorder"):
    angle[cl] = sum(angle[c] for c in cl.clades) / len(cl.clades)

depth = tree.depths()
rmax = max(depth.values())
scale = 1.0 / rmax

def polar(cl, r=None):
    rr = (depth[cl] if r is None else r) * scale
    return rr * math.cos(angle[cl]), rr * math.sin(angle[cl])

def ufboot(cl):
    lab = cl.name if cl.name else (str(cl.confidence) if cl.confidence is not None else "")
    m = re.match(r"^\s*([\d.]+)/([\d.]+)\s*$", lab)
    if m:
        return float(m.group(2))
    try:
        return float(lab)
    except ValueError:
        return None

fig = plt.figure(figsize=(13, 13))
ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
ax.set_aspect("equal"); ax.axis("off")

# branches
for cl in tree.get_nonterminals():
    r_par = depth[cl] * scale
    angs = [angle[c] for c in cl.clades]
    a0, a1 = min(angs), max(angs)
    arc = [(r_par * math.cos(a), r_par * math.sin(a)) for a in [a0 + (a1 - a0) * k / 40.0 for k in range(41)]]
    ax.plot([p[0] for p in arc], [p[1] for p in arc], color="#444444", lw=0.6, solid_capstyle="round")
    for c in cl.clades:
        x0, y0 = polar(c, depth[cl]); x1, y1 = polar(c)
        ax.plot([x0, x1], [y0, y1], color="#444444", lw=0.6, solid_capstyle="round")
    b = ufboot(cl)
    if b is not None and b >= 95 and cl is not tree.root:
        x, y = polar(cl)
        ax.plot(x, y, "o", ms=2.2, color="black", zorder=5)

# tips
label_r = 1.02
for t in tips:
    sp, acc = t.name.split("_", 1)
    x, y = polar(t)
    if sp == "Ha":
        desc = products.get(acc, "")
        sub = sub_ha.get(acc, L.subfamily_from_name(desc))
        lab = "%s %s" % (L.short_name(desc), acc)
        weight = "bold"
    else:
        desc = ref_names.get(acc, ("", ""))[0]
        sub = L.subfamily_from_name(desc)
        lab = "%s %s" % (L.short_name(desc), acc) if desc else acc
        weight = "normal"
    col = L.SUBFAM_COLORS.get(sub, "#7f7f7f")
    ax.plot(x, y, L.SPECIES_MARKERS[sp], ms=(4.6 if sp == "Ha" else 3.2), mfc=L.SPECIES_COLORS[sp], mec="white", mew=0.3, zorder=6)
    # connector from tip to label ring
    xr, yr = polar(t, rmax * label_r)
    ax.plot([x, xr], [y, yr], color="#cccccc", lw=0.3, zorder=1)
    a = angle[t]; deg = math.degrees(a)
    ha = "left"
    if 90 < deg % 360 < 270:
        deg += 180; ha = "right"
    ax.text(xr * 1.01, yr * 1.01, lab, rotation=deg, rotation_mode="anchor", ha=ha, va="center",
            fontsize=3.4, color=col, fontweight=weight)

# scale bar
ax.plot([0.95, 0.95 + 0.2 * scale], [-1.30, -1.30], color="black", lw=1)
ax.text(0.95 + 0.1 * scale, -1.32, "0.2 substitutions/site", fontsize=6, va="top", ha="center")
ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.35)

# legends
h1 = [Line2D([0], [0], marker=L.SPECIES_MARKERS[s], color="none", mfc=L.SPECIES_COLORS[s], mec="none", ms=6,
             label=L.SPECIES[s]) for s in ["Ha", "At", "Os", "Sl", "Ls"]]
h2 = [Line2D([0], [0], color=L.SUBFAM_COLORS[s], lw=4, label=s) for s in ["PIP", "TIP", "NIP", "SIP", "XIP"]]
h2.append(Line2D([0], [0], color="#7f7f7f", lw=4, label="no subfamily in name"))
h3 = [Line2D([0], [0], marker="o", color="none", mfc="black", mec="none", ms=4, label="UFBoot >= 95")]
leg1 = ax.legend(handles=h1, loc="upper left", fontsize=7, frameon=False, title="Species", title_fontsize=7)
ax.add_artist(leg1)
ax.legend(handles=h2 + h3, loc="lower left", fontsize=7, frameon=False, title="Subfamily (label color)", title_fontsize=7)

fig.savefig(out_prefix + ".pdf")
fig.savefig(out_prefix + ".png", dpi=300)
from collections import Counter
print("tips:", n, dict(Counter(t.name.split("_")[0] for t in tips)))
print("wrote", out_prefix + ".pdf/.png")
