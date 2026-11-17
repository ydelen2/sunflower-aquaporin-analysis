# -*- coding: utf-8 -*-
"""
Maximum-likelihood tree of aquaporins from five species.

Two layouts:
  fan     midpoint-rooted circular tree with one label per tip
  radial  unrooted equal-angle tree without tip labels; subfamily names are
          written next to the clades

In both layouts the tip marker shows the species, branches are colored by
subfamily, and black dots mark nodes with ultrafast bootstrap support >= 95.

The fan tree is Fig. 1 and the radial tree is Fig. S1.

usage: python fig_tree.py [treefile] [out_prefix] [fan|radial]
(an empty string keeps the default of an argument, for example: python fig_tree.py "" "" radial)
"""
import math, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
from Bio import Phylo

import figlib as L

layout = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] else "fan"
tree_path = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] else os.path.join(L.RESULTS, "05_phylogenetics", "results_v2", "aquaporin_tree.treefile")
out_prefix = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] else L.fig_path("FigS1_radial_tree" if layout == "radial" else "Fig1_phylogenetic_tree")

products = L.load_products()
sub_ha = L.load_subfamily_by_protein()
ref_names = L.load_ref_names()
GRAY = "#7f7f7f"
BRANCH = "#444444"
FIVE = ("PIP", "TIP", "NIP", "SIP", "XIP")


def tip_info(t):
    """species tag, subfamily, label and font weight of a tip"""
    sp, acc = t.name.split("_", 1)
    if sp == "Ha":
        desc = products.get(acc, "")
        sub = sub_ha.get(acc, L.subfamily_from_name(desc))
        return sp, sub, "%s %s" % (L.short_name(desc), acc), "bold"
    desc = ref_names.get(acc, ("", ""))[0]
    return sp, L.subfamily_from_name(desc), ("%s %s" % (L.short_name(desc), acc) if desc else acc), "normal"


def ufboot(cl):
    lab = cl.name if cl.name else (str(cl.confidence) if cl.confidence is not None else "")
    m = re.match(r"^\s*([\d.]+)/([\d.]+)\s*$", lab)
    if m:
        return float(m.group(2))
    try:
        return float(lab)
    except ValueError:
        return None


def branch_supports():
    """Ultrafast bootstrap of every internal branch of the unrooted tree, keyed by
    the set of tips on the side of the branch that does not hold a fixed tip."""
    t = Phylo.read(tree_path, "newick")
    names = frozenset(x.name for x in t.get_terminals())
    ref = min(names)
    out = {}
    for cl in t.get_nonterminals():
        s = frozenset(x.name for x in cl.get_terminals())
        if 1 < len(s) < len(names) - 1:
            out[s if ref not in s else names - s] = ufboot(cl)
    return out, names, ref


def node_supports(tree):
    """Support of the branch above each internal node of a rerooted tree. Rerooting
    moves the node labels of the Newick file away from their branches, so the values
    are looked up by the split that the branch defines."""
    by_split, names, ref = branch_supports()
    out = {}
    seen = set()
    for cl in sorted(tree.get_nonterminals(), key=lambda c: len(c.get_terminals())):
        if cl is tree.root:
            continue
        s = frozenset(x.name for x in cl.get_terminals())
        key = s if ref not in s else names - s
        if key in by_split and key not in seen:      # the two children of the root share one branch
            out[cl] = by_split[key]
            seen.add(key)
    return out


def clade_colors(tree):
    """Color of the branch leading to each clade: the subfamily color when all
    tips of the clade that carry a subfamily agree, dark gray otherwise."""
    subs = {}
    for cl in tree.find_clades(order="postorder"):
        if cl.is_terminal():
            s = tip_info(cl)[1]
            subs[cl] = {s} if s in FIVE else set()
        else:
            subs[cl] = set().union(*[subs[c] for c in cl.clades])
    return {cl: (L.SUBFAM_COLORS[next(iter(s))] if len(s) == 1 else BRANCH) for cl, s in subs.items()}


def legends(ax, fs, sub_title, with_gray):
    italic = FontProperties(style="italic", size=fs)
    h1 = [Line2D([0], [0], marker=L.SPECIES_MARKERS[s], color="none", mfc=L.SPECIES_COLORS[s], mec="none", ms=fs * 0.9,
                 label=L.SPECIES[s]) for s in ["Ha", "At", "Os", "Sl", "Ls"]]
    h2 = [Line2D([0], [0], color=L.SUBFAM_COLORS[s], lw=4, label=s) for s in FIVE]
    if with_gray:
        h2.append(Line2D([0], [0], color=GRAY, lw=4, label="no subfamily in name"))
    h3 = [Line2D([0], [0], marker="o", color="none", mfc="black", mec="none", ms=fs * 0.55, label="UFBoot ≥ 95")]
    leg1 = ax.legend(handles=h1, loc="upper left", prop=italic, frameon=False, title="Species", title_fontsize=fs + 1)
    leg1._legend_box.align = "left"
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=h2 + h3, loc="lower left", fontsize=fs, frameon=False, title=sub_title, title_fontsize=fs + 1)
    leg2._legend_box.align = "left"


# ============================================================== fan layout
def draw_fan():
    tree = Phylo.read(tree_path, "newick")
    tree.root_at_midpoint()
    tree.ladderize()
    tips = tree.get_terminals()
    n = len(tips)
    angle = {}
    for i, t in enumerate(tips):
        angle[t] = 2 * math.pi * i / n
    for cl in tree.get_nonterminals(order="postorder"):
        angle[cl] = sum(angle[c] for c in cl.clades) / len(cl.clades)
    depth = tree.depths()
    rmax = max(depth.values())
    scale = 1.0 / rmax
    color = clade_colors(tree)
    support = node_supports(tree)

    def polar(cl, r=None):
        rr = (depth[cl] if r is None else r) * scale
        return rr * math.cos(angle[cl]), rr * math.sin(angle[cl])

    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
    ax.set_aspect("equal"); ax.axis("off")
    for cl in tree.get_nonterminals():
        r_par = depth[cl] * scale
        angs = [angle[c] for c in cl.clades]
        a0, a1 = min(angs), max(angs)
        arc = [(r_par * math.cos(a), r_par * math.sin(a)) for a in [a0 + (a1 - a0) * k / 40.0 for k in range(41)]]
        ax.plot([p[0] for p in arc], [p[1] for p in arc], color=color[cl], lw=0.8, solid_capstyle="round")
        for c in cl.clades:
            x0, y0 = polar(c, depth[cl]); x1, y1 = polar(c)
            ax.plot([x0, x1], [y0, y1], color=color[c], lw=0.8, solid_capstyle="round")
        b = support.get(cl)
        if b is not None and b >= 95:
            x, y = polar(cl)
            ax.plot(x, y, "o", ms=2.4, color="black", zorder=5)
    label_r = 1.02
    for t in tips:
        sp, sub, lab, weight = tip_info(t)
        x, y = polar(t)
        col = L.SUBFAM_COLORS.get(sub, GRAY)
        ax.plot(x, y, L.SPECIES_MARKERS[sp], ms=(5.2 if sp == "Ha" else 3.6), mfc=L.SPECIES_COLORS[sp], mec="white", mew=0.3, zorder=6)
        xr, yr = polar(t, rmax * label_r)
        ax.plot([x, xr], [y, yr], color="#cccccc", lw=0.3, zorder=1)
        deg = math.degrees(angle[t]); ha = "left"
        if 90 < deg % 360 < 270:
            deg += 180; ha = "right"
        ax.text(xr * 1.01, yr * 1.01, lab, rotation=deg, rotation_mode="anchor", ha=ha, va="center",
                fontsize=4.6, color=col, fontweight=weight)
    lim = 1.47
    ax.plot([lim - 0.42, lim - 0.42 + 0.2 * scale], [-lim + 0.07, -lim + 0.07], color="black", lw=1.2)
    ax.text(lim - 0.42 + 0.1 * scale, -lim + 0.05, "0.2 substitutions/site", fontsize=11, va="top", ha="center")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    legends(ax, 12, "Subfamily (branch and label color)", True)
    return fig, tips


# =========================================================== radial layout
def draw_radial():
    tree = Phylo.read(tree_path, "newick")
    tree.ladderize()
    ntip = {}
    for cl in tree.find_clades(order="postorder"):
        ntip[cl] = 1 if cl.is_terminal() else sum(ntip[c] for c in cl.clades)
    pos = {tree.root: (0.0, 0.0)}
    wedge = {tree.root: (0.0, 2 * math.pi)}
    for cl in tree.find_clades(order="preorder"):
        a0, a1 = wedge[cl]
        a = a0
        for c in cl.clades:
            w = (a1 - a0) * ntip[c] / ntip[cl]
            wedge[c] = (a, a + w)
            mid = a + w / 2.0
            bl = c.branch_length or 0.0
            pos[c] = (pos[cl][0] + bl * math.cos(mid), pos[cl][1] + bl * math.sin(mid))
            a += w
    color = clade_colors(tree)
    tips = tree.get_terminals()

    fig = plt.figure(figsize=(11, 11))
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
    ax.set_aspect("equal", adjustable="datalim"); ax.axis("off")
    for cl in tree.get_nonterminals():
        for c in cl.clades:
            ax.plot([pos[cl][0], pos[c][0]], [pos[cl][1], pos[c][1]], color=color[c], lw=0.9, solid_capstyle="round", zorder=2)
        b = ufboot(cl)
        if b is not None and b >= 95 and cl is not tree.root:
            ax.plot(pos[cl][0], pos[cl][1], "o", ms=2.6, color="black", zorder=5)
    cent = {}
    for t in tips:
        sp, sub, lab, weight = tip_info(t)
        x, y = pos[t]
        ax.plot(x, y, L.SPECIES_MARKERS[sp], ms=(7.5 if sp == "Ha" else 4.6), mfc=L.SPECIES_COLORS[sp], mec="white", mew=0.4, zorder=6)
        if sub in FIVE:
            cent.setdefault(sub, []).append((x, y))
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    span = max(max(xs) - min(xs), max(ys) - min(ys))
    for sub, pts in cent.items():
        cx = sum(p[0] for p in pts) / len(pts); cy = sum(p[1] for p in pts) / len(pts)
        d = math.hypot(cx, cy)
        far = max(p[0] * cx / d + p[1] * cy / d for p in pts)
        ax.text(cx / d * (far + 0.05 * span), cy / d * (far + 0.05 * span), sub, color=L.SUBFAM_COLORS[sub],
                fontsize=17, fontweight="bold", ha="center", va="center")
    pad = 0.10 * span
    ax.set_xlim(min(xs) - pad, max(xs) + pad); ax.set_ylim(min(ys) - pad, max(ys) + pad)
    fig.canvas.draw()                      # the limits are final after the equal-aspect adjustment
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    xb, yb = x1 - 0.07 * (x1 - x0), y0 + 0.05 * (y1 - y0)
    ax.plot([xb - 0.2, xb], [yb, yb], color="black", lw=1.2)
    ax.text(xb - 0.1, yb - 0.008 * (y1 - y0), "0.2 substitutions/site", fontsize=11, va="top", ha="center")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    legends(ax, 12, "Subfamily (branch color)", False)
    return fig, tips


fig, tips = draw_radial() if layout == "radial" else draw_fan()
fig.savefig(out_prefix + ".pdf")
fig.savefig(out_prefix + ".png", dpi=300)
from collections import Counter
print("layout:", layout, "| tips:", len(tips), dict(Counter(t.name.split("_")[0] for t in tips)))
print("wrote", out_prefix + ".pdf/.png")
