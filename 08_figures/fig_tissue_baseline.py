# -*- coding: utf-8 -*-
"""
Baseline (control-library) expression of the 87 sunflower aquaporin genes by
BioProject and tissue (Additional file figure).
  A  heatmap of log2(TPM + 1), mean of the control libraries of each BioProject
     and tissue (genes grouped by subfamily, ordered by maximum expression)
  B  root/leaf log2 ratio of the control means in the two BioProjects that
     sampled both tissues

Inputs: 04_expression/results/aquaporin_expression_v2/aquaporin_TPM_control_means.tsv
(09_aquaporin_tpm_v2.sh), Table S1 (subfamily, product) and Table S4 (number of
significant contrasts per gene).

usage: python fig_tissue_baseline.py [out_prefix]
"""
import os, sys, math, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import spearmanr

import figlib as L

out_prefix = sys.argv[1] if len(sys.argv) > 1 else L.fig_path("FigS7_tissue_baseline")
D = os.path.join(L.RESULTS, "04_expression", "results", "aquaporin_expression_v2")
s1 = L.load_table_s1()
SUB_ORDER = ["PIP", "TIP", "NIP", "SIP"]
GROUP_ORDER = ["PRJNA869183_leaf", "PRJNA797473_leaf", "PRJNA908908_leaf", "PRJNA1041959_leaf", "PRJNA492303_leaf",
               "PRJNA1041959_root", "PRJNA492303_root", "PRJNA850121_root"]

rows = L.read_tsv(os.path.join(D, "aquaporin_TPM_control_means.tsv"))
n_lib = rows[0]; rows = rows[1:]
tpm = {r["gene_id"]: {g: float(r[g]) for g in GROUP_ORDER} for r in rows}
genes = [g for g in tpm if g in s1]

# number of significant contrasts per gene (Table S4 layout)
s4 = L.read_tsv(os.path.join(L.RESULTS, "09_figures", "supplementary", "Table_S4_aqp_deg_all_contrasts.tsv"))
ncon = {}
for r in s4:
    n = 0
    for k, v in r.items():
        if k.endswith("_LFC") and v not in ("NA", ""):
            p = r.get(k[:-4] + "_padj", "NA")
            if p not in ("NA", "") and abs(float(v)) > 1 and float(p) < 0.05:
                n += 1
    ncon[r["gene_id"]] = n

leaf_groups = [g for g in GROUP_ORDER if g.endswith("leaf")]
root_groups = [g for g in GROUP_ORDER if g.endswith("root")]
summary = []
for g in genes:
    r = s1[g]
    d = dict(gene_id=g, product=L.short_name(r.get("product", "")), subfamily=r["subfamily"],
             mean_leaf=np.mean([tpm[g][x] for x in leaf_groups]), mean_root=np.mean([tpm[g][x] for x in root_groups]),
             r1=math.log2((tpm[g]["PRJNA1041959_root"] + 0.1) / (tpm[g]["PRJNA1041959_leaf"] + 0.1)),
             r2=math.log2((tpm[g]["PRJNA492303_root"] + 0.1) / (tpm[g]["PRJNA492303_leaf"] + 0.1)),
             mx=max(tpm[g].values()), n_sig=ncon.get(g, 0))
    summary.append(d)
summary.sort(key=lambda d: (SUB_ORDER.index(d["subfamily"]) if d["subfamily"] in SUB_ORDER else 9, -d["mx"]))

expressed = [d for d in summary if d["mx"] >= 1]
root_pref = [d for d in expressed if d["r1"] >= 2 and d["r2"] >= 2]
leaf_pref = [d for d in expressed if d["r1"] <= -2 and d["r2"] <= -2]
rho, p = spearmanr([d["r1"] for d in summary], [d["r2"] for d in summary])
print("genes:", len(summary), "| expressed (max TPM >= 1):", len(expressed), "| >= 10:", sum(1 for d in summary if d["mx"] >= 10),
      "| >= 100:", sum(1 for d in summary if d["mx"] >= 100))
print("silent in every control group:", [d["gene_id"] for d in summary if d["mx"] < 1])
for sub in SUB_ORDER:
    ds = [d for d in summary if d["subfamily"] == sub]
    print("%s n=%d median leaf %.1f root %.1f" % (sub, len(ds), np.median([d["mean_leaf"] for d in ds]), np.median([d["mean_root"] for d in ds])))
print("root/leaf ratio consistency: rho = %.2f, P = %.1e" % (rho, p))
print("root-preferential (both datasets):", len(root_pref), dict(collections.Counter(d["subfamily"] for d in root_pref)))
print("leaf-preferential (both datasets):", len(leaf_pref), [d["gene_id"] for d in leaf_pref])
rho3, p3 = spearmanr([math.log2(d["mx"] + 1) for d in summary], [d["n_sig"] for d in summary])
print("log2(max TPM + 1) vs significant contrasts: rho = %.2f, P = %.1e" % (rho3, p3))

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(10, 12.5))
gs = fig.add_gridspec(1, 2, width_ratios=[2.3, 2.0], wspace=0.55)
ax = fig.add_subplot(gs[0, 0])
mat = np.array([[math.log2(tpm[d["gene_id"]][x] + 1) for x in GROUP_ORDER] for d in summary])
im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0, vmax=max(8, mat.max()))
ax.set_xticks(range(len(GROUP_ORDER)))
ax.set_xticklabels(["%s\n%s (n = %s)" % (x.split("_")[0], x.split("_")[1], n_lib[x]) for x in GROUP_ORDER], fontsize=5.5, rotation=90)
ax.set_yticks(range(len(summary)))
ax.set_yticklabels(["%s (%s)" % (d["product"], d["gene_id"]) for d in summary], fontsize=4.6)
for lab, d in zip(ax.get_yticklabels(), summary):
    lab.set_color(L.SUBFAM_COLORS.get(d["subfamily"], "#7f7f7f"))
ax.axvline(4.5, color="white", lw=1.5)
ax.text(2, -1.0, "leaf", ha="center", va="bottom", fontsize=7); ax.text(6, -1.0, "root", ha="center", va="bottom", fontsize=7)
pos = 0
for sub in SUB_ORDER:
    n = sum(1 for d in summary if d["subfamily"] == sub)
    if pos:
        ax.axhline(pos - 0.5, color="white", lw=1.2)
    ax.text(len(GROUP_ORDER) - 0.4, pos + n / 2.0, sub, color=L.SUBFAM_COLORS[sub], fontsize=7, va="center", ha="left", fontweight="bold")
    pos += n
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.12, shrink=0.5)
cb.set_label("log2(TPM + 1), mean of control libraries", fontsize=6.5); cb.ax.tick_params(labelsize=6)
ax.set_title("A", loc="left", fontweight="bold")

ax2 = fig.add_subplot(gs[0, 1])
ax2.set_box_aspect(1)
for d in summary:
    ax2.scatter(d["r1"], d["r2"], s=16, color=L.SUBFAM_COLORS.get(d["subfamily"], "#7f7f7f"),
                alpha=0.9 if d["mx"] >= 1 else 0.25, edgecolor="none", zorder=3)
vals = [d["r1"] for d in summary] + [d["r2"] for d in summary]
lim = max(abs(v) for v in vals) + 1.0
ax2.plot([-lim, lim], [-lim, lim], color="gray", lw=0.6, ls="--", zorder=1)
ax2.axhline(0, color="gray", lw=0.4, zorder=1); ax2.axvline(0, color="gray", lw=0.4, zorder=1)
ax2.set_xlim(-lim, lim); ax2.set_ylim(-lim, lim)
ax2.set_xlabel("log2(root / leaf TPM), PRJNA1041959 controls", fontsize=7)
ax2.set_ylabel("log2(root / leaf TPM), PRJNA492303 controls", fontsize=7)
ax2.tick_params(labelsize=6)
ax2.text(0.03, 0.97, "Spearman rho = %.2f\nP = %.1e\nn = %d genes" % (rho, p, len(summary)), transform=ax2.transAxes, fontsize=6.5, va="top")
# labels: the two leaf-preferential genes and the two multi-stress candidates (the
# root-preferential genes are too many to label here; they are listed in panel A)
LABEL_OFFSET = {"LOC110898014": (10, -4), "LOC110908794": (-8, -12), "LOC110904457": (12, -10), "LOC110916097": (12, 6)}
for d in summary:
    if d["gene_id"] in LABEL_OFFSET:
        dx, dy = LABEL_OFFSET[d["gene_id"]]
        ax2.annotate("%s\n%s" % (d["product"], d["gene_id"]), (d["r1"], d["r2"]), fontsize=5.5, xytext=(dx, dy),
                     textcoords="offset points", ha="left" if dx > 0 else "right", va="bottom" if dy > 0 else "top",
                     arrowprops=dict(arrowstyle="-", lw=0.5, color="gray", shrinkA=0, shrinkB=2), zorder=4)
ax2.text(0.97, 0.03, "%d genes root-preferential in both datasets\n(log2 ratio >= 2); %d leaf-preferential" % (len(root_pref), len(leaf_pref)),
         transform=ax2.transAxes, fontsize=6, ha="right", va="bottom")
ax2.legend(handles=[Patch(color=L.SUBFAM_COLORS[s], label=s) for s in SUB_ORDER], fontsize=6, loc="upper left", bbox_to_anchor=(0.02, 0.80), frameon=False)
ax2.set_title("B", loc="left", fontweight="bold")
fig.savefig(out_prefix + ".pdf", bbox_inches="tight"); fig.savefig(out_prefix + ".png", dpi=300, bbox_inches="tight")
print("wrote", out_prefix + ".pdf/.png")
