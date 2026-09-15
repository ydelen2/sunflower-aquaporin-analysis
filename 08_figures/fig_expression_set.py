# -*- coding: utf-8 -*-
"""
Expression figures from the aquaporin DEG table (Table S8) and the WGCNA tables:
  Fig3   up/down aquaporin DEG counts per stress contrast (main text)
  FigS3  log2 fold-change heatmap, genes x contrasts (significant cells marked)
  FigS4  overlap of aquaporin DEG sets across the seven stress types (UpSet-style)
  FigS5  broadly responsive genes: dot plot of log2FC across the 16 contrasts
  FigS6  WGCNA module-trait correlation heatmap for the aquaporin-containing modules

usage: python fig_expression_set.py
"""
import os, re, sys, itertools, collections, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import figlib as L

s1 = L.load_table_s1()
s4 = {r["gene_id"]: r for r in L.read_tsv(os.path.join(L.RESULTS, "09_figures", "supplementary", "Table_S4_aqp_deg_all_contrasts.tsv"))}
ORDER = ["Cold", "Heat", "Drought_869", "Salt_NaCl", "Rehydration", "Drought_7d", "Drought_14d", "Drought_21d",
         "PEG_72h", "Flooding", "Sclerotinia", "Orobanche_A", "Orobanche_B", "Orobanche_C", "Orobanche_D", "Orobanche_E"]
LABEL = {"Drought_869": "Drought (PRJNA869183)", "Salt_NaCl": "Salt (NaCl)", "PEG_72h": "PEG 72 h", "Drought_7d": "Drought 7 d",
         "Drought_14d": "Drought 14 d", "Drought_21d": "Drought 21 d"}
lab = lambda c: LABEL.get(c, c.replace("_", " "))
STRESS_TYPE = {"Cold": "Cold", "Heat": "Heat", "Drought_869": "Drought", "Drought_7d": "Drought", "Drought_14d": "Drought",
               "Drought_21d": "Drought", "PEG_72h": "Drought", "Salt_NaCl": "Salt", "Rehydration": "Rehydration",
               "Flooding": "Flooding", "Sclerotinia": "Sclerotinia", "Orobanche_A": "Orobanche", "Orobanche_B": "Orobanche",
               "Orobanche_C": "Orobanche", "Orobanche_D": "Orobanche", "Orobanche_E": "Orobanche"}
LFC_CUT, PADJ_CUT = 1.0, 0.05

def val(g, c):
    r = s4[g]; l, p = r.get(c + "_LFC", "NA"), r.get(c + "_padj", "NA")
    if l in ("NA", "") or p in ("NA", ""):
        return None, None
    return float(l), float(p)

def sig(g, c):
    l, p = val(g, c)
    return l is not None and abs(l) > LFC_CUT and p < PADJ_CUT

genes = [g for g in s4]
order = {"PIP": 0, "TIP": 1, "NIP": 2, "SIP": 3}
genes.sort(key=lambda g: (order.get(s1[g]["subfamily"], 9), L.short_name(s1[g]["product"]), g))
name = lambda g: "%s %s" % (L.short_name(s1[g]["product"]), g)

# ------------------------------------------------------------------ Fig 3
up = [sum(1 for g in genes if sig(g, c) and val(g, c)[0] > 0) for c in ORDER]
dn = [sum(1 for g in genes if sig(g, c) and val(g, c)[0] < 0) for c in ORDER]
fig, ax = plt.subplots(figsize=(9, 4.2))
x = np.arange(len(ORDER))
ax.bar(x, up, color="#c0392b", label="Upregulated")
ax.bar(x, [-v for v in dn], color="#2e86c1", label="Downregulated")
for i, (u, d_) in enumerate(zip(up, dn)):
    if u: ax.text(i, u + 0.6, str(u), ha="center", va="bottom", fontsize=7)
    if d_: ax.text(i, -d_ - 0.6, str(d_), ha="center", va="top", fontsize=7)
ax.axhline(0, color="black", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels([lab(c) for c in ORDER], rotation=55, ha="right", fontsize=8)
ax.set_ylabel("Aquaporin DEGs (n)", fontsize=9)
ax.set_ylim(-max(dn) - 6, max(up) + 6)
ax.set_yticks(ax.get_yticks()); ax.set_yticklabels([str(abs(int(t))) for t in ax.get_yticks()])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(L.fig_path("Fig3_deg_bar.pdf")); fig.savefig(L.fig_path("Fig3_deg_bar.png"), dpi=300)
print("Fig3: up", up, "down", dn)

# ------------------------------------------------------------------ Fig S3 heatmap
mat = np.full((len(genes), len(ORDER)), np.nan)
for i, g in enumerate(genes):
    for j, c in enumerate(ORDER):
        l, p = val(g, c)
        if l is not None:
            mat[i, j] = l
fig, ax = plt.subplots(figsize=(8.5, 0.16 * len(genes) + 1.8))
v = 4.0
im = ax.imshow(np.clip(mat, -v, v), cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
for i, g in enumerate(genes):
    for j, c in enumerate(ORDER):
        if np.isnan(mat[i, j]):
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#e8e8e8", edgecolor="none"))
        elif sig(g, c):
            ax.text(j, i, "•", ha="center", va="center", fontsize=6, color="black")
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([lab(c) for c in ORDER], rotation=60, ha="right", fontsize=7)
ax.set_yticks(range(len(genes))); ax.set_yticklabels([name(g) for g in genes], fontsize=4.8)
for t, g in zip(ax.get_yticklabels(), genes):
    t.set_color(L.SUBFAM_COLORS.get(s1[g]["subfamily"], "#7f7f7f"))
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02); cb.set_label("log2 fold change (clipped at ±4)", fontsize=7); cb.ax.tick_params(labelsize=6)
# subfamily separators
prev = None
for i, g in enumerate(genes):
    if prev is not None and s1[g]["subfamily"] != prev:
        ax.axhline(i - 0.5, color="black", lw=0.6)
    prev = s1[g]["subfamily"]
fig.tight_layout(); fig.savefig(L.fig_path("FigS3_expression_heatmap.pdf")); fig.savefig(L.fig_path("FigS3_expression_heatmap.png"), dpi=300)

# ------------------------------------------------------------------ Fig S4 UpSet-style overlap by stress type
types = ["Cold", "Heat", "Drought", "Salt", "Flooding", "Sclerotinia", "Orobanche"]
TYPE_COL = {"Cold": "#1f77b4", "Heat": "#d62728", "Drought": "#ff7f0e", "Salt": "#9467bd", "Flooding": "#17becf", "Sclerotinia": "#2ca02c", "Orobanche": "#8c564b"}
sets = {t: {g for g in genes for c in ORDER if STRESS_TYPE[c] == t and sig(g, c)} for t in types}
member = {g: frozenset(t for t in types if g in sets[t]) for g in genes}
combos = collections.Counter(m for m in member.values() if m)
combos = sorted(combos.items(), key=lambda kv: (-kv[1], -len(kv[0])))[:25]
fig = plt.figure(figsize=(11, 5.2))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 5], height_ratios=[2, 1.6], hspace=0.05, wspace=0.12)
axb = fig.add_subplot(gs[0, 1]); axm = fig.add_subplot(gs[1, 1], sharex=axb); axs = fig.add_subplot(gs[1, 0], sharey=axm)
xs = np.arange(len(combos))
axb.bar(xs, [n for _, n in combos], color=[TYPE_COL[next(iter(m))] if len(m) == 1 else "#555555" for m, _ in combos])
for i, (_, n) in enumerate(combos):
    axb.text(i, n + 0.3, str(n), ha="center", va="bottom", fontsize=7)
axb.set_ylabel("Genes", fontsize=8); axb.tick_params(labelbottom=False, labelsize=7)
for s in ("top", "right"):
    axb.spines[s].set_visible(False)
for j, t in enumerate(types):
    for i, (m, _) in enumerate(combos):
        axm.plot(i, j, "o", color=TYPE_COL[t] if t in m else "#dddddd", ms=6)
    for i, (m, _) in enumerate(combos):
        if t in m:
            rows = [k for k, tt in enumerate(types) if tt in m]
            axm.plot([i, i], [min(rows), max(rows)], color="#555555", lw=1.5, zorder=0)
axm.set_yticks(range(len(types))); axm.tick_params(labelleft=False, left=False); axm.set_xticks([]); axm.invert_yaxis()
for s in ("top", "right", "bottom"):
    axm.spines[s].set_visible(False)
axs.barh(range(len(types)), [len(sets[t]) for t in types], color=[TYPE_COL[t] for t in types])
for j, t in enumerate(types):
    axs.text(len(sets[t]) + 0.5, j, str(len(sets[t])), va="center", fontsize=7)
axs.invert_xaxis(); axs.set_xlabel("Set size", fontsize=8); axs.set_yticks(range(len(types))); axs.set_yticklabels(types, fontsize=8); axs.tick_params(labelsize=7)
for s in ("top", "left"):
    axs.spines[s].set_visible(False)
fig.savefig(L.fig_path("FigS4_upset.pdf"), bbox_inches="tight"); fig.savefig(L.fig_path("FigS4_upset.png"), dpi=300, bbox_inches="tight")

# ------------------------------------------------------------------ Fig S5 dot plot of broadly responsive genes
nsig = {g: sum(1 for c in ORDER if sig(g, c)) for g in genes}
top = [g for g in sorted(genes, key=lambda g: (-nsig[g], g)) if nsig[g] >= 5]
fig, ax = plt.subplots(figsize=(8.5, 0.28 * len(top) + 1.6))
for i, g in enumerate(top):
    for j, c in enumerate(ORDER):
        l, p = val(g, c)
        if l is None:
            continue
        s = sig(g, c)
        ax.scatter(j, i, s=min(abs(l), 6) * 22 + 6, c=[l], cmap="RdBu_r", vmin=-4, vmax=4, edgecolors="black" if s else "none", linewidths=0.6, alpha=1 if s else 0.35)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([lab(c) for c in ORDER], rotation=60, ha="right", fontsize=7)
ax.set_yticks(range(len(top))); ax.set_yticklabels(["%s (%d)" % (name(g), nsig[g]) for g in top], fontsize=6.5)
for t, g in zip(ax.get_yticklabels(), top):
    t.set_color(L.SUBFAM_COLORS.get(s1[g]["subfamily"], "#7f7f7f"))
ax.invert_yaxis(); ax.set_xlim(-0.6, len(ORDER) - 0.4)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-4, 4)); cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02); cb.set_label("log2 fold change", fontsize=7); cb.ax.tick_params(labelsize=6)
fig.tight_layout(); fig.savefig(L.fig_path("FigS5_broadly_responsive.pdf")); fig.savefig(L.fig_path("FigS5_broadly_responsive.png"), dpi=300)

# ------------------------------------------------------------------ Fig S6 WGCNA module-trait heatmap
mt = L.read_tsv(os.path.join(L.RESULTS, "04_expression", "results", "wgcna", "module_trait_table.tsv"))
mm = L.read_tsv(os.path.join(L.RESULTS, "04_expression", "results", "wgcna", "aquaporin_module_membership.tsv"))
aqp_per = collections.Counter("ME" + r["module_number"] for r in mm if r["module_number"] != "0")
mods = sorted(aqp_per, key=lambda m: int(m[2:]))
traits = [t for t in ["cold", "heat", "drought", "salt"] if any(r["trait"] == t for r in mt)]
M = np.array([[float(next(r["correlation"] for r in mt if r["module"] == m and r["trait"] == t)) for t in traits] for m in mods])
Pm = np.array([[float(next(r["pvalue"] for r in mt if r["module"] == m and r["trait"] == t)) for t in traits] for m in mods])
fig, ax = plt.subplots(figsize=(4.8, 0.32 * len(mods) + 1.5))
im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
for i in range(len(mods)):
    for j in range(len(traits)):
        star = "***" if Pm[i, j] < 0.001 else "**" if Pm[i, j] < 0.01 else "*" if Pm[i, j] < 0.05 else ""
        ax.text(j, i, "%.2f%s" % (M[i, j], star), ha="center", va="center", fontsize=6.5, color="white" if abs(M[i, j]) > 0.6 else "black")
ax.set_xticks(range(len(traits))); ax.set_xticklabels([t.capitalize() for t in traits], fontsize=8)
ax.set_yticks(range(len(mods))); ax.set_yticklabels(["%s (%d AQP)" % (m, aqp_per[m]) for m in mods], fontsize=7)
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03); cb.set_label("Module–trait correlation", fontsize=7); cb.ax.tick_params(labelsize=6)
fig.tight_layout(); fig.savefig(L.fig_path("FigS6_wgcna_module_trait.pdf")); fig.savefig(L.fig_path("FigS6_wgcna_module_trait.png"), dpi=300)
print("FigS4 combos:", len(combos), "| FigS5 genes:", len(top), "| FigS6 modules:", len(mods))
print("wrote Fig3 to", L.FIG_DIR, "and FigS3, FigS4, FigS5, FigS6 to", L.SUPP_FIG_DIR)
