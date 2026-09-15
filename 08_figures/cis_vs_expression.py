# -*- coding: utf-8 -*-
"""
Promoter cis-elements versus stress-responsive expression of the aquaporin
genes.

For each a-priori element/contrast pairing (e.g. ABRE with the drought and
salt contrasts, LTR with cold, HSE with heat, ARE with flooding, W-box with
the two biotic stresses), the genes are split by element presence in the 2 kb
promoter and by differential expression in that contrast (|log2FC| > 1,
adjusted P < 0.05). A Fisher exact test gives the odds ratio; a Spearman
correlation relates element copy number to |log2FC|. P values are adjusted
with Benjamini-Hochberg over all tests. The figure shows one row per tested
pairing: DE fractions with and without the element, the log2 odds ratio with
its 95% CI, and the Spearman rho.

Inputs
  cis counts : 06_cis_elements/results/cis_elements_counts.tsv (one row per
               protein isoform promoter; gene value = max over isoforms)
  DEG matrix : a TSV with gene_id and <Contrast>_LFC / <Contrast>_padj
               columns (Table S8 layout). Default: 09_figures/supplementary/
               Table_S4_aqp_deg_all_contrasts.tsv under results/, else sheet
               S8_AQP_DEGs of Supplementary_Tables.xlsx.
  gene map   : 02_gene_family/05_gene_structure/protein_to_gene_map.tsv

usage: python cis_vs_expression.py [deg_matrix.tsv] [out_prefix]
"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, spearmanr

import figlib as L

deg_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(L.RESULTS, "09_figures", "supplementary", "Table_S4_aqp_deg_all_contrasts.tsv")
if not os.path.exists(deg_path):
    deg_path = None
out_prefix = sys.argv[2] if len(sys.argv) > 2 else L.fig_path("Fig5_cis_vs_expression")
tab_prefix = os.path.join(L.ROOT, "tables_and_figures", "tables", "cis_vs_expression")

LFC_CUT, PADJ_CUT = 1.0, 0.05

# element -> list of contrasts it is expected to act in
PAIRINGS = {
    "ABRE":     ["Drought_869", "Drought_7d", "Drought_14d", "Drought_21d", "PEG_72h", "Salt_NaCl"],
    "DRE_CRT":  ["Cold", "Drought_869", "Drought_7d", "Drought_14d", "Drought_21d", "PEG_72h"],
    "LTR":      ["Cold"],
    "MBS":      ["Drought_869", "Drought_7d", "Drought_14d", "Drought_21d", "PEG_72h"],
    "HSE":      ["Heat"],
    "ARE":      ["Flooding"],
    "W_box":    ["Sclerotinia", "Orobanche_A", "Orobanche_B", "Orobanche_C", "Orobanche_D", "Orobanche_E"],
    "STRE":     ["Cold", "Heat", "Drought_14d", "Salt_NaCl", "Sclerotinia"],
}
ELEMENT_LABEL = {"DRE_CRT": "DRE/CRT", "W_box": "W-box", "G_box": "G-box"}

# ------------------------------------------------------------ cis counts
# The eight stress elements are scanned here directly from the 2 kb promoter
# FASTA with the consensus sequences given in the Methods (exact match, both
# strands), so that the table and the tests share one definition.
MOTIFS = {
    "ABRE":    ["ACGTG"],
    "DRE_CRT": ["CCGAC"],
    "W_box":   ["TTGAC"],
    "MBS":     ["CAACTG", "TAACTG"],
    "HSE":     ["AAAAAATTTC"],
    "LTR":     ["CCGAAA"],
    "ARE":     ["AAACCA"],
    "STRE":    ["AGGGG", "CCCCT"],
    "G_box":   ["CACGTG"],          # light-responsive core element, kept for Table 2
}
COMP = str.maketrans("ACGT", "TGCA")

def count_motif(seq, patterns):
    rc = seq.translate(COMP)[::-1]
    n = 0
    for pat in patterns:
        for s in (seq, rc):
            start = 0
            while True:
                k = s.find(pat, start)
                if k < 0:
                    break
                n += 1; start = k + 1
    return n

p2g = L.load_protein_gene()
prom_path = None
for cand in (os.path.join(L.RESULTS, "02_gene_family", "results", "promoters_2kb.fa"),
             os.path.join(L.RESULTS, "02_gene_family", "05_gene_structure", "promoter_sequences_2kb.fasta")):
    if os.path.exists(cand):
        prom_path = cand; break
seqs, sid = {}, None
for line in open(prom_path):
    line = line.strip()
    if line.startswith(">"):
        sid = line[1:].split("::")[0].split()[0]; seqs[sid] = []
    elif sid:
        seqs[sid].append(line.upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
elements = list(MOTIFS)
iso_counts = {p: {e: count_motif(s, MOTIFS[e]) for e in elements} for p, s in seqs.items()}
cis_gene = {}
for prot, d in iso_counts.items():
    g = p2g.get(prot, {}).get("gene_id")
    if not g:
        continue
    dg = cis_gene.setdefault(g, {e: 0 for e in elements})
    for e in elements:
        dg[e] = max(dg[e], d[e])
print("promoter sequences:", len(seqs), "| genes with promoter data:", len(cis_gene))
with open(tab_prefix + "_gene_counts.tsv", "w") as fh:
    fh.write("gene_id\t" + "\t".join(elements) + "\n")
    for g in sorted(cis_gene):
        fh.write(g + "\t" + "\t".join(str(cis_gene[g][e]) for e in elements) + "\n")
ng = len(cis_gene)
with open(tab_prefix + "_table2.tsv", "w") as fh:
    fh.write("element\tconsensus\tgenes_with_element\tpercent_of_genes\ttotal_hits_all_promoters\n")
    print("\nTable 2 style summary (this inventory):")
    for e in elements:
        k = sum(1 for g in cis_gene if cis_gene[g][e] > 0)
        hits = sum(d[e] for d in iso_counts.values())
        fh.write("%s\t%s\t%d\t%.1f\t%d\n" % (e, "/".join(MOTIFS[e]), k, 100.0 * k / ng, hits))
        print("  %-8s %-13s genes %2d (%3.0f%%)  hits %d" % (e, "/".join(MOTIFS[e]), k, 100.0 * k / ng, hits))

# ------------------------------------------------------------ DEG matrix
def load_deg_tsv(path):
    return {r["gene_id"]: r for r in L.read_tsv(path)}

def load_deg_xlsx():
    import openpyxl
    p = os.path.join(L.ROOT, "tables_and_figures", "supplementary_tables", "Supplementary_Tables.xlsx")
    rows = list(openpyxl.load_workbook(p, read_only=True)["S8_AQP_DEGs"].iter_rows(values_only=True))
    hi = next(i for i, r in enumerate(rows) if r and r[0] == "gene_id")
    hdr = rows[hi]
    return {r[0]: {hdr[k]: r[k] for k in range(len(hdr))} for r in rows[hi + 1:] if r and r[0]}

deg = load_deg_tsv(deg_path) if deg_path else load_deg_xlsx()
contrasts = sorted({re.sub(r"_LFC$", "", k) for k in next(iter(deg.values())) if k.endswith("_LFC")})
print("genes with expression data:", len(deg), "| contrasts:", len(contrasts))

def fnum(x):
    try:
        v = float(x); return None if np.isnan(v) else v
    except (TypeError, ValueError):
        return None

genes = sorted(set(cis_gene) & set(deg))
print("genes in both:", len(genes))

# ------------------------------------------------------------ tests
results = []
for elem, cons in PAIRINGS.items():
    if elem not in elements:
        print("  element not in counts table, skipped:", elem); continue
    for con in cons:
        if con not in contrasts:
            continue
        present, de, lfc, cnt = [], [], [], []
        for g in genes:
            l = fnum(deg[g].get(con + "_LFC")); p = fnum(deg[g].get(con + "_padj"))
            if l is None or p is None:
                continue
            present.append(cis_gene[g][elem] > 0)
            de.append(abs(l) > LFC_CUT and p < PADJ_CUT)
            lfc.append(abs(l)); cnt.append(cis_gene[g][elem])
        present, de = np.array(present), np.array(de)
        a = int(np.sum(present & de)); b = int(np.sum(present & ~de))
        c = int(np.sum(~present & de)); d = int(np.sum(~present & ~de))
        odds, pf = fisher_exact([[a, b], [c, d]])
        rho, ps = (spearmanr(cnt, lfc) if len(set(cnt)) > 1 else (np.nan, np.nan))
        results.append({"element": elem, "contrast": con, "n": len(de),
                        "present_DE": a, "present_notDE": b, "absent_DE": c, "absent_notDE": d,
                        "frac_DE_present": a / (a + b) if a + b else np.nan,
                        "frac_DE_absent": c / (c + d) if c + d else np.nan,
                        "odds_ratio": odds, "fisher_p": pf, "spearman_rho": rho, "spearman_p": ps})

def bh(pvals):
    p = np.array(pvals, dtype=float); n = len(p); order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(q, 1.0); return out

fq = bh([r["fisher_p"] for r in results])
sp = [r["spearman_p"] if not np.isnan(r["spearman_p"]) else 1.0 for r in results]
sq = bh(sp)
for r, q1, q2 in zip(results, fq, sq):
    r["fisher_q"] = q1; r["spearman_q"] = q2

cols = ["element", "contrast", "n", "present_DE", "present_notDE", "absent_DE", "absent_notDE",
        "frac_DE_present", "frac_DE_absent", "odds_ratio", "fisher_p", "fisher_q", "spearman_rho", "spearman_p", "spearman_q"]
with open(tab_prefix + "_results.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in results:
        fh.write("\t".join("%.4g" % r[c] if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

print("\n%-9s %-13s %4s %6s %6s %8s %8s %8s %7s %8s" % ("element", "contrast", "n", "DE|pr", "DE|ab", "OR", "Fisher p", "q", "rho", "rho q"))
for r in results:
    print("%-9s %-13s %4d %6.2f %6.2f %8.2f %8.3f %8.3f %7.2f %8.3f" % (
        r["element"], r["contrast"], r["n"], r["frac_DE_present"], r["frac_DE_absent"], r["odds_ratio"],
        r["fisher_p"], r["fisher_q"], r["spearman_rho"], r["spearman_q"]))
print("tests:", len(results), "| Fisher q < 0.05:", sum(1 for r in results if r["fisher_q"] < 0.05),
      "| Spearman q < 0.05:", sum(1 for r in results if r["spearman_q"] < 0.05))

# ------------------------------------------------------------ figure
# one row per tested pairing, grouped by element
#   A  fraction of DE genes among genes with / without the element
#   B  log2 odds ratio of DE given the element, 95% CI (Woolf, 0.5 continuity correction)
#   C  Spearman rho between element copy number and |log2FC|
import math
elems = [e for e in PAIRINGS if e in elements]
CON_ORDER = ["Cold", "Heat", "Drought_869", "Drought_7d", "Drought_14d", "Drought_21d", "PEG_72h", "Salt_NaCl", "Flooding",
             "Sclerotinia", "Orobanche_A", "Orobanche_B", "Orobanche_C", "Orobanche_D", "Orobanche_E"]
XLAB = {"Drought_869": "Drought (PRJNA869183)", "Salt_NaCl": "Salt (NaCl)", "PEG_72h": "PEG 72 h", "Drought_7d": "Drought 7 d",
        "Drought_14d": "Drought 14 d", "Drought_21d": "Drought 21 d"}
COL = dict(zip(elems, plt.cm.Dark2.colors))
rows = sorted(results, key=lambda r: (elems.index(r["element"]), CON_ORDER.index(r["contrast"])))

def woolf(a, b, c, d):
    a, b, c, d = [v + 0.5 for v in (a, b, c, d)]
    return math.log2(a * d / (b * c)), 1.96 * math.sqrt(1 / a + 1 / b + 1 / c + 1 / d) / math.log(2)

n = len(rows); y = np.arange(n)[::-1]
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 9.5), sharey=True, gridspec_kw=dict(width_ratios=[2.2, 1.6, 1.2], wspace=0.08))
for i, r in enumerate(rows):
    c = COL[r["element"]]; yy = y[i]
    a, b, cc, d = r["present_DE"], r["present_notDE"], r["absent_DE"], r["absent_notDE"]
    fp, fa = r["frac_DE_present"], r["frac_DE_absent"]
    if not np.isnan(fp) and not np.isnan(fa):
        ax1.plot([fa, fp], [yy, yy], color=c, lw=1.2, zorder=2)
    if not np.isnan(fp):
        ax1.scatter(fp, yy, s=34, color=c, edgecolor="black", linewidths=0.5, zorder=3)
    if not np.isnan(fa):
        ax1.scatter(fa, yy, s=34, facecolor="white", edgecolor=c, linewidths=1.2, zorder=3)
    if np.isnan(r["odds_ratio"]) or (a + b) == 0 or (cc + d) == 0:
        ax2.text(0, yy, "not estimable", ha="center", va="center", fontsize=6, color="#777777")
    else:
        lo, half = woolf(a, b, cc, d)
        ax2.errorbar(lo, yy, xerr=half, fmt="o", color=c, ms=4, capsize=2, lw=1, zorder=3)
        if r["fisher_q"] < 0.05:
            ax2.text(lo, yy + 0.35, "*", ha="center", fontsize=8)
    if not np.isnan(r["spearman_rho"]):
        ax3.barh(yy, r["spearman_rho"], color=c, height=0.6, zorder=2)
ax1.set_yticks(y)
ax1.set_yticklabels(["%s  (%d | %d)" % (XLAB.get(r["contrast"], r["contrast"].replace("_", " ")), r["present_DE"] + r["present_notDE"],
                                         r["absent_DE"] + r["absent_notDE"]) for r in rows], fontsize=7)
start = 0
for k_e, e in enumerate(elems):
    k = sum(1 for r in rows if r["element"] == e)
    if k == 0:
        continue
    top, bot = y[start] + 0.5, y[start + k - 1] - 0.5
    if k_e % 2 == 0:
        for ax in (ax1, ax2, ax3):
            ax.axhspan(bot, top, color="#f2f2f2", zorder=0)
    ax1.text(-0.66, (top + bot) / 2, ELEMENT_LABEL.get(e, e), transform=ax1.get_yaxis_transform(), ha="left", va="center",
             fontsize=7.5, fontweight="bold", color=COL[e])
    start += k
ax1.tick_params(axis="y", pad=4, length=0)
ax1.set_xlim(-0.02, 1.02); ax1.set_xlabel("Fraction of genes differentially expressed", fontsize=8)
ax1.scatter([], [], s=34, color="#555555", edgecolor="black", label="genes with the element")
ax1.scatter([], [], s=34, facecolor="white", edgecolor="#555555", linewidths=1.2, label="genes without the element")
ax1.legend(fontsize=6.5, loc="lower right", frameon=False)
ax1.set_title("A", loc="left", fontweight="bold")
ax2.axvline(0, color="gray", lw=0.6, ls="--"); ax2.set_xlim(-5, 5)
ax2.set_xlabel("log2 odds ratio of DE\ngiven the element (95% CI)", fontsize=8); ax2.set_title("B", loc="left", fontweight="bold")
ax3.axvline(0, color="gray", lw=0.6, ls="--"); ax3.set_xlim(-0.5, 0.5)
ax3.set_xlabel("Spearman rho\n(element copies vs |log2FC|)", fontsize=8); ax3.set_title("C", loc="left", fontweight="bold")
for ax in (ax1, ax2, ax3):
    ax.tick_params(axis="x", labelsize=7); ax.set_ylim(-0.7, n - 0.3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.savefig(out_prefix + ".pdf", bbox_inches="tight"); fig.savefig(out_prefix + ".png", dpi=300, bbox_inches="tight")
print("wrote", out_prefix + "_results.tsv and figure")
