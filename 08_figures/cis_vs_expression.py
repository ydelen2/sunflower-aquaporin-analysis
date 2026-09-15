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
with Benjamini-Hochberg over all tests.

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
elems = [e for e in PAIRINGS if e in elements]
cons_all = [c for c in ["Cold", "Heat", "Drought_869", "Drought_7d", "Drought_14d", "Drought_21d", "PEG_72h",
                        "Salt_NaCl", "Flooding", "Sclerotinia", "Orobanche_A", "Orobanche_B", "Orobanche_C",
                        "Orobanche_D", "Orobanche_E"] if c in contrasts]
mat = np.full((len(elems), len(cons_all)), np.nan)
lab = {}
for r in results:
    i, j = elems.index(r["element"]), cons_all.index(r["contrast"])
    o = r["odds_ratio"]
    if np.isnan(o):
        mat[i, j] = 0.0            # no contrast possible (element in all or no genes, or no DE genes)
    elif o == 0:
        mat[i, j] = -3.0
    elif np.isinf(o):
        mat[i, j] = 3.0
    else:
        mat[i, j] = max(-3.0, min(3.0, np.log2(o)))
    lab[(i, j)] = r
fig, ax = plt.subplots(figsize=(0.62 * len(cons_all) + 2.6, 0.5 * len(elems) + 1.6))
vmax = 3.0
im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
for (i, j), r in lab.items():
    txt = "%d/%d" % (r["present_DE"], r["present_DE"] + r["present_notDE"])
    star = "*" if r["fisher_q"] < 0.05 else ("+" if r["fisher_p"] < 0.05 else "")
    ax.text(j, i, txt + star, ha="center", va="center", fontsize=6.5)
for i in range(len(elems)):
    for j in range(len(cons_all)):
        if np.isnan(mat[i, j]) and (i, j) not in lab:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#f0f0f0", edgecolor="none"))
ax.set_xticks(range(len(cons_all))); XLAB = {"Drought_869": "Drought (PRJNA869183)", "Salt_NaCl": "Salt (NaCl)", "PEG_72h": "PEG 72 h", "Drought_7d": "Drought 7 d", "Drought_14d": "Drought 14 d", "Drought_21d": "Drought 21 d"}
ax.set_xticklabels([XLAB.get(c, c.replace("_", " ")) for c in cons_all], rotation=60, ha="right", fontsize=7)
ax.set_yticks(range(len(elems))); ax.set_yticklabels([ELEMENT_LABEL.get(e, e) for e in elems], fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("log2 odds ratio, DE given element present (clipped at +/-3; 0 where not estimable)", fontsize=7); cb.ax.tick_params(labelsize=6)
fig.savefig(out_prefix + ".pdf", bbox_inches="tight"); fig.savefig(out_prefix + ".png", dpi=300, bbox_inches="tight")
print("wrote", out_prefix + "_results.tsv and figure")
