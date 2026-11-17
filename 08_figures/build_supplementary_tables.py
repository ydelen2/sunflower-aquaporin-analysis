# -*- coding: utf-8 -*-
"""
Build Supplementary_Tables.xlsx (Additional file 1, Tables S1 to S11 and S13) from
the pipeline outputs under results/.

usage: python build_supplementary_tables.py
"""
import os, re, sys, csv, collections
import openpyxl
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter
import figlib as L

OUT = os.path.join(L.ROOT, "tables_and_figures", "supplementary_tables", "Supplementary_Tables.xlsx")
P = lambda *p: os.path.join(L.RESULTS, *p)
s1 = L.load_table_s1()
sub = L.load_subfamily_by_protein()
p2g = L.load_protein_gene()
prods = L.load_products()

def tsv(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], rows[1:]

def num(x):
    try:
        if x in ("NA", "", None):
            return x
        f = float(x)
        if f == 0.0 and re.search(r"[1-9]", str(x).lower().split("e")[0]):
            return x          # too small for a floating-point number (for example 2.4e-1574): keep the text
        return int(f) if f.is_integer() and "." not in str(x) else f
    except ValueError:
        return x

wb = openpyxl.Workbook(); wb.remove(wb.active)

# Labels of the result files that are written differently in the tables. PRJNA1041959 is a soil-drought
# experiment (irrigation withheld, sampled after 72 h) and the five groups A to E of PRJNA850121 are sampling
# times after inoculation with Orobanche cumana; the result files keep the group labels of the count matrices.
HOURS = {"A": "4h", "B": "12h", "C": "24h", "D": "48h", "E": "72h"}
RELABEL = [("drought_PEG_72h", "drought_72h"), ("PEG_72h", "Drought_72h"), ("PEG 72 h", "Drought 72 h"), ("1041:PEG", "1041:d72h"),
           ("PEG +14%", "Drought 72 h +14%"), ("~condition_per_stage", "~condition_per_time_point"), ("~genotype_fixed+condition", "~genotype+condition"),
           ("Stage A +52%, B +42%, C +57%, D +127%, E +35%", "4 h +52%, 12 h +42%, 24 h +57%, 48 h +127%, 72 h +35%")]
for k, v in HOURS.items():
    RELABEL += [("parasitized_" + k, "orobanche_" + v), ("Orobanche_" + k, "Orobanche_" + v), ("850:p" + k, "850:" + v),
                ("Orobanche cumana, stage " + k, "Orobanche cumana, %s h after inoculation" % v[:-1])]


def relabel(v):
    if isinstance(v, str):
        for a, b in RELABEL:
            if a in v:
                v = v.replace(a, b)
    return v


def sheet(name, title, header, rows, note=None):
    header = [relabel(x) for x in header]
    rows = [[relabel(x) for x in r] for r in rows]
    ws = wb.create_sheet(name)
    ws.cell(1, 1, title).font = Font(bold=True)
    if note:
        ws.cell(2, 1, note).font = Font(italic=True)
    for j, h in enumerate(header, 1):
        ws.cell(3, j, h).font = Font(bold=True)
    for i, r in enumerate(rows, 4):
        for j, v in enumerate(r, 1):
            ws.cell(i, j, num(v))
    # column widths from the longest entry of each column; long text is wrapped
    for j, h in enumerate(header, 1):
        w = max([len(str(h))] + [len(str(r[j - 1])) for r in rows[:1000] if j - 1 < len(r) and r[j - 1] is not None])
        ws.column_dimensions[get_column_letter(j)].width = min(w + 2, 60)
        if w > 60:
            for i in range(4, len(rows) + 4):
                ws.cell(i, j).alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "A4"
    return ws

# S1 gene inventory
h, r = tsv(P("09_figures", "supplementary", "Table_S1_aquaporin_genes.tsv"))
sheet("S1_Aquaporin_Genes", "Table S1. The %d sunflower aquaporin genes: subfamily by phylogenetic placement (subfamily) and by RefSeq product name (subfamily_from_name), tree support of the placement, position, gene structure of the representative mRNA and encoded isoforms." % len(r), h, r)

# S2 rejected candidates
h, r = tsv(P("02_gene_family", "03_domain_verification_v2", "failed_verification.tsv"))
rej = [x for x in r if x[1] == "yes"]
h2 = ["protein_id", "product", "length_aa", "tm_helices_DeepTMHMM", "reason_for_exclusion"]
r2 = []
tm = {x[0]: x for x in tsv(P("02_gene_family", "03_domain_verification_v2", "tm_helix_predictions.tsv"))[1]}
for x in rej:
    pid = x[0]
    r2.append([pid, prods.get(pid, "").replace(" [Helianthus annuus]", ""), tm.get(pid, ["", "", "", ""])[3], x[2], "fewer than four transmembrane helices predicted"])
sheet("S2a_Excluded_Candidates", "Table S2a. Candidates with a PF00230 domain that were excluded because DeepTMHMM predicted fewer than four transmembrane helices.", h2, r2)

# S2b loci found by tblastn of the aquaporin proteins against the genome sequence
h, r = tsv(P("02_gene_family", "11_tblastn_v3", "tblastn_loci.tsv"))
ci = h.index("sequence")
for x in r:
    x[ci] = L.HA_CHR.get(x[ci], x[ci])
h[ci] = "chromosome"
n_cls = collections.Counter(x[h.index("class")] for x in r)
sheet("S2b_tblastn_Loci", "Table S2b. Loci found by the tblastn search of the genome sequence with the %d aquaporin proteins (E-value <= 1e-10; hits on the same strand within 10 kb merged): position, best-matching protein, fraction of that protein covered, and overlap with annotated genes. Class inventory: locus of one of the 87 aquaporin genes (%d loci); annotated: overlaps another annotated gene (%d); unannotated: no gene model (%d)." % (len(sub), n_cls["inventory"], n_cls["annotated"], n_cls["unannotated"]), h, r)

# S5 physicochemical and structural features
h, r = tsv(P("02_gene_family", "04_protein_properties_v3", "protein_properties.tsv"))
mh, mr = tsv(P("02_gene_family", "12_npa_arR_v3", "npa_arR_by_alignment.tsv"))
motif_cols = ["NPA1", "NPA1_position", "NPA2", "NPA2_position", "NPA_spacing", "H2", "H5", "LE1", "LE2", "arR_filter"]
motif = {x[mh.index("protein_id")]: [x[mh.index(c)] for c in motif_cols] for x in mr if x[mh.index("species")] == "Ha"}
k = h.index("aromaticity") + 1
h = h[:1] + ["gene_id", "subfamily_tree"] + h[1:k] + motif_cols + h[k:]
r = [[x[0], p2g.get(x[0], {}).get("gene_id", ""), sub.get(x[0], "")] + x[1:k] + motif[x[0]] + x[k:] for x in r]
r.sort(key=lambda x: (x[1], x[0]))
sheet("S5_Physicochemical", "Table S5. Physicochemical properties, NPA motifs and ar/R selectivity filter residues of the %d aquaporin protein isoforms. The two NPA motifs and the four filter positions (H2, H5, LE1, LE2) were read from the columns of the multiple sequence alignment that hold the two motifs and Phe87, His216, Thr225 and Arg231 of Arabidopsis PIP2;1 (NP_190910.1). NPA1_position and NPA2_position give the residue number of the first motif residue, NPA_spacing the distance between the two motifs; a dash marks a position that is absent from the protein." % len(r), h, r)

# S7a duplication, S7c Ka/Ks (one source table)
h, r = tsv(P("02_gene_family", "10_duplication_synteny_v3", "duplication_synteny_class.tsv"))
CLASS = {"tandem": "tandem", "proximal": "proximal", "segmental_anchor": "segmental (collinear anchor pair)",
         "segmental_in_block": "segmental (matching positions in one collinear block)", "dispersed": "dispersed"}
ti, si = h.index("type"), h.index("synteny_class")
for x in r:
    x[ti] = CLASS[x[si]]
h[ti] = "duplication_class"
drop = ("aligned_codons", "Ka", "Ks", "Ka_Ks", "selection", "note", "synteny_class")
sheet("S7a_Duplication", "Table S7a. Paralogous aquaporin gene pairs (protein identity >= 70%%, query coverage >= 70%%; one pair per gene pair) and their duplication class (tandem <= 200 kb, proximal 200 kb to 1 Mb on the same chromosome; other pairs segmental when supported by a collinear block, as an anchor pair or at matching positions inside one block, and dispersed otherwise; block_id and block_anchors give the MCScanX block and its number of anchor pairs); n = %d." % len(r),
      [c for c in h if c not in drop],
      [[v for c, v in zip(h, x) if c not in drop] for x in r])
sheet("S7c_Ka_Ks", "Table S7c. Ka, Ks and Ka/Ks (Nei-Gojobori with Jukes-Cantor correction) for the same pairs; NA where fewer than 30 codons aligned, synonymous sites were saturated or Ks was zero.",
      ["gene1", "gene2", "protein1", "protein2", "duplication_class", "aligned_codons", "Ka", "Ks", "Ka_Ks", "selection", "note"],
      [[x[h.index(c)] for c in ["gene1", "gene2", "protein1", "protein2", "duplication_class", "aligned_codons", "Ka", "Ks", "Ka_Ks", "selection", "note"]] for x in r])

# S7b synteny
rows = []
for comp, fn, label in (("Ha_self", "Ha_self_intra_aqp_synteny.tsv", "Intra-genomic"), ("Ha_Ls", "Ha_Ls_inter_aqp_synteny.tsv", "Interspecies"), ("Ha_At", "Ha_At_inter_aqp_synteny.tsv", "Interspecies")):
    h, r = tsv(P("07_synteny", "results_v3", fn))
    for x in r:
        d = dict(zip(h, x))
        if label == "Interspecies" and d["species1"] == d["species2"]:
            continue
        rows.append([comp, label, L.chr_label(d["chr1"]) or d["chr1"], L.chr_label(d["chr2"]) or d["chr2"], d["species1"], d["gene1"], d["prot1"], d["species2"], d["gene2"], d["prot2"]])
sheet("S7b_Synteny", "Table S7b. Collinear aquaporin gene pairs from MCScanX: intra-genomic pairs (both genes in the inventory) and interspecies pairs with lettuce and Arabidopsis (sunflower member in the inventory); n = %d." % len(rows),
      ["Comparison", "Type", "Chr1", "Chr2", "Species1", "Gene1", "Protein1", "Species2", "Gene2", "Protein2"], rows)

# S8 DEG matrix
h, r = tsv(P("09_figures", "supplementary", "Table_S4_aqp_deg_all_contrasts.tsv"))
sheet("S8_AQP_DEGs", "Table S8. DESeq2 log2 fold change and adjusted P value of every aquaporin gene in the 16 stress contrasts (NA: no estimate).", h, r)

# S9, S10 WGCNA
h, r = tsv(P("04_expression", "results", "wgcna_v3", "aquaporin_module_membership.tsv"))
sheet("S9_WGCNA_Membership", "Table S9. WGCNA module assignment of the aquaporin genes present in the PRJNA869183 expression matrix (module 0 = unassigned) and module membership (kME) values.", h, r)
h, r = tsv(P("04_expression", "results", "wgcna_v3", "module_trait_table.tsv"))
traits = [t for t in ["cold", "heat", "drought", "rewatering", "salt", "SA"] if any(x[1] == t for x in r)]; mods = sorted({x[0] for x in r}, key=lambda m: int(m[2:]))
table = {(x[0], x[1]): (x[2], x[3]) for x in r}
sheet("S10_Module_Trait", "Table S10. Module-trait correlations (Pearson r and P value) between the eigengene of every WGCNA module of PRJNA869183 and the six treatment indicators (cold, heat, drought, rewatering, salt = NaCl, SA = the groups labeled SA1 to SA4).",
      ["module"] + ["r_%s" % t for t in traits] + ["p_%s" % t for t in traits],
      [[m] + [table.get((m, t), ("", ""))[0] for t in traits] + [table.get((m, t), ("", ""))[1] for t in traits] for m in mods])

# S11 GO enrichment
h, r = tsv(P("04_expression", "results", "enrichment", "GO_all_contrasts_v2.tsv"))
n_uni = sum(1 for line in open(P("04_expression", "results", "enrichment", "GO_all_contrasts_v2_universe.txt")) if line.strip())
sheet("S11_GO_Enrichment", "Table S11. GO terms enriched in the DEG set of each of the 16 stress contrasts (over-representation analysis with clusterProfiler against one common universe of {:,} genes; {:,} term-by-contrast combinations, {:,} distinct terms).".format(n_uni, len(r), len({x[h.index("ID")] for x in r})), h, r)

# S6 MEME motifs (protein)
h, r = tsv(P("05_phylogenetics", "results", "meme_motifs", "motif_summary.tsv"))
sheet("S6_MEME_Motifs", "Table S6. Conserved protein motifs identified by MEME (zoops model, 20 motifs searched, ten retained) in the %d aquaporin isoforms." % len(sub), h, r[:10])

# S4 LOC-TAIR mapping
h, r = tsv(P("09_figures", "supplementary", "Table_S9_loc_tair_mapping.tsv"))
sheet("S4_LOC_TAIR_Mapping", "Table S4. LOC-to-TAIR ortholog mapping used for GO enrichment (BLASTp best hit, E-value <= 1e-5, identity > 30%, query coverage > 50%; {:,} pairs).".format(len(r)), h, r)

# S3 DESeq2 design rationale; the PRJNA492303 row reflects the 96-run analysis
h, r = tsv(P("09_figures", "supplementary", "Table_S10_design_rationale.tsv"))
s = tsv(P("04_expression", "results", "deseq2", "PRJNA492303_fixed", "DEG_summary.tsv"))[1][0]
for x in r:
    if x[0] == "PRJNA492303":
        x[1] = "~tissue + condition (46 of the 96 runs)"
        x[2] = "~genotype + tissue + age + condition (all 96 runs)"
        x[3] = "The first analysis used an arbitrary block of 46 consecutive runs in which HA351 occurred only among controls; the experiment is balanced (2 genotypes x 2 treatments x 2 tissues x 3 ages x 4 replicates)"
        x[4] = "Flooding 1,671 -> %s genome-wide DEGs (%s up, %s down); aquaporin DEGs 5 -> 3" % (s[3], s[4], s[5])
    if x[0] == "PRJNA869183":
        x[2] = "~group, fitted separately for each contrast on the libraries of its two groups (Table S3b)"
sheet("S3a_Design_Rationale", "Table S3a. DESeq2 design rationale and correction summary for all six BioProjects.", h, r)

# S3b contrasts, groups and replicates; DEG numbers are those of the fits in 10_deseq2_designs_v3.sh
V3 = lambda f: P("04_expression", "results", "deseq2_designs_v3", f)
h, r = tsv(os.path.join(L.ROOT, "04_expression", "contrast_definitions.tsv"))
rh, rr = tsv(V3("reproduction_report.tsv"))
rep = {}
for x in rr:
    d = dict(zip(rh, x))
    if d["deg_new"] == d["deg_published"] and float(d["lfc_correlation"]) > 0.999999 and float(d["lfc_max_abs_diff"]) < 0.001:
        rep[(d["dataset"], d["contrast"])] = d
ah, ar = tsv(P("04_expression", "results", "aquaporin_deg_v2", "aquaporin_deg_counts.tsv"))
aq = {x[0]: x for x in ar}
KEY = {"Cold": "cold", "Heat": "heat", "Drought_869": "drought", "Salt_NaCl": "salt", "Rehydration": "rehydration", "Drought_7d": "drought_7d",
       "Drought_14d": "drought_14d", "Drought_21d": "drought_21d", "PEG_72h": "drought_PEG_72h", "Flooding": "flooding", "Sclerotinia": "sclerotinia"}
rows = []
for x in r:
    d = dict(zip(h, x))
    k = KEY.get(d["label"], "parasitized_" + d["label"][-1])
    g = rep[(d["bioproject"], k)]
    rows.append(x + [g["tested_new"], g["deg_new"], g["up_new"], g["down_new"], aq[d["label"]][3], aq[d["label"]][1], aq[d["label"]][2]])
sheet("S3b_Contrasts_Replicates", "Table S3b. The %d contrasts: model, groups compared, number of libraries per group and biological replicates, with the genome-wide and aquaporin DEG numbers (|log2FC| > 1, adjusted P < 0.05). The fits were re-run from the count matrix and reproduce the DEG tables used in the paper. In the flooding contrast of PRJNA492303 the treated libraries of the third sampling time (d30) come from plants on the sixth day of recovery after flooding." % len(rows),
      h + ["genes_tested", "DEGs", "DEGs_up", "DEGs_down", "aquaporin_DEGs", "aquaporin_up", "aquaporin_down"], rows)

# S3c interaction tests, S3d treatment contrasts within genotype, tissue or age
h, r = tsv(V3("interaction_tests_summary.tsv"))
r = [[re.sub(r" {2,}", " ", v) for v in x] for x in r]      # the formula text holds runs of blanks from deparse()
sheet("S3c_Interaction_Tests", "Table S3c. Likelihood-ratio tests of the interaction of treatment with genotype, tissue or plant age (full model against the model without the tested term; Benjamini-Hochberg adjusted P < 0.05), genome-wide and for the aquaporin genes. In PRJNA492303 the ages d16, d24 and d30 are the three sampling times: the first and the eighth day of flooding and the sixth day of recovery after flooding.", h, r)
h, r = tsv(V3("stratified_contrasts_summary.tsv"))
sheet("S3d_Stratified_Contrasts", "Table S3d. Treatment contrasts estimated within each genotype, tissue or plant age (Wald test in a model in which treatment is combined with that factor; |log2FC| > 1, adjusted P < 0.05), genome-wide and for the aquaporin genes. In PRJNA492303 the ages d16, d24 and d30 are the three sampling times: the first and the eighth day of flooding and the sixth day of recovery after flooding.", h, r)

# S3e per-contrast and joint multiple-testing correction for the aquaporin genes
h, r = tsv(V3("joint_fdr_aquaporin.tsv"))
jh, jr = tsv(V3("joint_fdr_summary.tsv"))
js = {x[0]: dict(zip(jh, x)) for x in jr}
r.sort(key=lambda x: (x[0], x[1]))
sheet("S3e_Multiple_Testing", "Table S3e. Aquaporin results under three corrections for multiple testing: padj, Benjamini-Hochberg within the contrast (de_within); padj_joint, one Benjamini-Hochberg correction over all tests of the 16 contrasts (de_joint); de_bonferroni_contrasts, padj < 0.05/16. A call also requires |log2FC| > 1. Genes with at least one call: %s, %s and %s of %s." % (
          js["per_contrast_BH"]["genes_ge1"], js["joint_BH_all_tests"]["genes_ge1"], js["bonferroni_over_contrasts"]["genes_ge1"], js["per_contrast_BH"]["aquaporin_genes_total"]), h, r)

# S13 baseline TPM of the control libraries by BioProject and tissue
h, r = tsv(P("04_expression", "results", "aquaporin_expression_v2", "aquaporin_TPM_control_means.tsv"))
n_lib = dict(zip(h, r[0])); r = r[1:]
ratio = {(x[0], x[1]): x[4] for x in tsv(P("04_expression", "results", "aquaporin_expression_v2", "tissue_summary.tsv"))[1]}
order = ["PRJNA869183_leaf", "PRJNA797473_leaf", "PRJNA908908_leaf", "PRJNA1041959_leaf", "PRJNA492303_leaf", "PRJNA1041959_root", "PRJNA492303_root", "PRJNA850121_root"]
grp = [g for g in order if g in h] + [g for g in h[1:] if g not in order]
idx = [h.index(g) for g in grp]
hdr13 = ["gene_id", "product", "subfamily"] + ["%s (n=%s)" % (g, n_lib[g]) for g in grp] + ["log2_root_leaf_PRJNA1041959", "log2_root_leaf_PRJNA492303"]
rows13 = [[x[0], L.short_name(s1.get(x[0], {}).get("product", "")), s1.get(x[0], {}).get("subfamily", "")] + [x[i] for i in idx] + [ratio.get((x[0], "PRJNA1041959"), ""), ratio.get((x[0], "PRJNA492303"), "")] for x in r]
sheet("S13_Baseline_TPM", "Table S13. Mean TPM of the %d aquaporin genes in the control libraries (untreated and pre-treatment samples) of each BioProject and tissue, and the root/leaf log2 ratio of the control means in the two BioProjects that sampled both tissues (0.1 added to each mean before the ratio)." % len(rows13), hdr13, rows13)

wb._sheets.sort(key=lambda ws: (int(re.match(r"S(\d+)", ws.title).group(1)), ws.title))
wb.properties.creator = "Yavuz Delen"; wb.properties.lastModifiedBy = "Yavuz Delen"
os.makedirs(os.path.dirname(OUT), exist_ok=True)
wb.save(OUT)
print("saved", OUT)
for ws in wb.worksheets:
    print("  %-24s rows %d" % (ws.title, ws.max_row))
