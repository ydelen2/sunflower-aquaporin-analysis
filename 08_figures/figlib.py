# -*- coding: utf-8 -*-
"""Shared loaders, paths and colours for the figure scripts."""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.environ.get("AQP_RESULTS", os.path.join(ROOT, "results"))
FIG_DIR = os.path.join(ROOT, "tables_and_figures", "figures")
SUPP_FIG_DIR = os.path.join(ROOT, "tables_and_figures", "supplementary_figures")


def fig_path(name):
    """Output path of a figure file: FigS* go to supplementary_figures."""
    d = SUPP_FIG_DIR if name.startswith("FigS") else FIG_DIR
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)

SUBFAM_COLORS = {"PIP": "#1f77b4", "TIP": "#2ca02c", "NIP": "#ff7f0e", "SIP": "#9467bd",
                 "XIP": "#d62728", "unassigned": "#7f7f7f", "other": "#7f7f7f"}
SPECIES = {"Ha": "Helianthus annuus", "At": "Arabidopsis thaliana", "Os": "Oryza sativa",
           "Sl": "Solanum lycopersicum", "Ls": "Lactuca sativa"}
SPECIES_COLORS = {"Ha": "#000000", "At": "#e377c2", "Os": "#17becf", "Sl": "#bcbd22", "Ls": "#8c564b"}
SPECIES_MARKERS = {"Ha": "*", "At": "s", "Os": "^", "Sl": "D", "Ls": "v"}

HA_CHR = {"NC_0354%d.2" % (32 + i): "Ha%d" % i for i in range(1, 18)}
AT_CHR = {"NC_003070.9": "At1", "NC_003071.7": "At2", "NC_003074.8": "At3", "NC_003075.7": "At4", "NC_003076.8": "At5"}
LS_CHR = {"NC_0566%d.2" % (22 + i): "Ls%d" % i for i in range(1, 10)}

NAME_RULES = [
    (re.compile(r"\bPIP\s?\d", re.I), "PIP"), (re.compile(r"\bTIP\s?\d", re.I), "TIP"),
    (re.compile(r"\bNIP\s?\d", re.I), "NIP"), (re.compile(r"\bSIP\s?\d", re.I), "SIP"),
    (re.compile(r"\bXIP\s?\d|X[- ]intrinsic", re.I), "XIP"), (re.compile(r"nodulin[- ]?26|NOD26", re.I), "NIP"),
    (re.compile(r"plasma membrane intrinsic", re.I), "PIP"), (re.compile(r"tonoplast intrinsic", re.I), "TIP"),
    (re.compile(r"small basic intrinsic|small and basic", re.I), "SIP"),
    (re.compile(r"\bPIP\b|\bTIP\b|\bNIP\b|\bSIP\b"), None),
]


def subfamily_from_name(desc):
    for rx, sub in NAME_RULES:
        m = rx.search(desc or "")
        if m:
            return sub if sub else m.group(0).upper()
    return "unassigned"


def short_name(desc):
    """'aquaporin PIP2-4 isoform X1 [Helianthus annuus]' -> 'PIP2-4'"""
    d = re.sub(r"\s*\[.*?\]\s*$", "", desc or "")
    d = re.sub(r"^LOW QUALITY PROTEIN:\s*", "", d)
    d = re.sub(r"\s+isoform X\d+$", "", d)
    d = re.sub(r"%2C", ",", d)
    d = re.sub(r",?\s*transcript variant X\d+$", "", d)
    if re.match(r"^(probable |putative )?aquaporin[- ]?(\d+|like)\b", d, flags=re.I):
        return re.sub(r"^(probable |putative )", "", d, flags=re.I)   # keep "aquaporin-5", "aquaporin-like"
    d = re.sub(r"^(probable |putative )?aquaporin[- ]?", "", d, flags=re.I)
    d = re.sub(r"^(probable |putative )", "", d, flags=re.I)
    d = d.replace("plasma membrane intrinsic protein", "PIP").replace("tonoplast intrinsic protein", "TIP")
    d = d.replace("NOD26-like intrinsic protein", "NIP").replace("small and basic intrinsic protein", "SIP")
    d = d.replace(" family protein", "").strip()
    return d or "aquaporin"


def read_fasta_headers(path):
    out = {}
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                p = line[1:].rstrip("\n").split(None, 1)
                out[p[0]] = p[1] if len(p) > 1 else ""
    return out


def read_tsv(path):
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return rows


def load_products(results=RESULTS):
    """XP accession -> product description (sunflower)"""
    for p in (os.path.join(results, "02_gene_family", "candidate_aquaporins.faa"),
              os.path.join(results, "02_gene_family", "02_blast_validation", "candidate_aquaporins.faa")):
        if os.path.exists(p):
            return read_fasta_headers(p)
    raise FileNotFoundError("candidate_aquaporins.faa not found under " + results)


def load_protein_gene(results=RESULTS):
    """XP -> dict(gene_id LOC..., chromosome, gene_start, gene_end, strand)"""
    out = {}
    for r in read_tsv(os.path.join(results, "02_gene_family", "05_gene_structure", "protein_to_gene_map.tsv")):
        r["gene_id"] = re.sub(r"^gene-", "", r["gene_id"])
        out[r["protein_id"]] = r
    return out


def load_table_s1(results=RESULTS):
    """loc_id -> row of Table S1 (gene inventory)"""
    rows = read_tsv(os.path.join(results, "09_figures", "supplementary", "Table_S1_aquaporin_genes.tsv"))
    out = {}
    for r in rows:
        sub = r.get("subfamily") or r.get("subfamily.x") or r.get("subfamily.y") or "unassigned"
        r["subfamily"] = sub
        out[r["loc_id"]] = r
    return out


def load_subfamily_by_protein(results=RESULTS):
    """XP -> subfamily. Uses 08_classification_v2/subfamily_by_tree.tsv when
    present (re-analysis), otherwise Table S1 via the protein-gene map."""
    cls = os.path.join(results, "02_gene_family", "08_classification_v2", "subfamily_by_tree.tsv")
    if os.path.exists(cls):
        return {r["protein_id"]: r["subfamily_tree"] for r in read_tsv(cls)}
    s1 = load_table_s1(results); p2g = load_protein_gene(results)
    return {xp: s1.get(g["gene_id"], {}).get("subfamily", "unassigned") for xp, g in p2g.items()}


def load_ref_names():
    """reference accession -> (title, organism), cached from NCBI"""
    out = {}
    p = os.path.join(HERE, "ref_tip_names.tsv")
    if os.path.exists(p):
        for r in read_tsv(p):
            out[r["accession"]] = (r["title"], r["organism"])
    return out


def load_chrom_lengths(results=RESULTS):
    """{'Ha1': len, ..., 'At1': ..., 'Ls1': ...}"""
    out = {}
    with open(os.path.join(results, "02_gene_family", "05_gene_structure", "chrom_sizes.txt")) as fh:
        for line in fh:
            acc, ln = line.split()[:2]
            if acc in HA_CHR:
                out[HA_CHR[acc]] = int(ln)
    p = os.path.join(HERE, "chromosome_names.tsv")
    if os.path.exists(p):
        for r in read_tsv(p):
            acc = r["accession"]
            if acc in AT_CHR: out[AT_CHR[acc]] = int(r["length"])
            if acc in LS_CHR: out[LS_CHR[acc]] = int(r["length"])
    return out


def chr_label(tagged):
    """'Ha_NC_035439.2' -> 'Ha7'; 'Ls_NC_056623.2' -> 'Ls1'; else None"""
    acc = tagged.split("_", 1)[1] if "_" in tagged else tagged
    return HA_CHR.get(acc) or AT_CHR.get(acc) or LS_CHR.get(acc)
