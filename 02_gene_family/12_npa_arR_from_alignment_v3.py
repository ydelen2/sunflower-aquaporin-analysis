# -*- coding: utf-8 -*-
"""
NPA motifs and aromatic/arginine (ar/R) selectivity-filter residues read from the
multiple sequence alignment.

The four ar/R positions (H2, H5, LE1, LE2) and the two NPA motifs are located once,
in a reference protein of the alignment, and then read in every sequence from the
same alignment columns. The reference is an Arabidopsis PIP2 protein; in a PIP2 the
positions are, relative to the Asn of the two NPA motifs (N1, N2):

    H2  = N1 - 20   (Phe;  Phe81 of spinach SoPIP2;1, NPA at 101-103)
    H5  = N2 - 12   (His;  His210, NPA at 222-224)
    LE1 = N2 - 3    (Thr;  Thr219)
    LE2 = N2 + 3    (Arg;  Arg225)

The script stops if the reference does not carry Phe, His, Thr and Arg at these
positions. Sequences with a gap in one of the columns get "-" for that position.

Inputs
  all_aquaporins_mafft.fa     untrimmed MAFFT alignment (05_phylogenetics/results_v2)
  subfamily_by_tree.tsv       subfamily of the sunflower proteins
Output
  npa_arR_by_alignment.tsv    one row per protein of the alignment

usage: python 12_npa_arR_from_alignment_v3.py [alignment] [subfamily table] [out.tsv] [reference tip]
"""
import collections, csv, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.environ.get("AQP_RESULTS", os.path.join(ROOT, "results"))
aln_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RES, "05_phylogenetics", "results_v2", "all_aquaporins_mafft.fa")
cls_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(RES, "02_gene_family", "08_classification_v2", "subfamily_by_tree.tsv")
out_path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(RES, "02_gene_family", "12_npa_arR_v3", "npa_arR_by_alignment.tsv")
REF = sys.argv[4] if len(sys.argv) > 4 else "At_NP_190910.1"      # Arabidopsis PIP2;1 (PIP2A)

aln = collections.OrderedDict()
name = None
for line in open(aln_path):
    if line.startswith(">"):
        name = line[1:].split()[0]; aln[name] = []
    else:
        aln[name].append(line.strip())
aln = collections.OrderedDict((k, "".join(v).upper()) for k, v in aln.items())
width = {len(v) for v in aln.values()}
assert len(width) == 1, "sequences of unequal length in the alignment"
assert REF in aln, "reference %s is not in the alignment" % REF


def columns(seq):
    """alignment column of every residue of a sequence"""
    return [i for i, c in enumerate(seq) if c != "-"]


ref = aln[REF]
ref_res = ref.replace("-", "")
ref_col = columns(ref)
npa = [m.start() for m in re.finditer("NPA", ref_res)]
assert len(npa) == 2, "the reference must contain exactly two NPA motifs, found %d" % len(npa)
n1, n2 = npa
pos = collections.OrderedDict([("H2", n1 - 20), ("H5", n2 - 12), ("LE1", n2 - 3), ("LE2", n2 + 3)])
expected = {"H2": "F", "H5": "H", "LE1": "T", "LE2": "R"}
for k, p in pos.items():
    assert ref_res[p] == expected[k], "reference %s: %s is %s, expected %s" % (REF, k, ref_res[p], expected[k])
col = {k: ref_col[p] for k, p in pos.items()}
npa_cols = {"NPA1": [ref_col[n1 + i] for i in range(3)], "NPA2": [ref_col[n2 + i] for i in range(3)]}
print("reference %s (%d residues): NPA at %d and %d; H2 %s%d, H5 %s%d, LE1 %s%d, LE2 %s%d" % (
    REF, len(ref_res), n1 + 1, n2 + 1, ref_res[pos["H2"]], pos["H2"] + 1, ref_res[pos["H5"]], pos["H5"] + 1,
    ref_res[pos["LE1"]], pos["LE1"] + 1, ref_res[pos["LE2"]], pos["LE2"] + 1))
print("alignment columns (1-based): " + ", ".join("%s %d" % (k, v + 1) for k, v in col.items()) +
      "; NPA1 %d-%d, NPA2 %d-%d" % (npa_cols["NPA1"][0] + 1, npa_cols["NPA1"][2] + 1, npa_cols["NPA2"][0] + 1, npa_cols["NPA2"][2] + 1))

sub = {}
if os.path.exists(cls_path):
    for r in csv.DictReader(open(cls_path, encoding="utf-8"), delimiter="\t"):
        sub["Ha_" + r["protein_id"]] = r["subfamily_tree"]

os.makedirs(os.path.dirname(out_path), exist_ok=True)
rows = []
for k, s in aln.items():
    res_index = {c: i for i, c in enumerate(columns(s))}          # column -> residue number (0-based)
    m1 = "".join(s[c] for c in npa_cols["NPA1"]); m2 = "".join(s[c] for c in npa_cols["NPA2"])
    f = [s[col[x]] for x in ("H2", "H5", "LE1", "LE2")]
    p1 = res_index.get(npa_cols["NPA1"][0]); p2 = res_index.get(npa_cols["NPA2"][0])
    rows.append([k.split("_", 1)[0], k.split("_", 1)[1], sub.get(k, ""), m1, "" if p1 is None else p1 + 1, m2, "" if p2 is None else p2 + 1,
                 "" if (p1 is None or p2 is None) else p2 - p1] + f + ["/".join(f)])
with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("species\tprotein_id\tsubfamily\tNPA1\tNPA1_position\tNPA2\tNPA2_position\tNPA_spacing\tH2\tH5\tLE1\tLE2\tarR_filter\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
ha = [r for r in rows if r[0] == "Ha"]
print("proteins: %d (sunflower %d)" % (len(rows), len(ha)))
for s_ in ("PIP", "TIP", "NIP", "SIP", "XIP"):
    rs = [r for r in ha if r[2] == s_]
    if not rs:
        continue
    filt = collections.Counter(r[12] for r in rs)
    mot = collections.Counter("%s/%s" % (r[3], r[5]) for r in rs)
    print("%s n=%d  filters: %s" % (s_, len(rs), ", ".join("%s %d" % kv for kv in filt.most_common())))
    print("        NPA1/NPA2: %s" % ", ".join("%s %d" % kv for kv in mot.most_common()))
print("distinct filters among the sunflower proteins: %d (without gaps: %d)" % (len({r[12] for r in ha}), len({r[12] for r in ha if "-" not in r[12]})))
print("written", out_path)
