# -*- coding: utf-8 -*-
"""
Physicochemical properties of the aquaporin proteins.

For every protein of verified_aquaporins.faa: length, molecular weight, isoelectric
point, GRAVY, instability index (stable when below 40), aromaticity, percentage of
charged (D, E, K, R, H) and of hydrophobic (A, V, I, L, M, F, Y, W) residues, and the
count of each amino acid. The values are computed with the ProtParam module of
Biopython.

The NPA motifs and the ar/R filter are read from the alignment by
12_npa_arR_from_alignment_v3.py, and the subfamily is assigned from the tree by
classify_by_tree.py.

Input
  verified_aquaporins.faa
Output
  04_protein_properties_v3/protein_properties.tsv

usage: python 04_protein_properties_v3.py [proteins.faa] [out.tsv]
"""
import os, re, sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.environ.get("AQP_RESULTS", os.path.join(ROOT, "results"))
faa = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RES, "02_gene_family", "verified_aquaporins.faa")
out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(RES, "02_gene_family", "04_protein_properties_v3", "protein_properties.tsv")

AA = [("Ala", "A"), ("Arg", "R"), ("Asn", "N"), ("Asp", "D"), ("Cys", "C"), ("Gln", "Q"), ("Glu", "E"), ("Gly", "G"),
      ("His", "H"), ("Ile", "I"), ("Leu", "L"), ("Lys", "K"), ("Met", "M"), ("Phe", "F"), ("Pro", "P"), ("Ser", "S"),
      ("Thr", "T"), ("Trp", "W"), ("Tyr", "Y"), ("Val", "V")]
CHARGED, HYDROPHOBIC = "DEKRH", "AVILMFYW"

seqs, name = [], None
for line in open(faa):
    if line.startswith(">"):
        name = line[1:].split()[0]
        seqs.append([name, ""])
    elif name:
        seqs[-1][1] += line.strip().rstrip("*")

header = ["protein_id", "length_aa", "molecular_weight_Da", "isoelectric_point", "gravy", "instability_index", "stability",
          "aromaticity", "charged_aa_pct", "hydrophobic_aa_pct"] + [a for a, _ in AA]
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, "w", encoding="utf-8", newline="\n") as out:
    out.write("\t".join(header) + "\n")
    for pid, seq in seqs:
        clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())
        pa = ProteinAnalysis(clean)
        instability = pa.instability_index()
        n = float(len(seq))
        row = [pid, len(seq), "%.2f" % pa.molecular_weight(), "%.2f" % pa.isoelectric_point(), "%.4f" % pa.gravy(),
               "%.2f" % instability, "stable" if instability < 40 else "unstable", "%.4f" % pa.aromaticity(),
               "%.1f" % (100.0 * sum(seq.count(a) for a in CHARGED) / n), "%.1f" % (100.0 * sum(seq.count(a) for a in HYDROPHOBIC) / n)]
        row += [seq.count(c) for _, c in AA]
        out.write("\t".join(str(x) for x in row) + "\n")
print("proteins:", len(seqs))
print("written", out_path)
