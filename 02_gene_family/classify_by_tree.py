#!/usr/bin/env python3
"""
Assign sunflower aquaporins to subfamilies by phylogenetic placement.

For every sunflower tip (label Ha_<protein accession>) in the IQ-TREE tree,
the K reference tips (At_, Os_, Sl_, Ls_) with the smallest patristic
distance are collected and the subfamily that the majority of them carry is
assigned. Reference subfamilies are read from the RefSeq product names in the
reference proteome FASTA headers (e.g. "aquaporin PIP2-1", "probable
aquaporin NIP5-1", "aquaporin SIP1-2", "nodulin-26"). The same name-based
call is also made for the sunflower protein itself, from the sunflower
proteome headers, and written next to the tree-based call so the two can be
compared.

Usage:
  classify_by_tree.py --tree aquaporin_tree.treefile \
      --ref arabidopsis/protein.faa --ref rice/protein.faa \
      --ref tomato/protein.faa --ref lettuce/protein.faa \
      --sunflower sunflower/protein.faa --k 5 --out subfamily_by_tree.tsv
"""
import argparse, re, sys
from collections import Counter

from Bio import Phylo

NAME_RULES = [
    (re.compile(r"\bPIP\s?\d", re.I), "PIP"),
    (re.compile(r"\bTIP\s?\d", re.I), "TIP"),
    (re.compile(r"\bNIP\s?\d", re.I), "NIP"),
    (re.compile(r"\bSIP\s?\d", re.I), "SIP"),
    (re.compile(r"\bXIP\s?\d", re.I), "XIP"),
    (re.compile(r"nodulin[- ]?26|NOD26", re.I), "NIP"),
    (re.compile(r"plasma membrane intrinsic", re.I), "PIP"),
    (re.compile(r"tonoplast intrinsic", re.I), "TIP"),
    (re.compile(r"small basic intrinsic|small and basic", re.I), "SIP"),
]


def subfamily_from_name(desc):
    for rx, sub in NAME_RULES:
        if rx.search(desc):
            return sub
    return "unassigned"


def read_headers(fasta):
    """accession -> description (header without the accession)"""
    out = {}
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].rstrip("\n").split(None, 1)
                out[parts[0]] = parts[1] if len(parts) > 1 else ""
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--ref", action="append", required=True, help="reference proteome FASTA (repeatable)")
    ap.add_argument("--sunflower", required=True, help="sunflower proteome FASTA (for product names)")
    ap.add_argument("--k", type=int, default=5, help="number of nearest reference tips (default 5)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    names = {}
    for f in a.ref:
        names.update(read_headers(f))
    ha_names = read_headers(a.sunflower)

    tree = Phylo.read(a.tree, "newick")
    tips = tree.get_terminals()
    ha_tips = [t for t in tips if t.name.startswith("Ha_")]
    ref_tips = [t for t in tips if not t.name.startswith("Ha_")]
    if not ha_tips or not ref_tips:
        sys.exit("ERROR: could not find both Ha_ tips and reference tips in the tree")

    ref_sub = {}
    for t in ref_tips:
        acc = t.name.split("_", 1)[1]
        ref_sub[t.name] = subfamily_from_name(names.get(acc, ""))
    n_unk = sum(1 for v in ref_sub.values() if v == "unassigned")
    print("reference tips: %d (%d without a subfamily in their product name)" % (len(ref_tips), n_unk), file=sys.stderr)

    # patristic distances: precompute depths to the root for speed
    depths = tree.depths()
    parent = {}
    for cl in tree.find_clades(order="level"):
        for ch in cl.clades:
            parent[ch] = cl

    def path_to_root(c):
        p = [c]
        while c in parent:
            c = parent[c]; p.append(c)
        return p

    root_paths = {t: path_to_root(t) for t in tips}

    def dist(x, y):
        px, py = root_paths[x], root_paths[y]
        sy = set(py)
        for c in px:
            if c in sy:
                return depths[x] + depths[y] - 2 * depths[c]
        return depths[x] + depths[y]

    labelled_refs = [t for t in ref_tips if ref_sub[t.name] != "unassigned"]

    with open(a.out, "w") as out:
        out.write("protein_id\tsubfamily_tree\tsupport_fraction\tnearest_reference\tnearest_distance\t"
                  "k_reference_subfamilies\tsubfamily_from_ncbi_name\tncbi_product\n")
        summary = Counter()
        for t in ha_tips:
            acc = t.name.split("_", 1)[1]
            nearest = sorted(labelled_refs, key=lambda r: dist(t, r))[: a.k]
            subs = [ref_sub[r.name] for r in nearest]
            call, votes = Counter(subs).most_common(1)[0]
            frac = votes / float(len(subs))
            desc = ha_names.get(acc, "")
            out.write("\t".join([
                acc, call, "%.2f" % frac, nearest[0].name, "%.4f" % dist(t, nearest[0]),
                ";".join(subs), subfamily_from_name(desc), desc,
            ]) + "\n")
            summary[call] += 1
        print("sunflower tips classified: %d" % len(ha_tips), file=sys.stderr)
        for sub, n in sorted(summary.items()):
            print("  %s: %d" % (sub, n), file=sys.stderr)


if __name__ == "__main__":
    main()
