#!/bin/bash
#SBATCH --job-name=aqp_dup_kaks
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/09_dup_kaks_v2_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/09_dup_kaks_v2_%j.err
#SBATCH --mail-type=END,FAIL

# ============================================================================
# 09_duplication_kaks_v2.sh
# Gene-level duplication classification and Ka/Ks for the aquaporin family.
#
# In the first analysis the gene-level tables (duplication_analysis.tsv and
# kaks_analysis.tsv, the sources of Tables S7a and S7c) were assembled from the
# 06_chromosomal_location.sh outputs in an unscripted step. This script makes
# that step reproducible and applies one consistent rule set:
#
#   pairs      : self-BLASTP hits from 06 (E <= 1e-10), protein identity >= 70 %
#                and query coverage >= 70 %, collapsed to one pair per gene pair
#                (highest bit score); same-gene isoform pairs are dropped
#   tandem     : same chromosome, gene starts <= 200 kb apart
#   proximal   : same chromosome, 200 kb < distance <= 1 Mb
#   segmental  : all other pairs
#   Ka/Ks      : MAFFT pairwise protein alignment, back-translation against the
#                RefSeq CDS, Nei-Gojobori (1986) with Jukes-Cantor correction,
#                pairs with fewer than 30 aligned codons excluded, Ks reported
#                as NA when the synonymous p-distance is saturated (>= 0.75)
#
# Inputs (from 02_gene_family/06_chromosomal_location, after 05 and 06 v2 runs):
#   self_blastp.outfmt6, chromosome_locations.tsv
#   02_gene_family/verified_aquaporins.faa
#   01_references/GCF_002127325.2_HanXRQr2.0_cds_from_genomic.fna
# Outputs:
#   02_gene_family/09_duplication_kaks_v2/duplication_kaks.tsv (one row per gene pair)
#   09_figures/supplementary/duplication_analysis.tsv and kaks_analysis.tsv (same
#   layout as the first version, regenerated from the table above)
# ============================================================================

set -eo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
source "${PROJ_DIR}/config.sh"

GF="${PROJ_DIR}/02_gene_family"
LOC_DIR="${GF}/06_chromosomal_location"
OUT_DIR="${GF}/09_duplication_kaks_v2"
SUPP="${PROJ_DIR}/09_figures/supplementary"
BLAST="${LOC_DIR}/self_blastp.outfmt6"
CHROM="${LOC_DIR}/chromosome_locations.tsv"
PROT="${GF}/verified_aquaporins.faa"
CDS="${PROJ_DIR}/01_references/GCF_002127325.2_HanXRQr2.0_cds_from_genomic.fna"

MIN_IDENT=70
MIN_QCOV=70
TANDEM_BP=200000
PROXIMAL_BP=1000000
MIN_CODONS=30

mkdir -p "${OUT_DIR}/tmp" "${SUPP}" "${PROJ_DIR}/logs"

module purge
try_modules() { local m; for m in "$@"; do module load "${m}" 2>/dev/null && { echo "  loaded module: ${m}"; return 0; }; done; echo "  note: no module loaded from: $*"; return 0; }
try_modules "mafft/7.526" "mafft" "MAFFT"
try_modules "miniforge/24.5" "miniforge"
if ! command -v mafft &>/dev/null; then
    for env in "${PROJ_DIR}/envs/mafft_env"; do [[ -x "${env}/bin/mafft" ]] && export PATH="${env}/bin:${PATH}"; done
fi
command -v mafft &>/dev/null || { echo "ERROR: mafft not found" >&2; exit 1; }

for f in "${BLAST}" "${CHROM}" "${PROT}" "${CDS}"; do
    [[ -s "${f}" ]] || { echo "ERROR: missing input ${f}" >&2; exit 1; }
done

python3 - "${BLAST}" "${CHROM}" "${PROT}" "${CDS}" "${OUT_DIR}" "${SUPP}" \
          "${MIN_IDENT}" "${MIN_QCOV}" "${TANDEM_BP}" "${PROXIMAL_BP}" "${MIN_CODONS}" << 'PYEOF'
import sys, os, re, math, subprocess, collections
blast_f, chrom_f, prot_f, cds_f, out_dir, supp_dir = sys.argv[1:7]
MIN_IDENT, MIN_QCOV = float(sys.argv[7]), float(sys.argv[8])
TANDEM_BP, PROXIMAL_BP, MIN_CODONS = int(sys.argv[9]), int(sys.argv[10]), int(sys.argv[11])
tmp = os.path.join(out_dir, 'tmp')

def read_fasta(path, key=lambda h: h.split()[0]):
    seqs, k = {}, None
    for line in open(path):
        if line.startswith('>'):
            k = key(line[1:].strip()); seqs[k] = []
        elif k is not None:
            seqs[k].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}

# protein -> gene / position
p2g, gpos = {}, {}
for i, line in enumerate(open(chrom_f)):
    f = line.rstrip('\n').split('\t')
    if i == 0:
        h = f; continue
    r = dict(zip(h, f))
    g = re.sub(r'^gene-', '', r['gene_id'])
    p2g[r['protein_id']] = g
    gpos[g] = (r['chromosome'], int(r['gene_start']), int(r['gene_end']), r['strand'])
chr_map = {'NC_0354%d.2' % (32 + i): 'Ha%d' % i for i in range(1, 18)}

# self-BLAST -> best protein pair per gene pair
def clean_id(x):
    # makeblastdb -parse_seqids writes subject ids as ref|XP_...|; strip the wrapper
    m = re.match(r'^(?:\w+\|)?([A-Za-z]{2}_\d+\.\d+)\|?$', x)
    return m.group(1) if m else x
best = {}
for line in open(blast_f):
    f = line.rstrip('\n').split('\t')
    if len(f) < 12 or f[0] == f[1]:
        continue
    q, s, pid, bits = clean_id(f[0]), clean_id(f[1]), float(f[2]), float(f[11])
    qcov = float(f[12]) if len(f) > 12 else 100.0
    if pid < MIN_IDENT or qcov < MIN_QCOV:
        continue
    g1, g2 = p2g.get(q), p2g.get(s)
    if not g1 or not g2 or g1 == g2:
        continue
    key = tuple(sorted([g1, g2]))
    if key not in best or bits > best[key]['bits']:
        best[key] = {'p1': q if key[0] == g1 else s, 'p2': s if key[0] == g1 else q,
                     'ident': pid, 'qcov': qcov, 'bits': bits}
print('gene pairs after identity/coverage filter:', len(best))

def classify(g1, g2):
    c1, s1, e1, _ = gpos[g1]; c2, s2, e2, _ = gpos[g2]
    if c1 != c2:
        return 'segmental', 'NA'
    d = abs(s1 - s2)
    if d <= TANDEM_BP:
        return 'tandem', d
    if d <= PROXIMAL_BP:
        return 'proximal', d
    return 'segmental', d

# Ka/Ks machinery (Nei-Gojobori with Jukes-Cantor), same as 06_chromosomal_location.sh
BASES = 'TCAG'
AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}

def syn_sites(codon):
    aa = CODON[codon]; s = 0.0
    for pos in range(3):
        tot = syn = 0
        for b in BASES:
            if b == codon[pos]:
                continue
            new = codon[:pos] + b + codon[pos + 1:]
            if CODON[new] == '*':
                continue
            tot += 1
            if CODON[new] == aa:
                syn += 1
        if tot:
            s += syn / tot
    return s, 3.0 - s

def nei_gojobori(c1s, c2s):
    S = N = Sd = Nd = 0.0; n = 0
    for i in range(len(c1s) // 3):
        c1, c2 = c1s[3 * i:3 * i + 3], c2s[3 * i:3 * i + 3]
        if '-' in c1 or '-' in c2 or c1 not in CODON or c2 not in CODON or CODON[c1] == '*' or CODON[c2] == '*':
            continue
        n += 1
        s1, n1 = syn_sites(c1); s2, n2 = syn_sites(c2)
        S += (s1 + s2) / 2; N += (n1 + n2) / 2
        if c1 == c2:
            continue
        diff = [p for p in range(3) if c1[p] != c2[p]]
        if len(diff) == 1:
            if CODON[c1] == CODON[c2]: Sd += 1
            else: Nd += 1
        else:
            # average over all substitution pathways that avoid stop codons
            import itertools
            paths = []
            for order in itertools.permutations(diff):
                cur = c1; sd = nd = 0.0; ok = True
                for p in order:
                    nxt = cur[:p] + c2[p] + cur[p + 1:]
                    if CODON[nxt] == '*':
                        ok = False; break
                    if CODON[cur] == CODON[nxt]: sd += 1
                    else: nd += 1
                    cur = nxt
                if ok:
                    paths.append((sd, nd))
            if paths:
                Sd += sum(p[0] for p in paths) / len(paths); Nd += sum(p[1] for p in paths) / len(paths)
            else:
                Nd += len(diff)
    if n < MIN_CODONS or S == 0 or N == 0:
        return n, None, None, None, 'fewer_than_%d_codons' % MIN_CODONS if n < MIN_CODONS else 'no_sites'
    pS, pN = Sd / S, Nd / N
    Ks = max(0.0, -0.75 * math.log(1 - 4.0 / 3.0 * pS)) if pS < 0.75 else None
    Ka = max(0.0, -0.75 * math.log(1 - 4.0 / 3.0 * pN)) if pN < 0.75 else None
    note = ''
    if Ks is None: note = 'Ks_saturated'
    if Ka is None: note = (note + ';' if note else '') + 'Ka_saturated'
    kaks = (Ka / Ks) if (Ka is not None and Ks not in (None, 0)) else None
    if Ks == 0 and Ka is not None:
        note = (note + ';' if note else '') + 'Ks_zero'
    return n, Ka, Ks, kaks, note

prot = read_fasta(prot_f)
def cds_key(h):
    m = re.search(r'\[protein_id=([^\]]+)\]', h)
    if m: return m.group(1)
    m = re.search(r'_cds_([A-Z]P_\d+\.\d+)_', h)
    return m.group(1) if m else h.split()[0]
cds = read_fasta(cds_f, key=cds_key)
print('proteins:', len(prot), '| CDS records mapped by protein_id:', sum(1 for p in prot if p in cds))

def align(p1, p2):
    inp = os.path.join(tmp, 'pair.faa'); out = os.path.join(tmp, 'pair.aln')
    with open(inp, 'w') as fh:
        fh.write('>%s\n%s\n>%s\n%s\n' % (p1, prot[p1], p2, prot[p2]))
    with open(out, 'w') as fo:
        subprocess.run(['mafft', '--auto', '--quiet', inp], stdout=fo, stderr=subprocess.DEVNULL, check=True, timeout=120)
    a = read_fasta(out)
    return a[p1], a[p2]

def back_translate(aln, seq):
    seq = seq.upper()
    core = aln.replace('-', '')
    if len(seq) < 3 * len(core):
        return None
    out, k = [], 0
    for ch in aln:
        if ch == '-':
            out.append('---')
        else:
            out.append(seq[3 * k:3 * k + 3]); k += 1
    return ''.join(out)

rows = []
for key in sorted(best):
    g1, g2 = key; b = best[key]
    typ, dist = classify(g1, g2)
    p1, p2 = b['p1'], b['p2']
    n = 0; Ka = Ks = kaks = None; note = ''
    if p1 in cds and p2 in cds:
        a1, a2 = align(p1, p2)
        c1, c2 = back_translate(a1, cds[p1]), back_translate(a2, cds[p2])
        if c1 and c2:
            n, Ka, Ks, kaks, note = nei_gojobori(c1, c2)
        else:
            note = 'cds_protein_length_mismatch'
    else:
        note = 'cds_missing'
    if kaks is None: sel = 'NA'
    elif kaks < 1: sel = 'purifying'
    elif kaks > 1: sel = 'positive'
    else: sel = 'neutral'
    c1n, c2n = gpos[g1][0], gpos[g2][0]
    rows.append([g1, g2, p1, p2, '%.3f' % b['ident'], '%.1f' % b['qcov'], c1n, c2n, chr_map.get(c1n, c1n), chr_map.get(c2n, c2n),
                 str(dist), typ, str(n),
                 '%.4f' % Ka if Ka is not None else 'NA', '%.4f' % Ks if Ks is not None else 'NA',
                 '%.4f' % kaks if kaks is not None else 'NA', sel, note])

hdr = ['gene1', 'gene2', 'protein1', 'protein2', 'identity', 'qcov', 'chr1', 'chr2', 'chr1_name', 'chr2_name',
       'distance_bp', 'type', 'aligned_codons', 'Ka', 'Ks', 'Ka_Ks', 'selection', 'note']
with open(os.path.join(out_dir, 'duplication_kaks.tsv'), 'w') as fh:
    fh.write('\t'.join(hdr) + '\n')
    for r in rows: fh.write('\t'.join(r) + '\n')
# compatibility copies in the first-version layouts
with open(os.path.join(supp_dir, 'duplication_analysis.tsv'), 'w') as fh:
    fh.write('gene1\tgene2\tprotein1\tprotein2\tidentity\tchr1\tchr2\tdistance\ttype\tchr1_name\tchr2_name\n')
    for r in rows: fh.write('\t'.join([r[0], r[1], r[2], r[3], r[4], r[6], r[7], r[10] if r[10] != 'NA' else '', r[11], r[8], r[9]]) + '\n')
with open(os.path.join(supp_dir, 'kaks_analysis.tsv'), 'w') as fh:
    fh.write('gene1\tgene2\tprotein1\tprotein2\tidentity\tKa\tKs\tKa_Ks\tselection\n')
    for r in rows: fh.write('\t'.join([r[0], r[1], r[2], r[3], r[4], r[13], r[14], r[15], r[16]]) + '\n')

types = collections.Counter(r[11] for r in rows)
valid = [r for r in rows if r[15] != 'NA']
print('pairs:', len(rows), '| ' + ', '.join('%s %d (%.1f%%)' % (t, n, 100.0 * n / len(rows)) for t, n in sorted(types.items())))
print('Ka/Ks estimable:', len(valid), '| Ks saturated or NA:', len(rows) - len(valid))
if valid:
    v = sorted(float(r[15]) for r in valid)
    print('mean Ka/Ks %.3f | median %.3f | <1: %d (%.1f%%) | >1: %d' % (sum(v) / len(v), v[len(v) // 2], sum(1 for x in v if x < 1), 100.0 * sum(1 for x in v if x < 1) / len(v), sum(1 for x in v if x > 1)))
for t in ('tandem', 'proximal', 'segmental'):
    ks = sorted(float(r[14]) for r in rows if r[11] == t and r[14] != 'NA')
    if ks: print('  %-9s Ks median %.3f (n=%d, range %.3f-%.3f)' % (t, ks[len(ks) // 2], len(ks), ks[0], ks[-1]))
genes_in_pairs = {g for r in rows for g in (r[0], r[1])}
print('genes involved in at least one pair:', len(genes_in_pairs))
PYEOF

echo "Outputs: ${OUT_DIR}/duplication_kaks.tsv, ${SUPP}/duplication_analysis.tsv, ${SUPP}/kaks_analysis.tsv"
echo "Finished: $(date)"
