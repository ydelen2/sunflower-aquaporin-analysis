#!/bin/bash
#SBATCH --job-name=aqp_tables
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/07_build_tables_v2_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/07_build_tables_v2_%j.err
#SBATCH --mail-type=END,FAIL

# ============================================================================
# 07_build_gene_tables_v2.sh
# Builds the gene-level tables that the expression, cis-element and synteny
# steps read.
#
# Inputs (all produced by earlier steps):
#   02_gene_family/verified_aquaporin_ids.txt                 (03 v2)
#   02_gene_family/candidate_aquaporins.faa                   (02; product names)
#   02_gene_family/05_gene_structure/protein_to_gene_map.tsv  (05)
#   02_gene_family/05_gene_structure/aqp_gene_structure_summary.tsv (05)
#   02_gene_family/05_gene_structure/aqp_genes.gff3           (05; mRNA, exon and CDS records)
#   02_gene_family/05_gene_structure/promoter_sequences_2kb.fasta   (05)
#   02_gene_family/08_classification_v2/subfamily_by_tree.tsv (08)
#
# Outputs:
#   02_gene_family/results/verified_aquaporins.txt   gene table (gene_id, product, subfamily, protein_ids)
#   02_gene_family/results/aquaporin_loc_mapping.tsv loc_id, product, subfamily
#   02_gene_family/results/promoters_2kb.fa          copy of the 05 promoter FASTA
#   09_figures/supplementary/Table_S1_aquaporin_genes.tsv
#
# Gene-level values: one representative mRNA per gene, the one with the
# longest coding sequence (the longest transcript where the coding sequences
# of several mRNAs are equal in length); exon and intron counts are those of
# the representative mRNA; mrna_variants is the number of verified isoforms
# of the gene.
# ============================================================================

set -eo pipefail

PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
source "${PROJ_DIR}/config.sh"

GF="${PROJ_DIR}/02_gene_family"
RES="${GF}/results"
SUPP="${PROJ_DIR}/09_figures/supplementary"
mkdir -p "${RES}" "${SUPP}" "${PROJ_DIR}/logs"

IDS="${GF}/verified_aquaporin_ids.txt"
FAA="${GF}/candidate_aquaporins.faa"
P2G="${GF}/05_gene_structure/protein_to_gene_map.tsv"
STRUCT="${GF}/05_gene_structure/aqp_gene_structure_summary.tsv"
GENEGFF="${GF}/05_gene_structure/aqp_genes.gff3"
PROM="${GF}/05_gene_structure/promoter_sequences_2kb.fasta"
CLASS="${GF}/08_classification_v2/subfamily_by_tree.tsv"

for f in "${IDS}" "${FAA}" "${P2G}" "${STRUCT}" "${GENEGFF}" "${PROM}" "${CLASS}"; do
    [[ -s "${f}" ]] || { echo "ERROR: missing input ${f}" >&2; exit 1; }
done

# keep first-version copies once
for f in "${RES}/verified_aquaporins.txt" "${RES}/aquaporin_loc_mapping.tsv" "${RES}/promoters_2kb.fa" "${SUPP}/Table_S1_aquaporin_genes.tsv"; do
    if [[ -e "${f}" || -L "${f}" ]] && [[ ! -e "${f}.v1_backup" ]]; then
        cp -L "${f}" "${f}.v1_backup" 2>/dev/null || true
    fi
done
rm -f "${RES}/verified_aquaporins.txt"

python3 - "${IDS}" "${FAA}" "${P2G}" "${STRUCT}" "${CLASS}" "${RES}" "${SUPP}" "${GENEGFF}" << 'PYEOF'
import sys, re, collections
ids_f, faa_f, p2g_f, struct_f, class_f, res_dir, supp_dir, gff_f = sys.argv[1:9]

ids = [l.strip() for l in open(ids_f) if l.strip()]
idset = set(ids)

# product names from the candidate FASTA headers
product = {}
for line in open(faa_f):
    if line.startswith('>'):
        parts = line[1:].rstrip('\n').split(None, 1)
        desc = parts[1] if len(parts) > 1 else ''
        desc = re.sub(r'\s*\[Helianthus annuus\]\s*$', '', desc)
        product[parts[0]] = desc

# protein -> gene / coordinates
p2g = {}
hdr = None
for i, line in enumerate(open(p2g_f)):
    f = line.rstrip('\n').split('\t')
    if i == 0:
        hdr = f; continue
    row = dict(zip(hdr, f))
    row['gene_id'] = re.sub(r'^gene-', '', row['gene_id'])
    p2g[row['protein_id']] = row

# gene length per gene (05 summary)
gene_length = {}
for i, line in enumerate(open(struct_f)):
    f = line.rstrip('\n').split('\t')
    if i == 0:
        shdr = f; continue
    row = dict(zip(shdr, f))
    gene_length[re.sub(r'^gene-', '', row['gene_id'])] = row['gene_length']

# structure of every mRNA: exons and coding sequence, mRNA -> protein through the CDS records
exons, cds, mrna_prot = collections.defaultdict(list), collections.defaultdict(list), {}
for line in open(gff_f):
    if line.startswith('#') or not line.strip():
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('exon', 'CDS'):
        continue
    a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
    if f[2] == 'exon':
        exons[a['Parent']].append((int(f[3]), int(f[4])))
    else:
        cds[a['Parent']].append((int(f[3]), int(f[4])))
        mrna_prot[a['Parent']] = re.sub(r'^cds-', '', a['ID'])
struct = {}
for m, p in mrna_prot.items():
    ex = sorted(exons[m])
    introns = [ex[i + 1][0] - ex[i][1] - 1 for i in range(len(ex) - 1)]
    st = {'cds_length': sum(e - s + 1 for s, e in cds[m]), 'transcript_length': sum(e - s + 1 for s, e in ex),
          'num_exons': str(len(ex)), 'num_introns': str(len(introns)),
          'avg_intron_length': '%.1f' % (sum(introns) / float(len(introns)) if introns else 0)}
    if p not in struct or (st['cds_length'], st['transcript_length']) > (struct[p]['cds_length'], struct[p]['transcript_length']):
        struct[p] = st

# subfamily by tree per protein
sub = {}
for i, line in enumerate(open(class_f)):
    f = line.rstrip('\n').split('\t')
    if i == 0:
        chdr = f; continue
    row = dict(zip(chdr, f))
    sub[row['protein_id']] = row

chr_map = {'NC_0354%d.2' % (32 + i): 'Ha%d' % i for i in range(1, 18)}

missing = [p for p in ids if p not in p2g]
if missing:
    sys.stderr.write('WARNING: %d verified proteins have no gene mapping (05_gene_structure): %s\n' % (len(missing), ', '.join(missing[:10])))

genes = collections.OrderedDict()
for p in ids:
    if p not in p2g:
        continue
    g = p2g[p]['gene_id']
    genes.setdefault(g, []).append(p)

def length_of(p):
    # longest coding sequence, then longest transcript
    s = struct.get(p)
    return (s['cds_length'], s['transcript_length']) if s else (0, 0)

rows = []
sub_counts = collections.Counter()
for g, prots in genes.items():
    rep = max(prots, key=length_of)
    info = p2g[rep]
    st = struct.get(rep, {})
    calls = [sub[p]['subfamily_tree'] for p in prots if p in sub]
    subfam = collections.Counter(calls).most_common(1)[0][0] if calls else 'unassigned'
    name_calls = [sub[p]['subfamily_from_ncbi_name'] for p in prots if p in sub]
    name_sub = collections.Counter(name_calls).most_common(1)[0][0] if name_calls else 'unassigned'
    conf = min(float(sub[p]['support_fraction']) for p in prots if p in sub) if calls else 0.0
    chrom = info.get('chromosome', '')
    rows.append({
        'loc_id': g,
        'product': re.sub(r' isoform X\d+$', '', product.get(rep, '')),
        'subfamily': subfam,
        'subfamily_from_name': name_sub,
        'tree_support': '%.2f' % conf,
        'chromosome': chrom,
        'chr_name': chr_map.get(chrom, chrom),
        'start': info.get('gene_start', ''),
        'end': info.get('gene_end', ''),
        'strand': info.get('strand', ''),
        'gene_length_bp': gene_length.get(g, ''),
        'exon_count': st.get('num_exons', ''),
        'intron_count': st.get('num_introns', ''),
        'avg_intron_length': st.get('avg_intron_length', ''),
        'mrna_variants': str(len(prots)),
        'representative_protein': rep,
        'protein_ids': ';'.join(prots),
    })
    sub_counts[subfam] += 1

def sort_key(r):
    m = re.match(r'Ha(\d+)', r['chr_name'])
    return (int(m.group(1)) if m else 99, int(r['start'] or 0))
rows.sort(key=sort_key)

cols = ['loc_id', 'product', 'subfamily', 'subfamily_from_name', 'tree_support', 'chromosome', 'chr_name',
        'start', 'end', 'strand', 'gene_length_bp', 'exon_count', 'intron_count', 'avg_intron_length',
        'mrna_variants', 'representative_protein', 'protein_ids']
with open(supp_dir + '/Table_S1_aquaporin_genes.tsv', 'w') as out:
    out.write('\t'.join(cols) + '\n')
    for r in rows:
        out.write('\t'.join(r[c] for c in cols) + '\n')

with open(res_dir + '/aquaporin_loc_mapping.tsv', 'w') as out:
    out.write('loc_id\tproduct\tsubfamily\n')
    for r in rows:
        out.write('%s\t%s\t%s\n' % (r['loc_id'], r['product'], r['subfamily']))

with open(res_dir + '/verified_aquaporins.txt', 'w') as out:
    out.write('gene_id\tgene_name\tsubfamily\tchromosome\tprotein_ids\n')
    for r in rows:
        out.write('%s\t%s\t%s\t%s\t%s\n' % (r['loc_id'], r['product'], r['subfamily'], r['chr_name'], r['protein_ids']))

print('verified proteins: %d | genes: %d' % (len(ids), len(rows)))
print('subfamilies (genes): ' + ', '.join('%s %d' % kv for kv in sorted(sub_counts.items())))
chrom_counts = collections.Counter(r['chr_name'] for r in rows)
print('chromosomes: ' + ', '.join('%s %d' % (c, chrom_counts[c]) for c in sorted(chrom_counts, key=lambda x: int(re.sub(r'\D', '', x) or 99))))
dis = [r for r in rows if r['subfamily'] != r['subfamily_from_name'] and r['subfamily_from_name'] != 'unassigned']
print('genes where tree call differs from the NCBI product name: %d' % len(dis))
for r in dis:
    print('   %s\t%s (tree)\t%s (name)\t%s' % (r['loc_id'], r['subfamily'], r['subfamily_from_name'], r['product']))
PYEOF

[[ "${PROM}" -ef "${RES}/promoters_2kb.fa" ]] || cp "${PROM}" "${RES}/promoters_2kb.fa"
echo "promoter FASTA copied: $(grep -c '^>' "${RES}/promoters_2kb.fa") sequences"
echo ""
echo "Outputs:"
echo "  ${RES}/verified_aquaporins.txt"
echo "  ${RES}/aquaporin_loc_mapping.tsv"
echo "  ${RES}/promoters_2kb.fa"
echo "  ${SUPP}/Table_S1_aquaporin_genes.tsv"
echo "Finished: $(date)"
