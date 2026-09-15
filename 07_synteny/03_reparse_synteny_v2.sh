#!/usr/bin/env bash
#SBATCH --job-name=syn_parse_v2
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/syn_parse_v2_%j.log
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/syn_parse_v2_%j.err

# ============================================================================
# 03_reparse_synteny_v2.sh
# Re-parse MCScanX collinearity files for aquaporin synteny.
#
# Difference from the first version: an intra-genomic pair is counted as an
# aquaporin syntenic pair only when BOTH genes are in the aquaporin inventory.
# The first version kept a pair when either gene was an aquaporin, which is
# how partners outside the inventory ended up in Table S7b. For interspecies
# pairs the sunflower member must be in the inventory; the partner's status
# is recorded in the new partner_in_reference_set column when a reference
# aquaporin list is available (see REF_AQP_LISTS below), otherwise NA.
#
# Output file names and columns match the first version so that
# 04_synteny_figure.sh keeps working; the previous results_v2 folder is kept
# as results_v2_before_inventory_fix.
# ============================================================================

set -euo pipefail

PROJ_DIR="/work/dweikat/ydelen2/aquaporin_study"
SYNTENY_DIR="${PROJ_DIR}/07_synteny"
AQP_LIST="${PROJ_DIR}/09_figures/supplementary/Table_S1_aquaporin_genes.tsv"

# Optional: one accession per line, reference-species aquaporin proteins
# (for example the sequence_ids written by 05_phylogenetics/results_v2).
# Leave empty strings if not available.
REF_AQP_AT="${PROJ_DIR}/05_phylogenetics/results_v2/At_hits.txt"
REF_AQP_LS="${PROJ_DIR}/05_phylogenetics/results_v2/Ls_hits.txt"

module load miniforge/24.5
conda activate "${PROJ_DIR}/conda_envs/aquaporin_env"

if [[ -d "${SYNTENY_DIR}/results_v2" && ! -d "${SYNTENY_DIR}/results_v2_before_inventory_fix" ]]; then
    cp -r "${SYNTENY_DIR}/results_v2" "${SYNTENY_DIR}/results_v2_before_inventory_fix"
    echo "kept previous results as results_v2_before_inventory_fix"
fi

python3 - "${SYNTENY_DIR}" "${AQP_LIST}" "${REF_AQP_AT}" "${REF_AQP_LS}" <<'PYFIX'
import os, sys

synteny_dir, aqp_file, ref_at, ref_ls = sys.argv[1:5]
data_dir = os.path.join(synteny_dir, 'data')
results_dir = os.path.join(synteny_dir, 'results_v2')
os.makedirs(results_dir, exist_ok=True)

all_prot_to_gene = {}
for sp in ['Ha', 'At', 'Ls']:
    map_file = os.path.join(data_dir, sp + '_prot_to_gene.tsv')
    if os.path.exists(map_file):
        n = 0
        with open(map_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    all_prot_to_gene[parts[0]] = (sp, parts[1]); n += 1
        print('  %s mappings: %d' % (sp, n))

aqp_locs = set()
with open(aqp_file) as f:
    next(f)
    for line in f:
        if line.strip():
            aqp_locs.add(line.strip().split('\t')[0])
aqp_prots = {p for p, (sp, loc) in all_prot_to_gene.items() if sp == 'Ha' and loc in aqp_locs}
print('AQP genes: %d, AQP proteins: %d' % (len(aqp_locs), len(aqp_prots)))

ref_sets = {}
for sp, path in (('At', ref_at), ('Ls', ref_ls)):
    if path and os.path.exists(path):
        ref_sets[sp] = {l.strip() for l in open(path) if l.strip()}
        print('  reference aquaporin list for %s: %d accessions' % (sp, len(ref_sets[sp])))

chr_map = {'NC_0354%d.2' % (32 + i): 'Ha%d' % i for i in range(1, 18)}

def parse_collinearity(path):
    blocks = []
    header = None; chr1 = chr2 = None; inter = False; genes = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith('## Alignment'):
                if header is not None:
                    blocks.append((header, chr1, chr2, inter, genes))
                header = s; genes = []
                parts = s.split()
                loc = [p for p in parts if '&' in p]
                if loc:
                    chr1, chr2 = loc[0].split('&')
                inter = (chr1.split('_')[0] != chr2.split('_')[0]) if chr1 and chr2 else False
            elif s and not s.startswith('#'):
                parts = s.split()
                prots = [p for p in parts if p.startswith(('XP_', 'NP_'))]
                if len(prots) >= 2:
                    genes.append((prots[0], prots[1]))
    if header is not None:
        blocks.append((header, chr1, chr2, inter, genes))
    return blocks

for comp in ['Ha_self_v2', 'Ha_At_v2', 'Ha_Ls_v2']:
    col_file = os.path.join(synteny_dir, 'mcscanx', comp, comp + '.collinearity')
    print('\n=== %s ===' % comp)
    if not os.path.exists(col_file):
        print('  collinearity file not found: ' + col_file); continue
    blocks = parse_collinearity(col_file)
    print('  total blocks: %d' % len(blocks))
    other_sp = None if comp == 'Ha_self_v2' else comp.split('_')[1]

    intra_pairs, inter_pairs = [], []
    dropped_one_sided = 0
    for header, chr1, chr2, inter, gene_pairs in blocks:
        for g1, g2 in gene_pairs:
            sp1, gene1 = all_prot_to_gene.get(g1, ('?', g1))
            sp2, gene2 = all_prot_to_gene.get(g2, ('?', g2))
            if other_sp is None:
                # intra-genomic: both members must be inventory aquaporins
                if not (g1 in aqp_prots and g2 in aqp_prots):
                    if g1 in aqp_prots or g2 in aqp_prots:
                        dropped_one_sided += 1
                    continue
                partner_flag = 'NA'
            else:
                ha_first = (sp1 == 'Ha')
                ha_prot, other_prot = (g1, g2) if ha_first else (g2, g1)
                if ha_prot not in aqp_prots:
                    continue
                ref = ref_sets.get(other_sp)
                partner_flag = 'NA' if ref is None else ('yes' if other_prot in ref else 'no')
            pair = {'chr1': chr1, 'chr2': chr2, 'sp1': sp1, 'gene1': gene1, 'prot1': g1,
                    'sp2': sp2, 'gene2': gene2, 'prot2': g2, 'partner': partner_flag}
            (inter_pairs if inter else intra_pairs).append(pair)   # a block within one genome is intra even in a two-species run

    if other_sp is None:
        print('  intra-genomic AQP-AQP pairs: %d (one-sided pairs excluded: %d)' % (len(intra_pairs), dropped_one_sided))
        genes = {p['gene1'] for p in intra_pairs} | {p['gene2'] for p in intra_pairs}
        print('  unique Ha AQP genes in AQP-AQP pairs: %d of %d' % (len(genes), len(aqp_locs)))
        chr_pairs = {}
        for p in intra_pairs:
            c1 = chr_map.get(p['chr1'].replace('Ha_', ''), p['chr1'])
            c2 = chr_map.get(p['chr2'].replace('Ha_', ''), p['chr2'])
            chr_pairs.setdefault('%s-%s' % (c1, c2), []).append((p['gene1'], p['gene2']))
        for cp in sorted(chr_pairs):
            print('    %s: %d pairs' % (cp, len(chr_pairs[cp])))
    else:
        ha_genes = {p['gene1'] if p['sp1'] == 'Ha' else p['gene2'] for p in inter_pairs}
        other_genes = {p['gene2'] if p['sp1'] == 'Ha' else p['gene1'] for p in inter_pairs}
        print('  interspecies pairs with a sunflower AQP: %d' % len(inter_pairs))
        print('  unique Ha AQP genes: %d of %d | unique %s partners: %d' % (len(ha_genes), len(aqp_locs), other_sp, len(other_genes)))
        if other_sp in ref_sets:
            n_yes = sum(1 for p in inter_pairs if p['partner'] == 'yes')
            print('  partners that are reference aquaporins: %d of %d' % (n_yes, len(inter_pairs)))

    tag = comp.replace('_v2', '')
    for pairs, suffix in ((intra_pairs, 'intra'), (inter_pairs, 'inter')):
        out_file = os.path.join(results_dir, '%s_%s_aqp_synteny.tsv' % (tag, suffix))
        with open(out_file, 'w') as f:
            f.write('chr1\tchr2\tspecies1\tgene1\tprot1\tspecies2\tgene2\tprot2\tpartner_in_reference_set\n')
            for p in pairs:
                f.write('\t'.join([p['chr1'], p['chr2'], p['sp1'], p['gene1'], p['prot1'],
                                   p['sp2'], p['gene2'], p['prot2'], p['partner']]) + '\n')
        print('  wrote %s (%d rows)' % (os.path.basename(out_file), len(pairs)))

print('\nDone. Results in results_v2/')
PYFIX
