#!/usr/bin/env bash
#SBATCH --job-name=syn_parse
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/syn_parse_%j.log
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/syn_parse_%j.err

set -euo pipefail

PROJ_DIR="/work/dweikat/ydelen2/aquaporin_study"
SYNTENY_DIR="${PROJ_DIR}/07_synteny"
AQP_LIST="${PROJ_DIR}/09_figures/supplementary/Table_S1_aquaporin_genes.tsv"

module load miniforge/24.5
conda activate "${PROJ_DIR}/conda_envs/aquaporin_env"

python3 - "${SYNTENY_DIR}" "${AQP_LIST}" <<'PYFIX'
import os, sys, re

synteny_dir = sys.argv[1]
aqp_file = sys.argv[2]
data_dir = os.path.join(synteny_dir, 'data')
results_dir = os.path.join(synteny_dir, 'results_v2')
os.makedirs(results_dir, exist_ok=True)

# Load protein -> gene mappings
all_prot_to_gene = {}
for sp in ['Ha', 'At', 'Ls']:
    map_file = os.path.join(data_dir, f'{sp}_prot_to_gene.tsv')
    if os.path.exists(map_file):
        count = 0
        with open(map_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    all_prot_to_gene[parts[0]] = (sp, parts[1])
                    count += 1
        print(f"  {sp} mappings: {count}")

# Load AQP genes
aqp_locs = set()
with open(aqp_file) as f:
    next(f)
    for line in f:
        aqp_locs.add(line.strip().split('\t')[0])

# AQP protein IDs
aqp_prots = set()
for prot, (sp, loc) in all_prot_to_gene.items():
    if sp == 'Ha' and loc in aqp_locs:
        aqp_prots.add(prot)
print(f"\nAQP genes: {len(aqp_locs)}, AQP proteins: {len(aqp_prots)}")

# NC accession -> chromosome name
chr_map = {}
for i in range(1, 18):
    nc = f'NC_0354{32+i}.2'
    chr_map[nc] = f'Ha{i}'

for comp in ['Ha_self_v2', 'Ha_At_v2', 'Ha_Ls_v2']:
    col_file = os.path.join(synteny_dir, 'mcscanx', comp, f'{comp}.collinearity')
    if not os.path.exists(col_file):
        print(f"\n{comp}: file not found")
        continue

    print(f"\n{'='*60}")
    print(f"  {comp}")
    print(f"{'='*60}")

    # Parse collinearity file into blocks
    blocks = []  # list of (header, chr1, chr2, is_inter, gene_pairs)
    current_header = ''
    current_chr1 = ''
    current_chr2 = ''
    current_is_inter = False
    current_genes = []

    with open(col_file) as f:
        for line in f:
            line_s = line.strip()
            if line_s.startswith('## Alignment'):
                if current_header:
                    blocks.append((current_header, current_chr1, current_chr2, current_is_inter, current_genes))
                current_header = line_s
                current_genes = []
                chr_match = re.findall(r'(\w+_NC_\d+\.\d+)', line_s)
                if len(chr_match) >= 2:
                    current_chr1 = chr_match[0]
                    current_chr2 = chr_match[1]
                    sp1_prefix = current_chr1.split('_')[0]
                    sp2_prefix = current_chr2.split('_')[0]
                    current_is_inter = (sp1_prefix != sp2_prefix)
                else:
                    current_chr1 = ''
                    current_chr2 = ''
                    current_is_inter = False
            elif line_s and not line_s.startswith('#'):
                parts = line_s.split()
                prots = [p for p in parts if p.startswith(('XP_', 'NP_'))]
                if len(prots) >= 2:
                    current_genes.append((prots[0], prots[1]))
    if current_header:
        blocks.append((current_header, current_chr1, current_chr2, current_is_inter, current_genes))

    print(f"  Total blocks: {len(blocks)}")

    intra_pairs = []
    inter_pairs = []

    for header, chr1, chr2, is_inter, gene_pairs in blocks:
        for g1, g2 in gene_pairs:
            if g1 not in aqp_prots and g2 not in aqp_prots:
                continue
            sp1_info = all_prot_to_gene.get(g1, ('?', g1))
            sp2_info = all_prot_to_gene.get(g2, ('?', g2))
            pair = {
                'chr1': chr1, 'chr2': chr2,
                'sp1': sp1_info[0], 'gene1': sp1_info[1], 'prot1': g1,
                'sp2': sp2_info[0], 'gene2': sp2_info[1], 'prot2': g2,
            }
            if is_inter:
                inter_pairs.append(pair)
            else:
                intra_pairs.append(pair)

    print(f"  Intra-species AQP pairs: {len(intra_pairs)}")
    print(f"  Inter-species AQP pairs: {len(inter_pairs)}")

    # Unique genes
    intra_genes = set()
    for p in intra_pairs:
        if p['gene1'] in aqp_locs: intra_genes.add(p['gene1'])
        if p['gene2'] in aqp_locs: intra_genes.add(p['gene2'])

    inter_genes_ha = set()
    inter_genes_other = set()
    for p in inter_pairs:
        if p['sp1'] == 'Ha' and p['gene1'] in aqp_locs: inter_genes_ha.add(p['gene1'])
        if p['sp2'] == 'Ha' and p['gene2'] in aqp_locs: inter_genes_ha.add(p['gene2'])
        if p['sp1'] != 'Ha' and p['sp1'] != '?': inter_genes_other.add(f"{p['sp1']}:{p['gene1']}")
        if p['sp2'] != 'Ha' and p['sp2'] != '?': inter_genes_other.add(f"{p['sp2']}:{p['gene2']}")

    print(f"  Unique Ha AQP genes (intra): {len(intra_genes)}")
    print(f"  Unique Ha AQP genes (inter): {len(inter_genes_ha)}")
    print(f"  Unique other-species orthologs: {len(inter_genes_other)}")

    # Chromosome pair summary for intra
    if intra_pairs:
        print(f"\n  INTRA-SPECIES (Ha-Ha) pairs:")
        chr_pairs = {}
        for p in intra_pairs:
            nc1 = p['chr1'].replace('Ha_', '')
            nc2 = p['chr2'].replace('Ha_', '')
            c1 = chr_map.get(nc1, nc1)
            c2 = chr_map.get(nc2, nc2)
            key = f"{c1}-{c2}"
            if key not in chr_pairs:
                chr_pairs[key] = []
            g = p['gene1'] if p['gene1'] in aqp_locs else p['gene2']
            chr_pairs[key].append(g)
        for cp in sorted(chr_pairs):
            genes = sorted(set(chr_pairs[cp]))
            print(f"    {cp}: {len(chr_pairs[cp])} pairs ({', '.join(genes)})")

    # Inter-species pairs detail
    if inter_pairs:
        print(f"\n  INTER-SPECIES pairs:")
        for p in inter_pairs:
            nc1 = p['chr1'].replace('Ha_', '').replace('At_', '').replace('Ls_', '')
            nc2 = p['chr2'].replace('Ha_', '').replace('At_', '').replace('Ls_', '')
            c1 = chr_map.get(nc1, nc1)
            c2 = chr_map.get(nc2, nc2)
            ha_gene = '?'
            other_gene = '?'
            other_sp = '?'
            if p['sp1'] == 'Ha' and p['gene1'] in aqp_locs:
                ha_gene = p['gene1']
                other_sp = p['sp2']
                other_gene = p['gene2']
            elif p['sp2'] == 'Ha' and p['gene2'] in aqp_locs:
                ha_gene = p['gene2']
                other_sp = p['sp1']
                other_gene = p['gene1']
            elif p['sp1'] == 'Ha':
                ha_gene = p['gene1']
                other_sp = p['sp2']
                other_gene = p['gene2']
            else:
                ha_gene = p['gene2']
                other_sp = p['sp1']
                other_gene = p['gene1']
            print(f"    {c1}-{c2}: Ha:{ha_gene} <-> {other_sp}:{other_gene}")

    # Save clean TSVs
    tag = comp.replace('_v2', '')
    for pairs_list, suffix in [(intra_pairs, 'intra'), (inter_pairs, 'inter')]:
        out_file = os.path.join(results_dir, f'{tag}_{suffix}_aqp_synteny.tsv')
        with open(out_file, 'w') as f:
            f.write('chr1\tchr2\tspecies1\tgene1\tprot1\tspecies2\tgene2\tprot2\n')
            for p in pairs_list:
                f.write(f"{p['chr1']}\t{p['chr2']}\t{p['sp1']}\t{p['gene1']}\t{p['prot1']}\t{p['sp2']}\t{p['gene2']}\t{p['prot2']}\n")

print("\n\nDone. Results in results_v2/")
PYFIX
