#!/usr/bin/env bash
#SBATCH --job-name=syn_fig
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/syn_fig_%j.log
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/syn_fig_%j.err

set -euo pipefail

PROJ_DIR="/work/dweikat/ydelen2/aquaporin_study"
RESULTS_DIR="${PROJ_DIR}/07_synteny/results_v2"
PLOTS_DIR="${PROJ_DIR}/09_figures/plots"

module load miniforge/24.5
conda activate "${PROJ_DIR}/conda_envs/aquaporin_env"

python3 - "${RESULTS_DIR}" "${PLOTS_DIR}" <<'PYFIG'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv
import sys
import os

results_dir = sys.argv[1]
plots_dir = sys.argv[2]

# Load data
def load_tsv(fname):
    rows = []
    fpath = os.path.join(results_dir, fname)
    if os.path.exists(fpath):
        with open(fpath) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                rows.append(row)
    return rows

intra = load_tsv('Ha_self_intra_aqp_synteny.tsv')
ha_at_inter = load_tsv('Ha_At_inter_aqp_synteny.tsv')
ha_ls_inter = load_tsv('Ha_Ls_inter_aqp_synteny.tsv')

print(f"Loaded: {len(intra)} intra, {len(ha_at_inter)} Ha-At, {len(ha_ls_inter)} Ha-Ls")

# NC -> Chr mapping
chr_map = {}
for i in range(1, 18):
    nc = f'NC_0354{32+i}.2'
    chr_map[f'Ha_{nc}'] = f'Ha{i}'
    chr_map[nc] = f'Ha{i}'

def get_chr(nc_str):
    return chr_map.get(nc_str, nc_str.split('_')[-1][:8])

# Count unique Ha AQP genes per comparison
def unique_ha_genes(pairs):
    genes = set()
    for p in pairs:
        if p['species1'] == 'Ha':
            genes.add(p['gene1'])
        if p['species2'] == 'Ha':
            genes.add(p['gene2'])
    return genes

intra_genes = unique_ha_genes(intra)
at_genes = unique_ha_genes(ha_at_inter)
ls_genes = unique_ha_genes(ha_ls_inter)

# Count pairs per chromosome for intra
chr_pair_counts = {}
for p in intra:
    c1 = get_chr(p['chr1'])
    c2 = get_chr(p['chr2'])
    key = tuple(sorted([c1, c2]))
    chr_pair_counts[key] = chr_pair_counts.get(key, 0) + 1

# Create figure with 2 panels
fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [1, 1.3]})

# Panel A: Bar chart comparing 3 comparisons
ax1 = axes[0]
comparisons = ['Ha self\n(intra-genomic)', 'Ha vs At\n(interspecies)', 'Ha vs Ls\n(interspecies)']
pair_counts = [len(intra), len(ha_at_inter), len(ha_ls_inter)]
gene_counts = [len(intra_genes), len(at_genes), len(ls_genes)]
colors_pairs = ['#2196F3', '#FF9800', '#4CAF50']
colors_genes = ['#1565C0', '#E65100', '#2E7D32']

x = range(len(comparisons))
width = 0.35
bars1 = ax1.bar([i - width/2 for i in x], pair_counts, width, label='Syntenic pairs',
                color=colors_pairs, edgecolor='white', linewidth=0.5)
bars2 = ax1.bar([i + width/2 for i in x], gene_counts, width, label='Unique Ha AQP genes',
                color=colors_genes, edgecolor='white', linewidth=0.5)

ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
ax1.set_title('A. Aquaporin Syntenic Relationships', fontsize=13, fontweight='bold', pad=10)
ax1.set_xticks(list(x))
ax1.set_xticklabels(comparisons, fontsize=10)
ax1.legend(fontsize=9, loc='upper left')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Add value labels
for bar in bars1:
    ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1,
             str(int(bar.get_height())), ha='center', va='bottom', fontsize=10, fontweight='bold')
for bar in bars2:
    ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1,
             str(int(bar.get_height())), ha='center', va='bottom', fontsize=10, fontweight='bold')

# Panel B: Chromosome pair heatmap for intra-genomic
ax2 = axes[1]
chromosomes = [f'Ha{i}' for i in range(1, 18)]
n = len(chromosomes)
matrix = [[0]*n for _ in range(n)]
for (c1, c2), count in chr_pair_counts.items():
    if c1 in chromosomes and c2 in chromosomes:
        i1 = chromosomes.index(c1)
        i2 = chromosomes.index(c2)
        matrix[i1][i2] = count
        matrix[i2][i1] = count

import numpy as np
mat = np.array(matrix, dtype=float)
mat[mat == 0] = np.nan

cmap = plt.cm.YlOrRd.copy()
cmap.set_bad('white')
im = ax2.imshow(mat, cmap=cmap, aspect='equal', interpolation='nearest')
ax2.set_xticks(range(n))
ax2.set_yticks(range(n))
ax2.set_xticklabels(chromosomes, fontsize=8, rotation=45, ha='right')
ax2.set_yticklabels(chromosomes, fontsize=8)
ax2.set_title('B. Intra-genomic AQP Syntenic Pairs\nby Chromosome Pair', fontsize=13, fontweight='bold', pad=10)

# Add text annotations
for i in range(n):
    for j in range(n):
        if matrix[i][j] > 0:
            ax2.text(j, i, str(matrix[i][j]), ha='center', va='center',
                    fontsize=8, fontweight='bold',
                    color='white' if matrix[i][j] >= 4 else 'black')

cbar = plt.colorbar(im, ax=ax2, shrink=0.8, label='Number of AQP pairs')

plt.tight_layout()

# Save
out_pdf = os.path.join(plots_dir, 'Fig9_synteny_summary.pdf')
out_png = os.path.join(plots_dir, 'Fig9_synteny_summary.png')
plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
plt.savefig(out_png, dpi=300, bbox_inches='tight')
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")

# Print summary
print(f"\nSummary:")
print(f"  Intra-genomic: {len(intra)} pairs, {len(intra_genes)} unique genes")
print(f"  Ha-At inter: {len(ha_at_inter)} pairs, {len(at_genes)} unique genes")
print(f"  Ha-Ls inter: {len(ha_ls_inter)} pairs, {len(ls_genes)} unique genes")
print(f"  Top chr pairs: {sorted(chr_pair_counts.items(), key=lambda x: -x[1])[:5]}")
PYFIG
