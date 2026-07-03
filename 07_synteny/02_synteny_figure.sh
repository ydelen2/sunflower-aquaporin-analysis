#!/usr/bin/env bash
#SBATCH --job-name=synteny_fig
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/synteny_fig_%j.log
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/synteny_fig_%j.err

# =============================================================================
# 02_synteny_figure.sh
# Generate synteny visualization (Circos-style + dual synteny plot)
# Run AFTER 01_synteny_analysis.sh completes
# =============================================================================

set -euo pipefail

source /work/dweikat/ydelen2/aquaporin_study/config.sh

module load miniforge/24.5
conda activate "${CONDA_ENV_PREFIX}"

SYNTENY_DIR="${PROJ_DIR}/07_synteny"
PLOT_DIR="${SYNTENY_DIR}/plots"
FIG_DIR="${PROJ_DIR}/09_figures/plots"
mkdir -p "${PLOT_DIR}" "${FIG_DIR}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "Generating synteny figures..."

python3 << 'FIGPY'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np
import os, csv

synteny_dir = '/work/dweikat/ydelen2/aquaporin_study/07_synteny'
results_dir = os.path.join(synteny_dir, 'results')
plot_dir = os.path.join(synteny_dir, 'plots')
fig_dir = '/work/dweikat/ydelen2/aquaporin_study/09_figures/plots'

# Load AQP gene list with chromosomes
aqp_genes = {}
s1_path = '/work/dweikat/ydelen2/aquaporin_study/09_figures/supplementary/Table_S1_aquaporin_genes.tsv'
if os.path.exists(s1_path):
    with open(s1_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('gene_id', row.get(reader.fieldnames[0], ''))
            chrom = row.get('chromosome', row.get('chr', ''))
            subfamily = row.get('subfamily', row.get('subfamily_predicted', ''))
            aqp_genes[gene] = {'chr': chrom, 'subfamily': subfamily}

# Color scheme for subfamilies
sf_colors = {
    'PIP': '#E63946',
    'TIP': '#457B9D', 
    'NIP': '#2A9D8F',
    'SIP': '#E9C46A',
    'unknown': '#999999'
}

# ─── Figure 9A: Intra-genomic synteny (Circos-like) ───
def parse_synteny_pairs(filepath):
    pairs = []
    if not os.path.exists(filepath):
        return pairs
    with open(filepath) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pairs.append((parts[1], parts[2]))
    return pairs

# ─── Figure 9B: Interspecies synteny summary ───
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

comparisons = [
    ('Ha_self', 'Sunflower\nintra-genomic'),
    ('Ha_At', 'Sunflower vs\nArabidopsis'), 
    ('Ha_Ls', 'Sunflower vs\nLettuce')
]

for idx, (comp, title) in enumerate(comparisons):
    ax = axes[idx]
    pairs = parse_synteny_pairs(os.path.join(results_dir, f'{comp}_aqp_synteny.tsv'))
    
    if not pairs:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14)
        ax.set_title(title, fontsize=13, fontweight='bold')
        continue
    
    # Count pairs by subfamily
    sf_counts = {}
    for g1, g2 in pairs:
        sf = aqp_genes.get(g1, {}).get('subfamily', 
             aqp_genes.get(g2, {}).get('subfamily', 'unknown'))
        sf_counts[sf] = sf_counts.get(sf, 0) + 1
    
    # Bar chart
    sfs = sorted(sf_counts.keys())
    counts = [sf_counts[s] for s in sfs]
    colors = [sf_colors.get(s, '#999999') for s in sfs]
    
    bars = ax.bar(sfs, counts, color=colors, edgecolor='white', linewidth=1.5)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_ylabel('Syntenic pairs', fontsize=11)
    ax.set_xlabel('Subfamily', fontsize=11)
    
    # Add count labels
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                str(count), ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Total annotation
    ax.text(0.95, 0.95, f'Total: {sum(counts)}', transform=ax.transAxes,
            ha='right', va='top', fontsize=11, 
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgray', alpha=0.7))

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'Fig9_synteny_summary.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(plot_dir, 'Fig9_synteny_summary.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(fig_dir, 'Fig9_synteny_summary.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(fig_dir, 'Fig9_synteny_summary.pdf'), bbox_inches='tight')
plt.close()

print("Figure 9 saved")

# ─── Detailed synteny table for supplementary ───
sup_path = '/work/dweikat/ydelen2/aquaporin_study/09_figures/supplementary'
with open(os.path.join(sup_path, 'Table_S3b_synteny.tsv'), 'w') as f:
    f.write('comparison\tgene1\tgene2\tsubfamily1\tsubfamily2\tblock\n')
    for comp, _ in comparisons:
        pairs_file = os.path.join(results_dir, f'{comp}_aqp_synteny.tsv')
        if not os.path.exists(pairs_file):
            continue
        with open(pairs_file) as rf:
            next(rf)  # header
            for line in rf:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    block = parts[0] if len(parts) > 0 else ''
                    g1 = parts[1]
                    g2 = parts[2]
                    sf1 = aqp_genes.get(g1, {}).get('subfamily', 'N/A')
                    sf2 = aqp_genes.get(g2, {}).get('subfamily', 'N/A')
                    f.write(f'{comp}\t{g1}\t{g2}\t{sf1}\t{sf2}\t{block}\n')

print("Table S3b (synteny) saved")
print("\nDone!")

FIGPY

# Copy results to supplementary
cp "${SYNTENY_DIR}/results/"*.tsv "${FIG_DIR}/../supplementary/" 2>/dev/null || true

log "Synteny figures complete!"
log "Output:"
log "  ${PLOT_DIR}/Fig9_synteny_summary.png"
log "  ${FIG_DIR}/Fig9_synteny_summary.png"
