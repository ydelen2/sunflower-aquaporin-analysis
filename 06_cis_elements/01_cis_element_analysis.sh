#!/bin/bash
#SBATCH --job-name=aqp_cis
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/06_cis_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/06_cis_%j.err
#SBATCH --mail-type=END,FAIL

# ============================================================================
# 01_cis_element_analysis.sh
# Cis-regulatory element analysis of sunflower aquaporin promoter regions.
#
# Strategy:
#   1. Local regex-based scan for known plant stress-responsive cis-elements
#      (ABRE, DRE/CRT, MBS, TC-rich, LTR, HSE, W-box, as-1, etc.)
#   2. FIMO scan against JASPAR plant TF motifs (if DB available)
#   3. Generate summary tables and R-based heatmap visualization
# ============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
source "${PROJ_DIR}/config.sh"

WORK_DIR="${PROJ_DIR}/06_cis_elements"
RESULTS_DIR="${WORK_DIR}/results"
SCRIPTS_DIR="${WORK_DIR}/scripts"
LOG_DIR="${PROJ_DIR}/logs"

mkdir -p "${RESULTS_DIR}" "${LOG_DIR}" "${SCRIPTS_DIR}"

PROMOTER_FA="${PROJ_DIR}/02_gene_family/results/promoters_2kb.fa"

if [[ ! -f "${PROMOTER_FA}" ]]; then
    echo "ERROR: Promoter sequences not found: ${PROMOTER_FA}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Modules
# ---------------------------------------------------------------------------
module purge
module load MEME/5.5
module load R/4.1

echo "================================================================"
echo "Pipeline: Cis-Regulatory Element Analysis"
echo "Started: $(date)"
echo "Job ID:  ${SLURM_JOB_ID}"
echo "================================================================"

SEQ_COUNT=$(grep -c '^>' "${PROMOTER_FA}")
echo "Promoter sequences: ${SEQ_COUNT}"

# ===================================================================
# PART 1: Local motif scanning — known plant cis-elements
# ===================================================================
echo "[$(date '+%H:%M:%S')] Running local cis-element scan..."

cat > "${SCRIPTS_DIR}/scan_cis_elements.py" << 'PYEOF'
#!/usr/bin/env python3
"""
scan_cis_elements.py
Scan promoter FASTA sequences for known plant cis-regulatory elements
using regex pattern matching on both strands.

Outputs:
  cis_elements_all.tsv       — every hit: gene, element, position, strand, sequence
  cis_elements_counts.tsv    — gene x element count matrix
  cis_elements_summary.tsv   — total counts per element type
"""

import sys
import re
import csv
import os
from collections import defaultdict, OrderedDict

# ---- Known plant cis-element motifs ----
# IUPAC: R=[AG], Y=[CT], W=[AT], S=[GC], M=[AC], K=[GT], N=[ACGT],
#         B=[CGT], D=[AGT], H=[ACT], V=[ACG]
IUPAC_MAP = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'W': '[AT]', 'S': '[GC]',
    'M': '[AC]', 'K': '[GT]', 'N': '[ACGT]',
    'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
}

def iupac_to_regex(motif):
    """Convert IUPAC motif string to regex pattern."""
    return ''.join(IUPAC_MAP.get(c.upper(), c) for c in motif)

def reverse_complement(seq):
    comp = str.maketrans('ACGTacgtRYWSMKBDHVNrywsmkbdhvn',
                         'TGCAtgcaYRWSKMVHDBNyrwskmvhdbn')
    return seq.translate(comp)[::-1]

# Curated plant cis-element database
# Format: (element_name, IUPAC_motif, description)
CIS_ELEMENTS = [
    # ABA responsive
    ('ABRE',           'ACGTG',       'ABA responsive element (core)'),
    ('ABRE_motifIIb',  'AGTACGTGGC',  'ABA responsive element motif IIb'),

    # Drought / dehydration / cold
    ('DRE_CRT',        'RCCGAC',      'Dehydration-responsive / C-repeat'),
    ('DRE_core',       'ACCGAC',      'DRE core'),
    ('LTRE_CRT',       'CCGAC',       'Low-temperature responsive / CRT core'),

    # MYB binding sites
    ('MBS',            'CAACTG',      'MYB binding site (drought)'),
    ('MBS_I',          'YAACKG',      'MYB binding site type I'),
    ('MYB_core',       'CNGTTR',      'MYB recognition core'),
    ('MYB_AACCA',      'AACCA',       'MYB binding motif AACCA'),

    # MYC binding
    ('MYC',            'CANNTG',      'MYC recognition site (E-box core)'),
    ('G_box',          'CACGTG',      'G-box (MYC/bHLH/light responsive)'),

    # Heat stress
    ('HSE',            'AAAAAATTTC',  'Heat shock element (perfect)'),
    ('HSE_core',       'NGAAN',       'HSE core repeat unit'),
    ('HSE_inverted',   'NTTCN',       'HSE inverted repeat unit'),

    # TC-rich repeats (defense/stress)
    ('TC_rich',        'ATTTTCTTCA',  'TC-rich repeats (defense/stress)'),
    ('TC_rich_var',    'GTTTTCTTAC',  'TC-rich repeat variant'),

    # LTR (low temperature responsive)
    ('LTR',            'CCGAAA',      'Low-temperature responsive element'),

    # WRKY binding (W-box)
    ('W_box',          'TTGACY',      'W-box (WRKY binding site)'),
    ('W_box_cluster',  'TTGACYTTGACY','W-box cluster (tandem)'),

    # as-1 element (oxidative stress, auxin)
    ('as1_element',    'TGACG',       'as-1 / ocs element core'),
    ('as1_full',       'TGACGTAA',    'as-1 element extended'),

    # Salicylic acid responsive
    ('TCA_element',    'CCATCTTTTT',  'TCA-element (SA responsive)'),
    ('SARE',           'TGACG',       'SA responsive element'),

    # Jasmonate responsive
    ('CGTCA_motif',    'CGTCA',       'MeJA-responsive CGTCA-motif'),
    ('TGACG_motif',    'TGACG',       'MeJA-responsive TGACG-motif'),

    # Ethylene responsive
    ('ERE',            'ATTTCAAA',    'Ethylene-responsive element'),
    ('GCC_box',        'AGCCGCC',     'GCC-box (ethylene/pathogen response)'),

    # Anaerobic / hypoxia
    ('ARE',            'AAACCA',      'Anaerobic responsive element'),

    # Gibberellin responsive
    ('GARE',           'TAACAAR',     'GA-responsive element'),
    ('P_box',          'CCTTTTG',     'GA-responsive P-box'),

    # Auxin responsive
    ('AuxRE',          'TGTCTC',      'Auxin-responsive element'),
    ('TGA_element',    'AACGAC',      'Auxin-responsive TGA element'),

    # Light responsive
    ('GT1_motif',      'GRWAAW',      'GT1 light-responsive motif'),
    ('I_box',          'GATAAG',      'I-box (light responsive)'),
    ('Box_4',          'ATTAAT',      'Box 4 (light responsive)'),
    ('ACE',            'ACGTGGA',     'ACE cis-acting element (light)'),
    ('Sp1',            'GGGCGG',      'Sp1 light responsive element'),

    # Circadian
    ('circadian',      'CAANNNNATC',  'Circadian control element'),

    # Wound responsive
    ('WUN_motif',      'AAATTTCCT',   'Wound-responsive element'),

    # Water stress
    ('MYB_drought',    'TAACTG',      'MYB drought-inducible'),

    # General stress
    ('STRE',           'AGGGG',       'Stress responsive element'),
    ('DRE_like',       'GCCGAC',      'DRE-like element'),
]

def scan_sequence(name, seq, elements):
    """Scan a sequence for all cis-elements on both strands."""
    hits = []
    seq_upper = seq.upper()
    rc_seq = reverse_complement(seq_upper)
    seq_len = len(seq_upper)

    for elem_name, iupac_motif, description in elements:
        pattern = iupac_to_regex(iupac_motif)
        regex = re.compile(pattern, re.IGNORECASE)

        # Forward strand
        for m in regex.finditer(seq_upper):
            hits.append({
                'gene': name,
                'element': elem_name,
                'motif_pattern': iupac_motif,
                'description': description,
                'position': m.start() + 1,  # 1-based
                'end': m.end(),
                'strand': '+',
                'sequence': m.group(),
                'distance_to_tss': -(seq_len - m.start()),  # negative = upstream
            })

        # Reverse strand
        for m in regex.finditer(rc_seq):
            # Convert RC position back to forward-strand coordinate
            fwd_pos = seq_len - m.end()
            hits.append({
                'gene': name,
                'element': elem_name,
                'motif_pattern': iupac_motif,
                'description': description,
                'position': fwd_pos + 1,
                'end': fwd_pos + len(m.group()),
                'strand': '-',
                'sequence': m.group(),
                'distance_to_tss': -(seq_len - fwd_pos),
            })

    return hits

def read_fasta(fasta_path):
    """Simple FASTA parser."""
    sequences = OrderedDict()
    name = None
    seqs = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    sequences[name] = ''.join(seqs)
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
    if name is not None:
        sequences[name] = ''.join(seqs)
    return sequences

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <promoters.fa> <output_dir>", file=sys.stderr)
        sys.exit(1)

    fasta_path = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)

    sequences = read_fasta(fasta_path)
    print(f"  Loaded {len(sequences)} promoter sequences")

    # Scan all sequences
    all_hits = []
    for name, seq in sequences.items():
        hits = scan_sequence(name, seq, CIS_ELEMENTS)
        all_hits.extend(hits)

    print(f"  Total cis-element hits: {len(all_hits)}")

    # ---- Output 1: All hits ----
    all_path = os.path.join(output_dir, 'cis_elements_all.tsv')
    fields = ['gene', 'element', 'motif_pattern', 'position', 'end', 'strand',
              'sequence', 'distance_to_tss', 'description']
    with open(all_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(sorted(all_hits, key=lambda x: (x['gene'], x['position'])))
    print(f"  Wrote all hits: {all_path}")

    # ---- Output 2: Gene x Element count matrix ----
    gene_elem_counts = defaultdict(lambda: defaultdict(int))
    elem_names_seen = set()
    for hit in all_hits:
        gene_elem_counts[hit['gene']][hit['element']] += 1
        elem_names_seen.add(hit['element'])

    # Sort element names by total frequency
    elem_totals = defaultdict(int)
    for hit in all_hits:
        elem_totals[hit['element']] += 1
    sorted_elems = sorted(elem_names_seen, key=lambda e: -elem_totals[e])

    count_path = os.path.join(output_dir, 'cis_elements_counts.tsv')
    with open(count_path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['gene'] + sorted_elems)
        for gene in sorted(gene_elem_counts.keys()):
            row = [gene] + [gene_elem_counts[gene].get(e, 0) for e in sorted_elems]
            writer.writerow(row)
    print(f"  Wrote count matrix: {count_path}")

    # ---- Output 3: Summary per element type ----
    summary_path = os.path.join(output_dir, 'cis_elements_summary.tsv')
    # Get description for each element
    elem_desc = {e[0]: e[2] for e in CIS_ELEMENTS}
    with open(summary_path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['element', 'total_hits', 'genes_with_element', 'description'])
        genes_with = defaultdict(set)
        for hit in all_hits:
            genes_with[hit['element']].add(hit['gene'])
        for elem in sorted_elems:
            writer.writerow([elem, elem_totals[elem], len(genes_with[elem]),
                             elem_desc.get(elem, '')])
    print(f"  Wrote summary: {summary_path}")

    # ---- Output 4: Stress-element focused table ----
    stress_elements = {
        'ABRE', 'ABRE_motifIIb',
        'DRE_CRT', 'DRE_core', 'LTRE_CRT',
        'MBS', 'MBS_I', 'MYB_core', 'MYB_drought',
        'TC_rich', 'TC_rich_var',
        'LTR',
        'HSE', 'HSE_core',
        'W_box', 'W_box_cluster',
        'as1_element', 'as1_full',
        'ARE', 'ERE', 'GCC_box',
        'STRE', 'DRE_like',
    }
    stress_path = os.path.join(output_dir, 'stress_elements_counts.tsv')
    stress_sorted = [e for e in sorted_elems if e in stress_elements]
    with open(stress_path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['gene'] + stress_sorted)
        for gene in sorted(gene_elem_counts.keys()):
            row = [gene] + [gene_elem_counts[gene].get(e, 0) for e in stress_sorted]
            writer.writerow(row)
    print(f"  Wrote stress element counts: {stress_path}")


if __name__ == '__main__':
    main()
PYEOF

python3 "${SCRIPTS_DIR}/scan_cis_elements.py" \
    "${PROMOTER_FA}" \
    "${RESULTS_DIR}"

# ===================================================================
# PART 2: FIMO scan against JASPAR plant motifs (optional)
# ===================================================================
echo ""
echo "[$(date '+%H:%M:%S')] Checking for JASPAR plant motif database..."

JASPAR_MEME="${PROJ_DIR}/01_references/JASPAR_plants.meme"

if [[ -f "${JASPAR_MEME}" ]]; then
    echo "  Found JASPAR database: ${JASPAR_MEME}"
    echo "  Running FIMO scan..."

    FIMO_OUT="${RESULTS_DIR}/fimo_output"

    fimo \
        --oc "${FIMO_OUT}" \
        --thresh 1e-4 \
        --max-stored-scores 500000 \
        --verbosity 1 \
        "${JASPAR_MEME}" \
        "${PROMOTER_FA}" \
        2>&1 | tee "${LOG_DIR}/fimo_${SLURM_JOB_ID}.log"

    if [[ -f "${FIMO_OUT}/fimo.tsv" ]]; then
        FIMO_HITS=$(tail -n +2 "${FIMO_OUT}/fimo.tsv" | grep -v '^#' | wc -l)
        echo "  FIMO hits: ${FIMO_HITS}"
    fi
else
    echo "  JASPAR database not found at: ${JASPAR_MEME}"
    echo "  Skipping FIMO scan."
    echo "  To enable: download JASPAR CORE plants (MEME format) from:"
    echo "    https://jaspar.elixir.no/downloads/"
    echo "  Save as: ${JASPAR_MEME}"
fi

# ===================================================================
# PART 3: Visualization — R heatmap
# ===================================================================
echo ""
echo "[$(date '+%H:%M:%S')] Generating heatmap visualizations..."

cat > "${SCRIPTS_DIR}/cis_element_heatmap.R" << 'REOF'
#!/usr/bin/env Rscript
# cis_element_heatmap.R
# Generate heatmaps: (1) all elements, (2) stress-related elements only
# Usage: Rscript cis_element_heatmap.R <counts_tsv> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript cis_element_heatmap.R <counts.tsv> <output_dir>")
}

counts_file <- args[1]
output_dir  <- args[2]

# Install/load packages
for (pkg in c("pheatmap", "RColorBrewer")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
}
library(pheatmap)
library(RColorBrewer)

# Read count matrix
df <- read.delim(counts_file, row.names = 1, check.names = FALSE)
mat <- as.matrix(df)

# Skip if empty
if (nrow(mat) == 0 || ncol(mat) == 0) {
    cat("No data to plot.\n")
    quit(status = 0)
}

# Remove columns (elements) with zero total hits
col_sums <- colSums(mat)
mat <- mat[, col_sums > 0, drop = FALSE]

# Remove rows (genes) with zero total hits
row_sums <- rowSums(mat)
mat <- mat[row_sums > 0, drop = FALSE, ]

if (nrow(mat) == 0 || ncol(mat) == 0) {
    cat("No non-zero data to plot after filtering.\n")
    quit(status = 0)
}

# ---- Heatmap 1: All elements ----
# Color palette
n_colors <- 50
colors <- colorRampPalette(c("white", "#FEE0D2", "#FC9272", "#DE2D26", "#67000D"))(n_colors)

# Determine figure size based on matrix dimensions
fig_height <- max(8, nrow(mat) * 0.25 + 3)
fig_width  <- max(10, ncol(mat) * 0.45 + 4)

out_pdf <- file.path(output_dir, "cis_elements_heatmap_all.pdf")
pdf(out_pdf, width = fig_width, height = fig_height)
pheatmap(
    mat,
    color = colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    display_numbers = TRUE,
    number_format = "%d",
    number_color = "grey30",
    fontsize = 8,
    fontsize_number = 6,
    fontsize_row = 7,
    fontsize_col = 8,
    angle_col = 45,
    main = "Cis-Regulatory Elements in Sunflower Aquaporin Promoters",
    border_color = "grey90"
)
dev.off()
cat("  Wrote:", out_pdf, "\n")

# PNG version
out_png <- file.path(output_dir, "cis_elements_heatmap_all.png")
png(out_png, width = fig_width * 100, height = fig_height * 100, res = 100)
pheatmap(
    mat,
    color = colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    display_numbers = TRUE,
    number_format = "%d",
    number_color = "grey30",
    fontsize = 8,
    fontsize_number = 6,
    fontsize_row = 7,
    fontsize_col = 8,
    angle_col = 45,
    main = "Cis-Regulatory Elements in Sunflower Aquaporin Promoters",
    border_color = "grey90"
)
dev.off()
cat("  Wrote:", out_png, "\n")

# ---- Heatmap 2: Stress elements only ----
stress_file <- file.path(output_dir, "stress_elements_counts.tsv")
if (file.exists(stress_file)) {
    df_stress <- read.delim(stress_file, row.names = 1, check.names = FALSE)
    mat_stress <- as.matrix(df_stress)

    col_sums_s <- colSums(mat_stress)
    mat_stress <- mat_stress[, col_sums_s > 0, drop = FALSE]
    row_sums_s <- rowSums(mat_stress)
    mat_stress <- mat_stress[row_sums_s > 0, drop = FALSE, ]

    if (nrow(mat_stress) > 0 && ncol(mat_stress) > 0) {
        fig_h2 <- max(8, nrow(mat_stress) * 0.25 + 3)
        fig_w2 <- max(10, ncol(mat_stress) * 0.5 + 4)

        out_stress_pdf <- file.path(output_dir, "stress_elements_heatmap.pdf")
        pdf(out_stress_pdf, width = fig_w2, height = fig_h2)
        pheatmap(
            mat_stress,
            color = colorRampPalette(c("white", "#DEEBF7", "#4292C6", "#08306B"))(n_colors),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            clustering_method = "ward.D2",
            display_numbers = TRUE,
            number_format = "%d",
            number_color = "grey30",
            fontsize = 8,
            fontsize_number = 7,
            fontsize_row = 7,
            fontsize_col = 8,
            angle_col = 45,
            main = "Stress-Responsive Cis-Elements in Aquaporin Promoters",
            border_color = "grey90"
        )
        dev.off()
        cat("  Wrote:", out_stress_pdf, "\n")

        out_stress_png <- file.path(output_dir, "stress_elements_heatmap.png")
        png(out_stress_png, width = fig_w2 * 100, height = fig_h2 * 100, res = 100)
        pheatmap(
            mat_stress,
            color = colorRampPalette(c("white", "#DEEBF7", "#4292C6", "#08306B"))(n_colors),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            clustering_method = "ward.D2",
            display_numbers = TRUE,
            number_format = "%d",
            number_color = "grey30",
            fontsize = 8,
            fontsize_number = 7,
            fontsize_row = 7,
            fontsize_col = 8,
            angle_col = 45,
            main = "Stress-Responsive Cis-Elements in Aquaporin Promoters",
            border_color = "grey90"
        )
        dev.off()
        cat("  Wrote:", out_stress_png, "\n")
    }
}

# ---- Element distribution bar chart ----
out_bar <- file.path(output_dir, "element_distribution.pdf")
totals <- sort(colSums(mat), decreasing = TRUE)
totals <- totals[totals > 0]
if (length(totals) > 0) {
    pdf(out_bar, width = max(8, length(totals) * 0.4 + 2), height = 6)
    par(mar = c(10, 5, 3, 1))
    barplot(
        totals,
        las = 2,
        col = brewer.pal(9, "Set1")[1],
        border = NA,
        ylab = "Total occurrences",
        main = "Cis-Element Distribution Across All Aquaporin Promoters",
        cex.names = 0.7
    )
    dev.off()
    cat("  Wrote:", out_bar, "\n")
}

cat("Visualization complete.\n")
REOF

Rscript "${SCRIPTS_DIR}/cis_element_heatmap.R" \
    "${RESULTS_DIR}/cis_elements_counts.tsv" \
    "${RESULTS_DIR}" \
    2>&1 | tee "${LOG_DIR}/cis_heatmap_${SLURM_JOB_ID}.log"

# ===================================================================
# Summary
# ===================================================================
echo ""
echo "================================================================"
echo "Cis-Element Analysis Results"
echo "================================================================"
echo ""
echo "Element type summary:"
column -t -s $'\t' "${RESULTS_DIR}/cis_elements_summary.tsv" | head -30
echo ""
echo "Output files:"
ls -lh "${RESULTS_DIR}/"
echo ""
echo "================================================================"
echo "Pipeline completed: $(date)"
echo "================================================================"
