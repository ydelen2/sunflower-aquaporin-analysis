#!/bin/bash
#SBATCH --job-name=aqp_motif
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/05_motif_%j.out
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/05_motif_%j.err
#SBATCH --mail-type=END,FAIL

# ============================================================================
# 02_motif_analysis.sh
# MEME motif discovery on sunflower aquaporin protein sequences.
# Identifies conserved protein motifs (up to 20), then parses output
# into a summary table mapping motifs to genes.
# ============================================================================

set -eo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJ_DIR='/work/dweikat/ydelen2/aquaporin_study'
source "${PROJ_DIR}/config.sh"

WORK_DIR="${PROJ_DIR}/05_phylogenetics"
RESULTS_DIR="${WORK_DIR}/results/meme_motifs"
LOG_DIR="${PROJ_DIR}/logs"
GENE_FAM_DIR="${PROJ_DIR}/02_gene_family/results"
SCRIPTS_DIR="${PROJ_DIR}/05_phylogenetics/scripts"

mkdir -p "${RESULTS_DIR}" "${LOG_DIR}" "${SCRIPTS_DIR}"

# Sunflower aquaporin protein sequences
INPUT_FA="${GENE_FAM_DIR}/aquaporin_proteins.fa"

if [[ ! -f "${INPUT_FA}" ]]; then
    echo "ERROR: Input file not found: ${INPUT_FA}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Modules
# ---------------------------------------------------------------------------
module purge
module load compiler/gcc/11 openmpi/4.1 MEME/5.5

echo "================================================================"
echo "Pipeline: MEME Motif Discovery"
echo "Started: $(date)"
echo "Job ID:  ${SLURM_JOB_ID}"
echo "================================================================"

SEQ_COUNT=$(grep -c '^>' "${INPUT_FA}")
echo "Input sequences: ${SEQ_COUNT}"

# ---------------------------------------------------------------------------
# Step 1: Run MEME motif discovery
# ---------------------------------------------------------------------------
MEME_OUT="${RESULTS_DIR}/meme_output"

echo "[$(date '+%H:%M:%S')] Running MEME..."

meme "${INPUT_FA}" \
    -protein \
    -oc "${MEME_OUT}" \
    -nmotifs 20 \
    -minw 6 \
    -maxw 50 \
    -mod zoops \
    -p "${SLURM_NTASKS_PER_NODE}" \
    -nostatus \
    2>&1 | tee "${LOG_DIR}/meme_${SLURM_JOB_ID}.log"

echo "[$(date '+%H:%M:%S')] MEME completed."

# Verify key output files exist
for f in "${MEME_OUT}/meme.xml" "${MEME_OUT}/meme.html" "${MEME_OUT}/meme.txt"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Expected MEME output not found: ${f}" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# Step 2: Parse MEME output to create motif summary tables
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Parsing MEME results..."

cat > "${SCRIPTS_DIR}/parse_meme.py" << 'PYEOF'
#!/usr/bin/env python3
"""
Parse MEME XML output and generate:
  1. motif_summary.tsv        - motif ID, width, sites, e-value, consensus
  2. gene_motif_table.tsv     - gene x motif presence/absence + positions
  3. motif_details.tsv        - per-site detail: gene, motif, start, end, p-value, sequence
"""

import sys
import xml.etree.ElementTree as ET
import csv
import os
from collections import defaultdict

def parse_meme_xml(xml_path, output_dir):
    tree = ET.parse(xml_path)
    root = tree.getroot()

    # ----- Extract motif metadata -----
    motifs_elem = root.find('.//motifs')
    if motifs_elem is None:
        print("ERROR: No <motifs> element found in MEME XML", file=sys.stderr)
        sys.exit(1)

    motif_info = []
    motif_ids = []
    for motif in motifs_elem.findall('motif'):
        mid = motif.get('id', motif.get('name', ''))
        motif_ids.append(mid)
        motif_info.append({
            'motif_id': mid,
            'name': motif.get('name', mid),
            'alt_name': motif.get('alt', ''),
            'width': motif.get('width', ''),
            'sites': motif.get('sites', ''),
            'evalue': motif.get('e_value', ''),
        })

    # Write motif summary
    summary_path = os.path.join(output_dir, 'motif_summary.tsv')
    with open(summary_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['motif_id', 'name', 'alt_name', 'width', 'sites', 'evalue'],
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(motif_info)
    print(f"  Wrote {len(motif_info)} motifs to {summary_path}")

    # ----- Extract per-site hits -----
    # gene -> motif_id -> list of (start, end, pvalue, site_seq)
    gene_motif_hits = defaultdict(lambda: defaultdict(list))
    detail_rows = []

    for motif in motifs_elem.findall('motif'):
        mid = motif.get('id', motif.get('name', ''))
        mwidth = int(motif.get('width', 0))
        sites = motif.find('contributing_sites')
        if sites is None:
            continue
        for site in sites.findall('contributing_site'):
            seq_id = site.get('sequence_id', '')
            position = int(site.get('position', 0))
            strand = site.get('strand', '+')
            pvalue = site.get('pvalue', '')

            # Resolve sequence_id to gene name via sequences section
            seq_name = seq_id  # fallback

            detail_rows.append({
                'sequence_id': seq_id,
                'motif_id': mid,
                'start': position + 1,  # 1-based
                'end': position + mwidth,
                'strand': strand,
                'pvalue': pvalue,
            })
            gene_motif_hits[seq_id][mid].append(position + 1)

    # Resolve sequence IDs to names
    seq_id_to_name = {}
    sequences = root.find('.//sequences')
    if sequences is not None:
        # Try training_set first (MEME 5.x)
        training = root.find('.//training_set')
        if training is not None:
            for seq in training.findall('sequence'):
                seq_id_to_name[seq.get('id', '')] = seq.get('name', seq.get('id', ''))
        else:
            for seq in sequences.findall('sequence'):
                seq_id_to_name[seq.get('id', '')] = seq.get('name', seq.get('id', ''))
    # Also try top-level training_set
    training = root.find('.//training_set')
    if training is not None:
        for seq in training.findall('sequence'):
            sid = seq.get('id', '')
            if sid not in seq_id_to_name:
                seq_id_to_name[sid] = seq.get('name', sid)

    # Map seq IDs to names in detail rows
    for row in detail_rows:
        row['gene'] = seq_id_to_name.get(row['sequence_id'], row['sequence_id'])

    # Write motif details
    details_path = os.path.join(output_dir, 'motif_details.tsv')
    with open(details_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['gene', 'motif_id', 'start', 'end', 'strand', 'pvalue'],
                                delimiter='\t')
        writer.writeheader()
        for row in detail_rows:
            writer.writerow({k: row[k] for k in ['gene', 'motif_id', 'start', 'end', 'strand', 'pvalue']})
    print(f"  Wrote {len(detail_rows)} motif sites to {details_path}")

    # ----- Gene x Motif presence/absence table -----
    all_genes = sorted(set(seq_id_to_name.values()))
    table_path = os.path.join(output_dir, 'gene_motif_table.tsv')
    with open(table_path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        header = ['gene'] + motif_ids
        writer.writerow(header)
        for sid, gname in sorted(seq_id_to_name.items(), key=lambda x: x[1]):
            row = [gname]
            for mid in motif_ids:
                hits = gene_motif_hits[sid].get(mid, [])
                if hits:
                    row.append(','.join(str(p) for p in hits))
                else:
                    row.append('-')
            writer.writerow(row)
    print(f"  Wrote gene x motif table for {len(all_genes)} genes to {table_path}")

    # ----- Motif count summary per gene -----
    count_path = os.path.join(output_dir, 'gene_motif_counts.tsv')
    with open(count_path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        header = ['gene', 'total_motifs'] + motif_ids
        writer.writerow(header)
        for sid, gname in sorted(seq_id_to_name.items(), key=lambda x: x[1]):
            counts = [len(gene_motif_hits[sid].get(mid, [])) for mid in motif_ids]
            row = [gname, sum(counts)] + counts
            writer.writerow(row)
    print(f"  Wrote motif counts to {count_path}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <meme.xml> <output_dir>", file=sys.stderr)
        sys.exit(1)
    parse_meme_xml(sys.argv[1], sys.argv[2])
PYEOF

python3 "${SCRIPTS_DIR}/parse_meme.py" \
    "${MEME_OUT}/meme.xml" \
    "${RESULTS_DIR}"

# ---------------------------------------------------------------------------
# Step 3: Print summary
# ---------------------------------------------------------------------------
echo ""
echo "================================================================"
echo "MEME Motif Discovery Results"
echo "================================================================"
echo ""
echo "Motifs found:"
column -t -s $'\t' "${RESULTS_DIR}/motif_summary.tsv" | head -25
echo ""
echo "Output directory: ${RESULTS_DIR}"
echo ""
ls -lh "${RESULTS_DIR}/"
echo ""
echo "MEME HTML report: ${MEME_OUT}/meme.html"
echo "================================================================"
echo "Pipeline completed: $(date)"
echo "================================================================"
