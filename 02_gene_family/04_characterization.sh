#!/bin/bash
#SBATCH --job-name=aqp_characterize
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/04_characterization_%j.out
#SBATCH --error=logs/04_characterization_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 04_characterization.sh
# Physicochemical characterization and subfamily classification of verified
# aquaporin proteins. Uses Biopython for protein analysis and custom logic
# for NPA motif / ar/R selectivity filter identification.
###############################################################################

set -euo pipefail

CONFIG="/work/dweikat/ydelen2/aquaporin_study/config.sh"
if [[ ! -f "${CONFIG}" ]]; then
    echo "FATAL: config file not found: ${CONFIG}" >&2
    exit 1
fi
source "${CONFIG}"

module purge
module load miniforge/24.5
module load mafft/7.526

# ── Directories ──────────────────────────────────────────────────────────────
WORK_DIR="${PROJ_DIR}/02_gene_family/04_characterization"
PARENT_DIR="${PROJ_DIR}/02_gene_family"
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${WORK_DIR}" "${LOG_DIR}"

LOGFILE="${LOG_DIR}/characterization_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "================================================================="
echo "  Aquaporin Characterization"
echo "  Started: $(date)"
echo "  Job ID:  ${SLURM_JOB_ID:-local}"
echo "================================================================="

# ── Input validation ─────────────────────────────────────────────────────────
VERIFIED_FASTA="${PARENT_DIR}/verified_aquaporins.faa"
VERIFIED_TABLE="${PARENT_DIR}/verified_aquaporins.tsv"

for f in "${VERIFIED_FASTA}" "${VERIFIED_TABLE}"; do
    if [[ ! -f "${f}" ]]; then
        echo "FATAL: Required input not found: ${f}" >&2
        echo "  Run 03_domain_verification.sh first." >&2
        exit 1
    fi
done

N_AQP=$(grep -c '^>' "${VERIFIED_FASTA}")
echo "Verified aquaporins: ${N_AQP}"

# ── Ensure Biopython is available ────────────────────────────────────────────
echo ""
echo "[Dependency check] Biopython..."
python3 -c "from Bio.SeqUtils.ProtParam import ProteinAnalysis; print('  Biopython available')" 2>/dev/null || {
    echo "  Installing Biopython..."
    pip install --user biopython
}

# ── Step 1: Create characterization Python script ────────────────────────────
echo ""
echo "[Step 1] Running physicochemical characterization..."

CHAR_SCRIPT="${WORK_DIR}/characterize_aquaporins.py"
cat > "${CHAR_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""
Comprehensive characterization of aquaporin proteins.

Calculates:
  - Molecular weight (MW)
  - Isoelectric point (pI)
  - GRAVY (Grand Average of Hydropathicity)
  - Instability index
  - Amino acid count and length
  - NPA motif identification (aquaporin signature)
  - ar/R selectivity filter residues (H2, H5, LE1, LE2)
  - Preliminary subfamily classification

Output: comprehensive TSV table
"""

import sys
import os
import re
from collections import OrderedDict

try:
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.SeqUtils import molecular_weight
except ImportError:
    print("FATAL: Biopython required. Install with: pip install biopython", file=sys.stderr)
    sys.exit(1)


def parse_fasta(fasta_file):
    """Parse FASTA using Biopython."""
    records = OrderedDict()
    for record in SeqIO.parse(fasta_file, "fasta"):
        records[record.id] = str(record.seq).upper().replace('*', '')
    return records


def find_npa_motifs(sequence):
    """
    Find NPA (Asn-Pro-Ala) motifs in aquaporin sequence.
    Canonical aquaporins have two NPA motifs in loops B and E.
    Also check for NPA variants: NPT, NPV, NPS, NPG, etc.
    Returns list of (position, motif_sequence) tuples.
    """
    npa_motifs = []
    # Strict NPA
    for m in re.finditer(r'NP[ASTVG]', sequence):
        npa_motifs.append((m.start() + 1, m.group()))  # 1-indexed

    return npa_motifs


def identify_ar_r_filter(sequence, npa_positions):
    """
    Identify ar/R selectivity filter residues.

    The ar/R filter consists of 4 residues at positions relative to the
    NPA motifs and specific TM helices:
      H2:  ~13 residues before the first NPA (in TM helix 2)
      H5:  ~9 residues after the second NPA (in TM helix 5)
      LE1: ~1 residue after the first NPA in loop E (actually in the second half)
      LE2: ~2 residues after LE1

    These positions are approximate and depend on alignment with reference
    aquaporins. For precise positions, use MSA with known aquaporin structures.

    Returns dict with H2, H5, LE1, LE2 residues.
    """
    ar_r = {'H2': 'NA', 'H5': 'NA', 'LE1': 'NA', 'LE2': 'NA'}

    if len(npa_positions) < 2:
        return ar_r

    # Sort NPA positions
    npa_pos = sorted([p[0] for p in npa_positions])[:2]

    npa1 = npa_pos[0] - 1  # 0-indexed
    npa2 = npa_pos[1] - 1  # 0-indexed

    # Approximate positions (these vary by protein, alignment-based is better)
    # H2: ~13 residues before first NPA
    h2_pos = npa1 - 13
    if 0 <= h2_pos < len(sequence):
        ar_r['H2'] = sequence[h2_pos]

    # H5: ~9 residues after second NPA
    h5_pos = npa2 + 9
    if 0 <= h5_pos < len(sequence):
        ar_r['H5'] = sequence[h5_pos]

    # LE1 and LE2: residues in loop E, after second NPA
    le1_pos = npa2 + 1
    le2_pos = npa2 + 2
    if 0 <= le1_pos < len(sequence):
        ar_r['LE1'] = sequence[le1_pos]
    if 0 <= le2_pos < len(sequence):
        ar_r['LE2'] = sequence[le2_pos]

    return ar_r


def classify_subfamily_by_features(sequence, npa_motifs, ar_r, seq_length):
    """
    Preliminary subfamily classification based on sequence features.

    Typical features:
      PIP: ~280-290 aa, NPA/NPA, ar/R: F-H-T-R, high water permeability
      TIP: ~240-260 aa, NPA/NPA or NPA/NPV, ar/R: H-I-A-V or H-I-G-R
      NIP: ~270-300 aa, NPA/NPA or NPA/NPS, ar/R: W-V-A-R or W-I-A-R
      SIP: ~230-250 aa, often NPT/NPA, ar/R: variable, ER-localized
      XIP: ~280-310 aa, NPA/NPA, ar/R: variable (found in dicots)

    This is a heuristic. Use phylogenetic classification for publication.
    """
    subfamily = "unclassified"

    # Get NPA motif strings
    npa_strings = [m[1] for m in npa_motifs[:2]]

    h2 = ar_r.get('H2', '')
    h5 = ar_r.get('H5', '')

    # PIP signatures: typically F at H2, H at H5
    if h2 == 'F' and h5 in ('H', 'R'):
        subfamily = "PIP"
    # TIP signatures: typically H at H2, I at H5
    elif h2 == 'H' and h5 in ('I', 'V', 'A'):
        subfamily = "TIP"
    # NIP signatures: typically W at H2, V or I at H5
    elif h2 == 'W' and h5 in ('V', 'I', 'A'):
        subfamily = "NIP"
    # SIP signatures: shorter sequences, often unusual NPA variants
    elif seq_length < 260:
        if any(m != 'NPA' for m in npa_strings):
            subfamily = "SIP"
        else:
            subfamily = "SIP_tentative"
    # XIP: dicot-specific, check for unique features
    elif h2 in ('V', 'I') and h5 in ('V', 'I', 'T'):
        subfamily = "XIP"
    # Length-based fallback
    elif 275 <= seq_length <= 295:
        subfamily = "PIP_tentative"
    elif 235 <= seq_length <= 265:
        subfamily = "TIP_tentative"

    return subfamily


def predict_localization(subfamily):
    """Predict subcellular localization based on subfamily."""
    loc_map = {
        'PIP': 'Plasma membrane',
        'PIP_tentative': 'Plasma membrane (tentative)',
        'TIP': 'Tonoplast (vacuolar membrane)',
        'TIP_tentative': 'Tonoplast (tentative)',
        'NIP': 'Plasma membrane / ER',
        'SIP': 'Endoplasmic reticulum',
        'SIP_tentative': 'Endoplasmic reticulum (tentative)',
        'XIP': 'Plasma membrane',
        'unclassified': 'Unknown'
    }
    return loc_map.get(subfamily, 'Unknown')


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <verified.faa> <output.tsv>", file=sys.stderr)
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_file = sys.argv[2]

    proteins = parse_fasta(input_fasta)
    print(f"Processing {len(proteins)} aquaporin sequences...", file=sys.stderr)

    header = [
        'protein_id', 'length_aa', 'molecular_weight_Da', 'isoelectric_point',
        'gravy', 'instability_index', 'stability', 'aromaticity',
        'npa_motif_count', 'npa_motif_1', 'npa_motif_1_pos',
        'npa_motif_2', 'npa_motif_2_pos',
        'ar_r_H2', 'ar_r_H5', 'ar_r_LE1', 'ar_r_LE2', 'ar_r_filter',
        'subfamily_predicted', 'subcellular_localization',
        'charged_aa_pct', 'hydrophobic_aa_pct',
        'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly',
        'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser',
        'Thr', 'Trp', 'Tyr', 'Val'
    ]

    results = []

    for prot_id, sequence in proteins.items():
        # Clean sequence (remove ambiguous residues for ProteinAnalysis)
        clean_seq = re.sub(r'[BXZJUO*]', '', sequence)
        if len(clean_seq) < 50:
            print(f"  WARNING: {prot_id} has very short clean sequence ({len(clean_seq)} aa), skipping.", file=sys.stderr)
            continue

        try:
            pa = ProteinAnalysis(clean_seq)

            # Basic properties
            mw = pa.molecular_weight()
            pi = pa.isoelectric_point()
            gravy = pa.gravy()
            instability = pa.instability_index()
            stability = "stable" if instability < 40 else "unstable"
            aromaticity = pa.aromaticity()

            # Amino acid composition
            aa_count = pa.count_amino_acids()
            aa_pct = pa.get_amino_acids_percent()

            # Charged and hydrophobic AA percentages
            charged = sum(aa_pct.get(aa, 0) for aa in 'DEKRH') * 100
            hydrophobic = sum(aa_pct.get(aa, 0) for aa in 'AVILMFYW') * 100

            # NPA motifs (use original sequence)
            npa_motifs = find_npa_motifs(sequence)
            npa_count = len(npa_motifs)
            npa1 = npa_motifs[0][1] if npa_count >= 1 else 'NA'
            npa1_pos = str(npa_motifs[0][0]) if npa_count >= 1 else 'NA'
            npa2 = npa_motifs[1][1] if npa_count >= 2 else 'NA'
            npa2_pos = str(npa_motifs[1][0]) if npa_count >= 2 else 'NA'

            # ar/R selectivity filter
            ar_r = identify_ar_r_filter(sequence, npa_motifs)
            ar_r_str = f"{ar_r['H2']}/{ar_r['H5']}/{ar_r['LE1']}/{ar_r['LE2']}"

            # Subfamily classification
            subfamily = classify_subfamily_by_features(
                sequence, npa_motifs, ar_r, len(sequence)
            )

            # Subcellular localization
            localization = predict_localization(subfamily)

            # Build row
            row = [
                prot_id, len(sequence), f"{mw:.2f}", f"{pi:.2f}",
                f"{gravy:.4f}", f"{instability:.2f}", stability, f"{aromaticity:.4f}",
                npa_count, npa1, npa1_pos, npa2, npa2_pos,
                ar_r['H2'], ar_r['H5'], ar_r['LE1'], ar_r['LE2'], ar_r_str,
                subfamily, localization,
                f"{charged:.1f}", f"{hydrophobic:.1f}",
                aa_count.get('A', 0), aa_count.get('R', 0), aa_count.get('N', 0),
                aa_count.get('D', 0), aa_count.get('C', 0), aa_count.get('Q', 0),
                aa_count.get('E', 0), aa_count.get('G', 0), aa_count.get('H', 0),
                aa_count.get('I', 0), aa_count.get('L', 0), aa_count.get('K', 0),
                aa_count.get('M', 0), aa_count.get('F', 0), aa_count.get('P', 0),
                aa_count.get('S', 0), aa_count.get('T', 0), aa_count.get('W', 0),
                aa_count.get('Y', 0), aa_count.get('V', 0)
            ]

            results.append(row)

        except Exception as e:
            print(f"  ERROR processing {prot_id}: {e}", file=sys.stderr)
            continue

    # Write output
    with open(output_file, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for row in results:
            out.write('\t'.join(str(x) for x in row) + '\n')

    # Print subfamily summary
    subfamilies = {}
    for row in results:
        sf = row[18]  # subfamily_predicted column
        subfamilies[sf] = subfamilies.get(sf, 0) + 1

    print(f"\nSubfamily distribution:", file=sys.stderr)
    for sf in sorted(subfamilies.keys()):
        print(f"  {sf}: {subfamilies[sf]}", file=sys.stderr)
    print(f"  Total: {len(results)}", file=sys.stderr)

    # NPA motif summary
    npa_counts = {}
    for row in results:
        nc = row[8]  # npa_motif_count
        npa_counts[nc] = npa_counts.get(nc, 0) + 1

    print(f"\nNPA motif count distribution:", file=sys.stderr)
    for nc in sorted(npa_counts.keys()):
        print(f"  {nc} NPA: {npa_counts[nc]} proteins", file=sys.stderr)


if __name__ == '__main__':
    main()
PYEOF

OUTPUT_TSV="${WORK_DIR}/aquaporin_characterization.tsv"
python3 "${CHAR_SCRIPT}" "${VERIFIED_FASTA}" "${OUTPUT_TSV}"

echo "  Characterization table: ${OUTPUT_TSV}"

# ── Step 2: Subfamily classification via phylogenetics (MSA-based) ───────────
echo ""
echo "[Step 2] Building MSA with reference aquaporins for subfamily classification..."

# Combine sunflower AQPs with reference sequences for alignment
BLAST_DIR="${PROJ_DIR}/02_gene_family/02_blast_validation"
REF_SEQS="${BLAST_DIR}/query_sequences/combined_aqp_query.faa"

if [[ -f "${REF_SEQS}" ]]; then
    COMBINED_FOR_MSA="${WORK_DIR}/combined_for_classification.faa"
    cat "${VERIFIED_FASTA}" "${REF_SEQS}" > "${COMBINED_FOR_MSA}"

    MSA_OUTPUT="${WORK_DIR}/classification_msa.faa"

    echo "  Running MAFFT alignment..."
    mafft \
        --auto \
        --thread "${SLURM_NTASKS_PER_NODE:-4}" \
        --reorder \
        "${COMBINED_FOR_MSA}" > "${MSA_OUTPUT}" 2>/dev/null

    echo "  MSA for classification: ${MSA_OUTPUT}"
    echo "  Use this MSA to build a neighbor-joining tree for definitive"
    echo "  subfamily classification based on phylogenetic clustering"
    echo "  with reference Arabidopsis/rice aquaporins."
else
    echo "  WARNING: Reference sequences not found at ${REF_SEQS}"
    echo "  Subfamily classification relies on heuristic only."
fi

# ── Step 3: Summary statistics ───────────────────────────────────────────────
echo ""
echo "[Step 3] Generating summary statistics..."

SUMMARY_SCRIPT="${WORK_DIR}/summary_stats.py"
cat > "${SUMMARY_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""Generate summary statistics from characterization table."""
import sys
import csv

infile = sys.argv[1]
outfile = sys.argv[2]

data = []
with open(infile) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        data.append(row)

with open(outfile, 'w') as out:
    out.write("Aquaporin Characterization Summary\n")
    out.write("=" * 50 + "\n\n")

    out.write(f"Total aquaporins: {len(data)}\n\n")

    # Subfamily counts
    subfamilies = {}
    for row in data:
        sf = row['subfamily_predicted']
        subfamilies[sf] = subfamilies.get(sf, 0) + 1

    out.write("Subfamily distribution:\n")
    for sf in sorted(subfamilies.keys()):
        out.write(f"  {sf:20s}: {subfamilies[sf]:3d}\n")

    # MW range
    mws = [float(row['molecular_weight_Da']) for row in data]
    out.write(f"\nMolecular weight range: {min(mws):.0f} - {max(mws):.0f} Da\n")

    # pI range
    pis = [float(row['isoelectric_point']) for row in data]
    out.write(f"Isoelectric point range: {min(pis):.2f} - {max(pis):.2f}\n")

    # GRAVY range
    gravys = [float(row['gravy']) for row in data]
    out.write(f"GRAVY range: {min(gravys):.4f} - {max(gravys):.4f}\n")

    # Length range
    lens = [int(row['length_aa']) for row in data]
    out.write(f"Length range: {min(lens)} - {max(lens)} aa\n")

    # Stability
    stable = sum(1 for row in data if row['stability'] == 'stable')
    out.write(f"\nStability: {stable} stable, {len(data)-stable} unstable\n")

    # NPA motifs
    two_npa = sum(1 for row in data if int(row['npa_motif_count']) == 2)
    out.write(f"Proteins with 2 NPA motifs: {two_npa}/{len(data)}\n")

    # ar/R filter diversity
    ar_r_filters = {}
    for row in data:
        f = row['ar_r_filter']
        ar_r_filters[f] = ar_r_filters.get(f, 0) + 1

    out.write(f"\nar/R selectivity filter diversity:\n")
    for f in sorted(ar_r_filters.keys(), key=lambda x: -ar_r_filters[x]):
        out.write(f"  {f:15s}: {ar_r_filters[f]:3d}\n")

    print(f"Summary written to {outfile}", file=sys.stderr)

PYEOF

SUMMARY_FILE="${WORK_DIR}/characterization_summary.txt"
python3 "${SUMMARY_SCRIPT}" "${OUTPUT_TSV}" "${SUMMARY_FILE}"

echo "  Summary: ${SUMMARY_FILE}"

# ── Copy outputs ─────────────────────────────────────────────────────────────
cp "${OUTPUT_TSV}" "${PARENT_DIR}/aquaporin_characterization.tsv"

echo ""
echo "================================================================="
echo "  Characterization Complete"
echo "  Main output: ${OUTPUT_TSV}"
echo ""
echo "  NOTE: Subfamily classification is heuristic-based."
echo "  For definitive classification, build a phylogenetic tree using"
echo "  the MSA at: ${WORK_DIR}/classification_msa.faa"
echo "  with IQ-TREE or RAxML, then classify by clade membership."
echo ""
echo "  ar/R filter positions are approximate. For precise identification,"
echo "  align with crystallographic reference (e.g., SoPIP2;1, PDB: 1Z98)"
echo "  and extract positions from the structural alignment."
echo ""
echo "  Completed: $(date)"
echo "================================================================="
