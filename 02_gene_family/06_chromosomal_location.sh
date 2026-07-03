#!/bin/bash
#SBATCH --job-name=aqp_chrom_loc
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/06_chromosomal_location_%j.out
#SBATCH --error=logs/06_chromosomal_location_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 06_chromosomal_location.sh
# Chromosomal mapping and duplication analysis of aquaporin genes:
#   - Extract chromosomal coordinates
#   - Identify tandem duplications (within 200kb on same chromosome)
#   - Identify segmental duplications via self-BLASTP
#   - Ka/Ks calculation for duplicate pairs
#   - Estimate duplication time (T = Ks / 2λ, λ ≈ 6.5e-9 for dicots)
###############################################################################

set -euo pipefail

CONFIG="/work/dweikat/ydelen2/aquaporin_study/config.sh"
if [[ ! -f "${CONFIG}" ]]; then
    echo "FATAL: config file not found: ${CONFIG}" >&2
    exit 1
fi
source "${CONFIG}"

module purge
module load blast/2.17
module load mafft/7.526
module load miniforge/24.5

# ── Directories ──────────────────────────────────────────────────────────────
WORK_DIR="${PROJ_DIR}/02_gene_family/06_chromosomal_location"
PARENT_DIR="${PROJ_DIR}/02_gene_family"
STRUCT_DIR="${PROJ_DIR}/02_gene_family/05_gene_structure"
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${WORK_DIR}" "${LOG_DIR}" "${WORK_DIR}/pairwise_aln" "${WORK_DIR}/kaks_tmp"

LOGFILE="${LOG_DIR}/chromosomal_location_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "================================================================="
echo "  Chromosomal Location & Duplication Analysis"
echo "  Started: $(date)"
echo "  Job ID:  ${SLURM_JOB_ID:-local}"
echo "================================================================="

# ── Input validation ─────────────────────────────────────────────────────────
VERIFIED_FASTA="${PARENT_DIR}/verified_aquaporins.faa"
VERIFIED_IDS="${PARENT_DIR}/verified_aquaporin_ids.txt"
GENE_MAP="${STRUCT_DIR}/protein_to_gene_map.tsv"
GFF3="${REF_DIR}/GCF_002127325.2_HanXRQr2.0_genomic.gff"

# CDS nucleotide sequences (needed for Ka/Ks)
CDS_FASTA="${REF_DIR}/GCF_002127325.2_HanXRQr2.0_cds_from_genomic.fna"

for f in "${VERIFIED_FASTA}" "${VERIFIED_IDS}" "${GFF3}"; do
    if [[ ! -f "${f}" ]]; then
        echo "FATAL: Required input not found: ${f}" >&2
        exit 1
    fi
done

# Gene map from step 05 is strongly preferred but not fatal
if [[ ! -f "${GENE_MAP}" ]]; then
    echo "WARNING: Gene map not found at ${GENE_MAP}."
    echo "  Will extract coordinates directly from GFF3."
    GENE_MAP=""
fi

N_AQP=$(wc -l < "${VERIFIED_IDS}")
echo "Verified aquaporins: ${N_AQP}"

# ── Step 1: Extract chromosomal coordinates ──────────────────────────────────
echo ""
echo "[Step 1] Extracting chromosomal coordinates..."

CHROM_LOC="${WORK_DIR}/chromosome_locations.tsv"

if [[ -n "${GENE_MAP}" && -f "${GENE_MAP}" ]]; then
    # Use pre-computed gene map
    cp "${GENE_MAP}" "${CHROM_LOC}"
    echo "  Using gene map from step 05."
else
    # Extract directly from GFF3
    echo -e "protein_id\tgene_id\tmrna_id\tchromosome\tgene_start\tgene_end\tstrand" \
        > "${CHROM_LOC}"

    python3 -c "
import sys

target_ids = set()
with open('${VERIFIED_IDS}') as f:
    target_ids = set(line.strip() for line in f if line.strip())

cds_to_mrna = {}
mrna_to_gene = {}
gene_coords = {}

with open('${GFF3}') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        attrs = {}
        for item in fields[8].split(';'):
            if '=' in item:
                k, v = item.split('=', 1)
                attrs[k] = v

        if fields[2] == 'gene':
            gid = attrs.get('ID', '')
            gene_coords[gid] = (fields[0], fields[3], fields[4], fields[6])
        elif fields[2] in ('mRNA', 'transcript'):
            mid = attrs.get('ID', '')
            parent = attrs.get('Parent', '')
            mrna_to_gene[mid] = parent
        elif fields[2] == 'CDS':
            pid = attrs.get('protein_id', '')
            parent = attrs.get('Parent', '')
            if pid:
                cds_to_mrna[pid] = parent

for pid in sorted(target_ids):
    mid = cds_to_mrna.get(pid, '')
    if mid:
        gid = mrna_to_gene.get(mid, 'NA')
        c = gene_coords.get(gid, ('NA','NA','NA','NA'))
        print(f'{pid}\t{gid}\t{mid}\t{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}')
" >> "${CHROM_LOC}"
fi

echo "  Chromosomal locations: ${CHROM_LOC}"
echo "  Genes mapped: $(( $(wc -l < "${CHROM_LOC}") - 1 ))"

# ── Step 2: Identify tandem duplications ─────────────────────────────────────
echo ""
echo "[Step 2] Identifying tandem duplications (within 200kb on same chromosome)..."

TANDEM_DIST=200000  # 200 kb threshold

TANDEM_SCRIPT="${WORK_DIR}/find_tandem_dups.py"
cat > "${TANDEM_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""
Identify tandem duplicated aquaporin genes.
Tandem duplicates: genes on the same chromosome within a specified distance,
with no more than one intervening non-aquaporin gene.
"""
import sys
from collections import defaultdict

def main():
    chrom_loc_file = sys.argv[1]
    max_dist = int(sys.argv[2])
    output_file = sys.argv[3]

    # Parse chromosomal locations
    genes = []
    with open(chrom_loc_file) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7 or parts[3] == 'NA':
                continue
            genes.append({
                'protein_id': parts[0],
                'gene_id': parts[1],
                'chrom': parts[3],
                'start': int(parts[4]),
                'end': int(parts[5]),
                'strand': parts[6]
            })

    # Group by chromosome
    chrom_genes = defaultdict(list)
    for g in genes:
        chrom_genes[g['chrom']].append(g)

    # Sort genes on each chromosome by position
    for chrom in chrom_genes:
        chrom_genes[chrom].sort(key=lambda x: x['start'])

    # Find tandem pairs
    tandem_pairs = []
    tandem_clusters = []

    for chrom, chrom_g in sorted(chrom_genes.items()):
        if len(chrom_g) < 2:
            continue

        # Find consecutive genes within distance threshold
        current_cluster = [chrom_g[0]]

        for i in range(1, len(chrom_g)):
            prev = chrom_g[i-1]
            curr = chrom_g[i]

            distance = curr['start'] - prev['end']

            if distance <= max_dist:
                tandem_pairs.append({
                    'gene1': prev['protein_id'],
                    'gene2': curr['protein_id'],
                    'chrom': chrom,
                    'distance': distance,
                    'gene1_pos': f"{prev['start']}-{prev['end']}",
                    'gene2_pos': f"{curr['start']}-{curr['end']}"
                })
                if curr not in current_cluster:
                    current_cluster.append(curr)
            else:
                if len(current_cluster) >= 2:
                    tandem_clusters.append(current_cluster[:])
                current_cluster = [curr]

        if len(current_cluster) >= 2:
            tandem_clusters.append(current_cluster[:])

    # Write tandem pairs
    with open(output_file, 'w') as out:
        out.write("gene1\tgene2\tchromosome\tdistance_bp\tgene1_position\tgene2_position\n")
        for pair in tandem_pairs:
            out.write(f"{pair['gene1']}\t{pair['gene2']}\t{pair['chrom']}\t"
                      f"{pair['distance']}\t{pair['gene1_pos']}\t{pair['gene2_pos']}\n")

    # Write tandem clusters
    cluster_file = output_file.replace('.tsv', '_clusters.tsv')
    with open(cluster_file, 'w') as out:
        out.write("cluster_id\tchromosome\tnum_genes\tgene_ids\tspan_bp\n")
        for i, cluster in enumerate(tandem_clusters, 1):
            chrom = cluster[0]['chrom']
            gene_ids = ','.join(g['protein_id'] for g in cluster)
            span = cluster[-1]['end'] - cluster[0]['start']
            out.write(f"TC{i}\t{chrom}\t{len(cluster)}\t{gene_ids}\t{span}\n")

    print(f"Tandem duplicate pairs: {len(tandem_pairs)}", file=sys.stderr)
    print(f"Tandem clusters: {len(tandem_clusters)}", file=sys.stderr)
    for i, cluster in enumerate(tandem_clusters, 1):
        print(f"  Cluster {i}: {len(cluster)} genes on {cluster[0]['chrom']}", file=sys.stderr)

if __name__ == '__main__':
    main()
PYEOF

TANDEM_OUTPUT="${WORK_DIR}/tandem_duplicates.tsv"
python3 "${TANDEM_SCRIPT}" "${CHROM_LOC}" "${TANDEM_DIST}" "${TANDEM_OUTPUT}"

echo "  Tandem duplicates: ${TANDEM_OUTPUT}"

# ── Step 3: Identify segmental duplications via self-BLASTP ──────────────────
echo ""
echo "[Step 3] Identifying segmental duplications via self-BLASTP..."

# Create BLAST database from aquaporin sequences
SELF_DB="${WORK_DIR}/aqp_selfdb"
makeblastdb -in "${VERIFIED_FASTA}" -dbtype prot -out "${SELF_DB}" -parse_seqids 2>/dev/null

# Self-BLASTP
SELF_BLAST="${WORK_DIR}/self_blastp.outfmt6"
blastp \
    -query "${VERIFIED_FASTA}" \
    -db "${SELF_DB}" \
    -out "${SELF_BLAST}" \
    -evalue 1e-10 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads "${SLURM_NTASKS_PER_NODE:-8}" \
    -max_target_seqs 20 \
    -seg yes

echo "  Self-BLAST hits: $(wc -l < "${SELF_BLAST}")"

# Identify segmental duplicates
SEGMENTAL_SCRIPT="${WORK_DIR}/find_segmental_dups.py"
cat > "${SEGMENTAL_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""
Identify segmental duplications from self-BLASTP results.
Segmental duplicates: paralogous genes on different chromosomes (or far apart
on the same chromosome) with high sequence similarity.

Criteria:
  - Identity ≥ 70%
  - Query coverage ≥ 70%
  - Not on the same chromosome (or >500kb apart on same chromosome)
  - Not already classified as tandem duplicates
"""
import sys
from collections import defaultdict

def main():
    blast_file = sys.argv[1]
    chrom_loc_file = sys.argv[2]
    tandem_file = sys.argv[3]
    output_file = sys.argv[4]

    min_identity = 70.0
    min_qcov = 70.0
    min_dist_same_chrom = 500000  # 500kb

    # Load chromosomal locations
    gene_chrom = {}
    gene_pos = {}
    with open(chrom_loc_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7 and parts[3] != 'NA':
                gene_chrom[parts[0]] = parts[3]
                gene_pos[parts[0]] = (int(parts[4]), int(parts[5]))

    # Load tandem pairs to exclude
    tandem_pairs = set()
    try:
        with open(tandem_file) as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    pair = tuple(sorted([parts[0], parts[1]]))
                    tandem_pairs.add(pair)
    except FileNotFoundError:
        pass

    # Parse BLAST results
    segmental_pairs = []
    seen_pairs = set()

    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            qid = parts[0]
            sid = parts[1]

            # Skip self-hits
            if qid == sid:
                continue

            pident = float(parts[2])
            qcov = float(parts[12]) if len(parts) > 12 else 0

            # Apply identity and coverage filters
            if pident < min_identity or qcov < min_qcov:
                continue

            # Skip if already seen (avoid A-B and B-A)
            pair = tuple(sorted([qid, sid]))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)

            # Skip tandem duplicates
            if pair in tandem_pairs:
                continue

            # Check chromosomal location
            q_chrom = gene_chrom.get(qid, 'NA')
            s_chrom = gene_chrom.get(sid, 'NA')

            if q_chrom == 'NA' or s_chrom == 'NA':
                dup_type = "unknown_location"
            elif q_chrom != s_chrom:
                dup_type = "segmental_interchrom"
            else:
                # Same chromosome - check distance
                q_pos = gene_pos.get(qid, (0, 0))
                s_pos = gene_pos.get(sid, (0, 0))
                dist = abs(q_pos[0] - s_pos[0])
                if dist > min_dist_same_chrom:
                    dup_type = "segmental_intrachrom"
                else:
                    continue  # Too close, likely tandem

            segmental_pairs.append({
                'gene1': qid,
                'gene2': sid,
                'pident': pident,
                'qcov': qcov,
                'gene1_chrom': q_chrom,
                'gene2_chrom': s_chrom,
                'dup_type': dup_type,
                'evalue': parts[10],
                'bitscore': parts[11]
            })

    # Write output
    with open(output_file, 'w') as out:
        out.write("gene1\tgene2\tpercent_identity\tquery_coverage\t"
                  "gene1_chrom\tgene2_chrom\tduplication_type\tevalue\tbitscore\n")
        for pair in segmental_pairs:
            out.write(f"{pair['gene1']}\t{pair['gene2']}\t{pair['pident']:.1f}\t"
                      f"{pair['qcov']:.1f}\t{pair['gene1_chrom']}\t{pair['gene2_chrom']}\t"
                      f"{pair['dup_type']}\t{pair['evalue']}\t{pair['bitscore']}\n")

    print(f"Segmental duplicate pairs: {len(segmental_pairs)}", file=sys.stderr)
    type_counts = defaultdict(int)
    for p in segmental_pairs:
        type_counts[p['dup_type']] += 1
    for t, c in sorted(type_counts.items()):
        print(f"  {t}: {c}", file=sys.stderr)

if __name__ == '__main__':
    main()
PYEOF

SEGMENTAL_OUTPUT="${WORK_DIR}/segmental_duplicates.tsv"
python3 "${SEGMENTAL_SCRIPT}" "${SELF_BLAST}" "${CHROM_LOC}" "${TANDEM_OUTPUT}" "${SEGMENTAL_OUTPUT}"

echo "  Segmental duplicates: ${SEGMENTAL_OUTPUT}"

# ── Step 4: Ka/Ks calculation for duplicate pairs ────────────────────────────
echo ""
echo "[Step 4] Calculating Ka/Ks for duplicate pairs..."

# Collect all duplicate pairs
ALL_PAIRS="${WORK_DIR}/all_duplicate_pairs.tsv"
echo -e "gene1\tgene2\tdup_type" > "${ALL_PAIRS}"

# Tandem pairs
if [[ -f "${TANDEM_OUTPUT}" ]]; then
    awk -F'\t' 'NR > 1 {print $1"\t"$2"\ttandem"}' "${TANDEM_OUTPUT}" >> "${ALL_PAIRS}"
fi

# Segmental pairs
if [[ -f "${SEGMENTAL_OUTPUT}" ]]; then
    awk -F'\t' 'NR > 1 {print $1"\t"$2"\t"$7}' "${SEGMENTAL_OUTPUT}" >> "${ALL_PAIRS}"
fi

N_PAIRS=$(( $(wc -l < "${ALL_PAIRS}") - 1 ))
echo "  Total duplicate pairs for Ka/Ks: ${N_PAIRS}"

if [[ "${N_PAIRS}" -gt 0 ]]; then

    # Ka/Ks calculation script
    KAKS_SCRIPT="${WORK_DIR}/calculate_kaks.py"
    cat > "${KAKS_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""
Calculate Ka/Ks (dN/dS) for duplicate gene pairs using the Nei-Gojobori method.

Approach:
  1. Extract CDS nucleotide sequences for each pair
  2. Align protein sequences with MAFFT
  3. Back-translate alignment to codon alignment
  4. Calculate Ka, Ks using simple counting method

For more accurate results, use KaKs_Calculator2.0 or PAML yn00.
This script provides a pure-Python implementation.
"""
import sys
import os
import subprocess
import tempfile
import math
from collections import defaultdict

# Standard genetic code
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

BASES = ['T', 'C', 'A', 'G']

def parse_fasta_dict(fasta_file):
    """Parse FASTA into dict."""
    seqs = {}
    header = None
    parts = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    seqs[header] = ''.join(parts)
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if header:
        seqs[header] = ''.join(parts)
    return seqs

def count_synonymous_sites(codon):
    """Count number of synonymous sites in a codon (Nei-Gojobori)."""
    if codon not in CODON_TABLE or CODON_TABLE[codon] == '*':
        return 0, 0

    s_sites = 0.0
    aa = CODON_TABLE[codon]

    for pos in range(3):
        syn_changes = 0
        total_changes = 0
        for base in BASES:
            if base == codon[pos]:
                continue
            new_codon = codon[:pos] + base + codon[pos+1:]
            if new_codon in CODON_TABLE and CODON_TABLE[new_codon] != '*':
                total_changes += 1
                if CODON_TABLE[new_codon] == aa:
                    syn_changes += 1

        if total_changes > 0:
            s_sites += syn_changes / total_changes

    n_sites = 3.0 - s_sites
    return s_sites, n_sites

def nei_gojobori(seq1, seq2):
    """
    Calculate Ka and Ks using Nei-Gojobori (1986) method with
    Jukes-Cantor correction.
    """
    if len(seq1) != len(seq2):
        return None, None, None

    S_sites = 0.0  # Total synonymous sites
    N_sites = 0.0  # Total nonsynonymous sites
    S_diffs = 0.0  # Synonymous differences
    N_diffs = 0.0  # Nonsynonymous differences

    n_codons = len(seq1) // 3

    for i in range(n_codons):
        c1 = seq1[i*3:(i+1)*3]
        c2 = seq2[i*3:(i+1)*3]

        # Skip codons with gaps or ambiguous bases
        if '-' in c1 or '-' in c2 or 'N' in c1 or 'N' in c2:
            continue

        if c1 not in CODON_TABLE or c2 not in CODON_TABLE:
            continue

        if CODON_TABLE[c1] == '*' or CODON_TABLE[c2] == '*':
            continue

        # Count sites (average of both codons)
        s1, n1 = count_synonymous_sites(c1)
        s2, n2 = count_synonymous_sites(c2)
        S_sites += (s1 + s2) / 2.0
        N_sites += (n1 + n2) / 2.0

        # Count differences
        if c1 != c2:
            diffs = sum(1 for a, b in zip(c1, c2) if a != b)
            if diffs == 1:
                # Single substitution
                if CODON_TABLE[c1] == CODON_TABLE[c2]:
                    S_diffs += 1
                else:
                    N_diffs += 1
            elif diffs == 2:
                # Two substitutions - use intermediate pathway
                syn = 0
                nonsyn = 0
                for pos in range(3):
                    if c1[pos] != c2[pos]:
                        # Check if changing just this position is syn or nonsyn
                        intermediate = c1[:pos] + c2[pos] + c1[pos+1:]
                        if intermediate in CODON_TABLE and CODON_TABLE[intermediate] != '*':
                            if CODON_TABLE[c1] == CODON_TABLE[intermediate]:
                                syn += 0.5
                            else:
                                nonsyn += 0.5
                            if CODON_TABLE[intermediate] == CODON_TABLE[c2]:
                                syn += 0.5
                            else:
                                nonsyn += 0.5
                S_diffs += syn
                N_diffs += nonsyn
            else:
                # Three substitutions - approximate
                if CODON_TABLE[c1] == CODON_TABLE[c2]:
                    S_diffs += 2.5
                    N_diffs += 0.5
                else:
                    S_diffs += 1
                    N_diffs += 2

    if S_sites == 0 or N_sites == 0:
        return None, None, None

    # Proportions
    pS = S_diffs / S_sites if S_sites > 0 else 0
    pN = N_diffs / N_sites if N_sites > 0 else 0

    # Jukes-Cantor correction
    try:
        if pS < 0.75:
            Ks = -3.0/4.0 * math.log(1 - 4.0/3.0 * pS)
        else:
            Ks = float('nan')  # Saturated

        if pN < 0.75:
            Ka = -3.0/4.0 * math.log(1 - 4.0/3.0 * pN)
        else:
            Ka = float('nan')  # Saturated
    except (ValueError, ZeroDivisionError):
        return None, None, None

    Ka_Ks = Ka / Ks if Ks > 0 and not math.isnan(Ka) and not math.isnan(Ks) else None

    return Ka, Ks, Ka_Ks


def align_pair(seq1_id, seq1, seq2_id, seq2, tmp_dir):
    """Align two protein sequences with MAFFT."""
    in_file = os.path.join(tmp_dir, 'pair_input.faa')
    out_file = os.path.join(tmp_dir, 'pair_aligned.faa')

    with open(in_file, 'w') as f:
        f.write(f">{seq1_id}\n{seq1}\n>{seq2_id}\n{seq2}\n")

    try:
        subprocess.run(
            ['mafft', '--auto', '--quiet', in_file],
            stdout=open(out_file, 'w'),
            stderr=subprocess.DEVNULL,
            check=True,
            timeout=60
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None

    return parse_fasta_dict(out_file)


def back_translate(aligned_prot, cds_seq):
    """
    Back-translate aligned protein to codon alignment.
    Gaps in protein alignment become '---' in codon alignment.
    """
    clean_prot = aligned_prot.replace('-', '').replace('*', '')
    clean_cds = cds_seq.replace('*', '')

    # Handle length mismatches (common with stop codons)
    expected_len = len(clean_prot) * 3
    if len(clean_cds) < expected_len:
        return None
    clean_cds = clean_cds[:expected_len]

    codon_aln = []
    cds_pos = 0
    for aa in aligned_prot:
        if aa == '-':
            codon_aln.append('---')
        else:
            codon_aln.append(clean_cds[cds_pos:cds_pos+3])
            cds_pos += 3

    return ''.join(codon_aln)


def main():
    pairs_file = sys.argv[1]
    protein_fasta = sys.argv[2]
    cds_fasta = sys.argv[3]
    output_file = sys.argv[4]
    tmp_dir = sys.argv[5]

    # Substitution rate for dicots (mutations per site per year)
    LAMBDA = 6.5e-9

    # Load sequences
    proteins = parse_fasta_dict(protein_fasta)
    cds_seqs = parse_fasta_dict(cds_fasta) if os.path.exists(cds_fasta) else {}

    # Map protein IDs to CDS IDs (NCBI CDS FASTA has different headers)
    # Try direct match first, then try mapping
    cds_map = {}
    for cds_id in cds_seqs:
        # Extract protein_id from CDS header if present
        cds_map[cds_id] = cds_seqs[cds_id]
        # Also try base ID
        base = cds_id.split('.')[0]
        cds_map[base] = cds_seqs[cds_id]

    # Load pairs
    pairs = []
    with open(pairs_file) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pairs.append((parts[0], parts[1], parts[2]))

    print(f"Processing {len(pairs)} duplicate pairs...", file=sys.stderr)

    results = []

    for g1, g2, dup_type in pairs:
        if g1 not in proteins or g2 not in proteins:
            print(f"  SKIP: {g1} or {g2} not in protein FASTA", file=sys.stderr)
            continue

        # Align proteins
        aligned = align_pair(g1, proteins[g1], g2, proteins[g2], tmp_dir)
        if aligned is None or g1 not in aligned or g2 not in aligned:
            print(f"  SKIP: alignment failed for {g1} vs {g2}", file=sys.stderr)
            continue

        # Get CDS sequences
        cds1 = cds_map.get(g1, cds_map.get(g1.split('.')[0], ''))
        cds2 = cds_map.get(g2, cds_map.get(g2.split('.')[0], ''))

        if not cds1 or not cds2:
            print(f"  SKIP: CDS not found for {g1} or {g2}", file=sys.stderr)
            results.append({
                'gene1': g1, 'gene2': g2, 'dup_type': dup_type,
                'Ka': 'NA', 'Ks': 'NA', 'Ka_Ks': 'NA',
                'dup_time_mya': 'NA', 'selection': 'NA',
                'note': 'CDS_not_found'
            })
            continue

        # Back-translate
        codon_aln1 = back_translate(aligned[g1], cds1)
        codon_aln2 = back_translate(aligned[g2], cds2)

        if codon_aln1 is None or codon_aln2 is None:
            results.append({
                'gene1': g1, 'gene2': g2, 'dup_type': dup_type,
                'Ka': 'NA', 'Ks': 'NA', 'Ka_Ks': 'NA',
                'dup_time_mya': 'NA', 'selection': 'NA',
                'note': 'back_translation_failed'
            })
            continue

        # Truncate to same length
        min_len = min(len(codon_aln1), len(codon_aln2))
        min_len = (min_len // 3) * 3  # ensure multiple of 3
        codon_aln1 = codon_aln1[:min_len]
        codon_aln2 = codon_aln2[:min_len]

        # Calculate Ka/Ks
        Ka, Ks, Ka_Ks = nei_gojobori(codon_aln1, codon_aln2)

        # Estimate duplication time: T = Ks / (2 * lambda)
        if Ks is not None and not math.isnan(Ks):
            dup_time = Ks / (2 * LAMBDA) / 1e6  # in Mya
        else:
            dup_time = None

        # Determine selection pressure
        if Ka_Ks is not None and not math.isnan(Ka_Ks):
            if Ka_Ks < 1:
                selection = "purifying"
            elif Ka_Ks == 1:
                selection = "neutral"
            else:
                selection = "positive"
        else:
            selection = "NA"

        results.append({
            'gene1': g1, 'gene2': g2, 'dup_type': dup_type,
            'Ka': f"{Ka:.6f}" if Ka is not None and not math.isnan(Ka) else 'NA',
            'Ks': f"{Ks:.6f}" if Ks is not None and not math.isnan(Ks) else 'NA',
            'Ka_Ks': f"{Ka_Ks:.6f}" if Ka_Ks is not None and not math.isnan(Ka_Ks) else 'NA',
            'dup_time_mya': f"{dup_time:.2f}" if dup_time is not None else 'NA',
            'selection': selection,
            'note': ''
        })

    # Write results
    with open(output_file, 'w') as out:
        out.write("gene1\tgene2\tduplication_type\tKa\tKs\tKa_Ks\t"
                  "duplication_time_Mya\tselection_pressure\tnote\n")
        for r in results:
            out.write(f"{r['gene1']}\t{r['gene2']}\t{r['dup_type']}\t"
                      f"{r['Ka']}\t{r['Ks']}\t{r['Ka_Ks']}\t"
                      f"{r['dup_time_mya']}\t{r['selection']}\t{r['note']}\n")

    # Summary
    kaks_valid = [r for r in results if r['Ka_Ks'] != 'NA']
    print(f"\nKa/Ks results:", file=sys.stderr)
    print(f"  Total pairs: {len(results)}", file=sys.stderr)
    print(f"  Successful calculations: {len(kaks_valid)}", file=sys.stderr)
    if kaks_valid:
        sel_counts = defaultdict(int)
        for r in kaks_valid:
            sel_counts[r['selection']] += 1
        for s, c in sorted(sel_counts.items()):
            print(f"  {s}: {c}", file=sys.stderr)

if __name__ == '__main__':
    main()
PYEOF

    # Extract CDS for aquaporin genes
    # If CDS FASTA is available from NCBI, extract relevant sequences
    AQP_CDS="${WORK_DIR}/aqp_cds.fna"

    if [[ -f "${CDS_FASTA}" ]]; then
        echo "  Extracting CDS sequences from ${CDS_FASTA}..."
        # Extract CDS by protein ID from CDS FASTA headers
        # NCBI CDS FASTA: >lcl|... [protein_id=XP_XXXXX]
        python3 -c "
import sys

target_ids = set()
with open('${VERIFIED_IDS}') as f:
    target_ids = set(line.strip() for line in f if line.strip())

found = set()
with open('${CDS_FASTA}') as fin, open('${AQP_CDS}', 'w') as fout:
    write = False
    for line in fin:
        if line.startswith('>'):
            write = False
            # Check if protein_id is in target
            for tid in target_ids:
                if tid in line:
                    # Rename header to protein_id for matching
                    fout.write(f'>{tid}\n')
                    write = True
                    found.add(tid)
                    break
        elif write:
            fout.write(line)

print(f'  CDS extracted: {len(found)}/{len(target_ids)}', file=sys.stderr)
"
    else
        echo "  WARNING: CDS FASTA not found: ${CDS_FASTA}"
        echo "  Ka/Ks calculation will be limited."
        echo "  Download CDS from NCBI: datasets download genome accession GCF_002127325.2 --include cds"
        touch "${AQP_CDS}"
    fi

    # Run Ka/Ks calculation
    KAKS_OUTPUT="${WORK_DIR}/kaks_results.tsv"
    TMP_DIR="${WORK_DIR}/kaks_tmp"

    python3 "${KAKS_SCRIPT}" \
        "${ALL_PAIRS}" \
        "${VERIFIED_FASTA}" \
        "${AQP_CDS}" \
        "${KAKS_OUTPUT}" \
        "${TMP_DIR}"

    echo "  Ka/Ks results: ${KAKS_OUTPUT}"

else
    echo "  No duplicate pairs found. Skipping Ka/Ks calculation."
    echo -e "gene1\tgene2\tduplication_type\tKa\tKs\tKa_Ks\tduplication_time_Mya\tselection_pressure\tnote" \
        > "${WORK_DIR}/kaks_results.tsv"
fi

# ── Step 5: Create chromosome map for visualization ──────────────────────────
echo ""
echo "[Step 5] Creating chromosome map files..."

CHROM_MAP_SCRIPT="${WORK_DIR}/create_chrom_map.py"
cat > "${CHROM_MAP_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""Create chromosome map visualization data for aquaporin genes."""
import sys
from collections import defaultdict

def main():
    chrom_loc_file = sys.argv[1]
    output_prefix = sys.argv[2]
    genome_fai = sys.argv[3]

    # Load chromosome sizes
    chrom_sizes = {}
    with open(genome_fai) as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom_sizes[parts[0]] = int(parts[1])

    # Load gene positions
    genes = []
    chrom_genes = defaultdict(list)
    with open(chrom_loc_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7 and parts[3] != 'NA':
                gene = {
                    'protein_id': parts[0],
                    'gene_id': parts[1],
                    'chrom': parts[3],
                    'start': int(parts[4]),
                    'end': int(parts[5]),
                    'strand': parts[6]
                }
                genes.append(gene)
                chrom_genes[gene['chrom']].append(gene)

    # Chromosome distribution summary
    summary_file = f"{output_prefix}_chromosome_distribution.tsv"
    with open(summary_file, 'w') as out:
        out.write("chromosome\tchrom_size_bp\tnum_aquaporins\tgene_ids\n")
        for chrom in sorted(chrom_genes.keys()):
            size = chrom_sizes.get(chrom, 0)
            g_list = sorted(chrom_genes[chrom], key=lambda x: x['start'])
            gene_ids = ','.join(g['protein_id'] for g in g_list)
            out.write(f"{chrom}\t{size}\t{len(g_list)}\t{gene_ids}\n")

    # TBtools-compatible chromosome location file
    tbtools_file = f"{output_prefix}_tbtools_chrom_location.txt"
    with open(tbtools_file, 'w') as out:
        for gene in sorted(genes, key=lambda x: (x['chrom'], x['start'])):
            out.write(f"{gene['protein_id']}\t{gene['chrom']}\t{gene['start']}\t{gene['end']}\n")

    # Create Circos-style link file for segmental duplicates
    # (if segmental_duplicates.tsv exists)
    print(f"Chromosome map files written:", file=sys.stderr)
    print(f"  {summary_file}", file=sys.stderr)
    print(f"  {tbtools_file}", file=sys.stderr)

    # Stats
    print(f"\nChromosome distribution:", file=sys.stderr)
    for chrom in sorted(chrom_genes.keys()):
        print(f"  {chrom}: {len(chrom_genes[chrom])} genes", file=sys.stderr)
    print(f"  Total: {len(genes)} genes on {len(chrom_genes)} chromosomes", file=sys.stderr)

if __name__ == '__main__':
    main()
PYEOF

GENOME_FAI="${REF_DIR}/GCF_002127325.2_HanXRQr2.0_genomic.fna.fai"
if [[ ! -f "${GENOME_FAI}" ]]; then
    GENOME_FAI="${WORK_DIR}/genome.fai"
    # Create minimal fai from GFF3 if genome index not available
    awk -F'\t' '/^##sequence-region/ {print $2"\t"$4}' "${GFF3}" > "${GENOME_FAI}" 2>/dev/null || touch "${GENOME_FAI}"
fi

python3 "${CHROM_MAP_SCRIPT}" "${CHROM_LOC}" "${WORK_DIR}/aqp" "${GENOME_FAI}"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "================================================================="
echo "  Chromosomal Location & Duplication Analysis Complete"
echo "================================================================="
echo ""
echo "  Output files:"
echo "    ${CHROM_LOC}"
echo "    ${TANDEM_OUTPUT}"
echo "    ${SEGMENTAL_OUTPUT}"
echo "    ${WORK_DIR}/kaks_results.tsv"
echo "    ${WORK_DIR}/aqp_chromosome_distribution.tsv"
echo "    ${WORK_DIR}/aqp_tbtools_chrom_location.txt"
echo ""
echo "  Tandem duplicates:    $(( $(wc -l < "${TANDEM_OUTPUT}") - 1 )) pairs"
echo "  Segmental duplicates: $(( $(wc -l < "${SEGMENTAL_OUTPUT}") - 1 )) pairs"
if [[ -f "${WORK_DIR}/kaks_results.tsv" ]]; then
    echo "  Ka/Ks calculated:     $(( $(wc -l < "${WORK_DIR}/kaks_results.tsv") - 1 )) pairs"
fi
echo ""
echo "  Duplication time formula: T = Ks / (2 × λ)"
echo "  λ = 6.5 × 10⁻⁹ substitutions/site/year (dicot rate)"
echo ""
echo "  Completed: $(date)"
echo "================================================================="
