#!/bin/bash
#SBATCH --job-name=aqp_gene_struct
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/05_gene_structure_%j.out
#SBATCH --error=logs/05_gene_structure_%j.err
#SBATCH --mail-type=END,FAIL

###############################################################################
# 05_gene_structure.sh
# Gene structure analysis of verified aquaporin genes:
#   - Exon-intron structure from GFF3 annotation
#   - Exon counts per gene
#   - Intron lengths
#   - 2kb upstream promoter sequences
#   - Output formatted for TBtools visualization
###############################################################################

set -euo pipefail

CONFIG="/work/dweikat/ydelen2/aquaporin_study/config.sh"
if [[ ! -f "${CONFIG}" ]]; then
    echo "FATAL: config file not found: ${CONFIG}" >&2
    exit 1
fi
source "${CONFIG}"

module purge
module load bedtools/2.31
module load samtools/1.23

# ── Directories ──────────────────────────────────────────────────────────────
WORK_DIR="${PROJ_DIR}/02_gene_family/05_gene_structure"
PARENT_DIR="${PROJ_DIR}/02_gene_family"
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${WORK_DIR}" "${LOG_DIR}"

LOGFILE="${LOG_DIR}/gene_structure_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "================================================================="
echo "  Gene Structure Analysis of Aquaporins"
echo "  Started: $(date)"
echo "  Job ID:  ${SLURM_JOB_ID:-local}"
echo "================================================================="

# ── Input validation ─────────────────────────────────────────────────────────
VERIFIED_IDS="${PARENT_DIR}/verified_aquaporin_ids.txt"
GFF3="${REF_DIR}/GCF_002127325.2_HanXRQr2.0_genomic.gff"
GENOME="${REF_DIR}/GCF_002127325.2_HanXRQr2.0_genomic.fna"

for f in "${VERIFIED_IDS}" "${GFF3}" "${GENOME}"; do
    if [[ ! -f "${f}" ]]; then
        echo "FATAL: Required input not found: ${f}" >&2
        exit 1
    fi
done

echo "Verified aquaporins: $(wc -l < "${VERIFIED_IDS}")"
echo "GFF3: ${GFF3}"
echo "Genome: ${GENOME}"

# Index genome if needed
if [[ ! -f "${GENOME}.fai" ]]; then
    echo "  Indexing genome..."
    samtools faidx "${GENOME}"
fi

# ── Step 1: Map protein IDs to gene IDs in GFF3 ─────────────────────────────
echo ""
echo "[Step 1] Mapping protein IDs to gene/mRNA features in GFF3..."

# NCBI GFF3 uses protein_id attribute on CDS features
# We need to find gene → mRNA → CDS linkage

GENE_MAP="${WORK_DIR}/protein_to_gene_map.tsv"
echo -e "protein_id\tgene_id\tmrna_id\tchromosome\tgene_start\tgene_end\tstrand" \
    > "${GENE_MAP}"

# Build protein → gene mapping from GFF3
MAPPING_SCRIPT="${WORK_DIR}/map_protein_to_gene.py"
cat > "${MAPPING_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""
Map protein IDs from NCBI GFF3 to gene features.
Handles GFF3 parent-child relationships: gene → mRNA → CDS (protein_id)
"""
import sys
import re

def parse_attributes(attr_str):
    """Parse GFF3 attributes into dict."""
    attrs = {}
    for item in attr_str.strip().split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key] = val
    return attrs

def main():
    gff_file = sys.argv[1]
    protein_ids_file = sys.argv[2]
    output_file = sys.argv[3]

    # Load target protein IDs
    with open(protein_ids_file) as f:
        target_ids = set(line.strip() for line in f if line.strip())

    # Parse GFF3 in two passes:
    # Pass 1: Build CDS protein_id → mRNA parent mapping
    # Pass 2: Build mRNA → gene parent mapping

    cds_to_mrna = {}      # protein_id → mRNA_id
    mrna_to_gene = {}     # mRNA_id → gene_id
    gene_coords = {}      # gene_id → (chr, start, end, strand)
    mrna_coords = {}      # mRNA_id → (chr, start, end, strand)

    print(f"Parsing GFF3: {gff_file}", file=sys.stderr)

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feat_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_attributes(fields[8])

            if feat_type == 'gene':
                gid = attrs.get('ID', '')
                if gid:
                    gene_coords[gid] = (chrom, start, end, strand)

            elif feat_type in ('mRNA', 'transcript'):
                mid = attrs.get('ID', '')
                parent = attrs.get('Parent', '')
                if mid and parent:
                    mrna_to_gene[mid] = parent
                    mrna_coords[mid] = (chrom, start, end, strand)

            elif feat_type == 'CDS':
                protein_id = attrs.get('protein_id', '')
                parent = attrs.get('Parent', '')
                if protein_id and parent:
                    cds_to_mrna[protein_id] = parent

    print(f"  Genes: {len(gene_coords)}", file=sys.stderr)
    print(f"  mRNAs: {len(mrna_to_gene)}", file=sys.stderr)
    print(f"  CDS→protein mappings: {len(cds_to_mrna)}", file=sys.stderr)

    # Map protein IDs to genes
    found = 0
    not_found = []
    with open(output_file, 'w') as out:
        out.write("protein_id\tgene_id\tmrna_id\tchromosome\tgene_start\tgene_end\tstrand\n")

        for pid in sorted(target_ids):
            mrna_id = cds_to_mrna.get(pid, '')
            if not mrna_id:
                # Try without version number
                pid_base = pid.split('.')[0]
                for k, v in cds_to_mrna.items():
                    if k.split('.')[0] == pid_base:
                        mrna_id = v
                        break

            if mrna_id:
                gene_id = mrna_to_gene.get(mrna_id, 'NA')
                if gene_id in gene_coords:
                    c = gene_coords[gene_id]
                    out.write(f"{pid}\t{gene_id}\t{mrna_id}\t{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\n")
                    found += 1
                elif mrna_id in mrna_coords:
                    c = mrna_coords[mrna_id]
                    out.write(f"{pid}\t{gene_id}\t{mrna_id}\t{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\n")
                    found += 1
                else:
                    not_found.append(pid)
            else:
                not_found.append(pid)

    print(f"  Mapped: {found}/{len(target_ids)}", file=sys.stderr)
    if not_found:
        print(f"  Not found: {len(not_found)}", file=sys.stderr)
        for nf in not_found[:10]:
            print(f"    {nf}", file=sys.stderr)

if __name__ == '__main__':
    main()
PYEOF

python3 "${MAPPING_SCRIPT}" "${GFF3}" "${VERIFIED_IDS}" "${GENE_MAP}"

# Extract gene IDs for downstream use
GENE_IDS="${WORK_DIR}/aqp_gene_ids.txt"
awk -F'\t' 'NR > 1 {print $2}' "${GENE_MAP}" | sort -u > "${GENE_IDS}"
echo "  Gene IDs: $(wc -l < "${GENE_IDS}")"

# ── Step 2: Extract exon-intron structure ────────────────────────────────────
echo ""
echo "[Step 2] Extracting exon-intron structure from GFF3..."

EXON_INTRON_SCRIPT="${WORK_DIR}/extract_gene_structure.py"
cat > "${EXON_INTRON_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""
Extract exon-intron structure for aquaporin genes from GFF3.
Outputs:
  - exon_structure.tsv: per-exon coordinates
  - intron_lengths.tsv: intron coordinates and lengths
  - gene_structure_summary.tsv: summary per gene
  - TBtools gene structure input files
"""
import sys
from collections import defaultdict

def parse_attributes(attr_str):
    attrs = {}
    for item in attr_str.strip().split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key] = val
    return attrs

def main():
    gff_file = sys.argv[1]
    gene_map_file = sys.argv[2]
    output_prefix = sys.argv[3]

    # Load gene → protein mapping
    gene_to_protein = {}
    mrna_ids = set()
    with open(gene_map_file) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                gene_to_protein[parts[1]] = parts[0]  # gene_id → protein_id
                mrna_ids.add(parts[2])  # mRNA IDs

    target_genes = set(gene_to_protein.keys())

    # Parse GFF3: collect exon and CDS features for target genes
    # Structure: gene → mRNA → exon/CDS
    mrna_parent = {}    # mRNA_id → gene_id
    exons = defaultdict(list)    # mRNA_id → [(start, end, strand), ...]
    cds_features = defaultdict(list)
    gene_features = {}   # gene_id → (chr, start, end, strand)
    utr5 = defaultdict(list)
    utr3 = defaultdict(list)

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feat_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_attributes(fields[8])

            if feat_type == 'gene':
                gid = attrs.get('ID', '')
                if gid in target_genes:
                    gene_features[gid] = (chrom, start, end, strand)

            elif feat_type in ('mRNA', 'transcript'):
                mid = attrs.get('ID', '')
                parent = attrs.get('Parent', '')
                if parent in target_genes:
                    mrna_parent[mid] = parent

            elif feat_type == 'exon':
                parent = attrs.get('Parent', '')
                if parent in mrna_ids or parent in mrna_parent:
                    exons[parent].append((start, end, strand, chrom))

            elif feat_type == 'CDS':
                parent = attrs.get('Parent', '')
                if parent in mrna_ids or parent in mrna_parent:
                    cds_features[parent].append((start, end, strand, chrom))

            elif feat_type == 'five_prime_UTR':
                parent = attrs.get('Parent', '')
                if parent in mrna_ids or parent in mrna_parent:
                    utr5[parent].append((start, end))

            elif feat_type == 'three_prime_UTR':
                parent = attrs.get('Parent', '')
                if parent in mrna_ids or parent in mrna_parent:
                    utr3[parent].append((start, end))

    # ── Output exon structure ────────────────────────────────────────────
    exon_file = f"{output_prefix}_exon_structure.tsv"
    with open(exon_file, 'w') as out:
        out.write("gene_id\tprotein_id\tmrna_id\tchromosome\texon_num\texon_start\texon_end\texon_length\tstrand\n")

        for mrna_id in sorted(mrna_parent.keys()):
            gene_id = mrna_parent[mrna_id]
            protein_id = gene_to_protein.get(gene_id, 'NA')

            if mrna_id not in exons:
                continue

            # Sort exons by position
            sorted_exons = sorted(exons[mrna_id], key=lambda x: x[0])

            for i, (s, e, strand, chrom) in enumerate(sorted_exons, 1):
                exon_len = e - s + 1
                out.write(f"{gene_id}\t{protein_id}\t{mrna_id}\t{chrom}\t{i}\t{s}\t{e}\t{exon_len}\t{strand}\n")

    # ── Output intron lengths ────────────────────────────────────────────
    intron_file = f"{output_prefix}_intron_lengths.tsv"
    with open(intron_file, 'w') as out:
        out.write("gene_id\tprotein_id\tchromosome\tintron_num\tintron_start\tintron_end\tintron_length\n")

        for mrna_id in sorted(mrna_parent.keys()):
            gene_id = mrna_parent[mrna_id]
            protein_id = gene_to_protein.get(gene_id, 'NA')

            if mrna_id not in exons or len(exons[mrna_id]) < 2:
                continue

            sorted_exons = sorted(exons[mrna_id], key=lambda x: x[0])
            chrom = sorted_exons[0][3]

            for i in range(len(sorted_exons) - 1):
                intron_start = sorted_exons[i][1] + 1
                intron_end = sorted_exons[i+1][0] - 1
                intron_len = intron_end - intron_start + 1

                if intron_len > 0:
                    out.write(f"{gene_id}\t{protein_id}\t{chrom}\t{i+1}\t{intron_start}\t{intron_end}\t{intron_len}\n")

    # ── Summary per gene ─────────────────────────────────────────────────
    summary_file = f"{output_prefix}_gene_structure_summary.tsv"
    with open(summary_file, 'w') as out:
        out.write("gene_id\tprotein_id\tchromosome\tgene_start\tgene_end\tstrand\t"
                  "gene_length\tnum_exons\ttotal_exon_length\tnum_introns\t"
                  "total_intron_length\tavg_intron_length\tnum_utr5\tnum_utr3\n")

        for mrna_id in sorted(mrna_parent.keys()):
            gene_id = mrna_parent[mrna_id]
            protein_id = gene_to_protein.get(gene_id, 'NA')

            if gene_id not in gene_features:
                continue

            chrom, gstart, gend, strand = gene_features[gene_id]
            gene_len = gend - gstart + 1

            sorted_exons = sorted(exons.get(mrna_id, []), key=lambda x: x[0])
            num_exons = len(sorted_exons)
            total_exon_len = sum(e - s + 1 for s, e, _, _ in sorted_exons)

            # Calculate introns
            intron_lengths = []
            for i in range(len(sorted_exons) - 1):
                il = sorted_exons[i+1][0] - sorted_exons[i][1] - 1
                if il > 0:
                    intron_lengths.append(il)

            num_introns = len(intron_lengths)
            total_intron_len = sum(intron_lengths)
            avg_intron_len = total_intron_len / num_introns if num_introns > 0 else 0

            n_utr5 = len(utr5.get(mrna_id, []))
            n_utr3 = len(utr3.get(mrna_id, []))

            out.write(f"{gene_id}\t{protein_id}\t{chrom}\t{gstart}\t{gend}\t{strand}\t"
                      f"{gene_len}\t{num_exons}\t{total_exon_len}\t{num_introns}\t"
                      f"{total_intron_len}\t{avg_intron_len:.1f}\t{n_utr5}\t{n_utr3}\n")

    # ── TBtools gene structure format ────────────────────────────────────
    # TBtools Gene Structure View expects: gene_id, feature_type, start, end
    tbtools_file = f"{output_prefix}_tbtools_gene_structure.txt"
    with open(tbtools_file, 'w') as out:
        for mrna_id in sorted(mrna_parent.keys()):
            gene_id = mrna_parent[mrna_id]
            protein_id = gene_to_protein.get(gene_id, 'NA')

            if gene_id not in gene_features:
                continue

            chrom, gstart, gend, strand = gene_features[gene_id]

            # Gene line
            out.write(f"{protein_id}\t{gstart}\t{gend}\n")

            # Exons (CDS) - use CDS features for TBtools
            for s, e, _, _ in sorted(cds_features.get(mrna_id, exons.get(mrna_id, [])),
                                      key=lambda x: x[0]):
                out.write(f"{protein_id}\tCDS\t{s}\t{e}\n")

            # UTRs
            for s, e in sorted(utr5.get(mrna_id, [])):
                out.write(f"{protein_id}\t5UTR\t{s}\t{e}\n")
            for s, e in sorted(utr3.get(mrna_id, [])):
                out.write(f"{protein_id}\t3UTR\t{s}\t{e}\n")

    print(f"Gene structure files written:", file=sys.stderr)
    print(f"  {exon_file}", file=sys.stderr)
    print(f"  {intron_file}", file=sys.stderr)
    print(f"  {summary_file}", file=sys.stderr)
    print(f"  {tbtools_file}", file=sys.stderr)

    # Print summary stats
    print(f"\nExon count distribution:", file=sys.stderr)
    exon_counts = defaultdict(int)
    for mrna_id in mrna_parent:
        n = len(exons.get(mrna_id, []))
        exon_counts[n] += 1
    for n in sorted(exon_counts.keys()):
        print(f"  {n} exons: {exon_counts[n]} genes", file=sys.stderr)

if __name__ == '__main__':
    main()
PYEOF

OUTPUT_PREFIX="${WORK_DIR}/aqp"
python3 "${EXON_INTRON_SCRIPT}" "${GFF3}" "${GENE_MAP}" "${OUTPUT_PREFIX}"

# ── Step 3: Extract 2kb upstream promoter sequences ─────────────────────────
echo ""
echo "[Step 3] Extracting 2kb upstream promoter sequences..."

# Create BED file for promoter regions
PROMOTER_BED="${WORK_DIR}/promoter_regions_2kb.bed"
PROMOTER_FASTA="${WORK_DIR}/promoter_sequences_2kb.fasta"

awk -F'\t' 'NR > 1 {
    chrom = $4
    start = $5
    end = $6
    strand = $7
    gene_id = $2
    prot_id = $1

    if (strand == "+") {
        # Promoter is upstream of TSS on + strand
        prom_end = start - 1
        prom_start = prom_end - 2000
        if (prom_start < 0) prom_start = 0
    } else {
        # Promoter is downstream of gene end on - strand
        prom_start = end
        prom_end = prom_start + 2000
    }

    # BED is 0-based half-open
    printf "%s\t%d\t%d\t%s\t0\t%s\n", chrom, prom_start, prom_end, prot_id, strand
}' "${GENE_MAP}" > "${PROMOTER_BED}"

echo "  Promoter BED: ${PROMOTER_BED}"
echo "  Promoter regions: $(wc -l < "${PROMOTER_BED}")"

# Get chromosome sizes for BED validation
CHROM_SIZES="${WORK_DIR}/chrom_sizes.txt"
awk -F'\t' '{print $1"\t"$2}' "${GENOME}.fai" > "${CHROM_SIZES}"

# Clip to chromosome boundaries and extract sequences
PROMOTER_BED_CLIPPED="${WORK_DIR}/promoter_regions_2kb_clipped.bed"
bedtools slop -i "${PROMOTER_BED}" -g "${CHROM_SIZES}" -b 0 \
    | awk -F'\t' '{if ($3 > $2) print}' \
    > "${PROMOTER_BED_CLIPPED}"

bedtools getfasta \
    -fi "${GENOME}" \
    -bed "${PROMOTER_BED_CLIPPED}" \
    -fo "${PROMOTER_FASTA}" \
    -name \
    -s

echo "  Promoter sequences extracted: $(grep -c '^>' "${PROMOTER_FASTA}" || echo 0)"
echo "  Output: ${PROMOTER_FASTA}"

# ── Step 4: Create GFF3 subset for TBtools ───────────────────────────────────
echo ""
echo "[Step 4] Creating GFF3 subset for TBtools visualization..."

AQP_GFF3="${WORK_DIR}/aqp_genes.gff3"

# Extract relevant lines from GFF3 for aquaporin genes
python3 -c "
import sys

gene_ids = set()
with open('${GENE_IDS}') as f:
    gene_ids = set(line.strip() for line in f if line.strip())

# Also collect all child IDs
child_ids = set()
keep_parents = set()

# Two-pass: first identify all relevant feature IDs
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

        feat_id = attrs.get('ID', '')
        parent = attrs.get('Parent', '')

        if feat_id in gene_ids:
            child_ids.add(feat_id)
        if parent in gene_ids or parent in child_ids:
            child_ids.add(feat_id)

# Second pass: extract all features belonging to target genes
with open('${AQP_GFF3}', 'w') as out:
    out.write('##gff-version 3\n')
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

            feat_id = attrs.get('ID', '')
            parent = attrs.get('Parent', '')

            if feat_id in gene_ids or feat_id in child_ids or parent in gene_ids or parent in child_ids:
                out.write(line)

print(f'  AQP GFF3 features extracted', file=sys.stderr)
"

echo "  AQP GFF3: ${AQP_GFF3}"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "================================================================="
echo "  Gene Structure Analysis Complete"
echo "================================================================="
echo "  Output files:"
echo "    ${OUTPUT_PREFIX}_exon_structure.tsv"
echo "    ${OUTPUT_PREFIX}_intron_lengths.tsv"
echo "    ${OUTPUT_PREFIX}_gene_structure_summary.tsv"
echo "    ${OUTPUT_PREFIX}_tbtools_gene_structure.txt"
echo "    ${PROMOTER_FASTA}"
echo "    ${AQP_GFF3}"
echo ""
echo "  For TBtools visualization:"
echo "    1. Load aqp_genes.gff3 in Gene Structure View"
echo "    2. Or use tbtools_gene_structure.txt for custom view"
echo "    3. Add phylogenetic tree on the left for publication figure"
echo ""
echo "  Completed: $(date)"
echo "================================================================="
