#!/usr/bin/env bash
#SBATCH --job-name=synteny_v2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/work/dweikat/ydelen2/aquaporin_study/logs/synteny_v2_%j.log
#SBATCH --error=/work/dweikat/ydelen2/aquaporin_study/logs/synteny_v2_%j.err

set -euo pipefail

PROJ_DIR="/work/dweikat/ydelen2/aquaporin_study"
SYNTENY_DIR="${PROJ_DIR}/07_synteny"
REF_DIR="${PROJ_DIR}/01_references"

HA_GFF="${REF_DIR}/sunflower/genomic.gff"
HA_PROT="${REF_DIR}/sunflower/protein.faa"
AT_PROT="${REF_DIR}/arabidopsis/protein.faa"
LS_PROT="${REF_DIR}/lettuce/protein.faa"
AQP_LIST="${PROJ_DIR}/09_figures/supplementary/Table_S1_aquaporin_genes.tsv"

module load blast/2.17
module load miniforge/24.5
conda activate "${PROJ_DIR}/conda_envs/aquaporin_env"

THREADS=8
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

mkdir -p "${SYNTENY_DIR}"/{data,blast,mcscanx,results}
mkdir -p "${REF_DIR}/arabidopsis" "${REF_DIR}/lettuce"

# =============================================================================
# STEP 1: Download REAL GFF files from NCBI
# =============================================================================
log "STEP 1: Downloading real GFF files from NCBI"

# Arabidopsis thaliana GCF_000001735.4 (TAIR10.1)
AT_GFF="${REF_DIR}/arabidopsis/genomic.gff"
if [[ ! -f "${AT_GFF}" ]]; then
    log "  Downloading Arabidopsis GFF..."
    wget -q -O "${AT_GFF}.gz" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz"
    gunzip "${AT_GFF}.gz"
    log "  Arabidopsis GFF: $(wc -l < ${AT_GFF}) lines"
else
    log "  Arabidopsis GFF already exists"
fi

# Also download Arabidopsis proteome from NCBI (to match GFF protein_ids)
AT_PROT_NCBI="${REF_DIR}/arabidopsis/protein_ncbi.faa"
if [[ ! -f "${AT_PROT_NCBI}" ]]; then
    log "  Downloading Arabidopsis proteome from NCBI..."
    wget -q -O "${AT_PROT_NCBI}.gz" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_protein.faa.gz"
    gunzip "${AT_PROT_NCBI}.gz"
    log "  Arabidopsis proteome: $(grep -c '^>' ${AT_PROT_NCBI}) sequences"
else
    log "  Arabidopsis NCBI proteome already exists"
fi

# Lactuca sativa GCF_002870075.4 (Lsat_Salinas_v11)
LS_GFF="${REF_DIR}/lettuce/genomic.gff"
if [[ ! -f "${LS_GFF}" ]]; then
    log "  Downloading lettuce GFF..."
    wget -q -O "${LS_GFF}.gz" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/870/075/GCF_002870075.4_Lsat_Salinas_v11/GCF_002870075.4_Lsat_Salinas_v11_genomic.gff.gz"
    gunzip "${LS_GFF}.gz"
    log "  Lettuce GFF: $(wc -l < ${LS_GFF}) lines"
else
    log "  Lettuce GFF already exists"
fi

LS_PROT_NCBI="${REF_DIR}/lettuce/protein_ncbi.faa"
if [[ ! -f "${LS_PROT_NCBI}" ]]; then
    log "  Downloading lettuce proteome from NCBI..."
    wget -q -O "${LS_PROT_NCBI}.gz" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/870/075/GCF_002870075.4_Lsat_Salinas_v11/GCF_002870075.4_Lsat_Salinas_v11_protein.faa.gz"
    gunzip "${LS_PROT_NCBI}.gz"
    log "  Lettuce proteome: $(grep -c '^>' ${LS_PROT_NCBI}) sequences"
else
    log "  Lettuce NCBI proteome already exists"
fi

# =============================================================================
# STEP 2: Parse GFFs to create MCScanX-compatible files with protein IDs
# =============================================================================
log "STEP 2: Creating MCScanX GFF files from REAL coordinates"

python3 - "${HA_GFF}" "${AT_GFF}" "${LS_GFF}" "${SYNTENY_DIR}/data" <<'PYGFF'
import re, sys, os

ha_gff, at_gff, ls_gff, data_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
os.makedirs(data_dir, exist_ok=True)

def parse_ncbi_gff(gff_path, species_prefix):
    """Parse NCBI GFF to extract: protein_id -> (chromosome, start, end)
    and gene_id (LOC) -> first protein_id mapping."""
    
    # Pass 1: gene coordinates (gene_id -> chr, start, end)
    genes = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.strip().split('\t')
            if len(p) < 9 or p[2] != 'gene':
                continue
            m = re.search(r'GeneID:(\d+)', p[8])
            if m:
                loc = f'LOC{m.group(1)}'
                genes[loc] = (p[0], int(p[3]), int(p[4]))
    
    # Pass 2: protein_id -> gene mapping from CDS
    prot_to_gene = {}
    gene_to_prot = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.strip().split('\t')
            if len(p) < 9 or p[2] != 'CDS':
                continue
            prot_m = re.search(r'protein_id=([^;]+)', p[8])
            gene_m = re.search(r'GeneID:(\d+)', p[8])
            if prot_m and gene_m:
                prot_id = prot_m.group(1)
                loc = f'LOC{gene_m.group(1)}'
                prot_to_gene[prot_id] = loc
                if loc not in gene_to_prot:
                    gene_to_prot[loc] = prot_id
    
    # Write GFF with protein IDs and REAL coordinates
    entries = []
    for loc, (chrom, start, end) in genes.items():
        if loc in gene_to_prot:
            prot_id = gene_to_prot[loc]
            entries.append((chrom, prot_id, start, end))
    
    # Sort by chromosome then position
    entries.sort(key=lambda x: (x[0], x[2]))
    
    out_path = os.path.join(data_dir, f'{species_prefix}.gff')
    with open(out_path, 'w') as f:
        for chrom, prot_id, start, end in entries:
            f.write(f'{species_prefix}_{chrom}\t{prot_id}\t{start}\t{end}\n')
    
    # Save mapping
    map_path = os.path.join(data_dir, f'{species_prefix}_prot_to_gene.tsv')
    with open(map_path, 'w') as f:
        for prot, gene in sorted(prot_to_gene.items()):
            f.write(f'{prot}\t{gene}\n')
    
    return len(entries), len(prot_to_gene)

# Parse all three species
n_ha, m_ha = parse_ncbi_gff(ha_gff, 'Ha')
print(f'Sunflower:   {n_ha} genes, {m_ha} protein mappings')

n_at, m_at = parse_ncbi_gff(at_gff, 'At')
print(f'Arabidopsis: {n_at} genes, {m_at} protein mappings')

n_ls, m_ls = parse_ncbi_gff(ls_gff, 'Ls')
print(f'Lettuce:     {n_ls} genes, {m_ls} protein mappings')
PYGFF

# =============================================================================
# STEP 3: BLASTp with NCBI proteomes (real protein IDs)
# =============================================================================
log "STEP 3: BLASTp with NCBI proteomes"

SCRATCH="/scratch/${USER}/synteny_v2_$$"
mkdir -p "${SCRATCH}"

# Use NCBI proteomes for At and Ls to match GFF protein_ids
AT_PROT_USE="${REF_DIR}/arabidopsis/protein_ncbi.faa"
LS_PROT_USE="${REF_DIR}/lettuce/protein_ncbi.faa"

# If NCBI proteomes don't exist, fall back to existing ones
[[ ! -f "${AT_PROT_USE}" ]] && AT_PROT_USE="${AT_PROT}"
[[ ! -f "${LS_PROT_USE}" ]] && LS_PROT_USE="${LS_PROT}"

# Ha vs At
log "  BLASTp: Ha vs At..."
cat "${HA_PROT}" "${AT_PROT_USE}" > "${SCRATCH}/Ha_At.faa"
makeblastdb -in "${SCRATCH}/Ha_At.faa" -dbtype prot -out "${SCRATCH}/Ha_At_db"
blastp -query "${SCRATCH}/Ha_At.faa" -db "${SCRATCH}/Ha_At_db" \
       -out "${SYNTENY_DIR}/blast/Ha_At_v2.blast" -evalue 1e-10 \
       -num_threads ${THREADS} -outfmt 6 -max_target_seqs 5
log "  Ha vs At: $(wc -l < ${SYNTENY_DIR}/blast/Ha_At_v2.blast) hits"

# Ha vs Ls
log "  BLASTp: Ha vs Ls..."
cat "${HA_PROT}" "${LS_PROT_USE}" > "${SCRATCH}/Ha_Ls.faa"
makeblastdb -in "${SCRATCH}/Ha_Ls.faa" -dbtype prot -out "${SCRATCH}/Ha_Ls_db"
blastp -query "${SCRATCH}/Ha_Ls.faa" -db "${SCRATCH}/Ha_Ls_db" \
       -out "${SYNTENY_DIR}/blast/Ha_Ls_v2.blast" -evalue 1e-10 \
       -num_threads ${THREADS} -outfmt 6 -max_target_seqs 5
log "  Ha vs Ls: $(wc -l < ${SYNTENY_DIR}/blast/Ha_Ls_v2.blast) hits"

# Ha self (reuse from v1)
if [[ ! -f "${SYNTENY_DIR}/blast/Ha_self.blast" ]]; then
    log "  BLASTp: Ha self..."
    makeblastdb -in "${HA_PROT}" -dbtype prot -out "${SCRATCH}/Ha_self_db"
    blastp -query "${HA_PROT}" -db "${SCRATCH}/Ha_self_db" \
           -out "${SYNTENY_DIR}/blast/Ha_self.blast" -evalue 1e-10 \
           -num_threads ${THREADS} -outfmt 6 -max_target_seqs 5
else
    log "  Ha self BLAST reused from v1"
fi

rm -rf "${SCRATCH}"

# =============================================================================
# STEP 4: Run MCScanX
# =============================================================================
log "STEP 4: MCScanX"

# Ha self (intra-genomic)
log "  MCScanX: Ha self..."
mkdir -p "${SYNTENY_DIR}/mcscanx/Ha_self_v2"
cp "${SYNTENY_DIR}/data/Ha.gff" "${SYNTENY_DIR}/mcscanx/Ha_self_v2/Ha_self_v2.gff"
cp "${SYNTENY_DIR}/blast/Ha_self.blast" "${SYNTENY_DIR}/mcscanx/Ha_self_v2/Ha_self_v2.blast"
cd "${SYNTENY_DIR}/mcscanx/Ha_self_v2"
MCScanX Ha_self_v2 -s 5 2>&1 || true
if [[ -f Ha_self_v2.collinearity ]]; then
    NGENES=$(grep "Number of collinear genes" Ha_self_v2.collinearity | grep -oP '\d+' | head -1)
    log "  Ha self: ${NGENES} collinear genes"
fi

# Ha vs At (interspecies)
log "  MCScanX: Ha vs At..."
mkdir -p "${SYNTENY_DIR}/mcscanx/Ha_At_v2"
cat "${SYNTENY_DIR}/data/Ha.gff" "${SYNTENY_DIR}/data/At.gff" \
    > "${SYNTENY_DIR}/mcscanx/Ha_At_v2/Ha_At_v2.gff"
cp "${SYNTENY_DIR}/blast/Ha_At_v2.blast" "${SYNTENY_DIR}/mcscanx/Ha_At_v2/Ha_At_v2.blast"
cd "${SYNTENY_DIR}/mcscanx/Ha_At_v2"
MCScanX Ha_At_v2 -s 5 2>&1 || true
if [[ -f Ha_At_v2.collinearity ]]; then
    NGENES=$(grep "Number of collinear genes" Ha_At_v2.collinearity | grep -oP '\d+' | head -1)
    log "  Ha vs At: ${NGENES} collinear genes"
fi

# Ha vs Ls (interspecies)
log "  MCScanX: Ha vs Ls..."
mkdir -p "${SYNTENY_DIR}/mcscanx/Ha_Ls_v2"
cat "${SYNTENY_DIR}/data/Ha.gff" "${SYNTENY_DIR}/data/Ls.gff" \
    > "${SYNTENY_DIR}/mcscanx/Ha_Ls_v2/Ha_Ls_v2.gff"
cp "${SYNTENY_DIR}/blast/Ha_Ls_v2.blast" "${SYNTENY_DIR}/mcscanx/Ha_Ls_v2/Ha_Ls_v2.blast"
cd "${SYNTENY_DIR}/mcscanx/Ha_Ls_v2"
MCScanX Ha_Ls_v2 -s 5 2>&1 || true
if [[ -f Ha_Ls_v2.collinearity ]]; then
    NGENES=$(grep "Number of collinear genes" Ha_Ls_v2.collinearity | grep -oP '\d+' | head -1)
    log "  Ha vs Ls: ${NGENES} collinear genes"
fi

# =============================================================================
# STEP 5: Extract AQP syntenic pairs with species distinction
# =============================================================================
log "STEP 5: Extracting aquaporin syntenic pairs"

python3 - "${SYNTENY_DIR}" "${AQP_LIST}" <<'PY5'
import os, sys, re

synteny_dir = sys.argv[1]
aqp_file = sys.argv[2]
data_dir = os.path.join(synteny_dir, 'data')
results_dir = os.path.join(synteny_dir, 'results')
os.makedirs(results_dir, exist_ok=True)

# Load all protein -> gene mappings (all species)
all_prot_to_gene = {}
for sp in ['Ha', 'At', 'Ls']:
    map_file = os.path.join(data_dir, f'{sp}_prot_to_gene.tsv')
    if os.path.exists(map_file):
        with open(map_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    all_prot_to_gene[parts[0]] = (sp, parts[1])

# Load sunflower AQP gene list
aqp_locs = set()
with open(aqp_file) as f:
    next(f)
    for line in f:
        aqp_locs.add(line.strip().split('\t')[0])

# Build AQP protein set (sunflower only)
aqp_prots = set()
for prot, (sp, loc) in all_prot_to_gene.items():
    if sp == 'Ha' and loc in aqp_locs:
        aqp_prots.add(prot)

print(f"Loaded: {len(aqp_locs)} AQP genes, {len(aqp_prots)} AQP proteins")
print(f"Mappings: Ha={sum(1 for _,(s,_) in all_prot_to_gene.items() if s=='Ha')}, "
      f"At={sum(1 for _,(s,_) in all_prot_to_gene.items() if s=='At')}, "
      f"Ls={sum(1 for _,(s,_) in all_prot_to_gene.items() if s=='Ls')}")

for comp in ['Ha_self_v2', 'Ha_At_v2', 'Ha_Ls_v2']:
    col_file = os.path.join(synteny_dir, 'mcscanx', comp, f'{comp}.collinearity')
    if not os.path.exists(col_file):
        print(f"\n  {comp}: collinearity file not found")
        continue
    
    total_blocks = 0
    aqp_blocks = 0
    aqp_pairs_intra = []   # Ha-Ha pairs
    aqp_pairs_inter = []   # Ha-At or Ha-Ls pairs
    block_id = ''
    block_genes = []
    block_chrs = ('', '')
    
    with open(col_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('## Alignment'):
                # Process previous block
                if block_id:
                    has_aqp = any(g in aqp_prots for pair in block_genes for g in pair)
                    if has_aqp:
                        aqp_blocks += 1
                        for g1, g2 in block_genes:
                            if g1 in aqp_prots or g2 in aqp_prots:
                                sp1 = all_prot_to_gene.get(g1, ('?','?'))
                                sp2 = all_prot_to_gene.get(g2, ('?','?'))
                                loc1 = sp1[1] if sp1[0] != '?' else g1
                                loc2 = sp2[1] if sp2[0] != '?' else g2
                                
                                if sp1[0] == sp2[0]:  # same species = intra
                                    aqp_pairs_intra.append((block_id, sp1[0], loc1, sp2[0], loc2, g1, g2))
                                else:  # different species = inter
                                    aqp_pairs_inter.append((block_id, sp1[0], loc1, sp2[0], loc2, g1, g2))
                
                total_blocks += 1
                block_id = line
                block_genes = []
                
                # Extract chromosomes
                chr_m = re.findall(r'(\w+_\w+\.\d+)', line)
                block_chrs = (chr_m[0] if len(chr_m) > 0 else '', 
                              chr_m[1] if len(chr_m) > 1 else '')
                
            elif line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 3:
                    block_genes.append((parts[1], parts[2]))
    
    # Last block
    if block_id:
        has_aqp = any(g in aqp_prots for pair in block_genes for g in pair)
        if has_aqp:
            aqp_blocks += 1
            for g1, g2 in block_genes:
                if g1 in aqp_prots or g2 in aqp_prots:
                    sp1 = all_prot_to_gene.get(g1, ('?','?'))
                    sp2 = all_prot_to_gene.get(g2, ('?','?'))
                    loc1 = sp1[1] if sp1[0] != '?' else g1
                    loc2 = sp2[1] if sp2[0] != '?' else g2
                    if sp1[0] == sp2[0]:
                        aqp_pairs_intra.append((block_id, sp1[0], loc1, sp2[0], loc2, g1, g2))
                    else:
                        aqp_pairs_inter.append((block_id, sp1[0], loc1, sp2[0], loc2, g1, g2))
    
    print(f"\n  {comp}:")
    print(f"    Total blocks: {total_blocks}")
    print(f"    Blocks with AQP: {aqp_blocks}")
    print(f"    Intra-species AQP pairs: {len(aqp_pairs_intra)}")
    print(f"    Inter-species AQP pairs: {len(aqp_pairs_inter)}")
    
    # Save results
    tag = comp.replace('_v2', '')
    
    # Intra
    out_intra = os.path.join(results_dir, f'{tag}_intra_aqp_synteny.tsv')
    with open(out_intra, 'w') as f:
        f.write('block\tspecies1\tgene1\tspecies2\tgene2\tprot1\tprot2\n')
        for row in aqp_pairs_intra:
            f.write('\t'.join(str(x) for x in row) + '\n')
    
    # Inter
    out_inter = os.path.join(results_dir, f'{tag}_inter_aqp_synteny.tsv')
    with open(out_inter, 'w') as f:
        f.write('block\tspecies1\tgene1\tspecies2\tgene2\tprot1\tprot2\n')
        for row in aqp_pairs_inter:
            f.write('\t'.join(str(x) for x in row) + '\n')
    
    # Print inter-species pairs detail
    if aqp_pairs_inter:
        print(f"    Inter-species pairs:")
        for row in aqp_pairs_inter:
            print(f"      {row[1]}:{row[2]} <-> {row[3]}:{row[4]}")

PY5

log "DONE."
log "Results in: ${SYNTENY_DIR}/results/"
