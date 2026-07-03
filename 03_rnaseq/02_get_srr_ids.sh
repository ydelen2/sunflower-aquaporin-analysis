#!/bin/bash
#SBATCH --job-name=get_srr_ids
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/02_get_srr_ids_%j.out
#SBATCH --error=logs/02_get_srr_ids_%j.err

# =============================================================================
# 02_get_srr_ids.sh
# Fetch SRR accession IDs and BioSample metadata for all BioProjects
# Uses NCBI Entrez Direct (esearch/efetch/elink)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

module purge
module load "$MOD_ENTREZ"

init_dirs

log_msg "Starting SRR accession retrieval for ${#BIOPROJECTS[@]} BioProjects"

# ---------------------------------------------------------------------------
# Function: fetch SRR accessions for a single BioProject
# ---------------------------------------------------------------------------
fetch_srr_ids() {
    local bioproject="$1"
    local outfile="${ACCESSION_DIR}/${bioproject}_srr.txt"
    local metafile="${ACCESSION_DIR}/${bioproject}_metadata.tsv"

    log_msg "Fetching SRR IDs for ${bioproject}..."

    # Query SRA database for all runs in this BioProject
    esearch -db sra -query "${bioproject}[BioProject]" \
        | efetch -format runinfo \
        | grep -v "^Run" \
        | cut -d',' -f1 \
        | sort -u \
        | grep -E "^[SED]RR" \
        > "$outfile" || true

    local count=$(wc -l < "$outfile" 2>/dev/null || echo 0)

    if [[ "$count" -eq 0 ]]; then
        log_msg "WARNING: No SRR IDs found for ${bioproject}"
        return 1
    fi

    log_msg "  Found ${count} runs for ${bioproject}"

    # Fetch detailed metadata (runinfo CSV -> TSV with selected columns)
    log_msg "  Fetching sample metadata for ${bioproject}..."

    esearch -db sra -query "${bioproject}[BioProject]" \
        | efetch -format runinfo \
        > "${ACCESSION_DIR}/${bioproject}_runinfo.csv" || true

    # Extract key columns: Run, Sample, BioSample, Experiment, LibraryLayout,
    # spots, bases, avgLength, SampleName
    if [[ -s "${ACCESSION_DIR}/${bioproject}_runinfo.csv" ]]; then
        head -1 "${ACCESSION_DIR}/${bioproject}_runinfo.csv" \
            | tr ',' '\t' > "$metafile"

        tail -n +2 "${ACCESSION_DIR}/${bioproject}_runinfo.csv" \
            | grep -v "^$" \
            | tr ',' '\t' >> "$metafile"

        log_msg "  Metadata saved to ${metafile}"
    fi

    # Fetch BioSample attributes (treatment, tissue, genotype)
    log_msg "  Fetching BioSample attributes for ${bioproject}..."
    local biosample_file="${ACCESSION_DIR}/${bioproject}_biosample_attrs.tsv"

    # Get unique BioSample IDs from runinfo
    local biosamples=$(tail -n +2 "${ACCESSION_DIR}/${bioproject}_runinfo.csv" \
        | cut -d',' -f26 \
        | sort -u \
        | grep -v "^$" || true)

    if [[ -n "$biosamples" ]]; then
        echo -e "BioSample\tAttribute\tValue" > "$biosample_file"

        for bs in $biosamples; do
            # Rate limit to avoid NCBI throttling
            sleep 0.4

            efetch -db biosample -id "$bs" -format docsum \
                | xtract -pattern DocumentSummary \
                    -block Attribute \
                    -sep "\t" \
                    -element "&ACCN" "@attribute_name" "Attribute" \
                2>/dev/null >> "$biosample_file" || true
        done

        log_msg "  BioSample attributes saved to ${biosample_file}"
    fi

    return 0
}

# ---------------------------------------------------------------------------
# Main: iterate over all BioProjects
# ---------------------------------------------------------------------------
SUCCESS=0
FAIL=0

for entry in "${BIOPROJECTS[@]}"; do
    parse_bioproject "$entry"
    log_msg "Processing: ${BP_ID} (${BP_NAME} - ${BP_DESC})"

    if fetch_srr_ids "$BP_ID"; then
        ((SUCCESS++))
    else
        ((FAIL++))
        log_msg "FAILED: ${BP_ID}"
    fi

    # Rate limit between BioProjects
    sleep 2
done

log_msg "Completed: ${SUCCESS} succeeded, ${FAIL} failed out of ${#BIOPROJECTS[@]} BioProjects"

# Build master sample list
build_master_sample_list

# Summary
log_msg "=== Accession Summary ==="
for entry in "${BIOPROJECTS[@]}"; do
    parse_bioproject "$entry"
    local acc_file="${ACCESSION_DIR}/${BP_ID}_srr.txt"
    if [[ -f "$acc_file" ]]; then
        local n=$(wc -l < "$acc_file")
        log_msg "  ${BP_ID} (${BP_NAME}): ${n} samples"
    else
        log_msg "  ${BP_ID} (${BP_NAME}): MISSING"
    fi
done

log_msg "Done."
