#!/bin/bash
#SBATCH --job-name=align_stats
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/07_alignment_stats_%j.out
#SBATCH --error=logs/07_alignment_stats_%j.err

# =============================================================================
# 07_alignment_stats.sh
# Parse HISAT2 alignment summaries and featureCounts summary
# Generate a consolidated table and flag low-quality samples
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

init_dirs

SUMMARY_FILE="${STATS_DIR}/pipeline_summary.tsv"
FLAGGED_FILE="${STATS_DIR}/flagged_samples.tsv"
FC_SUMMARY="${COUNTS_DIR}/gene_counts_raw.txt.summary"

log_msg "Collecting alignment and counting statistics..."

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
echo -e "sample\tbioproject\ttotal_reads\tmapped_reads\tmapping_rate\tuniquely_mapped\tassigned_reads\tassigned_pct\tflag" \
    > "$SUMMARY_FILE"

echo -e "sample\tbioproject\tissue\treason" \
    > "$FLAGGED_FILE"

# ---------------------------------------------------------------------------
# Parse featureCounts summary into associative array
# ---------------------------------------------------------------------------
declare -A FC_ASSIGNED

if [[ -f "$FC_SUMMARY" ]]; then
    log_msg "Parsing featureCounts summary..."

    # featureCounts summary has columns: Status, sample1, sample2, ...
    # We need the "Assigned" row
    python3 -c "
import sys, re

with open('${FC_SUMMARY}') as f:
    lines = f.readlines()

# Header line has BAM paths
header = lines[0].strip().split('\t')[1:]
samples = []
for h in header:
    m = re.search(r'([SED]RR\d+)\.sorted\.bam', h)
    samples.append(m.group(1) if m else h)

# Find Assigned row
for line in lines[1:]:
    parts = line.strip().split('\t')
    if parts[0] == 'Assigned':
        for s, v in zip(samples, parts[1:]):
            print(f'{s}\t{v}')
        break
" > "${STATS_DIR}/_fc_assigned.tmp"

    while IFS=$'\t' read -r srr count; do
        FC_ASSIGNED["$srr"]="$count"
    done < "${STATS_DIR}/_fc_assigned.tmp"

    rm -f "${STATS_DIR}/_fc_assigned.tmp"
fi

# ---------------------------------------------------------------------------
# Iterate over all samples, parse HISAT2 logs
# ---------------------------------------------------------------------------
TOTAL_PROCESSED=0
TOTAL_FLAGGED=0

while IFS= read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue

    BP_ID=$(find_bioproject_for_srr "$SAMPLE")
    ALIGN_LOG="${BAM_DIR}/${BP_ID}/${SAMPLE}.hisat2.log"

    # Defaults
    TOTAL_READS=0
    MAPPED=0
    RATE=0
    UNIQ=0
    ASSIGNED=${FC_ASSIGNED[$SAMPLE]:-0}
    FLAGS=""

    if [[ -f "$ALIGN_LOG" ]]; then
        # Parse HISAT2 new-summary format
        eval $(python3 -c "
import re
with open('${ALIGN_LOG}') as f:
    text = f.read()

total_m = re.search(r'Total reads:\s+(\d+)', text) or re.search(r'(\d+) reads', text)
total = int(total_m.group(1)) if total_m else 0

rate_m = re.search(r'Overall alignment rate:\s+([\d.]+)%', text) or re.search(r'([\d.]+)% overall', text)
rate = float(rate_m.group(1)) if rate_m else 0.0

uniq_m = re.search(r'Aligned 1 time:\s+(\d+)', text) or re.search(r'(\d+).*aligned exactly 1 time', text)
uniq = int(uniq_m.group(1)) if uniq_m else 0

mapped = int(total * rate / 100)

print(f'TOTAL_READS={total}')
print(f'MAPPED={mapped}')
print(f'RATE={rate}')
print(f'UNIQ={uniq}')
")
    else
        FLAGS="missing_log"
    fi

    # Compute assigned percentage
    ASSIGNED_PCT=0
    if [[ "$TOTAL_READS" -gt 0 && "$ASSIGNED" -gt 0 ]]; then
        ASSIGNED_PCT=$(python3 -c "print(round(100.0 * ${ASSIGNED} / ${TOTAL_READS}, 2))")
    fi

    # Flag checks
    if (( $(echo "$RATE < $MIN_MAPPING_RATE" | bc -l) )); then
        FLAGS="${FLAGS:+${FLAGS},}low_mapping_rate"
    fi
    if [[ "$TOTAL_READS" -lt "$MIN_TOTAL_READS" ]]; then
        FLAGS="${FLAGS:+${FLAGS},}low_read_count"
    fi

    # Write to summary
    echo -e "${SAMPLE}\t${BP_ID}\t${TOTAL_READS}\t${MAPPED}\t${RATE}\t${UNIQ}\t${ASSIGNED}\t${ASSIGNED_PCT}\t${FLAGS:-PASS}" \
        >> "$SUMMARY_FILE"

    # Write flagged samples
    if [[ -n "$FLAGS" && "$FLAGS" != "PASS" ]]; then
        echo -e "${SAMPLE}\t${BP_ID}\t-\t${FLAGS}" >> "$FLAGGED_FILE"
        ((TOTAL_FLAGGED++))
    fi

    ((TOTAL_PROCESSED++))

done < "${RNASEQ_DIR}/all_samples.txt"

# ---------------------------------------------------------------------------
# Print summary to stdout
# ---------------------------------------------------------------------------
log_msg "=== Pipeline Statistics Summary ==="
log_msg "Total samples processed: ${TOTAL_PROCESSED}"
log_msg "Flagged samples: ${TOTAL_FLAGGED}"
log_msg ""
log_msg "Summary table: ${SUMMARY_FILE}"
log_msg "Flagged samples: ${FLAGGED_FILE}"

if [[ "$TOTAL_FLAGGED" -gt 0 ]]; then
    log_msg ""
    log_msg "=== Flagged Samples ==="
    column -t -s $'\t' "$FLAGGED_FILE"
fi

# ---------------------------------------------------------------------------
# Per-BioProject summary
# ---------------------------------------------------------------------------
log_msg ""
log_msg "=== Per-BioProject Summary ==="

for entry in "${BIOPROJECTS[@]}"; do
    parse_bioproject "$entry"

    N_TOTAL=$(grep -c "$BP_ID" "$SUMMARY_FILE" 2>/dev/null || echo 0)
    N_FLAGGED=$(grep "$BP_ID" "$FLAGGED_FILE" 2>/dev/null | wc -l || echo 0)
    N_PASS=$((N_TOTAL - N_FLAGGED))

    # Average mapping rate
    AVG_RATE=$(awk -F'\t' -v bp="$BP_ID" '$2 == bp {sum+=$5; n++} END {if(n>0) printf "%.1f", sum/n; else print "NA"}' "$SUMMARY_FILE")

    log_msg "  ${BP_ID} (${BP_NAME}): ${N_TOTAL} samples, ${N_PASS} pass, ${N_FLAGGED} flagged, avg mapping rate: ${AVG_RATE}%"
done

log_msg "Done."
