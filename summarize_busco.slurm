#!/bin/bash
#SBATCH --job-name=busco_summary
#SBATCH --output=logs/busco_summary_%j.out
#SBATCH --error=logs/busco_summary_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=normal_q
#SBATCH --account=introtogds

set -euo pipefail

echo "=== BUSCO summary started at $(date) ==="

BUSCO_DIR="./busco_proteins"
OUTFILE="busco_summary.csv"

mkdir -p logs

echo "sample,complete_percent,single_copy,duplicated,fragmented,missing,total_buscos" > "$OUTFILE"

# Find BUSCO short summary files (v5 format)
SUMMARIES=$(find "$BUSCO_DIR" -mindepth 2 -maxdepth 3 -name "short_summary*.txt")

if [[ -z "$SUMMARIES" ]]; then
    echo "ERROR: No BUSCO summary files found."
    exit 1
fi

for SUMMARY in $SUMMARIES; do
    # Extract sample name from first directory
    SAMPLE=$(basename "$(dirname "$SUMMARY")")

    echo "Processing sample: $SAMPLE"

    # Extract values using flexible regex
    C=$(grep -oP 'C:\s*\K[\d.]+' "$SUMMARY" || echo "NA")
    S=$(grep -oP 'S:\s*\K[\d.]+' "$SUMMARY" || echo "NA")
    D=$(grep -oP 'D:\s*\K[\d.]+' "$SUMMARY" || echo "NA")
    F=$(grep -oP 'F:\s*\K[\d.]+' "$SUMMARY" || echo "NA")
    M=$(grep -oP 'M:\s*\K[\d.]+' "$SUMMARY" || echo "NA")
    N=$(grep -oP 'n:\s*\K[\d.]+' "$SUMMARY" || echo "NA")

    echo "${SAMPLE},${C},${S},${D},${F},${M},${N}" >> "$OUTFILE"
done

echo "=== BUSCO summary complete ==="
echo "CSV written to: ${OUTFILE}"
