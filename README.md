# I2GDS_xylella
This repo is for the I2GDS class individual project. This project runs BUSCO and OrthoFinder on Xylella protein sequences and assembles a phylogenetic tree in R.

```
module load Miniconda3
conda create -n 
busco-env -c conda-forge -c bioconda busco
```
```
#!/bin/bash
#SBATCH --job-name=busco_proteins
#SBATCH --output=logs/busco_%x_%j.out
#SBATCH --error=logs/busco_%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=normal_q
#SBATCH --account=introtogds
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== BUSCO protein mode batch started at $(date) ==="

# Load conda
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate busco-env

# Directories
ANNOT_DIR="./annotations"
OUT_BASE="./busco_proteins"

mkdir -p "$OUT_BASE" logs

# Loop over all .faa files (each inside its own directory)
for FAA in ${ANNOT_DIR}/*/*.faa; do
    SAMPLE=$(basename "${FAA}" .faa)
    OUT_DIR="${OUT_BASE}/${SAMPLE}"
    
    mkdir -p "${OUT_DIR}"

    echo "----------------------------------------------"
    echo "Running BUSCO for: ${SAMPLE}"
    echo "Input:   ${FAA}"
    echo "Output:  ${OUT_DIR}"
    echo "----------------------------------------------"

    busco \
        -i "${FAA}" \
        -l bacteria_odb10 \
        -o "${SAMPLE}" \
        -m proteins \
        --cpu "${SLURM_CPUS_PER_TASK}" \
        --out_path "${OUT_DIR}"

    echo "Finished BUSCO for ${SAMPLE}"
    echo
done

echo "=== BUSCO batch complete at $(date)! ==="
```

```
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
```


