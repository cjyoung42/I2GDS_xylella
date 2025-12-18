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
