#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH --output=logs/orthofinder_%j.out
#SBATCH --error=logs/orthofinder_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --partition=normal_q
#SBATCH --account=introtogds
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Starting OrthoFinder at $(date) ==="

# Load conda
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate orthofinder-env

# Input directory with cleaned FASTA files
INPUT_DIR="./orthofinder_input_clean"

# Output directory
OUTPUT_DIR="./orthofinder_results"
mkdir -p logs

echo "Input directory:  $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo

# --- Run OrthoFinder ---
orthofinder \
    -f "$INPUT_DIR" \
    -o "$OUTPUT_DIR" \
    -t "$SLURM_CPUS_PER_TASK" \
    -a "$SLURM_CPUS_PER_TASK" \
    -M msa \
    -S diamond_ultra_sens \
    -T fasttree

echo
echo "=== OrthoFinder completed at $(date) ==="
