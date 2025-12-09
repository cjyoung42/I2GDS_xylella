# I2GDS_xylella
This repo is for the I2GDS class individual project. This project runs BUSCO and OrthoFinder on Xylella protein sequences and assembles a phylogenetic tree in R.

```
module load Miniconda3
conda create -n busco-env python=3.10 -y

source activate busco-env

conda install -c conda-forge mamba -y

mamba install -c conda-forge -c bioconda busco=5.7.1

export AUGUSTUS_CONFIG_PATH=$(cond
a env list | grep busco-env | awk '{print $2}')/config/augustus

busco --download bacteria_odb10 --offline
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

# Loop over all .faa files from Prokka
for FAA in ${ANNOT_DIR}/*/*.faa; do
    SAMPLE=$(basename "${FAA}" .faa)
    OUT_DIR="${OUT_BASE}/${SAMPLE}"
    
    echo "----------------------------------------------"
    echo "Running BUSCO for: ${SAMPLE}"
    echo "Input:   ${FAA}"
    echo "Output:  ${OUT_DIR}"
    echo "----------------------------------------------"

    busco \
        -i "${FAA}" \
        -l bacteria_odb10 \
        -o "${SAMPLE}" \
        -m protein \
        --cpu "${SLURM_CPUS_PER_TASK}" \
        --out_path "${OUT_DIR}" \
        --offline

    echo "Finished BUSCO for ${SAMPLE}"
    echo
done
```

echo "=== BUSCO batch complete at $(date)! ==="
