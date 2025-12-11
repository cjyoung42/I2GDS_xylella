# I2GDS_xylella
This repo is for the I2GDS class individual project. This project runs BUSCO and OrthoFinder on Xylella protein sequences and assembles a phylogenetic tree in R.

```
module load Miniconda3
conda create -n 
busco-env -c conda-forge -c bioconda busco
source activate busco-env
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

```
#!/bin/bash
set -euo pipefail

# Where your annotation folders are (Prokka, NCBI, etc.)
SOURCE_DIR="./annotations"   # change if needed

# Where you want to put all FAA files for OrthoFinder
TARGET_DIR="./orthofinder_input"

mkdir -p "$TARGET_DIR"

echo "Collecting .faa files from: $SOURCE_DIR"
echo "Copying to: $TARGET_DIR"
echo

find "$SOURCE_DIR" -type f -name "*.faa" | while read FAA; do
    BASENAME=$(basename "$FAA")

    # If two files have the same name, prefix with parent folder
    PARENT=$(basename "$(dirname "$FAA")")
    NEWNAME="${PARENT}_${BASENAME}"

    echo "Copying:"
    echo "  $FAA"
    echo "  → $TARGET_DIR/$NEWNAME"
    echo

    cp "$FAA" "$TARGET_DIR/$NEWNAME"
done

echo "Done! All FASTA files are now in:"
echo "  $TARGET_DIR"
```
```
bash collect_faas.sh
```

```
#!/bin/bash
set -euo pipefail

# --- Configuration ---
IN_DIR="./orthofinder_input"           # original files
OUT_DIR="./orthofinder_input_clean2"   # cleaned files output
PREVIEW_LINES=6                        # how many lines to show from one cleaned file for quick check
# ----------------------

mkdir -p "$OUT_DIR"

# make sure patterns that don't match expand to nothing
shopt -s nullglob

echo "Input dir : $IN_DIR"
echo "Output dir: $OUT_DIR"
echo

found=0
for FAA in "$IN_DIR"/*.faa "$IN_DIR"/*.fa; do
    # if pattern didn't match, loop body is skipped (thanks to nullglob)
    found=1

    BASENAME=$(basename "$FAA")
    # remove file extensions and _prokka substring from base name
    NAME_NOEXT="${BASENAME%.*}"
    CLEAN_SAMPLE=$(echo "$NAME_NOEXT" | sed 's/_prokka//g')

    OUTFILE="$OUT_DIR/${CLEAN_SAMPLE}.faa"

    echo "Processing: $BASENAME"
    echo " → Writing: $OUTFILE"

    # AWK logic:
    # - For each header, try to extract the last "_" token if it's numeric (common pattern)
    # - If numeric token not found, use a sequential counter (padded) to guarantee uniqueness
    # - Strip any description after first whitespace
    awk -v SAMPLE="$CLEAN_SAMPLE" '
        BEGIN { seq=0 }
        /^>/ {
            # get first whitespace-delimited token of header (removes descriptions)
            header = $1
            # remove leading ">"
            if (substr(header,1,1) == ">") header = substr(header,2)

            # remove SAMPLE or SAMPLE_prokka prefix if present
            gsub("^"SAMPLE"_?prokka_?", "", header)   # remove SAMPLE_prokka_ or SAMPLE_
            gsub("^"SAMPLE"_?", "", header)

            # take last underscore-separated token
            n = split(header, parts, "_")
            candidate = parts[n]

            # if candidate is all digits, use it; else increment seq and use seq
            if (candidate ~ /^[0-9]+$/) {
                id = candidate
            } else {
                seq++
                id = sprintf("%06d", seq)
            }

            print ">" SAMPLE "|" id
            next
        }
        { print }
    ' "$FAA" > "$OUTFILE"

    # make sure output ended with a newline
    printf "\n" >> "$OUTFILE"

    echo
done

if [[ $found -eq 0 ]]; then
    echo "No .faa or .fa files found in $IN_DIR — nothing done."
    exit 1
fi

echo "Cleaning finished. Example preview of first cleaned file:"
first=$(ls -1 "$OUT_DIR"/*.faa | head -n1)
if [[ -n "$first" ]]; then
    echo "File: $first"
    echo "------"
    head -n "$PREVIEW_LINES" "$first"
else
    echo "No files found in $OUT_DIR (unexpected)."
fi

echo
echo "You can now run OrthoFinder on: $OUT_DIR"
```
```
bash clean_headers.sh
```

```
conda remove -n orthofinder-env --all -y

conda create -n orthofinder-env -y \
    -c bioconda -c conda-forge \
    orthofinder scikit-learn diamond mafft blast fasttree

source activate orthofinder-env
```

```
```
