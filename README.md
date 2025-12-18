# I2GDS Xylella Ortholog Analysis 
This repo is for the I2GDS class individual project. This project runs BUSCO and OrthoFinder on _Xylella_ protein sequences and assembles a phylogenetic tree using the R package ggtree. The input data for this pipeline is annotated sequences produced by Prokka, such as the outputs from the Linux portion of the Group 2 project; for previous steps in this pipeline, consult https://github.com/LiLabAtVT/I2GDS2025/tree/main/Group2. Input data for this project can be accessed on the arc at /projects/intro2gds/I2GDS2025/G2_PlantDisease/Carter/annotations. A relatively small dataset of 47 _Xylella_ proteomes was selected to run the following steps quickly. A known _Xylella_ species was selected to be an "unknown" to be identified using phylogenetic inference as a practice task. 

#### Linux workflow
##### Step 1: Quality Control.
Check the quality and completeness of Prokka proteome outputs (.faa files) using **BUSCO**
##### Step 2: Orthology Inference.
Identify orthologroups and generate a species tree using **OrthoFinder**
#### R workflow
##### Step 3: Generate Tree
Generate a phylogenetic tree figure with **ggtree** based on OrthoFinder's species tree. Manually root tree and edit aesthetics to make it visually pleasing.
Below is the final result of this pipeline.

<img width="445" height="669" alt="Screenshot 2025-12-16 at 10 34 22 PM" src="https://github.com/user-attachments/assets/718de075-eadb-4050-9755-2a90de5f2132" />

Scripts and image files are available in the files section of this repo.

## Step 1. Quality Control in BUSCO
**1.1**
First, an environment must be created to install the package BUSCO and all its dependencies. 
```
module load Miniconda3
conda create -n 
busco-env -c conda-forge -c bioconda busco
source activate busco-env
```
**1.2**
The following slurm batch script checks the annotated Prokka .faa files against the bacterial BUSCO lineage for proteome completeness.

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

#Load conda. This script won't run unless conda is specifically used.
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate busco-env

#Set directories, input annotations and output QC.
ANNOT_DIR="./annotations" #CHANGE if necessary
OUT_BASE="./busco_proteins"

mkdir -p "$OUT_BASE" logs

#Loop over all .faa files (each inside its own directory).
for FAA in ${ANNOT_DIR}/*/*.faa; do
    SAMPLE=$(basename "${FAA}" .faa)
    OUT_DIR="${OUT_BASE}/${SAMPLE}"
    
    mkdir -p "${OUT_DIR}"

    echo "----------------------------------------------"
    echo "Running BUSCO for: ${SAMPLE}"
    echo "Input:   ${FAA}"
    echo "Output:  ${OUT_DIR}"
    echo "----------------------------------------------"

#BUSCO should still run if BUSCO dataset bacteria_odb10 has not been manually downloaded.
#If this part gives an error, try busco --download bacteria_odb10 to download that dataset manually.

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


**1.3** Next, for ease of access, a csv file of BUSCO quality checks is made. This step is important if these aren't clean genomes/proteomes. Since they have already been filtered to some degree by Prokka, the results should all be good, which they are (completeness is nearly all 100%).
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

# Find BUSCO short summary files
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
## Step 2. OrthoFinder orthology inference
**2.1** To prepare .faa files for OrthoFinder analysis, move them to a single folder using the following script **collect_faas.sh**. This script creates a new directory for OrthoFinder inputs and moves the .faa files there.
```
#!/bin/bash
set -euo pipefail

# Where your annotation folders are (Prokka, NCBI, etc.)
SOURCE_DIR="./annotations"   #CHANGE if necessary

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
**2.2** Data cleaning is necessary to make the file names and contents more readable. Since the files contain unecessary _prokka endings in their names, and the headers of the protein sequences within the files also contain names that are too long, the following script renames things so that they are less confusing.
```
#!/bin/bash
set -euo pipefail

# --- Configuration ---
IN_DIR="./orthofinder_input"           # original files
OUT_DIR="./orthofinder_input_clean"    # cleaned files output
PREVIEW_LINES=6                        # how many lines to show from one cleaned file during quick check
# ----------------------

mkdir -p "$OUT_DIR"

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
**2.3** Install OrthoFinder and its dependencies using an environment. **Note**: Some packages, like OrthoFinder, need to have all their dependencies listed out during install on the ARC. Trying to create an OrthoFinder environment without listing all dependencies didn't work and gave errors. 
```
conda create -n orthofinder-env -y \
    -c bioconda -c conda-forge \
    orthofinder scikit-learn diamond mafft blast fasttree

source activate orthofinder-env
```
**2.4** Run OrthoFinder in a slurm batch script.
```
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
```
## Step 3. Figure generation in R. 
**3.1** First, load data data and create some preliminary trees to check if they are structured properly. Due to long branch attraction, the outgroup _Xylella taiwanensis_ places in between _X. fastidiosa fastidiosa_ and _X. f. muliplex_ clades. Therefore, the tree must be rooted manually on _Xylella taiwanensis_.
```
library(ape)
library(ggtree)
library(ggplot2)

results_dir <- "YOUR_PATH_HERE/orthofinder_results/Results_Example" #CHANGE to your path
tree_file <- file.path(results_dir, "Species_Tree", "SpeciesTree_rooted.txt")
t <- read.tree(tree_file)

#re-root tree on X. taiwanensis, the outgroup
t_rooted <- root(
  t,
  outgroup = "GCF_013177435.1_GCF_013177435.1",
  resolve.root = TRUE)

```
**3.2** Although the tree looks alright, it is difficult to interpret with only the accession codes as names for each sample. The labels on the plot must be renamed so that they display the correct species/subspecies and also the host, in order to infer lineages. A metadata file, **xylella_ids**, is needed to rename the labels.
```
library(ape)
library(ggtree)
library(dplyr)
library(stringr)

# set paths
tree_file <- "YOUR_PATH_HERE/orthofinder_results/Species_Tree/SpeciesTree_rooted.txt"
meta_file <- "YOUR_PATH_HERE/xylella_ids.csv"

# read tree
t <- read.tree(tree_file)

# read metadata
meta <- read.csv(meta_file, stringsAsFactors = FALSE)


strip_duplicate_id <- function(x) {
  n <- nchar(x)
  substr(x, 1, n / 2)
}
tip_df <- data.frame(label = t$tip.label) %>%
  mutate(sample_id = strip_duplicate_id(label))

#combine metadata to give samples unique labels
tip_df <- tip_df %>%
  left_join(meta, by = "sample_id") %>%
  mutate(
    pretty_label = paste(sample_id, subspecies, host, sep = "_")
  )

#visualize the default, uncolored tree
ggtree(t_rooted) %<+% tip_df +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.4))) +
  geom_tiplab(aes(label = pretty_label), size = 2.5) +
  theme_tree2()
```
**3.3** Next, I added colors to the tree to make it easier to read. A mutate function is necessary to standardize the labels so that colors can be assigned; previously, I wanted labels to be either lowercase or uppercase (to display the "UNKNOWN" sample) but this creates problems for adding colors. **Note**: ggplots/trees are very complicated to edit in R without breaking them. Oftentimes, adding a single feature or changing one color leads to multiple other issues that need correction. Adding a colored rectangle behind one label in order to highlight it is particularly problematic, and led to many entirely unreadable plots. Therefore, make sure to build the tree as an object first, and then add the colors and other aesthetics in a separate line of code.

```
tip_df <- tip_df %>%
  mutate(
    subspecies = toupper(trimws(as.character(subspecies))),
    subspecies = ifelse(is.na(subspecies), "UNKNOWN", subspecies)
  )

#build tree object first
p <- ggtree(t_rooted) %<+% tip_df

#add aesthetics next
print(
  p +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.4))) +

    geom_tiplab(
      data = subset(p$data, isTip & subspecies == "UNKNOWN"),
      aes(label = pretty_label),
      geom = "label",
      fill = "yellow",
      color = "black",
      label.size = 0,
      size = 2.5
    ) +

    geom_tiplab(
      aes(label = pretty_label, color = subspecies),
      size = 2.5
    ) +

    scale_color_manual(values = c(
      "TAIWANENSIS" = "forestgreen",
      "FASTIDIOSA" = "steelblue",
      "MULTIPLEX"  = "purple",
      "UNKNOWN"    = "black"
    )) +

    theme_tree2() +
    theme(legend.position = "none")
)
```

**3.4** Finally, save the plot as a .pdf using ggsave.
```
ggsave(
  filename = "xylella_tree.pdf",
  plot = last_plot(),
  width = 8,
  height = 12,
  units = "in",
  device = cairo_pdf
)
```
<img width="445" height="669" alt="Screenshot 2025-12-16 at 10 34 22 PM" src="https://github.com/user-attachments/assets/718de075-eadb-4050-9755-2a90de5f2132" />

## Results
This figure is a rooted phylogenetic tree showing relationships among sampled _Xylella_ spp. isolates, with tip labels colored by species/subspecies designation. Tips are labeled with sample identifiers and colored according to subspecies assignment: _X. taiwanensis_ (green), _X. f. fastidiosa_ (blue), and _X. f. multiplex_ (purple), while a sample lacking a confident subspecies designation is shown in black. Tips corresponding to unknown subspecies are additionally highlighted with a yellow background to emphasize their placement within the tree. The topology illustrates the phylogenetic structure of the dataset and the degree to which subspecies assignments correspond to distinct clades and allows identification of the unknown strain as _Xylella fastidiosa fastidiosa_. Additionally, host-specific clades within _fastidiosa_ form, with grape-infecting strains clustering apart from almond-infecting strains.
