#!/bin/bash

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set paths as variables
fastq_dir="/mnt/data/gtex-muscle-fastq/"
out_dir="../metadata"
output="${out_dir}/downloaded_gtex_muscle_ids.txt"

# List files, extract IDs (everything before first underscore), sort and keep unique IDs
for f in ${fastq_dir}/*.fastq.gz; do
# use basename to extract file name from the path and sed to keep everything before the first underscore (the accession ID)
    basename "$f" | sed -E 's/_.*//'
done | sort -u > "$output"
