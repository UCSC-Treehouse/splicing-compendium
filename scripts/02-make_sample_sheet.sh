#!/bin/bash

# a sample sheet containing a "group" (e.g. target, gtex) and a "sample" (accession ID) column is needed to run the snakemake workflow
# this scripts generates this sample sheet from files that are present in the results/ directory

# cause nonzero exit status and undefined variables to stop script
set -euo pipefail

# set working directory to directory location of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# array of sample groups
groups=("target" "gtex")

# define data directories
results_dir="../results"
config_dir="../config"

# define output file
sample_sheet="${config_dir}/sample_sheet.tsv"

# create column headers in blank sample sheet
# -e allows for special characters like \t
echo -e "samples\tgroup" > $sample_sheet

# sample results are organized for each group in the following structure:
# results/target|gtex/shiba/accession
# loop over sample results directories for each group
# and append the group + sample accession to the sample sheet
for group in "${groups[@]}"; do
    # find all paths to accession directories in results/
    # printf '%f\n' prints the basename of the file paths, followed by a new line
    find "${results_dir}/${group}/shiba" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
    | while read accession; do # loop through accession basenames and add them + group to sample sheet
        echo -e "${accession}\t${group}" >> $sample_sheet
    done
done
