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

# instantiate sample sheet
touch $sample_sheet

# add column headers to sample sheet
# -e allows for special characters like \t
echo -e "group\tsamples" > $sample_sheet

# sample results are organized for each group in the following structure:
# results/target|gtex/shiba/accession
# loop over sample results directories for each group
# and append the grooup + sample accession to the sample sheet
for group in "${groups[@]}"; do
    # find all paths to accession directories in results/
    find "${results_dir}/${group}/shiba" -mindepth 1 -maxdepth 1 -type d \
    # execute basename on each accession path to obtain accession id
    # the escaped \ (\;) marks end of exec command
    exec basename {} \; \
    | while read accession; do # loop through accession basenames and add them + group to sample sheet
        echo -e "${group}\t${accession}" >> $sample_sheet
    done
done
