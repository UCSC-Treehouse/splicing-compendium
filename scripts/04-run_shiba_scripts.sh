#!/bin/bash

# Create experiment.tsv file and run Shiba for 88 TARGET bulk RNA-seq samples

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# set sample group
threads=15

# Set paths as variables
git_path=$(git rev-parse --git-dir)
root_dir=$(dirname "$git_path")
reference_gtf="../references/gencode.v47.primary_assembly.annotation.gtf"
in_gtf="../results/target/shiba/merged_gtf.gtf"
log_dir="../logs"
exploration_dir="../exploration/merged_shiba_test/merged_method"
out_dir="${exploration_dir}/shiba_run"
# I needed to clone shiba v0.8.1 because I did not have permissions to run its scripts within the directory
shiba_dir="/home/ubuntu/Shiba/src"

# Define log file
log_file="${log_dir}/${current_datetime}_merge_shiba_run.txt"

# create directories if they do not already exist
mkdir -p $log_dir
mkdir -p $exploration_dir

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

# run gtf2event.py to extract splice events from gtf
# Do I add the path to gtf2event.py to an env variable?
# need to run chmod +x first on gtf2event.py
python $shiba_dir/gtf2event.py -i $in_gtf -r $reference_gtf -o OUTPUT -p $threads -v
