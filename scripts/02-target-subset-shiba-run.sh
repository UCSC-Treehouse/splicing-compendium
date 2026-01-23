#!/bin/bash

# Create experiment.tsv file and run Shiba for 88 TARGET bulk RNA-seq samples

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# Set paths as variables
sample_list="../shiba-run-data/target/fastq/target_md5sum.txt"
config_file="target_subset_shiba_config.yaml"
experiment_file="target_subset_pilot_shiba_experiment.tsv"
log_dir="../logs"
splicing_results_dir="../results/target-subset-shiba-output/results/splicing"

# set experiment variables
group="target"
threads=15 # number of threads on OpenStack

# Define log files
log_file="${log_dir}/${current_datetime}_TARGET_subset_shiba_run.txt"

# create directories if they do not already exist
mkdir -p $log_dir

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

# Create experiment.tsv for Shiba run
Rscript generate_experiment_file.R --input=$sample_list --output=$experiment_file --group=$group

# Run Shiba
time shiba.py -p $threads $config_file

# zip splicing results to save space
pigz "${splicing_results_dir}"/*.txt
