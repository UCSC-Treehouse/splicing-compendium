#!/bin/bash

# When OpenStack instances become full before data processing is done, we need to offload raw sequence files to clear space.
# This script makes a list of fastq files to remove from the accession IDs of sligned files present in a directory
# Then it removes these files

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# set working directory to be the directory of this file (scripts/)
cd "$(dirname "${BASH_SOURCE[0]}")"

# set time and date as variables
datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# set paths as variables
# git_path=$(git rev-parse --git-dir)
# root_dir=$(dirname "$git_path")
data_dir="../data"
group_dir="${data_dir}/target"
bam_dir="${group_dir}/star-output"
fastq_dir="${group_dir}/fastq"
log_dir="../logs"

# create directories if they don't already exist
mkdir -p $log_dir

# define output files
fastq_to_offload="${log_dir}/fastq_to_offload_${datetime}.txt"

# create list of subdirectories holding bam files from workflow
# use mapfile to read lines from find to put directory names into an array
mapfile -t dirs < <(find "${bam_dir}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')

# iterate through array to remove the fastq files corresponding to the accession ID
for accession in "${dirs[@]}"; do
    # construct pats to fastqs based on accessions in the bam dir
    fastq1_path="${fastq_dir}/${accession}_1.fastq.gz"
    fastq2_path="${fastq_dir}/${accession}_2.fastq.gz"

    # remove fastq files if they exist
    if [ -f "${fastq1_path}" ]; then
        rm ${fastq1_path}
    fi

        if [ -f "${fastq2_path}" ]; then
        rm ${fastq2_path}
    fi
done
