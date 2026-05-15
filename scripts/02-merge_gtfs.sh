#!/bin/bash

# Create experiment.tsv file and run Shiba for 88 TARGET bulk RNA-seq samples

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# set sample group
group="target"
threads=15

# Set paths as variables
reference_gtf="../references/gencode.v47.primary_assembly.annotation.gtf"
sample_list="../config/two_sample_merge_samples.tsv"
log_dir="../logs"
shiba_dir="../results/${group}/shiba"
gtf_list="../config/two_sample_gtf_list.txt"

# output gtf
out_gtf="${shiba_dir}/merged_gtf.gtf"

# Define log file
log_file="${log_dir}/${current_datetime}_gtf_list.txt"

# create directories if they do not already exist
mkdir -p $log_dir

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

# instantiate gtf list file
touch ${gtf_list}

# generate gtf list of files to merge from samples in gtf list
while read sample; do
    # check if file is not already unzipped
    if [[ ! -f "${shiba_dir}/${sample}/annotation/assembled_annotation.gtf" ]]; then
        # unzip the gtf for stringtie merge
        gunzip ${shiba_dir}/${sample}/annotation/assembled_annotation.gtf.gz
    fi

    # add gtf path to gtf list for stringtie merge
    echo "${shiba_dir}/${sample}/annotation/assembled_annotation.gtf" >> ${gtf_list}
done < ${sample_list}

# run stringtie merge on gtfs
stringtie --merge -p ${threads} -G ${reference_gtf} -o ${out_gtf} ${gtf_list}
