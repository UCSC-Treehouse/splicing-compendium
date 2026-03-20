#!/bin/bash

# Download script for TARGET bulk RNA-seq for pediatric cancer splicing analysis
# this script downloads TARGET cancer datasets obtained from filter_target_sra.qmd, which includes all TARGET paired-end RNA-seq samples, excluding ssRNA-seq samples, cell lines and xenografts, and samples known to be ribo-deplete

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# dbGaP key files
target_key="/home/ubuntu/prj_11732_D43588.ngc"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# Set paths as variables
storage_dir="/mnt"
git_path=$(git rev-parse --git-dir)
root_dir=$(dirname "$git_path")
data_dir="${root_dir}/data"
bulk_dir="${storage_dir}/bulk"
group_dir="${bulk_dir}/$1"
fastq_dir="${group_dir}/fastq"
log_dir="${data_dir}/logs"
metadata_dir="${data_dir}/metadata"
target_gtex_dir="${metadata_dir}/filter_target_gtex"

# check that sample group inputted is valid
# set dbgap key based on whether gtex or target files are to be downloaded
if [ $1 == "target" ]; then
    dbgap_key=$target_key
elif [ $1 == "gtex" ]; then
    dbgap_key=$gtex_key
else
    echo "please use valid sample group"
    # cause script to fail due to error
    exit 1
fi

# create directories if they do not already exist
mkdir -p $log_dir
mkdir -p $bulk_dir
mkdir -p $group_dir
mkdir -p $fastq_dir

# Define log files
log_file="${log_dir}/${current_datetime}_${1}_download.txt"
error_log_file="${log_dir}/${current_datetime}_${1}_download_errors.txt"

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

# Run filtering script to obtain accession IDs to download
quarto render filter_target_sra.qmd

# Use accession IDs in first column of accessions file to download fastqs for analysis
# iterate through first column of TARGET accession metadata, which has accession IDs
# skip the first line, which consists of column headers
time for ID in $(awk 'NR>1{print $1}' "${target_gtex_dir}/$1_accessions.tsv"); do
    # check that fastq or zipped fastq does not already exist
    if [ ! -f "${fastq_dir}/${ID}_1.fastq" ] && [ ! -f "${fastq_dir}/${ID}_1.fastq.gz" ]; then
        # prefetch file dependencies
        # set max size of prefetch file to 50GB based on max file size of samples in accessions lists
        # continue even if an error is thrown
        prefetch --ngc $dbgap_key --max-size 50000000 --output-directory $fastq_dir $ID || true
        # Check if prefetch file directory has been created
        # Lack of directory means prefetch has failed; skip to next accession
        if [ ! -d "${fastq_dir}/$ID" ]; then
            # record the accession ID not downloaded
            echo -e "$ID\tprefetch" >> $error_log_file
            continue
        fi

        # download fastq file for each accession
        fasterq-dump --temp $group_dir --ngc $dbgap_key $ID --threads 4 --outdir $fastq_dir || true

        # Check if fastq file exists
        # If file is missing, record error in fastq-dump step and move on to next accession
        if [ ! -f "${fastq_dir}"/${ID}_1.fastq || ! -f "${fastq_dir}"/${ID}_2.fastq]; then
            # record the accession ID not downloaded
            echo -e "$ID\tfasterq-dump" >> $error_log_file
            continue
        fi

        # zip fastqs to save space
        # use default number of processes, which is 8 or number of online processors
        pigz "${fastq_dir}"/${ID}_*.fastq

        # offload prefetch prerequisites
        rm -R "${fastq_dir}/${ID}"

    fi

done
