#!/bin/bash

# Download script of RNA-seq for pediatric cancer splicing analysis
# usage: scripts/01-target-download.sh [dataset] [batch number]
# example for downloading the first 1.2TB batch of target: scripts/01-target-download.sh target 1

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# dbGaP key files
target_key="/home/ubuntu/prj_11732_D43588.ngc"
gtex_key="/home/ubuntu/prj_9508_D43588.ngc"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# Set paths as variables
data_dir="../data"
group_dir="${data_dir}/${1}"
fastq_dir="${group_dir}/fastq"
log_dir="../logs/${1}_download"
metadata_dir="../metadata/filter_target_gtex"
notebooks_dir="../notebooks/data_filtering"
target_filtering="${notebooks_dir}/filter_target_sra.qmd"
gtex_filtering="${notebooks_dir}/filter_gtex_sra.qmd"
sample_sheet="${metadata_dir}/${1}_batch_${2}_samples.txt"

# check that sample group inputted is valid
# set dbgap key based on whether gtex or target files are to be downloaded
if [ "$1" == "target" ]; then
    dbgap_key=$target_key
    filtering=$target_filtering

elif [ "$1" == "gtex" ]; then
    dbgap_key=$gtex_key
    filtering=$gtex_filtering

else
    echo "please use valid sample group"
    # cause script to fail due to error
    exit 1
fi

# create directories if they do not already exist
mkdir -p $log_dir
mkdir -p $data_dir
mkdir -p $group_dir
mkdir -p $fastq_dir

# Define log files
log_file="${log_dir}/${current_datetime}_${1}_download.txt"
error_log_file="${log_dir}/${current_datetime}_${1}_download_errors.txt"

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

# Run filtering script to obtain accession IDs to download
quarto render ${filtering}

# Select 1.2TB batch of accessions to download based on user input
# pass batch ID from user input (2nd input) into awk and look for values matching batch in last column
# skip first line (header)
# make sure value in the last column equals the batch given by user input and print the first column (accession ID)
# print the set of accession IDs to be downloaded to spot-check
awk -v batch="${2}" 'NR>1 && $NF == batch {print $1}' "${metadata_dir}/${1}_accessions.tsv"

# iterate through first column of TARGET accession metadata after this filtering, which has accession IDs
time awk -v batch="${2}" 'NR>1 && $NF == batch {print $1}' "${metadata_dir}/${1}_accessions.tsv" | while read -r ID; do
    # check that fastq or zipped fastq does not already exist
    if [ ! -f "${fastq_dir}/${ID}_1.fastq" ] && [ ! -f "${fastq_dir}/${ID}_1.fastq.gz" ]; then

        # prefetch file dependencies if neither fastq nor zipped fastq exist
        # set max size of prefetch file to 40GB based on max file size of samples in accessions lists
        # continue even if an error is thrown
        prefetch --ngc $dbgap_key --max-size 40000000000 --output-directory $fastq_dir $ID || true

        # Check if prefetch file directory has been created
        # Lack of directory means prefetch has failed; skip to next accession
        if [ ! -d "${fastq_dir}/$ID" ]; then
            # record the accession ID not downloaded
            echo -e "$ID\tprefetch" >> $error_log_file
            continue
        fi

        # download fastq file for each accession
        fasterq-dump --temp $group_dir --ngc $dbgap_key $ID --threads 15 --outdir $fastq_dir || true

        # Check if fastq file exists
        # If file is missing, record error in fastq-dump step and move on to next accession
        if [[ ! -f "${fastq_dir}"/${ID}_1.fastq || ! -f "${fastq_dir}"/${ID}_2.fastq ]]; then
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

# write files in fastq dir into sample sheet for snakemake
find $fastq_dir -type f -name "*.fastq.gz" \
# strip path from each file so we just get the filename
    | xargs -n1 basename \
    # remove file suffix (everything after first underscore) to obtain sample ID
    | sed -E 's/_[12].fastq\.gz//' \
    # only keep unique accessions
    | sort -u > $sample_sheet
