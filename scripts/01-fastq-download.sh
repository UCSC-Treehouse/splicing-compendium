#!/bin/bash

# Download script of RNA-seq for pediatric cancer splicing analysis
## 2027_xx ##
# usage: bash scripts/01-fastq-download.sh target_osteo

## 2026_09 ##
# usage: scripts/01-fastq-download.sh [dataset] [batch number]
# example for downloading the first 1.2TB batch of target: scripts/01-target-download.sh target 1

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# dbGaP key files
target_key="/private/home/celiang/ngc_keys/prj_11732_D43588.ngc"
gtex_key="/private/home/celiang/ngc_keys/prj_9508_D43588.ngc"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# Set paths as variables
data_dir="../data"
group_dir="${data_dir}/${1}"
fastq_dir="${group_dir}/2027/fastq" # should I make a 2027 dir or keep these new osteo samples all together?
log_dir="../logs/2027/${1}_download"
metadata_dir="../metadata/filter_target_gtex"
config_dir="../config/2027"
notebooks_dir="../notebooks/data_filtering"
target_filtering="${notebooks_dir}/filter_target_sra.qmd"
gtex_filtering="${notebooks_dir}/filter_gtex_sra.qmd"
sample_sheet="${config_dir}/${1}_samples.tsv"

# check that sample group inputted is valid
# set dbgap key based on whether gtex or target files are to be downloaded
if [[ "$1" == "target" ]]; then
    dbgap_key=$target_key
    filtering=$target_filtering
    accessions_file="$metadata_dir/target_accessions.tsv"

elif [[ "$1" == "target_osteo" ]]; then
    dbgap_key=$target_key
    filtering=$target_filtering
    accessions_file="$metadata_dir/target_osteo_accessions.txt"

    echo "target osteo filtering"

elif [[ "$1" == "gtex" ]]; then
    dbgap_key=$gtex_key
    filtering=$gtex_filtering
    accessions_file="$metadata_dir/gtex_accessions.txt"

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
mkdir -p $config_dir

# Define log files
log_file="${log_dir}/${current_datetime}_${1}_download.txt"
error_log_file="${log_dir}/${current_datetime}_${1}_download_errors.txt"

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

# Run filtering script to obtain accession IDs to download
quarto render ${filtering}

echo "filtering notebook rendered, accessions generated for $1"

if [[ "$1" == "target_osteo" ]]; then

    echo "$1 read in for sra loop"

    # print the set of accession IDs to be downloaded to spot-check
    # set field separator as tab
    # skip first line (header)
    # print first column (accession IDs)
    awk -F '\t' 'NR > 1 {print $1}' ${accessions_file}

    # iterate through first column of TARGET accession metadata after this filtering, which has accession IDs
    # read each item in the column as ID object
    time for ID in $(awk -F '\t' 'NR>1{print $1}' "${accessions_file}"); do
        # check that fastq or zipped fastq does not already exist
        if [[ ! -f "${fastq_dir}/${ID}_1.fastq" ]] && [[ ! -f "${fastq_dir}/${ID}_1.fastq.gz" ]]; then

            # prefetch file dependencies if neither fastq nor zipped fastq exist
            # set max size of prefetch file to 40GB based on max file size of samples in accessions lists
            # continue even if an error is thrown
            prefetch --ngc $dbgap_key --max-size 40000000000 --output-directory $fastq_dir $ID || true

            # Check if prefetch file directory has been created
            # Lack of directory means prefetch has failed; skip to next accession
            if [[ ! -d "${fastq_dir}/$ID" ]]; then
                # record the accession ID not downloaded
                echo -e "$ID\tprefetch" >> $error_log_file
                continue
            fi

            # download fastq file for each accession
            fasterq-dump --temp $group_dir --ngc $dbgap_key $ID --threads 1 --outdir $fastq_dir || true

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

fi

if [[ "$1" == "target" || "$1" == "gtex" ]]; then

    echo "batched accession download beginning"

    # print accessions of the dataset

    awk -v batch="${2}" 'NR>1 && $NF == batch {print $1}' ${accessions_file}

    # iterate through first column of TARGET accession metadata after this filtering, which has accession IDs
    time awk -v batch="${2}" 'NR>1 && $NF == batch {print $1}' ${accessions_file} | while read -r ID; do
        # check that fastq or zipped fastq does not already exist
        if [[ ! -f "${fastq_dir}/${ID}_1.fastq" ]] && [[ ! -f "${fastq_dir}/${ID}_1.fastq.gz" ]]; then

            # prefetch file dependencies if neither fastq nor zipped fastq exist
            # set max size of prefetch file to 40GB based on max file size of samples in accessions lists
            # continue even if an error is thrown
            prefetch --ngc $dbgap_key --max-size 40000000000 --output-directory $fastq_dir $ID || true

            # Check if prefetch file directory has been created
            # Lack of directory means prefetch has failed; skip to next accession
            if [[ ! -d "${fastq_dir}/$ID" ]]; then
                # record the accession ID not downloaded
                echo -e "$ID\tprefetch" >> $error_log_file
                continue
            fi

            # download fastq file for each accession
            fasterq-dump --temp $group_dir --ngc $dbgap_key $ID --threads 1 --outdir $fastq_dir || true

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

fi

# create sample sheet file with "samples" column header
echo "samples" > $sample_sheet

# write files in fastq dir into sample sheet for snakemake
# strip path from each file so we just get the filename
# remove file suffix (everything after first underscore) to obtain sample ID
# only keep unique accessions
find $fastq_dir -type f -name "*.fastq.gz" \
    | xargs -n1 basename \
    | sed -E 's/_[12].fastq\.gz//' \
    | sort -u >> $sample_sheet
