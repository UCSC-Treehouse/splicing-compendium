
#!/bin/bash

# Download script for human reference and annotation files used fpr alignment and splicing analysis

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# Set paths as variables
data_dir=$(realpath ..)
log_dir="${data_dir}/logs"
ref_dir="${data_dir}/references"
metadata_dir="${data_dir}/metadata"
target_gtex_dir="${metadata_dir}/filter_target_gtex"

# make arrays of URLs to download
reference_urls=(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz"
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz"
)

gtex_metadata_urls=(
    "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
)

# Create directories
mkdir -p $log_dir
mkdir -p $ref_dir
mkdir -p $target_gtex_dir

# Redirect stdout and stderr to log file
exec > >(tee "${log_dir}/${current_datetime}_reference-file-download-log.txt") 2>&1

for url in ${reference_urls[@]}; do
    file="${ref_dir}/$(basename $url)"

    # check that file exists, zipped or unzipped
    # if it does not exist, download the file and unzip it
    if [[ ! -f $file && ! -f ${file%.gz} ]]; then
        wget -nv $url -O $file
    fi

    # Unzip gtf if it is not alrady unzipped
    if [[ $file =~ \.gtf\.gz$ && ! -f ${file%.gz} ]]; then
        gunzip -c $file > ${file%.gz}
    fi

    # Unzip barcode whitelists if they are not already unzipped
    if [[ $file =~ \.txt\.gz$ && ! -f ${file%.gz} ]]; then
        gunzip -c $file > ${file%.gz}
    fi

done

# Download metadata files
for url in ${gtex_metadata_urls[@]}; do
    file="${target_gtex_dir}/$(basename $url)"

    # check that file exists, zipped or unzipped
    # if it does not exist, download the file
    if [[ ! -f $file && ! -f ${file%.gz} ]]; then
        wget -nv $url -O $file
    fi

done
