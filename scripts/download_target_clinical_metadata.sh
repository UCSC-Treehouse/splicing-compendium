## description ##
# this script downloads clinical metadata for the TARGET cancer types in splice compendium v1.
# clinical metadata are obtained from the UCSC Xena browser: https://xenabrowser.net/datapages/?hub=https://gdc.xenahubs.net:443
# two sets of metadata under "phenotype" in each TARGET cancer type's page are downloaded:
# "phenotype" (or "clinical" in url) which includes non-survival clinical fields such as age at diagnosis, ethnicity, gender
# "survival data" which includes "OS" and "OS.time" fields

## script ##

#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

## data directories ##
metadata_dir="../metadata"
clinical_metadata_dir="${metadata_dir}/clinical"
target_clin_dir="${clinical_metadata_dir}/target"
log_dir="../logs"

## file paths ##

## links to metadata ##
# array of metadata urls to download
clinical_metadata_urls=(
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-ALL-P1.clinical.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-ALL-P1.survival.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-ALL-P2.clinical.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-ALL-P2.survival.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-AML.clinical.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-AML.survival.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-NBL.clinical.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-NBL.survival.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-RT.clinical.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-RT.survival.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-WT.clinical.tsv.gz"
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TARGET-WT.survival.tsv.gz"
)

# make directories if they don't exist
mkdir -p $target_clin_dir
mkdir -p $log_dir

# Redirect stdout and stderr to log file
exec > >(tee "${log_dir}/${current_datetime}_target-clinical-metadata-download-log.txt") 2>&1

# download metadata
for url in ${clinical_metadata_urls[@]}; do
    file="${target_clin_dir}/$(basename $url)"

    # check that file exists, zipped or unzipped
    # if it does not exist, download the file and unzip it
    if [[ ! -f $file && ! -f ${file%.gz} ]]; then
        wget -nv $url -O $file
        gunzip $file
    fi
done
