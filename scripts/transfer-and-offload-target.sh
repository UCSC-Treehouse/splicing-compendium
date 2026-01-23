#!/bin/bash

# There is not enough space to hold all raaw sequence files and alignments, so offload TARGET sequence files to another instance
# this script will either initiate file transfer, or check md5sums of transferred files, depending on user options
# usage for file transfer: bash transfer-and-offload-target.sh transfer
# usage for md5sum check: bash transfer-and-offload-target.sh checksums
# usage for offloading transferred files after checksums: bash transfer-and-offload-target.sh offload data/logs/2025-10-01T04:12:24_transfer_checksum.txt

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# Set paths as variables
# repo paths
git_path=$(git rev-parse --git-dir)
root_dir=$(dirname "$git_path")
data_dir="${root_dir}/data"
log_dir="${data_dir}/logs"
# storage paths outside of repo
storage_dir="/mnt"
bulk_dir="${storage_dir}/bulk"
group_dir="${bulk_dir}/target"
fastq_dir="${group_dir}/fastq"

# create directories if they do not already exist
mkdir -p $log_dir

# define output files
md5sums="${fastq_dir}/target_md5sum.txt"

# validate user input
if [ $1 == "transfer" ]; then
    # define log file
    log_file="${log_dir}/${current_datetime}_file-transfer.txt"
    echo "transferring files"
elif [ $1 == "checksums" ]; then
    # define log file
    log_file="${log_dir}/${current_datetime}_transfer_checksum.txt"
    echo "perform md5 checksum of files"
elif [ $1 == "offload" ]; then
    # define log file
    log_file="${log_dir}/${current_datetime}_offload_files.txt"
    echo "offload files"
else
    echo "please use valid option"
    # cause script to fail due to error
    exit 1
fi

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

### file transfer (assumes you are in openstack instance with all the raw sequences) ###
if [ $1 == "transfer" ]; then
    # generate md5 checksum file of raw sequence files to be transferred if it does not already exist
    if [ ! -f "${md5sums}" ]; then
        cd $fastq_dir
        md5sum *.fastq.gz > $md5sums
    fi

    # transfer sequence files to separate openstack environment
    rsync -avP ${fastq_dir} ubuntu@10.50.100.47:${fastq_dir}
fi

### md5sum check files that have been transferred (assumes you are in openstack instance with transferred files) ###
if [ $1 == "checksums" ]; then
    cd $fastq_dir
    md5sum -c $md5sums
fi

### offload successfully transferred files ###
if [ $1 == "offload" ]; then

    # be in fastq directory to use relative paths
    cd $fastq_dir

    # check what files have matching md5sums
    # md5sum logfile has lines like this if checksum succeeds: '/mnt/bulk/target/fastq/SRR2083188_2.fastq.gz: OK'
    # print first and second columns of each line in checksum log and only remove lines with ":OK" substring
    for file in $(awk -F': ' '$2 == "OK" {print $1}' $2); do
        rm $file
    done
fi
