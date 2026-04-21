#!/bin/bash

# To clear space on OpenStack for the continued processing of data, we need to offload processed data and results to /private/spinning/treehouse
# this directory has a 100TB quota but cannot perform with more than 3 threads
# this script transfers bam files and their indices to /private/spinning/treehouse

# Usage example
# bash scripts/transfer-and-offload-bams.sh target transfer

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# set floating IP of openstack with data as variable
ip="openstack"

# Set paths as variables
git_path=$(git rev-parse --git-dir)
# we need the root dir so that repo dirs can be accessed within data/ like from within star-output/
root_dir=$(dirname "$git_path")
data_dir="${root_dir}/data"
log_dir="${root_dir}/logs"
# storage paths outside of repo
group_dir="${data_dir}/$1"
# in each star_output sample dir, there is the bam, bam.bai, logs, ReadsPerGene.out.tab, and SJ.out.tab
bam_dir="${group_dir}/star-output"
# destination directory
destination_dir="/private/spinning/treehouse"

# create directories if they do not already exist
mkdir -p $log_dir

# define output files
md5sums="${log_dir}/${ip}_${1}_md5sum.txt"

# validate user input
if [ $2 == "transfer" ]; then
    # define log file
    log_file="${log_dir}/${current_datetime}_file-transfer.txt"
    echo "transferring files"
elif [ $2 == "checksums" ]; then
    # define log file
    log_file="${log_dir}/${current_datetime}_transfer_checksum.txt"
    echo "perform md5 checksum of files"
elif [ $2 == "offload" ]; then
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

### file transfer (assumes you are in openstack instance with all the files to be transferred) ###
if [ $2 == "transfer" ]; then
    # generate md5 checksum file of star output files to be transferred
    if [ ! -f "${md5sums}" ]; then
        cd $bam_dir
        # bam files are in subdirectories labeled by sample type
        # recursively make md5sums of everthing in subdirectories
        find -type f -exec md5sum '{}' \; > "${md5sums}"
    fi

    # transfer sequence files to Ceph storage
    rsync -avP ${bam_dir} celiang@mustard.prism:${destination_dir}
fi

### md5sum check files that have been transferred (assumes you are in mustard directory with transferred files) ###
if [ $2 == "checksums" ]; then
    cd $bam_dir
    md5sum -c $md5sums
fi

### offload successfully transferred files ###
if [ $2 == "offload" ]; then

    # be in fastq directory to use relative paths
    cd $bam_dir

    # check what files have matching md5sums
    # md5sum logfile has lines like this if checksum succeeds: '/mnt/bulk/target/fastq/SRR2083188_2.fastq.gz: OK'
    # print first and second columns of each line in checksum log and only remove files corresponding to lines with ":OK" substring
    for file in $(awk -F': ' '$2 == "OK" {print $1}' $2); do
        rm $file
    done
fi
