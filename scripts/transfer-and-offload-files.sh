#!/bin/bash

# To clear space on OpenStack for the continued processing of data, we need to offload processed data and results to /private/spinning/treehouse
# this directory has a 100TB quota but cannot perform with more than 3 threads
# this script transfers bam files and their indices to /private/spinning/treehouse

# Usage example for data transfer
# bash scripts/transfer-and-offload-bams.sh target transfer shiba

# Usage example for checksumming transferred files (must scp md5sum files to mustard first)
# bash scripts/transfer-and-offload-files.sh target checksums shiba

# Usage example for offloading files after transfer
# bash scripts/transfer-and-offload-bams.sh target checksums shiba 2026-04-24T18:33:08_openstack_target_transfer_checksum.txt

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# set floating IP of openstack that is the file source as variable
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
# splice results dir
results_dir="${root_dir}/results"
# shiba results dir - within here are annotation, events, junction, and splice/gene expression results files
shiba_dir="${results_dir}/$1/shiba"
# root destination directory
destination_root_dir="/private/spinning/treehouse"
# bam results destination dir
bam_dest_dir="${destination_root_dir}/data/$1"
# shiba results destination dir
shiba_dest_dir="${destination_root_dir}/results/$1"
# path to md5sum check results file
md5sum_checks="${log_dir}"/${4:-"onlyForOffloading"}

# create directories if they do not already exist
mkdir -p $log_dir

# define output files
md5sums=""${log_dir}"/"${ip}"_"${1}"_"${3}"_md5sum.txt"

# validate user input for data group
if [ $1 == "gtex" ]; then
    echo "gtex sample group"
elif [ $1 == "target" ]; then
    echo "target sample group"
else
    echo "please use valid option for sample group"
    # cause script to fail due to error
    exit 1
fi

# validate user input for script action
if [ $2 == "transfer" ]; then
    # define log file
    log_file="${log_dir}/${current_datetime}_${ip}_${1}_file-transfer.txt"
    echo "transferring files"
elif [ $2 == "checksums" ]; then
    # define log file
    log_file="${log_dir}/${current_datetime}_${ip}_${1}_transfer_checksum.txt"
    echo "perform md5 checksum of files"
elif [ $2 == "offload" ]; then
    # define log file
    log_file="${log_dir}/${current_datetime}_${ip}_${1}_offload_files.txt"
    echo "offload files"
else
    echo "please use valid option for what action to perform"
    # cause script to fail due to error
    exit 1
fi

# validate user input for what files to action on
if [ $3 == "bam" ]; then
    file_dir=$bam_dir
    destination_dir=$bam_dest_dir
    echo "star-output files"
elif [ $3 == "shiba" ]; then
    file_dir=$shiba_dir
    destination_dir=$shiba_dest_dir
    echo "shiba results files"
else
    echo "please use valid option for what files to act on"
    # cause script to fail due to error
    exit 1
fi

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

### file transfer (assumes you are in openstack instance with all the files to be transferred) ###
if [ $2 == "transfer" ]; then
    # generate md5 checksum file of star output files to be transferred
    if [ ! -f "${md5sums}" ]; then
        cd $file_dir
        # bam files are in subdirectories labeled by sample type
        # recursively make md5sums of everthing in subdirectories
        find -type f -exec md5sum '{}' \; > "${md5sums}"
    fi

    # transfer sequence files to Ceph storage
    rsync -avP ${file_dir} celiang@mustard.prism:${destination_dir}

fi

### md5sum check files that have been transferred (assumes you are in mustard directory with transferred files) ###
if [ $2 == "checksums" ]; then
    cd $file_dir
    md5sum -c $md5sums
fi

### offload successfully transferred files ###
# make sure this is done on the OpenStack instance, not mustard

if [ $2 == "offload" ]; then

    # be in fastq directory to use relative paths
    cd $file_dir

    # check what files have matching md5sums
    # md5sum logfile has lines like this if checksum succeeds: '/mnt/bulk/target/fastq/SRR2083188_2.fastq.gz: OK'
    # print first and second columns of each line in checksum log and only remove files corresponding to lines with ":OK" substring
    for file in $("awk -F': ' '$2 == "OK" {print $1}' $md5sum_checks"); do
        echo "deleting $file"
        rm $file
    done

    # check if any files failed md5sum check and list which files remain
    if [ -z "$(find $file_dir -type f -print -quit 2>/dev/null)" ]; then
        echo "$file_dir is empty, all files removed"
        # delete the directory (otherwise empty subdirectories will hang around and make it confusing to keep track of progress in new batches)
        rm -R $file_dir
    else
        echo "Directory is not empty. Files found:"
        find "$file_dir" -type f
    fi

fi
