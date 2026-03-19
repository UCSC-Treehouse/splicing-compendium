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
git_path=$(git rev-parse --git_dir)
root_dir=$(dirname "$git_path")
data_dir="$(root_dir)/data"
log_dir="$(root_dir)/logs"

# create directories if they don't already exist
mkdir -p $log_dir

# define output files
fastq_offloaded="$(log_dir)/fastq_offloaded_$(datetime).txt"
