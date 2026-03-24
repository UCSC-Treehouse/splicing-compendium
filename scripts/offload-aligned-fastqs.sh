#!/bin/bash

# When OpenStack instances become full before data processing is done, we need to offload raw sequence files to clear space.
# This script removes raw fastq files belonging to accession IDs of aligned files present in the star output directory

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# set working directory to be the directory of this file (scripts/)
cd "$(dirname "${BASH_SOURCE[0]}")"

# set time and date as variables
datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# set paths as variables
data_dir="../data"
group_dir="${data_dir}/target"
bam_dir="${group_dir}/star-output"
fastq_dir="${group_dir}/fastq"

# iterate through star output directories to obtain accessions of samples with alignments present
for accession_dir in "${bam_dir}"/*/; do
  if [ -f "${accession_dir}/Aligned.sortedByCoord.out.bam" ]; then
    accession=$(basename "$accession_dir")
  else
    # go to the next directory if the aligned file does not exist
    continue
  fi
  
    # construct paths to fastqs based on accessions in the bam dir
    fastq1_path="${fastq_dir}/${accession}_1.fastq.gz"
    fastq2_path="${fastq_dir}/${accession}_2.fastq.gz"

    # remove fastq files if they exist
   rm -f ${fastq1_path} ${fastq2_path}
done
