#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# read in project accession as variable
accession="PRJNA597870"

# Set paths as variables
data_dir=$(realpath ..)
ref_dir="${data_dir}/references"
genome_dir="${ref_dir}/gencode.v47.primary_assembly-STAR-database"
fa_file="${ref_dir}/GRCh38.primary_assembly.genome.fa.gz"
gtf_file="${ref_dir}/gencode.v47.primary_assembly.annotation.gtf"
log_dir="${data_dir}/logs"

# Make output directory for STAR database
mkdir -p $genome_dir

# Set time and date as a variable
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# Redirect stdout and stderr to log file
exec > ${log_dir}/${current_datetime}_star-database-log.txt 2>&1

# Unzip reference fasta for STAR indexing
if [ ! -f ${fa_file%.gz} ]; then
  gunzip -c $fa_file > ${fa_file%.gz}
fi

# Generate genome index if there are no files in genome_dir
if [ ! "$(ls -A $genome_dir)" ]; then
  time STAR \
    --genomeDir $genome_dir \
    --genomeFastaFiles ${fa_file%.gz} \
    --sjdbGTFfile ${gtf_file} \
    --runMode genomeGenerate \
    --runThreadN 20
fi

# Remove uncompressed reference fasta
rm ${fa_file%.gz}
