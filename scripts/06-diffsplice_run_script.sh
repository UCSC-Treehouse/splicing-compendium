#!/bin/bash
# shell script for running diffsplice script through slurm
# usage: sbatch scripts/06-diffsplice_run_script.sh

#SBATCH --mem=10gb
#SBATCH --partition=long
#SBATCH --job-name=diffsplice_compendium
#SBATCH --time=2:00:00
#SBATCH --account=standard
#SBATCH --mail-user=celiang@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output logs/diffsplice_%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# run the rscript
Rscript scripts/06-diffsplice_calculation.R --reference_group "Whole Blood" --query_group AML --metadata_file metadata/combined_compendium_metadata.tsv --psi_file test_files/subsetted_splice_matrix.txt --out_file test_diffsplice.tsv
