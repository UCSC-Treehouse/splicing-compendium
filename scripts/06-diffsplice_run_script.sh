#!/bin/bash
# shell script for running diffsplice script through slurm
# usage: sbatch scripts/06-diffsplice_run_script.sh

#SBATCH --mem=40gb
#SBATCH --partition=medium
#SBATCH --job-name=diffsplice_compendium
#SBATCH --time=8:00:00
#SBATCH --account=standard
#SBATCH --mail-user=celiang@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output logs/diffsplice_%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# run the rscript
Rscript scripts/06-diffsplice_calculation.R --reference_group "Whole Blood" --query_group AML --metadata_file combined_compendium_metadata.tsv --in_psi cleaned_psi_matrix.txt --out_file test_diffsplice.tsv
