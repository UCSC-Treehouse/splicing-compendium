#!/bin/bash
# shell script for running diffsplice script through slurm
# usage: sbatch scripts/06-diffsplice_run_script.sh

#SBATCH --mem=20gb
#SBATCH --job-name=diffsplice_compendium
#SBATCH --time=60
#SBATCH --account=standard
#SBATCH --mail-user=celiang@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output logs/diffsplice_%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# run the rscript
Rscript scripts/06-diffsplice_calculation.R --reference_group "Whole Blood" --query_group AML --metadata_file combined_compendium_metadata.tsv --in_psi cleaned_psi_matrix.txt --out_file test_diffsplice.tsv
