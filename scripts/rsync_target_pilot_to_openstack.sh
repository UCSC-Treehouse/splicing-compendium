#!/bin/bash
# usage: bash scripts/rsync_target_pilot_to_openstack.sh target

# cause nonzero exit status and undefined variables to stop the script
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# set floating IP of openstack that is the file source as variable
ip="10.50.100.156"

## define directories ##
data_dir="../data"
log_dir="../logs"
group_dir="${data_dir}/target"
# destination repo dir
mustard_repo_dir="/private/groups/treehouse/working-projects/celiang/splicing-compendium"
# splice results dir
mustard_results_dir="${mustard_repo_dir}/results"
# shiba results dir - within here are annotation, events, junction, and splice/gene expression results files
mustard_shiba_dir="${results_dir}/target/shiba"

## define input files ##
# experiment_tsv of pilot results
target_pilot_sample_file="../exploration/merging-psi-tables/target_pilot/experiment.tsv"

## define output files ##
# log file
log_file="${log_dir}/${current_datetime}_grab_target_pilot_results.txt"
target_pilot_manifest_file="${log_dir}/target_pilot_shiba_results_manifest.txt"

# create directories if they do not already exist
mkdir -p $log_dir

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

# instantiate .txt file to contain only TARGET pilot sample IDs
echo "" > $target_pilot_manifest_file

# print only TARGET pilot sample IDs to manifest file
awk 'NR > 1 {print $1}' $target_pilot_sample_file | while read sample; do
    echo $sample >> $target_pilot_manifest_file
done

# rsync all sample IDs from manifest sheet to openstack from prism
rsync -r -avP --files-from=$target_pilot_manifest_file celiang@mustard.prism:$mustard_shiba_dir ../results/target/shiba
