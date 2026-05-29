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

# we need the root dir so that repo dirs can be accessed within data/ like from within star-output/
data_dir="../data"
log_dir="../logs"
# storage paths outside of repo
group_dir="${data_dir}/$1"
# destination repo dir
mustard_repo_dir="/private/groups/treehouse/working-projects/celiang/splicing-compendium"
# splice results dir
results_dir="${mustard_repo_dir}/results"
# shiba results dir - within here are annotation, events, junction, and splice/gene expression results files
shiba_dir="${results_dir}/$1/shiba"
# log dir
log_file="${log_dir}/${current_datetime}_grab_target_pilot_results"
# experiment_tsv of pilot results
target_pilot_sample_file="../exploration/merging-psi-tables/target_pilot/experiment.tsv"

# create directories if they do not already exist
mkdir -p $log_dir

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

# Redirect stdout and stderr to log file and print to stdout
exec > >(tee $log_file) 2>&1

touch target_pilot_shiba_results_manifest.txt

awk 'NR > 1 {print $1}' $target_pilot_sample_file | while read sample; do
    echo $sample >> target_pilot_shiba_results_manifest.txt
done

rsync -r -avP --files-from=target_pilot_shiba_results_manifest.txt celiang@mustard.prism:$shiba_dir ../results/$1/shiba
