#!/usr/bin/env bash

# usage:
# bash exploration/merged-method-target-pilot-eval/check_gtf_merge_order.sh exploration/merged-method-target-pilot-eval/slurm_pilot_gtf_merge_order.txt exploration/merging-psi-tables/target_pilot/experiment.tsv

set -euo pipefail

slurm_gtf_merge_order_file="$1"
experiment_file="$2"

# Extract SRR IDs from slurm output message
# use mapfile to read SRR IDs grepped from slurm message output into array
# -t strips newlines from input
mapfile -t slurm_merge_ids < <(grep -oE 'SRR[0-9]+' "$slurm_gtf_merge_order_file")

# Extract the sample column (first column) into an array, skipping the header
mapfile -t experiment_ids < <(tail -n +2 "$experiment_file" | cut -f1)

# Compare lengths of each list of IDs
if [[ ${#slurm_merge_ids[@]} -ne ${#experiment_ids[@]} ]]; then
    echo "Number of IDs differs:"
    echo "  Output:      ${#slurm_merge_ids[@]}"
    echo "  Experiment:  ${#experiment_ids[@]}"
fi

# Compare order of IDs
# set max length to number of target pilot samples
max=88
# start at position 0 and increment position until position is equal to max number of samples
for ((i=0; i<max; i++)); do
# check if sample ID at each position is equal in both arrays of IDs
    if [[ "${slurm_merge_ids[$i]}" != "${experiment_ids[$i]}" ]]; then
        echo "Order differs at position $((i+1)):"
        echo "  Output:      ${slurm_merge_ids[$i]}"
        echo "  Experiment:  ${experiment_ids[$i]}"
        exit 1
    else
        echo "Sample ID order is identical between slurm merged method and experiment.tsv"
    fi
done
