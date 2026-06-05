#!/bin/bash

# quit if there are undefined variables or a nonzero exit status
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set time and date as variables
current_datetime=$(date +"%Y-%m-%dT%H:%M:%S")

# files and directories
# directories
exploration_dir="../exploration"
target_pilot_combined_results_dir="${exploration_dir}/merging-psi-tables/target_pilot/shiba_combined/"
target_pilot_combined_junctions_dir="${target_pilot_combined_results_dir}/junctions"

# files
# target pilot combined junctions bedfile
combined_junctions_bedfile="${target_pilot_combined_junctions_dir}/junctions.bed"
#output file
target_pilot_deduplicated_file="${target_pilot_combined_junctions_dir}/deduplicated_junctions.bed"


# Remove duplicated fields in bedfiles separated by ";" inside the tab-delimited bedfile
awk '
BEGIN{FS=OFS="\t"}

{ # set tab as delimiter
    for (i = 1; i <= NF; i++) {
        if ($i ~ /;/) { # only process fields containing ";"
            n = split($i, parts, ";") # split into parts on semicolons
            new = parts[1]

            for (j = 2; j <= n; j++) {
                if (parts[j] != parts[j-1]) { # skip consecutive duplicates
                    new = new ";" parts[j] # reassemble if values are different
                }
            }

            $i = new
        }
    }

    print
}
' $combined_junctions_bedfile > $target_pilot_deduplicated_file
