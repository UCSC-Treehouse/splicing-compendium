## description ##
# this script downloads clinical metadata for the TARGET cancer types in splice compendium v1.
# clinical metadata are

#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

## data directories ##
metadata_dir="../metadata"
clinical_metadata_dir="${metadata_dir}/clinical"
target_clin_dir="${clinical_metadata_dir}/target"
log_dir="../logs"

## file paths ##


# make directories if they don't exist
mkdir -p $target_clin_dir
mkdir -p $log_dir

# download metadata
