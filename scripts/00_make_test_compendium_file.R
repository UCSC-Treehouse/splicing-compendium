# make test compendium file by subsetting Treehouse Splice Compendium v1
# usage: rscript 00_make_test_compendium_file.R

## Setup ##
# set seed for random slicing
set.seed(0)

# files and directories
# directories
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
results_dir <- file.path(repo_root, "results")
test_dir <- file.path(repo_root, "test_files")

# files
splice_results <- file.path(results_dir, "merged_shiba", "splice_compendium_v1", "merged_persample_psi", "cleaned_psi_matrix.txt")
subsetted_results <- file.path(test_dir, "subsetted_splice_matrix.txt")

# create directories if they don't exist
if (!dir.exists(test_dir)) {
  dir.create(test_dir, recursive = TRUE)
}

## read in files ##

psi_matrix <- readr::read_tsv(splice_results)

# subset splice compendium by one GTEx and one TARGET disease type for testing downstream analysis scripts

subsetted_matrix <- psi_matrix |>
  # make a dataframe consisting of a random 10 position IDs per event type for testing
  dplyr::slice_sample(by = event_type, n = 10)

# write subsetted file to output
readr::write_tsv(subsetted_matrix, subsetted_results)
