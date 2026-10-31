# make test compendium file by subsetting Treehouse Splice Compendium v1

# test file consists of random 10 pos_ids per splice event type per gene

# usage: rscript 00_make_test_compendium_file.R

## Setup ##
# set seed for random slicing
set.seed(0)

# files and directories
# directories
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
results_dir <- file.path(repo_root, "results")
test_dir <- file.path(repo_root, "test_files")
metadata_dir <- file.path(repo_root, "metadata")

# files
splice_results <- file.path(results_dir, "merged_shiba", "splice_compendium_v1", "merged_persample_psi", "cleaned_psi_matrix.txt")
metadata_file <- file.path(metadata_dir, "combined_compendium_metadata.tsv")
subsetted_results <- file.path(test_dir, "subsetted_splice_matrix.txt")

# create directories if they don't exist
if (!dir.exists(test_dir)) {
  dir.create(test_dir, recursive = TRUE)
}

## read in files ##

psi_matrix <- readr::read_tsv(splice_results)

# only select run and tissue type cols to subset random accessions from
psi_metadata <- readr::read_tsv(
  metadata_file,
  col_select = c(
    Run,
    plot_tissue_type
  )
)

# subset metadata by random 10 samples of an abritrary tumor and healthy tissue type (ALL and EBV lymphocytes)
subsetted_metadata <- psi_metadata |>
  dplyr::filter(
    plot_tissue_type %in% c("ALL", "EBV lymphocytes")
  ) |>
  dplyr::slice_sample(by = plot_tissue_type, n = 10)

# pull accession vector from subsetted metadata
subsetted_accessions <- subsetted_metadata |>
  dplyr::pull(Run)

# subset splice compendium by one GTEx and one TARGET disease type for testing downstream analysis scripts

# create list of columns to select for the subsetted splice matrix
subset_columns <- c(
  subsetted_accessions,
  "event_type",
  "label",
  "gene_name",
  "gene_id",
  "pos_id"
)

# subset splice matrix by 10 samples of reference and query group and 10 pos_id per splice event type
subsetted_matrix <- psi_matrix |>
  # filter for subsetted accessions from arbitrary ref and query group
  dplyr::select(
    all_of(
      subset_columns
    )
  ) |>
  # make a dataframe consisting of a random 10 position IDs per event type for testing
  dplyr::slice_sample(by = event_type, n = 10)


# write subsetted file to output
readr::write_tsv(subsetted_matrix, subsetted_results)
