### merge per-sample PSI tables produced by merge_results.smk ###
#
# reads every PSI output within the sample_psi/ directory and combines them so that all samples' PSI results appear in each file.
# samples may have NA values in their PSI column if junction counts for the event < 10
#
# the sample name for each file is taken from the sample sheet passed to the workflow's configfile
#
# usage:
#   Rscript 04-merge_sample_psi.R --sample_sheet=sample_sheet.tsv --version_dir=v1.1.0_two_sample --output_dir=out_dir

# load packages
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(duckplyr)
})

# set up options to Rscript with optparse
option_list <- list(
  make_option(
    opt_str = "--sample_sheet",
    type = "character",
    action = "store",
    help = "path to sample sheet with list of sample IDs"
  ),
  make_option(
    opt_str = "--version_dir",
    type = "character",
    action = "store",
    help = "name of version directory of separate PSI results"
  ),
  make_option(
    opt_str = "--output_dir",
    type = "character",
    action = "store",
    help = "name of output directory of merged PSI results"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## Validate output options ##
# user must provide sample_sheet and version_dir to script
if ((is.null(opt$sample_sheet) || is.null(opt$version_dir)) || is.null(opt$output_dir)) {
  stop("Specify --sample_sheet, --version_dir, and --output_dir.")
}

### read in file paths ###

## Directories ##
# find project root to access separate results dir
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# result directories
results_dir <- file.path(repo_root, "results")
compendium_results_dir <- file.path(results_dir, "merged_shiba")
version_dir <- file.path(compendium_results_dir, opt$version_dir)
psi_dir <- file.path(version_dir, "sample_psi")

# merged psi table output dir
output_dir <- file.path(opt$output_dir)

# create output dir if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Files ##

# sample sheet file with sample names
samples_file <- file.path(opt$sample_sheet)

# read in sample sheet and create list of samples from samples column
samples <- readr::read_tsv(samples_file, col_types = list(.default = "c")) |>
  dplyr::pull(sample) |>
  as.list()

# make sample psi table paths
sample_paths <- file.path(
  psi_dir,
  samples
)

# output files
out_file_list <- c(
  se = "PSI_SE.txt",
  afe = "PSI_AFE.txt",
  ale = "PSI_ALE.txt",
  five = "PSI_FIVE.txt",
  three = "PSI_THREE.txt",
  mse = "PSI_MSE.txt",
  mxe = "PSI_MXE.txt",
  ri = "PSI_RI.txt",
  matrix = "PSI_matrix_sample.txt"
)

# construct output file paths
out_paths <- file.path(
  output_dir,
  out_file_list
)

names(out_paths) <- names(out_file_list)

### define functions ###

# function for reading per-sample sample PSI matrices ("PSI_matrix_sample.txt")
read_sample_psi_matrix <- function(psi_path) {
  # returns a PSI matrix data frame of all samples in sample sheet
  # construct paths to psi sample matrices
  file_paths <- file.path(psi_path, "PSI_matrix_sample.txt")
  # read psi matrix files and merge them together with duckplyr
  matrix_tables <- purrr::map(file_paths, \(file) {
    duckplyr::read_csv_duckdb(file, options = list(delim = "\t"))
  })

  purrr::reduce(matrix_tables, \(x, y) dplyr::left_join(x, y, by = c("event_id", "pos_id"))) |>
    dplyr::collect() # only materialize once, at the very end
}

# function for reading each event types' per-sample PSI tables (e.g. "PSI_SE.txt")
read_event_table <- function(sample_paths, event_table_name) {
  # returns one data frame of skipped exon event types with all samples' psi values as separate columns
  # construct paths to each psi output for each sample
  file_paths <- file.path(sample_paths, event_table_name)

  # read in all paths into tables and merge them together
  duckplyr::read_csv_duckdb(
    file_paths,
    options = list(delim = "\t")
  ) |>
    # move identifying columns to the front for spot-checking
    dplyr::relocate(
      "event_id",
      "pos_id",
      "gene_id"
    ) |>
    dplyr::collect()
}

### read in and merge psi tables ###
# merge each PSI table
# then write to output and remove it after writing

# make list of event types to loop through
event_types <- c("se", "afe", "ale", "five", "three", "mse", "mxe", "ri")

# loop through event types
for (event_type in event_types) {
  # print message for log
  message("Merging ", event_type, " PSI tables")

  # create merged table object
  merged_event_table <- read_event_table(sample_paths, out_file_list[event_type])
  readr::write_tsv(merged_event_table, out_paths[[event_type]])

  # remove table and clear from memory
  rm(merged_event_table)
  gc()

}

# print message when merging matrix
message("Merging PSI sample matrix")
# merge matrices and write merged matrix to output
merged_matrix <- read_sample_psi_matrix(sample_paths)
readr::write_tsv(merged_matrix, out_paths[["matrix"]])

# remove merged matrix from memory
rm(merged_matrix)
gc()
