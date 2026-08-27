### merge per-sample PSI tables produced by merge_results.smk ###
#
# reads every PSI output within the sample_psi/ directory and combines them so that all samples' PSI results appear in each file.
# samples may have NA values in their PSI column if junction counts for the event < 10
#
# the sample name for each file is taken from the sample sheet passed to the workflow's configfile
#
# usage:
#   Rscript 04-merge_sample_psi.R --sample_sheet=sample_sheet.tsv --version_dir=v1.1.0_two_sample

# load packages
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

# set up options to Rscript with optparse
option_list <- list(
  make_option(
    opt_str = "--sample_sheet",
    type = "character",
    action = "store",
    help = "name of sample sheet with list of sample IDs"
  ),
  make_option(
    opt_str = "--version_dir",
    type = "character",
    action = "store",
    help = "name of version directory of separate PSI results"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## Validate output options ##
# user must provide sample_sheet and version_dir to script
if ((is.null(opt$sample_sheet) || is.null(opt$version_dir))) {
  stop("Specify --sample_sheet and --version_dir.")
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

# sample sheet dir
config_dir <- file.path(repo_root, "config")

# merged psi table output dir
output_dir <- file.path(version_dir, "merged_psi")

# create output dir if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Files ##

# sample sheet file with sample names
samples_file <- file.path(config_dir, opt$sample_sheet)

# read in sample sheet and create list of samples from samples column
samples <- readr::read_tsv(samples_file, col_types = list(.default = "c")) |>
  dplyr::pull(sample) |>
  as.list()

# define list of PSI results files corresponding to event types quantified by Shiba bulk analysis
psi_file_list <- c(
  se = "PSI_SE.txt",
  afe = "PSI_AFE.txt",
  ale = "PSI_ALE.txt",
  five = "PSI_FIVE.txt",
  three = "PSI_THREE.txt",
  mse = "PSI_MSE.txt",
  mxe = "PSI_MXE.txt",
  ri = "PSI_RI.txt"
)

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
  # read psi matrix files
  matrix_list <- purrr::map(file_paths, \(file) {
    readr::read_tsv(file, col_types = readr::cols(.default = "c"))
  }) |>
    # Join them all by a common column name
    purrr::reduce(dplyr::left_join, by = c("event_id", "pos_id"))
}

# function for reading each event types' per-sample PSI tables (e.g. "PSI_SE.txt")
read_event_table <- function(sample_paths, event_table_name) {
  # returns one data frame of skipped exon event types with all samples' psi values as separate columns
  # construct paths to each psi output for each sample
  file_paths <- file.path(sample_paths, event_table_name)

  # loop over all samples in sample_paths
  purrr::map(file_paths, \(file) {
    # read each file into a table
    readr::read_tsv(file, col_names = TRUE, col_types = readr::cols(.default = "c"))
  }) |>
    # merge per-sample tables into one vertically
    dplyr::bind_rows() |>
    # move identifying columns to the front for spot-checking
    dplyr::relocate(
      "event_id",
      "pos_id",
      "gene_id"
    )
}

### read in and merge psi tables ###

merged_se_table <- read_event_table(sample_paths, out_file_list["se"])
merged_matrix <- read_sample_psi_matrix(sample_paths)

### write output ###

# roll outputs into a list for writing outputs with purrr
all_outputs <- list(
  se = merged_se_table,
  matrix = merged_matrix
)

# give purrr iwalk the list of named output files and write output using named output paths
purrr::iwalk(all_outputs, \(df, name) {
  readr::write_tsv(df, out_paths[[name]])
})
