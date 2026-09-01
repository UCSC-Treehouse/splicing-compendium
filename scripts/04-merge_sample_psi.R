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
samples <- duckplyr::read_csv_duckdb(samples_file,
options = list(
  delim = "\t",
  union_by_name = TRUE,
  header = TRUE
  )) |>
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
read_sample_psi_matrix <- function(psi_path, samples, out_path) {
  # returns a PSI matrix data frame of all samples in sample sheet
  # construct paths to psi sample matrices
  file_paths <- file.path(psi_path, "PSI_matrix_sample.txt")

  # obtain duckdb connection object so query won't open a new DuckDB session
  con <- duckplyr:::get_default_duckdb_connection()

  # increase global max expression depth in case of large samples
  DBI::dbExecute(con, "SET GLOBAL max_expression_depth TO 10000")

  # first reshape wide matrices into long format prior to combining
  # cols should be event_id, pos_id, sample, psi
  # build a SQL select statement for each sample to read each file, match columns by name, and rename PSI column with sample name for merging
  selects <- purrr::map2_chr(file_paths, samples, \(file_path, sample) {
    sprintf(
      "SELECT event_id, pos_id, '%s' AS sample, \"%s\" AS psi FROM read_csv('%s', delim='\t', header=true, union_by_name=true)",
      sample, sample, file_path
    )
  })

  # merge all sample tables vertically into a long table
  union_sql <- paste(selects, collapse = " UNION ALL ")

  # pivot long table back to wide form
  # then write merged table directly from stream to output
  merged_table <- sprintf(
    "COPY (
       PIVOT (%s)
       ON sample
       USING first(psi)
       GROUP BY event_id, pos_id
     ) TO '%s' (DELIMITER '\t', HEADER)",
    union_sql, out_path
  )

  # execute above queries in the connection
  DBI::dbExecute(con, merged_table)

}

# function for reading each event types' per-sample PSI tables (e.g. "PSI_SE.txt")
assemble_event_table <- function(sample_paths, event_table_name, out_path) {
  # returns one data frame; rows are events; all samples' psi values are in separate columns
  # construct vector of paths to each psi output for each sample
  file_paths <- file.path(sample_paths, event_table_name)

  # format paths into a SQL list by joining the paths into a string separated by , and enclosed by brackets
  files_sql <- paste0("['", paste(file_paths, collapse = "','"), "']")

  # obtain duckdb connection object so query won't open a new DuckDB session
  con <- duckplyr:::get_default_duckdb_connection()

  # build SQL query string for files in the file list
  # read all files in list of paths as one table with cols matched by header names
  # reorders columns so that event_id, pos_id, and gene_id are at teh front
  # then streams the result to out_path
  query <- sprintf(
    "COPY (
       SELECT event_id, pos_id, gene_id, * EXCLUDE (event_id, pos_id, gene_id)
       FROM read_csv(%s, delim = '\t', header=true, union_by_name = true)
     ) TO '%s' (DELIMITER '\t', HEADER)",
    files_sql, out_path
  )

  # execute query on the connection
  DBI::dbExecute(con, query)
}

### read in and merge psi tables ###

# print message when merging matrix
message("Merging PSI sample matrix")
# merge matrices and write merged matrix to output
read_sample_psi_matrix(sample_paths, samples, out_paths[["matrix"]])

# make list of event types to loop through
event_types <- c("se", "afe", "ale", "five", "three", "mse", "mxe", "ri")

# loop through event types
for (event_type in event_types) {
  # print message for log
  message("Merging ", event_type, " PSI tables")

  # create merged table object
  assemble_event_table(sample_paths, out_file_list[[event_type]], out_paths[[event_type]])

}
