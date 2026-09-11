# This code takes the following files as input:
# - merged sample PSI matrix
# - one sample's event tables (PSI_SE.txt, etc), excluding the retained intron table.
# Rows should correspond to all splice event coordinate IDs quantified in all samples and contain additional information for each splice event like gene ID and annotation status.
#
# This code performs the following:
# - Select for only gene_id, pos_id, and label (annotation status) rows in each event table
# - Merges one sample's event tables together to create one PSI table of all event types with event metadata columns for each event type
# - Creates an "event type" column in the merged table corresponding to each event type
# - Creates an "event type" column in the merged sample PSI matrix and removes RI events
# - rbinds the tables by pos_id

# usage: Rscript scripts/05-clean_sample_psi.R --sample_sheet={params.sample_sheet} --version_dir={params.version} --output_dir={params.out_dir}

### Load libraries ###
suppressPackageStartupMessages({
  library(rtracklayer)
  library(optparse)
})

### Read in options ###
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
    opt_str = "--in_matrix",
    type = "character",
    action = "store",
    help = "path to merged PSI results"
  ),
  make_option(
    opt_str = "--output",
    type = "character",
    action = "store",
    help = "path to output matrix"
  )
  make_option(
    opt_str = "--gtf",
    type = "character",
    action = "store",
    help = "path to reference annotation file"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## Validate output options ##
# user must provide sample_sheet and version_dir to script
if ((is.null(opt$sample_sheet) || is.null(opt$version_dir)) || is.null(opt$in_matrix) || is.null(opt$output) || is.null(opt$gtf)) {
  stop("Specify --sample_sheet, --version_dir, --in_matrix, --output, and --gtf.")
}

### Read in files and directories ###

## directories ##

# define the data directories
results_dir <- file.path(opt$version_dir)
sample_psi_dir <- file.path(results_dir, "sample_psi")

## files ##

# merged PSI results
merged_matrix <- file.path(opt$in_matrix)

# gtf file for converting ensg id to gene names
gtf_file <- file.path(opt$gtf)

# read in sample sheet and create list of samples from samples column
samples <- readr::read_tsv(
  opt$sample_sheet,
  col_types = readr::cols(.default = "c")
 ) |>
  dplyr::pull(sample)

# path to psi results of one sample from sample sheet
one_sample_results_dir <- file.path(sample_psi_dir, samples[[1]])

# define list of PSI event table results files corresponding to event types quantified by Shiba bulk analysis
event_files <- c(
  SE = "PSI_SE.txt",
  AFE = "PSI_AFE.txt",
  ALE = "PSI_ALE.txt",
  FIVE = "PSI_FIVE.txt",
  THREE = "PSI_THREE.txt",
  MSE = "PSI_MSE.txt",
  MXE = "PSI_MXE.txt"
)

# construct psi table paths of one sample's compendium Shiba results
event_psi_paths <- file.path(one_sample_results_dir, event_files)
# name file paths according to event type
names(event_psi_paths) <- names(event_files)

# output file path
out_matrix <- file.path(opt$output)

message("file paths loaded")

### Read in GTF and extract gene names ###
# import gtf
gtf <- rtracklayer::import(gtf_file, filter = list(type = "gene"))

# make a named vector of gene names to IDs
gene_names <- setNames(gtf$gene_name, gtf$gene_id)

# remove the gtf after names are extracted
rm(gtf)
# garbage collector to clear memory
gc()

message("gene names extracted from GTF")

### Read in PSI tables ###

# read in one sample's event tables
event_tables <- event_psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types = readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label)
  })

# combine event tables into one, with additional column labeling event type
all_events_table <- purrr::list_rbind(event_tables, names_to = "event_type")

message("event tables from one sample merged")

# read in merged sample matrix
sample_matrix <- readr::read_tsv(merged_matrix_file, col_types = readr::cols(.default = "c")) |>
  # split event_id (shiba-assigned event ID with values like SE_1, SE_2) into just event type acronym
  # split by underscore and keep first element (the event type)
  # remove retained intron events from matrix
  dplyr::filter(string
    stringr::str_starts(event_type, "RI_", negate = TRUE)
  )

message("RI events removed from merged matrix")

# join tables by pos_id column
annotated_matrix <- dplyr::left_join(sample_matrix,
  all_events_table,
  by = c("pos_id", "event_type")
) |>
  # make gene name column based off gene IDs
  dplyr::mutate(
    gene_name = gene_names[gene_id]
  ) |>
  # arrange descriptive columns to the front for ease of reading
  dplyr::relocate(
    event_type,
    label,
    gene_name,
    gene_id,
    pos_id
  )

message("cleaned matrix created")

## Write output ##
readr::write_tsv(annotated_matrix, out_matrix)
