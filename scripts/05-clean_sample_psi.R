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

### Load libraries ###
library(rtracklayer)

### Read in options ###
# options should consist of config version (directory where to grab the files from)

### Read in files and directories ###

## directories ##
# find the root-level repo directory so we can access the other files
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# define the data directories
references_dir <- file.path(repo_root, "references")
results_dir <- file.path(repo_root, "results")
merged_shiba_dir <- file.path(results_dir, "merged_shiba")
target_persample_pilot_results_dir <- file.path(merged_shiba_dir, "v1.1_reuse_table_target_pilot")
sample_psi_dir <- file.path(target_persample_pilot_results_dir, "sample_psi")
one_sample_results_dir <- file.path(sample_psi_dir, "SRR1559043")

# merged PSI results
merged_matrix_dir <- file.path(target_persample_pilot_results_dir, "merged_persample_psi")

## files ##

# gtf file for converting ensg id to gene names
gtf_file <- file.path(references_dir, "gencode.v47.primary_assembly.annotation.gtf")

# merged psi sample matrix
merged_matrix_file <- file.path(merged_matrix_dir, "PSI_matrix_sample.txt")

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
out_matrix <- file.path(merged_matrix_dir, "cleaned_psi_matrix.txt")

### Read in GTF and extract gene names ###
# import gtf
gtf <- rtracklayer::import(gtf_file, filter = list(type = "gene"))

# make a named vector of gene names to IDs
gene_names <- setNames(gtf$gene_name, gtf$gene_id)

# remove the gtf after names are extracted
rm(gtf)
# garbage collector to clear memory
gc()

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

# read in merged sample matrix
sample_matrix <- readr::read_tsv(merged_matrix_file, col_types = readr::cols(.default = "c")) |>
  # split event_id (shiba-assigned event ID with values like SE_1, SE_2) into just event type acronym
  # split by underscore and keep first element (the event type)
  # remove retained intron events from matrix
  dplyr::filter(string
    stringr::str_starts(event_type, "RI_", negate = TRUE)
  )

# join tables by pos_id column
annotated_matrix <- dplyr::left_join(sample_matrix,
  all_events_table,
  by = c("pos_id", "event_type")
) |>
  # make gene name column based off gene IDs
  dplyr::mutate(
    gene_name = gene_names[gene_id],
  ) |>
  # arrange descriptive columns to the front for ease of reading
  dplyr::relocate(
    event_type,
    label,
    gene_name,
    gene_id,
    pos_id
  )

## Write output ##
readr::write_tsv(annotated_matrix, out_matrix)
