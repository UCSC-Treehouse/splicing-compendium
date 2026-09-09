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

### Read in files and directories ###

## directories ##
# find the root-level repo directory so we can access the other files
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
# define the data directories
results_dir <- file.path(repo_root, "results")
merged_shiba_dir <- file.path(results_dir, "merged_shiba")
target_persample_pilot_results_dir <- file.path(merged_shiba_dir, "v1.1_reuse_table_target_pilot")
sample_psi_dir <- file.path(target_persample_pilot_results_dir, "sample_psi")
one_sample_results_dir <- file.path(sample_psi_dir, "SRR1559043")

# merged PSI results
merged_matrix_dir <- file.path(target_persample_pilot_results_dir, "merged_persample_psi")

## files ##

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

### Read in PSI tables ###

event_tables <- event_psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types = readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label)
  })

# combine psi tables into one, with additional column labeling event type
psi_table <- purrr::list_rbind(psi_results, names_to = "event_type")
