### Introduction ###

# This notebook performs differential analysis on splice event usage (PSI) values generated from the splice compendium workflow.
# The analyses included in this notebook consist of:
#
# 1. Statistical test for differential splicing: A Wilcoxon ranksum test is used to compare the PSI value distributions for splice events used in a reference and query group of samples.
# P-values are corrected for multiple testing using the Benjamini-Hochberg method.
#
# 2. Calculating the magnitude of differential splicing: The magnitude of differential splicing is calculated for each splice event by taking the difference between the median PSI value between the query and reference group.
#
# 3. Filtering splice event types by significance score and ranking by magnitude of differential splicing: splice event position IDs, their dPSI, and adjusted p-values are finally exported as a TSV based on user-defined adjusted p-value and dPSI filters.
# These splice events are ranked in decreasing dPSI magnitude.

### Setup ###

## load libraries
suppressPackageStartupMessages({
  library(optparse)
})

## read in options
option_list <- list(
  make_option(
    opt_str = "--reference_group",
    type = "character",
    action = "store",
    required = TRUE,
    help = "name of tissue type to serve as reference group for differential splicing analysis.
    see metadata file or notebooks/compendium_summary/compendium_v1_summary.md for all possible tissue types in compendium."
  ),
  make_option(
    opt_str = "--query_group",
    type = "character",
    action = "store",
    required = TRUE,
    help = "name of tissue type to serve as query group for differential splicing analysis.
    see metadata file or notebooks/compendium_summary/compendium_v1_summary.md for all possible tissue types in compendium."
  ),
  make_option(
    opt_str = "--min_samples",
    type = "character",
    action = "store",
    required = FALSE,
    default = 5,
    help = "minimum number of samples in reference or query group that must have quantified (non-NA) PSI values for differential splicing to be calculated for an event.
    default = 5"
  ),
  make_option(
    opt_str = "--metadata_file",
    type = "character",
    action = "store",
    required = TRUE,
    help = "file name of metadata with tissue/disease type information"
  ),
  make_option(
    opt_str = "--in_psi",
    type = "character",
    action = "store",
    required = TRUE,
    help = "file name of PSI matrix with gene names and sample accession IDs to calculate differential splicing"
  ),
  make_option(
    opt_str = "--out_file",
    type = "character",
    action = "store",
    required = TRUE,
    help = "file name of output dPSI table with p values"
  ),
  make_option(
    opt_str = "--padj",
    type = "character",
    action = "store",
    default = 0.05,
    required = FALSE,
    help = "adjusted p-value threshold for defining significantly differentially spliced events.
    default = 0.05"
  ),
  make_option(
    opt_str = "--dpsi",
    type = "character",
    action = "store",
    default = 0.15,
    required = FALSE,
    help = "dPSI threshold for defining differentially spliced events.
    default = 0.15"
  )
)

## directories and files ##
# however since users can download both metadata and psi input from zenodo, they may not be in metadata or results dirs
# should I instead have user supply the full path to their input files?
# or explain in script that the files have to be in specific directories?
# directories
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
metadata_dir <- file.path(repo_root, "metadata")
results_dir <- file.path(repo_root, "results")

# files
# metadata file
metadata_file <- file.path(metadata_dir, opt$metadata_file)
# target subset combined psi results
# maybe set version (splice_compendium_v1) as option?
in_psi_file <- file.path(results_dir, "merged_shiba", "splice_compendium_v1", "merged_persample_psi", opt$in_psi)

# output
diff_splice_results <- file.path(results_dir, opt$out_file)

## read in files ##
# read in compendium metadata
metadata <- readr::read_csv(metadata_file, col_types = c(.default = "c")) |>
  # filter for ref and query groups to save memory
  # to preserve the origin of the metadata info, original tissue type columns are labeled target_study_name
  # or body_site (for gtex)
  # the only shared column with tissue type info is the plot_tissue_type column
  dplyr::filter(
    plot_tissue_type %in% c(opt$reference_group, opt$query_group)
  ) |>
  # select only tissue type and accession ID columns
  dplyr::select(plot_tissue_type, Run)

# read in psi results
long_psi_matrix <- readr::read_csv(in_psi_file, col_types = c(.default = "c")) |>
  # pivot longer to allow accession IDs (PSI column names) to be associated with tissue type from metadata
  tidyr::pivot_longer(
    cols = contains("_PSI"),
    names_pattern = "(.*)_PSI",
    names_to = "Run",
    values_to = "PSI"
  ) |>
  # filter to positions where there are sufficient samples for calculating differential splicing
  dplyr::group_by(pos_id) |>
  # count number of samples in each group with non-NA PSI values for a given splice event
  dplyr::mutate(
    ref_count = sum(plot_tissue_type == opt$reference_group),
    query_count = sum(plot_tissue_type == opt$query_group)
  ) |>
  dplyr::filter(
    ref_count >= opt$min_samples,
    query_count >= opt$min_samples
  ) |>
  # filter out events where all PSI values are identical (nothing different to compare)
  dplyr::filter(any(PSI != PSI[1])) |>
  dplyr::ungroup()

# associate accession IDs to tissue type info
dplyr::inner_join(
  filtered_target_metadata,
  long_psi_matrix,
  by = c("Run")
)

## calculate significance of differential splicing ##
# run wilcoxon ranksum test
wilcox_test_events <- cancer_type_psi_table_filtered |>
  dplyr::group_by(pos_id) |>
  rstatix::wilcox_test(
    PSI ~ group,
    ref.group = "ref"
  ) |>
  # add multiple testing adjustment
  dplyr::mutate(
    p_adj = p.adjust(p, method = "BH")
  )

### Calculate magnitude of differential splicing (dPSI) values
# Calculate median PSIs for each splice event, per group
dPSI_table <- cancer_type_psi_table_filtered |>
  dplyr::group_by(pos_id, group) |>
  dplyr::mutate(
    median_PSI =
      median(PSI)
  ) |>
  # Calculate dPSI values by taking the difference between the ref and query group medians
  dplyr::group_by(pos_id) |>
  dplyr::mutate(
    dPSI = unique(median_PSI[group == "query"]) - unique(median_PSI[group == "ref"])
  )

# Merge dPSI values with p values from wilcox differential splicing test
p_dPSI_table <- wilcox_test_events |>
  # select for informative columns from wilcox test
  dplyr::select(pos_id, p, p_adj) |>
  # merge wilcox test p values with dPSI table
  dplyr::full_join(dPSI_table, by = "pos_id")

readr::write_tsv(p_dPSI_table, diff_splice_results)
