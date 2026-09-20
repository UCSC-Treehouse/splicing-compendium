### Introduction ###

# This notebook performs differential analysis on splice event usage (PSI) values generated from the splice compendium workflow.
# The analyses included in this notebook consist of:
#
#   1. Statistical test for differential splicing: A Wilcoxon ranksum test is used to compare the PSI value distributions for splice events used in a reference and query group of samples.
# P-values are corrected for multiple testing using the Benjamini-Hochberg method.
#
# 2. Calculating the magnitude of differential splicing: The magnitude of differential splicing is calculated for each splice event by taking the difference between the median PSI value between the query and reference group.
#
# 3. Visualizing splice events by significance and magnitude of differential splicing: Volcano plots are created to plot the adjusted p-value and dPSI values of each splice event type between the query and reference groups.
# One Volcano plot is created for all splice event types, and separate volcano plots are provided to split splice event types by category.
#
# 4. Filtering splice event types by significance score and ranking by magnitude of differential splicing: splice event position IDs, their dPSI, and adjusted p-values are finally exported as a TSV based on user-defined adjusted p-value and dPSI filters.
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
data_dir <- file.path("data")

# metadata file
target_metadata_file <- file.path(data_dir, "SraRunTable-TARGET.csv")
# target subset combined psi results
combined_psi_file <- file.path(data_dir, "combined_psi_results.rds")

# output
wilcox_results_file <- file.path(data_dir, "wilcox_test_results.tsv")
diff_splice_results <- file.path(data_dir, "target_pilot_diffsplice.tsv")
