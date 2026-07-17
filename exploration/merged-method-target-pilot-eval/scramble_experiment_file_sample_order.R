# title: "Scramble experiment.tsv sample order"
# author: "Cindy Liang (celiang@ucsc.edu)"

# Usage: Rscript exploration/merged-method-target-pilot-eval/scramble_experiment_file_sample_order.R --experiment_file experiment_scratch.tsv

# Background

# This R script reads in `experiment.tsv` file that is used as input to the shiba workflow and scrambles the order of sample information (all rows after the header rows).
# The results of this scrambled file will be used to create a merged GTF containing all samples' transcript models with the unaltered Shiba v0.8.1 `bam2gtf.py` and splice event coordinates with `gtf2events.py`.
# The position IDs of the event coordinates produced from the order-scrambled `experiment.tsv` will be compared with the position IDs created from an unscrambled experiment.tsv with the same number of threads, 
# to understand the impact of GTF merge order on splice event IDs. 

# --- Load libraries ---
library("optparse")

# --- Set up options to Rscript with optparse ---
option_list <-list(
  make_option(
    opt_str = "--experiment_file",
    type = "character",
    default = "experiment.tsv",
    action = "store",
    help = "Name of experiment.tsv file to scramble"),
  
  make_option(
    opt_str = "--output",
    type = "character",
    default = "scrambled_experiment.tsv",
    action = "store",
    help = "Name of output scrambled experiment file")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# --- Read in file paths ---
# find the root-level repo directory so files outside this script's dir can be accessed
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# -- Parent directories --
exploration_dir <- file.path(repo_root, "exploration")

method_tests_dir <- file.path(exploration_dir, "merged_shiba_target_pilot")
splice_event_test_dir <- file.path(method_tests_dir, "shiba_gtf_to_event_pilot_tests")
# -- Files --
# target pilot experiment.tsv
experiment_file <- file.path(splice_event_test_dir , opt$experiment_file)
# output experiment.tsv with scrambled row (sample) order
out_file <- file.path(splice_event_test_dir, opt$output)

# --- Read in files ---
# set seed for reproducibility
set.seed(1)

experiment_table <- readr::read_tsv(experiment_file, col_types = "c")

# scramble order of rows in dataframe
experiment_table <- experiment_table |> dplyr::slice_sample(n = Inf)

# --- Write output ---
readr::write_tsv(experiment_table, file = out_file)