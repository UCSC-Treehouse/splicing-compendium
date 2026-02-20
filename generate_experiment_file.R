# Run this with: Rscript generate_experiment_file.R --sample=sample_accession --bam=../path/to/bam --output=shiba_experiment.tsv --group=sample_group

# To run Shiba, a .tsv file grouping each bam file path into Reference or Alternative categories is required.
# This script generates one experiment.tsv per input sample .bam file

###########

# load library
library("optparse")

# Set up options to Rscript with optparse
option_list <-list(
  make_option(
    opt_str = "--group",
    type = "character",
    help = "Specify name of sample group for analysis (e.g. TARGET vs GTEx)"),
  make_option(
    opt_str = "--bam",
    type = "character",
    help = "Specify input path for sample .bam file to obtain accession IDs from"),
  make_option(
    opt_str = "--sample",
    type = "character",
    help = "Specify accession ID of bam"),
  make_option(
    opt_str = "--output",
    type = "character",
    help = "Specify output path for a tsv file which contain sample information for running shiba.")
  )

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# find the root-level repo directory
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# find the project directory (compendium-shiba-run)
base_dir <- here::here()

# define paths to files
input_file <- file.path(opt$bam)

# check that input bam exists
if(!file.exists(input_file)) {stop("Please enter valid input file for --input")}

# create output directory if it does not exist
out_dir <- dirname(opt$output)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Construct input bam paths using accession numbers
# TARGET sample path should look like: `shiba-run-data/target/star-output/{accession}/Aligned.sortedByCoord.out.bam`
# Note: shiba-run-data is a symlink to the "/mnt/bulk" directory used in the current Snakefile.
# However, this symlink must be redifined for all OpenStack instances.
# This is done with the command: `ln -s /mnt/bulk shiba-run-data`

# construct experiment.tsv for one sample at a time
experiment_df <- data.frame(
    sample = opt$sample,
    bam_path = opt$bam,
    group = opt$group,
    technology = "short",
    stringsAsFactors = FALSE
  )

# Write output files
readr::write_tsv(experiment_df, file = opt$output)
