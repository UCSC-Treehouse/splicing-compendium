# Run this with: Rscript generate_experiment_file.R --input=../shiba-run-data/target/fastq/target_md5sum.txt --output=target_subset_pilot_shiba_experiment.tsv --group=target

# To run Shiba, a .tsv file grouping each bam file path into Reference or Alternative categories is required.
# To obtain a list of the TARGET input files, this script uses the md5 checksum document obtained from running `bash transfer-and-offload-target.sh checksums`

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
    opt_str = "--input",
    type = "character",
    help = "Specify input path for md5sum .txt file to obtain accession IDs from. The columns contain md5sum hashes and the file name."),
  make_option(
    opt_str = "--star_dir",
    type = "character",
    default = "../shiba-run-data/target/star-output",
    help = "Specify output file path. Output is a tab-separated .tsv listing the samples to run shiba on [default %default]"),
  # allow for toggling of number of samples in experiment.tsv
  make_option(
    c("--n_samples", "-n"),
    type = "integer",
    help = "number of samples to include"
  ),
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
input_file <- file.path(opt$input)

# check that input file exists
if(!file.exists(input_file)) {stop("Please enter valid input file for --input")}

# read in files
target_accessions <- readr::read_table(input_file, col_names = c("md5", "fastq_file"))

# create output directory if it does not exist
out_dir <- dirname(opt$output)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Construct input bam paths using accession numbers
# TARGET sample path should look like: `shiba-run-data/target/star-output/{accession}/Aligned.sortedByCoord.out.bam`
# Note: shiba-run-data is a symlink to the "/mnt/bulk" directory used in the current Snakefile.
# However, this symlink must be redifined for all OpenStack instances.
# This is done with the command: `ln -s /mnt/bulk shiba-run-data`

# Construct bam file path for each accession ID in the md5sum file

bam_file <- "Aligned.sortedByCoord.out.bam"

# construct bam paths for TARGET accessions
target_experiment <- target_accessions |>
  # we only need accession numbers from each metadata file, which is in the fastq_file column
  dplyr::select(fastq_file) |>
  dplyr::mutate(
    # split fastq file names by underscore to obtain accession numbers
    accession = stringr::str_split_i(fastq_file, "_", 1)) |>
  # each accession will be duplicated due to the paired fastqs, so obtain unique accessions
  dplyr::distinct(accession) |>
  # construct paths
  dplyr::mutate(
    # Make bam path column
    bam_path = file.path(opt$star_dir, accession, bam_file),
    # Make group label column
    group = opt$group,
    # Make column for labeling sequencing technology (short or long-read)
    technology = "short",
    # Make a column that numbers the comparison group using row number
    sample = accession
  ) |>
  # Select for only fields needed for experiment file
  dplyr::select(c(sample, bam_path, group, technology))

# subset experiment file for testing
if (!is.null(opt$n_samples)){
  target_experiment <- head(target_experiment, opt$n_samples)
}

# Write output files
readr::write_tsv(target_experiment, file = opt$output)
