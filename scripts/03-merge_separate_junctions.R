### merge junctions.bed files produced by separate shiba runs ###
# usage: Rscript 03-merge_separate_junctions.R

# Load library
library("optparse")

# Set up options to Rscript with optparse
option_list <-list(
  make_option(
    opt_str = "--junctions",
    type = "character",
    action = "store",
    help = "Specify input path for sample junction .bed file to merge"),
    # junctions bedfile path should look like results/group/shiba/sample/junctions/junctions.bed
  make_option(
    opt_str = "--output",
    type = "character",
    help = "Specify output path for merged junction counts bedfile.")
    # output should look like results/merged_shiba/timestamp/merged_junctions.bed
  )

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## Directories and files ##
# find the root-level repo directory
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# find the project directory (compendium-shiba-run)
base_dir <- here::here()

# define paths to files
# create output directory if it does not exist
out_dir <- dirname(opt$output) # output is a file so this probably needs to be updated
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## make directories if they dont exist ##
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Split space-separated input files into vector
junction_paths <- strsplit(opt$junctions, " ")[[1]]

# check that input junction bedfile exists - needs to be updated to act on a vector
if(!file.exists(junction_paths)) {stop("Please enter valid input file for --junctions")}

# Use filenames as sample names
sample_names <- basename(dirname(dirname(junction_paths)))

names(junction_paths) <- sample_names

###### this section below needs to be reworded to use the junction file paths as input

# name the separate junction file paths
names(junction_paths) <- names(samples)


## read in files and merge##
# read in junctions.bed files created from separate Shiba runs
merged_junctions <- purrr::map(junction_paths, \(file) {
    read.table(file,
               header = TRUE,
               sep="\t",
               stringsAsFactors=FALSE,
               quote="",
               # make sure columns are all the same class for merging
               colClasses = "character")

  }) |>
    # merge junctions tables from multiple samples
    # the resulting table separates junction counts from each sample by columns with the sample ID
    purrr::reduce(\(x, y) dplyr::full_join(x, y, by = c("ID", "start", "end", "chr")))

# convert NAs to 0
merged_junctions[is.na(merged_junctions)] <- 0

## Save merged junction counts as output
readr::write_tsv(merged_junctions, out_junctions)
