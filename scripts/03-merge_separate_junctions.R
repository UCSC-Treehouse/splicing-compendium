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
    help = "Specify input path to manifest file of paths to sample junction .bed files to merge"),

  make_option(
    opt_str = "--output",
    type = "character",
    help = "Specify output path for merged junction counts bedfile.")
  )

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## Directories and files ##
# find the root-level repo directory
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# find the project directory (compendium-shiba-run)
base_dir <- here::here()

# Read in junctions bed manifest as dataframe
junctions_manifest_df <- read.delim(opt$junctions, sep="\t")

# extract junction bed paths from manifest into a vector
junction_paths <- junctions_manifest_df$junction_bed

# extract sample names from manifest into vector
sample_names <- junctions_manifest_df$sample

# assign sample names to junction paths
names(junction_paths) <- sample_names

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

# convert NAs from low junction counts to 0
merged_junctions[is.na(merged_junctions)] <- 0

## Save merged junction counts as output
readr::write_tsv(merged_junctions, out_junctions)
