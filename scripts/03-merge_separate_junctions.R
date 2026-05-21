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
    help = "Comma-separated list of input file paths to merge"),

  make_option(
    opt_str = "--output",
    type = "character",
    help = "Specify output path for merged junction counts bedfile.")
  )

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## File paths ##
# Read in junctions bed list into a vector
# each element of list looks like results/target/shiba/SRR4376025/junctions/junctions.bed
junction_paths <- strsplit(opt$junctions, ",")[[1]]

## read in files and merge ##
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
readr::write_tsv(merged_junctions, opt$output)
