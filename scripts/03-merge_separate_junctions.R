### merge junctions.bed files produced by separate shiba runs ###
# usage: Rscript 03-merge_separate_junctions.R

# avoid dplyr warnings
options(conflicts.policy = list(warn = FALSE))

# Load Packages
library(optparse)
library(dplyr)
library(duckplyr)

# Set up options to Rscript with optparse
option_list <- list(
  make_option(
    opt_str = "--junctions",
    type = "character",
    action = "store",
    help = "Input directory of deduplicated junction bedfiles"
  ),

  make_option(
    opt_str = "--output",
    type = "character",
    help = "Specify output path for merged junction counts bedfile."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## File paths ##
# Read in deduplicated junctions bed paths from temp dir into a vector
# each element of list looks like tempdir/1_junctions.bed
junction_paths <- list.files(
  path = opt$junctions,
  pattern = ".bed",
  full.names = TRUE
)

## read in files and merge ##
# have duckdb read all files into a single table
raw <- read_csv_duckdb(
  junction_paths,
  options = list(
    delim = "\t",
    union_by_name = TRUE,
    types = list(c(
      chr = "VARCHAR",
      start = "INTEGER",
      end = "INTEGER",
      ID = "VARCHAR"
    ))
  )
)

# get the value column names, which are already sample-specific
value_cols <- setdiff(names(raw), c("chr", "start", "end", "ID"))

# For each sample value, we will take the max for that column across all rows with the same junction ID.
# after removing NAs, this will be unique!
merged_junctions <- raw |>
  summarise(
    across(all_of(value_cols), \(x) max(x, na.rm = TRUE)),
    .by = c(chr, start, end, ID)
  ) |>
  # materialize duckdb query into a tibble
  as_tibble()

# Check if junction IDs are duplicated
if (any(duplicated(merged_junctions$ID))) {
  dup_rows <- merged_junctions[duplicated(merged_junctions$ID), ]
  # quit and send error message about duplicates
  stop(
    paste0("Found duplicate junction ID: ", "\n"),
    readr::format_tsv(dup_rows)
  )
}

# convert junctions not found in a sample from NA to 0
merged_junctions[is.na(merged_junctions)] <- 0

## Save merged junction counts as output
readr::write_tsv(merged_junctions, opt$output)
