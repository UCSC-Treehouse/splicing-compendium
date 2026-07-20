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
# read in junctions.bed files created from separate Shiba runs
merged_junctions <- junction_paths |>
  purrr::map(\(file) {
    # read in separate bedfiles for each sample
    junctions <- read_csv_duckdb(
      file,
      options = list(
        delim = "\t",
        types = list(c(
          chr = "VARCHAR",
          start = "INTEGER",
          end = "INTEGER",
          ID = "VARCHAR"
        ))
      )
    )

    # return junctions object as duckdplyr tibble
    junctions
  }) |>
  # merge junctions tables from multiple samples
  # the resulting table separates junction counts from each sample by columns with the sample ID
  purrr::reduce(
    \(x, y) {
      result <- full_join(
        x,
        y,
        by = c("chr", "start", "end", "ID")
      ) |>
        # compute results to make sure we don't get a huge query plan
        duckplyr::compute()

      result
    }
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
