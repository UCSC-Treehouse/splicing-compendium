### merge junctions.bed files produced by separate shiba runs ###
# usage: Rscript 03-merge_separate_junctions.R

# Load library
library("optparse")
library("duckplyr")

# Set up options to Rscript with optparse
option_list <-list(
  make_option(
    opt_str = "--junctions",
    type = "character",
    action = "store",
    help = "Input directory of deduplicated junction bedfiles"),

  make_option(
    opt_str = "--output",
    type = "character",
    help = "Specify output path for merged junction counts bedfile.")
  )

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## File paths ##
# Read in deduplicated junctions bed paths from temp dir into a vector
# each element of list looks like tempdir/1_junctions.bed
junction_paths <- list.files(path = opt$junctions, pattern = ".bed", full.names = TRUE)

## read in files and merge ##
# read in junctions.bed files created from separate Shiba runs
merged_junctions <- purrr::map(junction_paths, \(file) {
  # read in separate bedfiles for each sample
  junctions <- readr::read_tsv(
    file,
    # make sure columns are types that we expect
    col_types = readr::cols(
      .default = "d",
      ID = "c",
      chr = "c",
      start = "d",
      end = "d"
    )
  )
  # check if there are NAs in bedfile and print rows with NAs
  if (anyNA(junctions)) {
    # select rows in bedfile with NAs for printing
    na_bed_rows <- junctions[!complete.cases(junctions), ]
    stop(
      paste0("NAs found in bedfile", "\n"),
      readr::format_tsv(a_bed_rows)
    )
  }

  # return junctions object to merge
  junctions

  }) |>
  # merge junctions tables from multiple samples
  # the resulting table separates junction counts from each sample by columns with the sample ID
  purrr::reduce(
    \(x, y) duckplyr::full_join(
      x, y,
      by = c("chr", "start", "end", "ID")
    )
  )

# Check if junction IDs are duplicated
if (any(duplicated(merged_junctions$ID))) {
  dup_rows <- merged_junctions[duplicated(merged_junctions$ID),]
  # quit and send error message about duplicates
  stop(paste0("Found duplicate junction ID: ", "\n"),
       readr::format_tsv(dup_rows)
       )
    }

# convert junctions not found in a sample from NA to 0
merged_junctions[is.na(merged_junctions)] <- 0

## Save merged junction counts as output
readr::write_tsv(merged_junctions, opt$output)
