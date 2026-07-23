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
junction_paths <- list.files(
  path = opt$junctions,
  pattern = ".bed",
  full.names = TRUE
)

# Get sample names from junction file headers
# with file path as names
sample_df <- junction_paths |>
  purrr::set_names() |>
  purrr::map_chr(\(path) {
    colnames <- readr::read_tsv(
      path,
      n_max = 0,
      col_types = c(.default = "c")
    ) |>
      names()
    # get the last value
    colnames[length(colnames)]
  }) |>
  tibble::enframe(
    name = "path",
    value = "sample"
  )


## read in files and merge ##
# have duckdb read all files into a single long table
long_junctions <- read_csv_duckdb(
  junction_paths,
  options = list(
    delim = "\t",
    union_by_name = TRUE,
    header = TRUE,
    # add a column with the filename for later pivot
    filename = TRUE,
    # override sample column name to "count" for consistent structure
    names = list(c("chr", "start", "end", "ID", "count")),
    types = list(c(
      chr = "VARCHAR",
      start = "INTEGER",
      end = "INTEGER",
      ID = "VARCHAR",
      count = "INTEGER"
    ))
  )
) |>
  # replace filename with sample name
  left_join(sample_df, by = c("filename" = "path")) |>
  select(!filename)


merged_junctions <- long_junctions |>
  as_tibble() |>
  tidyr::pivot_wider(
    names_from = sample,
    values_from = count,
    values_fill = 0
  ) |>
  arrange(chr, start, end)

# Check if junction IDs are duplicated
if (any(duplicated(merged_junctions$ID))) {
  dup_rows <- merged_junctions[duplicated(merged_junctions$ID), ]
  # quit and send error message about duplicates
  stop(
    paste0("Found duplicate junction ID: ", "\n"),
    readr::format_tsv(dup_rows)
  )
}

## Save merged junction counts as output
readr::write_tsv(merged_junctions, opt$output)
