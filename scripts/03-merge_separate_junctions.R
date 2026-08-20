### create individual bed files with merged
# usage: Rscript 03-merge_separate_junctions.R

# Load Packages
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(duckplyr)
})

# Set up options to Rscript with optparse
option_list <- list(
  make_option(
    opt_str = "--junctions",
    type = "character",
    action = "store",
    help = "Input directory of deduplicated junction bedfiles"
  ),

  make_option(
    opt_str = "--output_dir",
    type = "character",
    help = paste(
      "Output directory for per-sample junction counts bedfiles.",
      "One <sample>.bed is written per sample, all sharing the same",
      "union set of junctions (missing counts filled with 0)."
    )
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## File paths ##
junction_paths <- list.files(
  path = opt$junctions,
  pattern = "\\.bed",
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


# the union of junctions seen in any sample
# (this replaces the pivot: every output file uses this same row set)
all_junctions <- long_junctions |>
  distinct(chr, start, end, ID)

# Check if junction IDs are duplicated
dup_rows <- all_junctions |>
  group_by(ID) |>
  filter(n() > 1) |>
  collect()

if (nrow(dup_rows) > 0) {
  # quit and send error message about duplicates
  stop(
    paste0("Found duplicate junction ID: ", "\n"),
    readr::format_tsv(dup_rows)
  )
}

## Save one junction counts bedfile per sample
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

purrr::walk(sample_df$sample, \(sample_name) {
  all_junctions |>
    left_join(
      long_junctions |>
        filter(sample == .env$sample_name) |>
        select(chr, start, end, ID, count),
      by = c("chr", "start", "end", "ID")
    ) |>
    # samples missing a junction get a zero count
    mutate(count = coalesce(count, 0L)) |>
    rename("{sample_name}" := count) |>
    arrange(chr, start, end) |>
    as_duckdb_tibble(prudence = "stingy") |>
    duckplyr::compute_csv(
      file.path(opt$output_dir, paste0(sample_name, ".bed")),
      options = list(delim = "\t", header = TRUE)
    )
})
