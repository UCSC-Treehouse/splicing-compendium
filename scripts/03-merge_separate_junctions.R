### merge junction bedfiles produced by separate Shiba runs ###
#
# Reads every *.bed junction count file in --junctions and combines them
# so that all samples use the same set of junctions.
# If a sample does not have a particular junction, it is treated as having
# zero reads for that junction.
#
# The sample name for each file is taken from the header of its final column
# (Shiba names that column after the sample), not from the filename.
#
# There are two output modes, depending on which output argument is given
#
# --output <file.bed>
#   - One merged wide table: chr, start, end, ID, then one count column per sample.
#     Must end in .bed.
#
# --output_dir <dir>
#  - A unified junctions.bed file with all junctions present in any sample.
#  - One <sample>_junction_counts.tsv per sample, all within the specified output directory.
#
#
# usage:
#   Rscript 03-merge_separate_junctions.R --junctions junction_dir/ --output merged.bed
#   Rscript 03-merge_separate_junctions.R --junctions junction_dir/ --output_dir split/

# Load Packages
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(duckplyr)
  library(dbplyr)
})

## Define functions

# define a function to read junction files with duckdb
# count column is always given the name `"count"`
read_junctions_duckdb <- function(path) {
  junctions <- read_csv_duckdb(
    path,
    options = list(
      delim = "\t",
      header = TRUE,
      # any additional columns will still be read with their headers
      names = list(c("chr", "start", "end", "ID", "count")),
      types = list(c(
        chr = "VARCHAR",
        start = "INTEGER",
        end = "INTEGER",
        ID = "VARCHAR",
        count = "INTEGER"
      ))
    )
  )
}

merge_junctions <- function(merged_bed_df, junction_path) {
  # function to take a data frame representing a merged bed file
  # and a path to a junction file, and produce a merged junction file
  # The input merged_bed_df should have chr, start, end, and ID columns.

  # read in the new file, select only location columns
  junction_df <- read_junctions_duckdb(junction_path) |>
    select(chr, start, end, ID)

  # combine it with the old, keeping only distinct rows
  merged_df <- bind_rows(merged_bed_df, junction_df) |>
    distinct()
}


# Set up options to Rscript with optparse
option_list <- list(
  make_option(
    opt_str = "--junctions_dir",
    help = "Input directory of deduplicated junction bedfiles"
  ),

  make_option(
    opt_str = "--output_dir",
    help = paste(
      "Specify output directory for per-sample junction counts bedfiles.",
      "One <sample>.bed is written per sample, all sharing the same",
      "union set of junctions (missing counts filled with 0)."
    )
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# output_dir must be a directory (created if it does not exist)
if (file.exists(opt$output_dir) && !dir.exists(opt$output_dir)) {
  stop("--output_dir exists but is not a directory: ", opt$output_dir)
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(opt$output_dir)) {
  stop("Could not create --output_dir: ", opt$output_dir)
}

## File paths ##
junction_paths <- list.files(
  path = opt$junctions,
  pattern = "\\.bed$",
  full.names = TRUE
)

# Get sample names from junction file headers
# with file path as names
sample_df <- junction_paths |>
  purrr::map(\(path) {
    colnames <- readr::read_tsv(
      path,
      n_max = 0,
      col_types = readr::cols(.default = "c")
    ) |>
      names()

    # return a 1 line data frame
    data.frame(
      # the last colname
      sample = colnames[length(colnames)],
      path = path
    )
  }) |>
  purrr::list_rbind()


## read in files and generate bedfile

# start with an empty data frame to hold all junctions
all_junctions <- data.frame(
  chr = character(),
  start = integer(),
  end = integer(),
  ID = character()
)

# combine all junction files into a single sorted junction table
all_junctions <- sample_df$path |>
  purrr::reduce(merge_junctions, .init = all_junctions) |>
  arrange(chr, start, end)

# Check if any junction IDs are duplicated
dup_rows <- all_junctions |>
  summarise(n = n(), .by = ID) |>
  filter(n > 1) |>
  collect()

if (nrow(dup_rows) > 0) {
  # quit and send error message about duplicates
  stop(
    paste0("Found duplicate junction ID: ", "\n"),
    readr::format_tsv(dup_rows)
  )
}

# output the all junction file
duckplyr::compute_csv(
  all_junctions,
  file.path(
    opt$output_dir,
    "all_junctions.bed"
  ),
  options = list(delim = "\t", header = TRUE)
) |>
  # don't print to stdout, just write to file
  invisible()

# create output files by merging inputs with all junctions
# to fill in zeros
sample_df |>
  purrr::pwalk(\(sample, path) {
    sample_counts <- left_join(
      all_junctions,
      read_junctions_duckdb(path),
      by = join_by(chr, start, end, ID)
    ) |>
      # sorting can be lost in the join, so re-sort
      arrange(chr, start, end) |>
      # fill in zero counts for the count column and select just that column
      # 0L to ensure integer type, not double
      mutate(count = coalesce(count, 0L)) |>
      # rename counts with sample name
      select("{sample}" := count)

    # write the output table for the sample
    duckplyr::compute_csv(
      sample_counts,
      file.path(opt$output_dir, paste0(sample, "_junction_counts.tsv")),
      options = list(delim = "\t", header = TRUE)
    )
  })
