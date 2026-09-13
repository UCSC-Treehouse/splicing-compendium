### merge junction bedfiles produced by separate Shiba runs ###
#
# Reads every *.bed junction count file in --junctions and combines them so that all samples use the same set of junctions.
# If a sample does not have a particular junction, it is treated as having zero reads for that junction.
#
# The sample name for each file is taken from the header of its final column
# (Shiba names that column after the sample), not from the filename.
#
# There are two output modes, depending on which output argument is given
#
# --output <file.bed>
#   - One merged wide table: chr, start, end, ID, then one count column per sample. Must end in .bed.
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


# Set up options to Rscript with optparse
option_list <- list(
  make_option(
    "--threads",
    type = "integer",
    default = parallel::detectCores(),
    help = "DuckDB worker threads"
    ),

  make_option(
    opt_str = "--junctions",
    type = "character",
    action = "store",
    help = "Input directory of deduplicated junction bedfiles"
  ),

  make_option(
    opt_str = "--output",
    dest = "output_file",
    type = "character",
    help = paste(
      "Specify output path for a single merged junction counts bedfile",
      "(must end in .bed). Mutually exclusive with --output_dir."
    )
  ),

  make_option(
    opt_str = "--output_dir",
    type = "character",
    help = paste(
      "Specify output directory for per-sample junction counts bedfiles.",
      "One <sample>.bed is written per sample, all sharing the same",
      "union set of junctions (missing counts filled with 0).",
      "Mutually exclusive with --output."
    )
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## Validate output options ##
# exactly one of --output / --output_dir
if (!xor(is.null(opt$output_file), is.null(opt$output_dir))) {
  stop("Specify exactly one of --output (a .bed file) or --output_dir.")
}

output_merged <- !is.null(opt$output_file)

if (output_merged) {
  # merged mode: must be a .bed file path, not an existing directory
  if (!grepl("\\.bed$", opt$output_file, ignore.case = TRUE)) {
    stop("--output must be a file path ending in .bed: ", opt$output_file)
  }
  if (dir.exists(opt$output_file)) {
    stop(
      "--output is an existing directory, expected a .bed file: ",
      opt$output_file
    )
  }
} else {
  # separate mode: must be a directory (created if it does not exist)
  if (file.exists(opt$output_dir) && !dir.exists(opt$output_dir)) {
    stop("--output_dir exists but is not a directory: ", opt$output_dir)
  }
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(opt$output_dir)) {
    stop("Could not create --output_dir: ", opt$output_dir)
  }
}

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
      col_types = readr::cols(.default = "c")
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

# turn long table into a 'tbl' duckdb table reference
lj_tbl <- duckplyr::as_tbl(long_junctions)
# extract duckdb connection object for raw sql statement execution
con <- dbplyr::remote_con(lj_tbl)

# use threads for the rule
DBI::dbExecute(con, sprintf("PRAGMA threads=%d", opt$threads))


# the union of junctions seen in any sample
# (this replaces the pivot: every output file uses this same row set)
all_junctions <- long_junctions |>
  distinct(chr, start, end, ID) |>
  arrange(chr, start, end)

# Check if junction IDs are duplicated
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

nm <- as.character(dbplyr::remote_name(lj_tbl))
pivot_sql <- glue::glue_sql(
  'PIVOT {`nm`} ON sample USING coalesce(first(count), 0) GROUP BY chr, "start", "end", ID',
  .con = con
)

# materialize the pivot ONCE, reuse for every output file
wide <- tbl(con, sql(pivot_sql)) |>
  arrange(chr, start, end) |>
  compute() |>
  as_duckdb_tibble(prudence = "stingy")

if (output_merged) {
  ## Merged mode: one wide table with a column per sample ##

   wide |>
    duckplyr::compute_csv(
      opt$output_file,
      options = list(delim = "\t", header = TRUE))

} else {
  ## Separate mode: one count file per sample + junction bed ##

  wide |>
    select(chr, start, end, ID) |>
    duckplyr::compute_csv(
      file.path(
        opt$output_dir, "all_junctions.bed"
      ),
      options = list(delim = "\t", header = TRUE)
    )

  # create one count file per sample
  purrr::walk(sample_df$sample, \(sample_name) {
    wide |>
      select(all_of(sample_name)) |>
      duckplyr::compute_csv(
        file.path(opt$output_dir, paste0(sample_name, "_junction_counts.tsv")),
        options = list(delim = "\t", header = TRUE)
      )
  })
}
