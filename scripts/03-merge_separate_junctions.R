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


# get the table name and connection for direct dbplyr SQL query
lj_tbl <- duckplyr::as_tbl(long_junctions)
con <- dbplyr::remote_con(lj_tbl) # duckplyr's DuckDB connection
nm <- as.character(dbplyr::remote_name(lj_tbl)) # the temp view's name

# create the SQL query to pivot the long table within DuckDB
pivot_sql <- glue::glue_sql(
  '
  PIVOT {`nm`}
  ON sample
  USING coalesce(first(count), 0)
  GROUP BY chr, "start", "end", ID
',
  .con = con
)

merged_junctions <- tbl(con, sql(pivot_sql)) |>
  arrange(chr, start, end) |>
  compute()

# Check if junction IDs are duplicated
dup_rows <- merged_junctions |>
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

## Save merged junction counts as output
merged_junctions |>
  as_duckdb_tibble(prudence = "stingy") |>
  duckplyr::compute_csv(
    opt$output,
    options = list(delim = "\t", header = TRUE)
  )
