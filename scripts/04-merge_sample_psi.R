### merge per-sample PSI tables produced by merge_results.smk ###
#
# reads every PSI output within the sample_psi/ directory and combines them so that all samples' PSI results appear in each file.
# samples may have NA values in their PSI column if junction counts for the event < 10
#
# the sample name for each file is taken from the sample sheet passed to the workflow's configfile
#
# usage:
#   Rscript 04-merge_sample_psi.R --sample_sheet=sample_sheet.tsv --version_dir=v1.1.0_two_sample --output_dir=out_dir

# load packages
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(duckplyr)
})

# environment variables
Sys.setenv(DUCKPLYR_TEMP_DIR = "/data/tmp/celiang/fs")

# set up options to Rscript with optparse
option_list <- list(
  make_option(
    opt_str = "--sample_sheet",
    type = "character",
    action = "store",
    help = "path to sample sheet with list of sample IDs"
  ),
  make_option(
    opt_str = "--version_dir",
    type = "character",
    action = "store",
    help = "name of version directory of separate PSI results"
  ),
  make_option(
    opt_str = "--output_dir",
    type = "character",
    action = "store",
    help = "name of output directory of merged PSI results"
  ),
    make_option(
    opt_str = "--mode",
    type = "character",
    action = "store",
    help = "specify which PSI tables to merge (options: matrix, se, afe, ale, five, three, mse, mxe, ri)"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## Validate output options ##
# user must provide sample_sheet and version_dir to script
if ((is.null(opt$sample_sheet) || is.null(opt$version_dir)) || is.null(opt$output_dir)) {
  stop("Specify --sample_sheet, --version_dir, and --output_dir.")
}

### read in file paths ###

## Directories ##
# find project root to access separate results dir
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# result directories
results_dir <- file.path(repo_root, "results")
compendium_results_dir <- file.path(results_dir, "merged_shiba")
version_dir <- file.path(compendium_results_dir, opt$version_dir)
psi_dir <- file.path(version_dir, "sample_psi")

# merged psi table output dir
output_dir <- file.path(opt$output_dir)

# create output dir if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Files ##

# sample sheet file with sample names
samples_file <- file.path(opt$sample_sheet)

# read in sample sheet and create list of samples from samples column
samples <- duckplyr::read_csv_duckdb(samples_file,
options = list(
  delim = "\t",
  union_by_name = TRUE,
  header = TRUE
  )) |>
  dplyr::pull(sample)

# output files
out_file_list <- c(
  se = "PSI_SE.txt",
  afe = "PSI_AFE.txt",
  ale = "PSI_ALE.txt",
  five = "PSI_FIVE.txt",
  three = "PSI_THREE.txt",
  mse = "PSI_MSE.txt",
  mxe = "PSI_MXE.txt",
  ri = "PSI_RI.txt",
  matrix = "PSI_matrix_sample.txt"
)

# construct output file paths
out_paths <- file.path(
  output_dir,
  out_file_list
)

names(out_paths) <- names(out_file_list)

### define functions ###

# function for reading per-sample sample PSI matrices ("PSI_matrix_sample.txt")
read_sample_psi_matrix <- function(psi_path, samples, out_path, chunk_size = 200) {
  # returns a PSI matrix data frame of all samples in sample sheet, built by
  # pivoting samples in batches of `chunk_size`, then
  # joining the batches back together by column name (event_id, pos_id) so
  # that each row corresponds to a unique pos_id
  # the default chunk size is 200 but the user can change this throught he function call
  # e.g. read_sample_psi_matrix(psi_dir, samples, out_paths[["matrix"]], chunk_size = 20)

  # obtain duckdb connection object so query won't open a new DuckDB session
  con <- duckplyr:::get_default_duckdb_connection()

  # set preserve row order to false to save memory
  DBI::dbExecute(con, "SET preserve_insertion_order = false;")

  # create temp dir to store intermediate batched tables
  tmp_dir <- tempfile("psi_matrix_chunks_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # split samples into chunks of specified sample size
  chunks <- split(samples, ceiling(seq_along(samples) / chunk_size))
  # make empty vector to fill with paths to chunked parquet tables later
  chunk_paths <- character(length(chunks))

  # loop through chunks and write chunks into a wide parquet file
  for (i in seq_along(chunks)) {
    message("writing chunk ", i, " into parquet")

    # make a list of samples in each chunk
    chunk_samples <- chunks[[i]]
    # construct paths to PSI files in the chunk
    chunk_file_paths <- file.path(psi_path, chunk_samples, "PSI_matrix_sample.txt")
    # construct output paths to parquet file for the chunk
    chunk_out <- file.path(tmp_dir, sprintf("chunk_%03d.parquet", i))
    # add output path to list of chunk parquet paths
    chunk_paths[i] <- chunk_out

    message("Pivoting chunk ", i, " of ", length(chunks),
            " (", length(chunk_samples), " samples)")

    # first reshape chunk's wide matrix into long format prior to combining
    # cols should be event_id, pos_id, sample, psi
    # build a SQL select statement for each sample to read each file, match columns by name, and rename PSI column with sample name for merging
    selects <- purrr::map2_chr(chunk_file_paths, chunk_samples, \(file_path, sample) {
      sprintf(
        "SELECT event_id, pos_id, '%s' AS sample, \"%s\" AS psi FROM read_csv('%s', delim='\t',
        header=true, union_by_name=true)",
        sample, sample, file_path
      )
    })

    # merge the chunk's sample tables vertically into a long table
    union_sql <- paste(selects, collapse = " UNION ALL ")

    # pivot chunk's long table back to wide form
    # then write merged table directly from stream to parquet file
    chunk_query <- sprintf(
      "COPY (
        PIVOT (%s)
        ON sample
        USING first(psi)
        GROUP BY event_id, pos_id
      ) TO '%s' (FORMAT PARQUET)",
      union_sql, chunk_out
    )

    # execute chunk pivot query on the connection
    DBI::dbExecute(con, chunk_query)

  }

  # join all chunk-level wide tables together on event_id, pos_id
  # start from the first chunk and left-join the rest in sequence
  message("Joining ", length(chunk_paths), " chunks into final matrix")

  # build a FROM clause chaining each chunk file with JOIN ... USING, which
  # DuckDB resolves by matching column names (event_id, pos_id) rather than
  # needing per-table aliases for the join condition
  # name each parquet "c1" "c2"
  chunk_refs <- sprintf("read_parquet('%s') AS c%d", chunk_paths, seq_along(chunk_paths))

  from_sql <- Reduce(function(acc, tbl) {
    paste(acc, "JOIN", tbl, "USING (event_id, pos_id)")
  }, chunk_refs[-1], init = chunk_refs[1])

  # select all columns; since joins are done with USING, event_id/pos_id
  # appear once each in the result rather than once per chunk table
  final_query <- sprintf(
    "COPY (
       SELECT event_id, pos_id, * EXCLUDE (event_id, pos_id)
       FROM %s
     ) TO '%s' (DELIMITER '\t', HEADER)",
    from_sql, out_path
  )

  # execute final join query on the connection
  DBI::dbExecute(con, final_query)

}

## functions for reading and merging event tables ##
# event tables have different columns depending on the event type,
# so these functions operate on event tables that share columns

# skipped exon event tables
read_se_table <- function(psi_path, samples, out_path, chunk_size = 200) {
  # returns one data frame of skipped exon event types with all samples' psi values as separate columns
  # each row should correspond to a unique position ID (pos_id)
  # columns include those shared across all samples (event_id, pos_id, gene_id, exon, intron_a, intron_b, intron_c, strand, gene_name, label)
  # as well as sample-specific columns (sample1_junction_a, sample2_junction_b, sample1_junction_c, sample1_PSI)

  # obtain duckdb connection object
  con <- duckplyr:::get_default_duckdb_connection()

  # set preserve row order to false to save memory
  DBI::dbExecute(con, "SET preserve_insertion_order = false;")

  # create temp dir to store intermediate batched tables
  tmp_dir <- tempfile("psi_event_chunks_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # split samples into chunks
  chunks <- split(samples, ceiling(seq_along(samples) / chunk_size))
  # make empty vector to fill with paths to chunked parquet tables later
  chunk_paths <- character(length(chunks))

  # loop through chunks, merge sample tables in chunks, and write into a wide parquet file
  for (i in seq_along(chunks)) {
    message("writing chunk ", i, " into parquet")

    # make list of samples in each chunk
    chunk_samples <- chunks [[i]]
    # construct SE PSI table paths to files in the chunk
    chunk_file_paths <- file.path(psi_path, chunk_samples, "PSI_SE.txt")
    # construct parquet output paths
    chunk_out <- file.path(tmp_dir, sprintf("chunk_%03d.parquet", i))
    # add output path to list of chunk parquet paths
    chunk_paths[i] <- chunk_out

    message("merging samples in chunk ", i. "of ", length(chunks),
    " (", length(chunk_samples), " samples)")

    # read in each sample in chunk and reshape into long format
    # cols should be
    # event_id, pos_id, gene_id, exon, intron_a, intron_b, intron_c, strand, gene_name, label,
    # junction_a, junction_b, junction_c, PSI, sample

    # build SQL statement to read each file, rename sample-specific cols to be shared across samples, and match cols by name,
    selects <- purrr::map2_chr(chunk_file_paths, chunk_samples, \(file_path, sample) {

      # query to pivot tables longer
      chunk_query <- sprintf(
        "SELECT
        event_id,
        pos_id,
        gene_id,
        exon,
        intron_a,
        intron_b,
        intron_c,
        strand,
        gene_name,
        label,
        regexp_extract(sample_col, '^(.+)_(junction_[abc]|PSI)$', 1) AS sample,
        regexp_extract(sample_col, '_(junction_[abc]|PSI)$', 1) AS measure,
        value
    FROM your_table
    UNPIVOT (
        value FOR sample_col IN (
            COLUMNS(
                * EXCLUDE (
                    event_id,
                    pos_id,
                    gene_id,
                    exon,
                    intron_a,
                    intron_b,
                    intron_c,
                    strand,
                    gene_name,
                    label
                )
            )
        )
    )",
      union_sql,
      chunk_out
    )
    })

    # merge the chunk's sample tables vertically into a long table
    union_sql <- paste(selects, collapse = " UNION ALL ")

    # execute chunk pivot query on the connection
    DBI::dbExecute(con, chunk_query)

  }

}

### read in and merge psi tables ###

# print message when merging matrix
message("Merging PSI sample matrix")
# merge matrices and write merged matrix to output
read_sample_psi_matrix(psi_dir, samples, out_paths[["matrix"]], chunk_size = 200)

# # make list of event types to loop through
# event_types <- c("se", "afe", "ale", "five", "three", "mse", "mxe", "ri")

# # loop through event types
# for (event_type in event_types) {
#   # print message for log
#   message("Merging ", event_type, " PSI tables")

#   # create merged table object
#   assemble_event_table(sample_paths, out_file_list[[event_type]], out_paths[[event_type]])

# }
