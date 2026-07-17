# Merge TARGET Pilot Shiba PSI tables
Cindy Liang (celiang@ucsc.edu)
2026-04-05

Shiba is run one sample at a time using our snakemake workflow. We then
merge PSI tables produced for each sample.

## Set up

## Directories and files

``` r
# load list of TARGET pilot accession IDs from config file of the run done with pilot samples together
samples_file <- file.path("experiment.tsv")

# find project root to access separate results dir
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# splice results directory of separate runs
group_dir <- file.path(repo_root, "results", "target", "shiba")

# shiba splice results path
splice_results_path <- "results/splicing"

# output directory and file of merged results
out_dir <- file.path("merged_results")
out_file <- file.path(out_dir, "merged_psi_table.tsv")

# create output directory if it does not exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
```

Construct paths to separate TARGET pilot shiba output files

``` r
# read in experiment.tsv and create list of samples from samples column
samples <- readr::read_tsv(samples_file) |>
  dplyr::pull(sample) |>
  as.list()
```

    Rows: 88 Columns: 4
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): sample, bam_path, group, technology

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# make sample paths
sample_paths <- file.path(
  group_dir,
  samples,
  splice_results_path
)

# define list of PSI results files corresponding to event types quantified by Shiba bulk analysis
psi_file_list <- c(
  se = "PSI_SE.txt.gz",
  afe = "PSI_AFE.txt.gz",
  ale = "PSI_ALE.txt.gz",
  five = "PSI_FIVE.txt.gz",
  three = "PSI_THREE.txt.gz",
  mse = "PSI_MSE.txt.gz",
  mxe = "PSI_MXE.txt.gz",
  ri = "PSI_RI.txt.gz"
)
```

Define functions

``` r
# Make function for reading event types for a single sample
read_sample <- function(sample_id, sample_path, psi_files = psi_file_list) {
  # construct file path for each sample
  file_paths <- file.path(sample_path, psi_files)
  # name each PSI table file path by event type
  names(file_paths) <- names(psi_files)

  # read in files
  purrr::map(file_paths, \(file) {
      readr::read_tsv(file, col_types = readr::cols(.default = "c")) |>
      # select for columns that will be used downstream in analysis so all dataframes have uniform column names
      dplyr::select(pos_id, gene_id, label, contains("_PSI"))
    }) |>
      # merge individual event type PSI tables to get one PSI table per sample
      dplyr::bind_rows(.id = "event_type")
}
```

Read in files

``` r
# read in one samples' PSI table
# take in list of samples and sample paths as input
sample_psi_tables <- purrr::map2(samples, sample_paths, read_sample) |>
  # name the list of data frames the sample names
  purrr::set_names(samples)
```

Combine PSI dataframes for all samples into one

``` r
# combine individual event types' PSI tables into one
merged_psi <- sample_psi_tables |>
  purrr::reduce( \(x, y) {
    dplyr::full_join(x, y, by = c("pos_id", "gene_id", "label", "event_type") )
  })
```

Save output as a table

``` r
readr::write_tsv(merged_psi, out_file)
```
