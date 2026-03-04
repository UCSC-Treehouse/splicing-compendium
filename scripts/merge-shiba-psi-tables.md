# Merge Shiba PSI tables
Cindy Liang (celiang@ucsc.edu)
2026-03-03

Shiba is run one sample at a time using our snakemake workflow. We then
merge PSI tables produced for each sample.

## Set up

## Directories and files

``` r
# samples list
# noting here that it will become terrible to paste infinite sample names in this list
# so maybe we want to ultimately have the script read in a file with all the samples for the big run
samples <- c("SRR601500", "SRR604528")

# find the project directory
base_dir <- here::here()

# define the data directories
# splice results directory
results_dir <- file.path(base_dir, "results")

# general shiba results paths
shiba_gtex_path <- "gtex/shiba"
splice_results_path <- "results/splicing"

# define list of PSI results files corresponding to event types quantified by Shiba bulk analysis
psi_files <- c(
  se = "PSI_SE.txt.gz",
  afe = "PSI_AFE.txt.gz",
  ale = "PSI_ALE.txt.gz",
  five = "PSI_FIVE.txt.gz",
  three = "PSI_THREE.txt.gz",
  mse = "PSI_MSE.txt.gz",
  mxe = "PSI_MXE.txt.gz",
  ri = "PSI_RI.txt.gz"
)

# Build a nested list of PSI result paths, where the outer list is the sample IDs and the inner list is the splice event PSI matrices for each sample
sample_psi_result_paths <- lapply(samples, function(sample) {
  
  # make path for splice results directory for each sample in samples list
  sample_splice_results_dir <- file.path(
    results_dir,
    shiba_gtex_path,
    sample,
    splice_results_path
  )
  
  # construct path to each splice event type PSI result
  paths <- file.path(sample_splice_results_dir, psi_files)
  
  # name PSI result paths by event type
  names(paths) <- names(psi_files)
  
  # save the result to variable
  return(paths)
})

# name the outer list of paths by sample ID
names(sample_psi_result_paths) <- samples

# output directory of results
out_dir <- file.path(results_dir, "merged_shiba_results_tables")

# output merged file for test samples
out_file <- file.path(out_dir, "test_merged_psi_table")

# create output directory if it does not already exist
dir.create(out_dir, showWarnings = FALSE)
```

Read in files

``` r
# Read in files in the nested list of samples' PSI values
# should return a nested list, where the outer list is each sample and the inner lists are data frames for each splice event type
splice_results <- sample_psi_result_paths |>
  # iterate over each sample''s PSI file paths
  purrr::map(\(sample_paths) {
    # iterate over each splice event type's file path corresponding to a sample
    purrr::map(sample_paths, \(path) {
      # read each event types' PSI results
      readr::read_tsv(path, col_types = readr::cols(.default = "c")) |>
        # select for columns that will be used downstream in analysis so all dataframes have uniform column names
        dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
        # turn PSI values numeric to ensure we can do math on them
        dplyr::mutate(across(contains("_PSI"), as.numeric)) 
    }
    )
  })
```

Combine PSI dataframes for all samples into one

``` r
# combine individual event types' PSI tables into one to obtain a named list of all splice events per sample
sample_splice_list <- purrr::map(splice_results, dplyr::bind_rows, .id = "event_type")

# combine each samples' PSI tables to obtain one PSI table for all samples
combined_splice_df <- dplyr::bind_rows(sample_splice_list, .id = "sample")
```

Low hanging fruit: Is there anything shared at all between the two
samples? If not, we might need to change the actual shiba run to combine
the junctions file upstram of PSI calculation

``` r
combined_splice_df |>
  dplyr::filter(!is.na(SRR601500_PSI) & !is.na(SRR604528_PSI))
```

| sample | event_type | pos_id | gene_id | label | SRR601500_PSI | SRR604528_PSI |
|:-------|:-----------|:-------|:--------|:------|--------------:|--------------:|

Looks like we need to go in and fiddle with the actual Shiba run to get
consistent position ID assignments across samples
