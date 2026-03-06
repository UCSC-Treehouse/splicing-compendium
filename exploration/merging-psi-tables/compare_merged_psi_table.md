# Compare merged PSI table vs. Shiba PSI tables
Cindy Liang (celiang@ucsc.edu)
2026-03-05

## Set up

## Directories and files

``` r
# find the project directory
base_dir <- here::here()

# define the data directories
# splice results directory
results_dir <- file.path(base_dir, "results")

# merged shiba PSI table of two samples
merged_psi_table_dir <- file.path(results_dir, "merged_shiba_results_tables")

# merged PSI file of two samples
merged_psi_file <- file.path(merged_psi_table_dir, "merged_psi_table.tsv")

# general shiba results paths
shiba_gtex_path <- file.path(results_dir, "gtex/shiba")

# path to splice results of shiba run done on two samples at once
merged_splice_results_path <- "two_samples_together/results/splicing"

# define list of PSI results files corresponding to event types quantified by Shiba bulk analysis
psi_files <- c(
  se = "PSI_SE.txt",
  afe = "PSI_AFE.txt",
  ale = "PSI_ALE.txt",
  five = "PSI_FIVE.txt",
  three = "PSI_THREE.txt",
  mse = "PSI_MSE.txt",
  mxe = "PSI_MXE.txt",
  ri = "PSI_RI.txt"
)

# construct psi table paths of shiba results of two samples run together
psi_paths <-file.path(shiba_gtex_path, merged_splice_results_path, psi_files)
names(psi_paths) <- names(psi_files)
```

Read in files

``` r
# read in manually combined splice table
merged_psi_table <- readr::read_tsv(merged_psi_file)
```

    Rows: 352888 Columns: 6
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): event_type, pos_id, gene_id, label
    dbl (2): SRR601500_PSI, SRR604528_PSI

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read in splice results to compare against manually combined results

splice_results <- psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
      dplyr::mutate(across(contains("_PSI"), as.numeric))
  })

# combine psi values of samples run together into one dataframe to compare against manually combined dataframe
combined_splice_results_ref <- purrr::list_rbind(splice_results, names_to = "event_type")
```

Use R setdiff For annotated sites are they all the same If differences
are in the novel sites, how different are they

``` r
# produce dataframe of elements in combined splice results reference dataframe but not in merged psi df
unique_reference <- probs::setdiff(combined_splice_results_ref, merged_psi_table)

# produce dataframe of elements in merged splice results reference dataframe but not in the reference df
unique_merged <- probs::setdiff(merged_psi_table, combined_splice_results_ref)
```

Start by counting the number of events that are annotated
vs. unannotated that are different

``` r
table(unique_merged$label)
```


      annotated unannotated 
           5752        3954 

``` r
table(unique_reference$label)
```


      annotated unannotated 
           4255        3870 

Weirdly, several events look like they have the same info for both
tables, e.g.  SE@chr10@119580471-119580846@119580010-119581922
