# compare merged sample psi matrix target pilot
Cindy Liang (celiang@ucsc.edu)
2026-08-26

Check if merged PSI sample matrix on target pilot are the same as
canonical PSI sample matrix hopefully this will give a quick view on
whether my changes to the separate junctions code altered PSi
calculation

``` r
## Directories ##
# find project root to access separate results dir
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# result directories
results_dir <- file.path(repo_root, "results")
compendium_results_dir <- file.path(results_dir, "merged_shiba")
version_dir <- file.path(compendium_results_dir, "v1.1_reuse_table_target_pilot")
psi_dir <- file.path(version_dir, "merged_psi")

# canonical shiba psi matrix dir
exploration_dir <- file.path(repo_root, "exploration")
merge_exploration_dir <- file.path(exploration_dir, "merging-psi-tables")
target_pilot_dir <- file.path(merge_exploration_dir, "target_pilot")
# directory of shiba results produced from the same shiba run
combined_dir <- file.path(target_pilot_dir, "shiba_combined", "results")
# psi table directory
combined_splice_results_dir <- file.path(combined_dir, "splicing")

## Files ##

# merged sample psi matrix
merged_matrix_file <- file.path(psi_dir, "PSI_matrix_sample.txt")
# canonical shiba run sample psi matrix
combined_matrix_file <- file.path(combined_splice_results_dir, "PSI_matrix_sample.txt")
```

Read in matrices to compare and filter out ri events

``` r
merged_matrix_df <- readr::read_tsv(merged_matrix_file, col_types = readr::cols(.default = "c")) |>
    dplyr::filter(!stringr::str_detect(event_id,"RI_")) |>
  # sort by pos id
  dplyr::arrange(pos_id)

combined_matrix_df <- readr::read_tsv(combined_matrix_file, col_types = readr::cols(.default = "c")) |>
    dplyr::filter(!stringr::str_detect(event_id,"RI_")) |>
  # sort by pos id
  dplyr::arrange(pos_id)
```

check if matrices are identical

``` r
all.equal(merged_matrix_df, combined_matrix_df)
```

    [1] "Component \"event_id\": 71280 string mismatches"

Check if, excluding event IDs (SE_1, SE_2, etc), matrices are identical:

``` r
merged_matrix_df_no_id <- merged_matrix_df |>
  dplyr::select(!event_id) |>
  # sort by pos id
  dplyr::arrange(pos_id)

combined_matrix_df_no_id <- combined_matrix_df |>
  dplyr::select(!event_id) |>
  # sort by pos id
  dplyr::arrange(pos_id)

all.equal(merged_matrix_df_no_id, combined_matrix_df_no_id)
```

    [1] TRUE

Merged PSI matrices, excluding RI events, are the same as the canonical
PSI sample matrix, except for shiba-assigned event IDs
