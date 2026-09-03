# compare merged sample psi matrix target pilot
Cindy Liang (celiang@ucsc.edu)
2026-09-02

## Introduction

The Treehouse Splice Compendium workflow has been modified to allow for
all 2581 samples of the compendium to be run given the Pheonix cluster’s
memory limits. These changes consist of taking the merged junctions
bedfile containing counts and junction positions of all samples,
splitting it into one file per sample so that each per-sample bedfile
contains junction positions from all junctions, but only quantifications
for one sample, and calculating per-sample PSI values using this
per-sample junction count with event coordinate positions defined from
the GTF made from all samples.

Next, to create the PSI tables for compendium users, these per-sample
PSI tables are merged into one table per Shiba output. The resulting
merged tables are:


    se = "PSI_SE.txt",
    afe = "PSI_AFE.txt",
    ale = "PSI_ALE.txt",
    five = "PSI_FIVE.txt", 
    three = "PSI_THREE.txt",
    mse = "PSI_MSE.txt",
    mxe = "PSI_MXE.txt",
    ri = "PSI_RI.txt",
    matrix = "PSI_matrix_sample.txt"

Where “PSI\_\[EVENT\]” refers to PSI quantifications for each event
type, with each row being a unique event’s position ID and each column
being event information such as gene ID assigned to the event, position
coordinates of features such as introns and exons used to define that
event, annotation status of that event, junction counts for each of
these features for each sample, and PSI values corresponding to that
event per sample.

“PSI_matrix_sample.txt” corresponds to the sample-level event matrix.
Splice events for all event types are present in this matrix. Each row
is a unique splice event. In this martix, columns include the position
ID (coordinates) of the splice event, a Shiba-assigned “event_ID”
consisting of the event type in acronyms, f ollowed by a number, and the
PSI values for each sample.

The goal of this notebook is to check if the merged tables generated
from the target pilot pilot samples are the same as the canonical PSI
tables.

## Setup

``` r
## Directories ##
# find project root to access separate results dir
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# result directories
results_dir <- file.path(repo_root, "results")
compendium_results_dir <- file.path(results_dir, "merged_shiba")
version_dir <- file.path(compendium_results_dir, "v1.1_reuse_table_target_pilot")
psi_dir <- file.path(version_dir, "merged_persample_psi")

# canonical shiba psi matrix dir
exploration_dir <- file.path(repo_root, "exploration")
merge_exploration_dir <- file.path(exploration_dir, "merging-psi-tables")
target_pilot_dir <- file.path(merge_exploration_dir, "target_pilot")
# directory of shiba results produced from the same shiba run
combined_dir <- file.path(target_pilot_dir, "shiba_combined", "results")
# psi table directory
combined_splice_results_dir <- file.path(combined_dir, "splicing")

## Files ##

merged_file_list <- c(
  se = "PSI_SE.txt.gz",
  afe = "PSI_AFE.txt.gz",
  ale = "PSI_ALE.txt.gz",
  five = "PSI_FIVE.tx.gz",
  three = "PSI_THREE.txt.gz",
  mse = "PSI_MSE.txt.gz",
  mxe = "PSI_MXE.txt.gz",
  ri = "PSI_RI.txt.gz",
  matrix = "PSI_matrix_sample.txt.gz"
)

combined_file_list <- c(
  se = "PSI_SE.txt",
  afe = "PSI_AFE.txt",
  ale = "PSI_ALE.txt",
  five = "PSI_FIVE.tx",
  three = "PSI_THREE.txt",
  mse = "PSI_MSE.txt",
  mxe = "PSI_MXE.txt",
  ri = "PSI_RI.txt",
  matrix = "PSI_matrix_sample.txt"
)

# Construct paths to merged and combined sets of files

merged_files <- file.path(psi_dir, merged_file_list)
names(merged_files) <- names(merged_file_list)
combined_files <- file.path(combined_splice_results_dir, combined_file_list)
names(combined_files) <- names(combined_file_list)
```

Based on previous comparison work done in `exploration`, we know that
retained intron events are different (different PSI values, even for
events with the same position ID between splice compendium workflow
method results and canonical Shiba results). So RI event types are
filtered out in all matrices prior to comparison.

## Compare PSI sample matrices

Read in matrices to compare and filter out ri events

``` r
merged_matrix_df <- readr::read_tsv(merged_files[["matrix"]], col_types = readr::cols(.default = "c")) |>
    dplyr::filter(!stringr::str_detect(event_id,"RI_")) |>
  # sort by pos id
  dplyr::arrange(pos_id)

combined_matrix_df <- readr::read_tsv(combined_files[["matrix"]], col_types = readr::cols(.default = "c")) |>
    dplyr::filter(!stringr::str_detect(event_id,"RI_")) |>
  # sort by pos id
  dplyr::arrange(pos_id)
```

Check if matrices are identical - matrices are not identical on
Shiba-assigned event IDs

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
PSI sample matrix, except for Shiba-assigned event IDs

Remove PSI matrices from memory prior to comparing other event tables

``` r
# removed matrices that won't be compared anymore
rm(merged_matrix_df)
rm(combined_matrix_df)
rm(merged_matrix_df_no_id)
rm(combined_matrix_df_no_id)

# garbage collector
gc()
```

              used (Mb) gc trigger   (Mb) limit (Mb)  max used   (Mb)
    Ncells  831063 44.4    3633622  194.1         NA   4255102  227.3
    Vcells 3530031 27.0  278327615 2123.5      16384 347841300 2653.9
