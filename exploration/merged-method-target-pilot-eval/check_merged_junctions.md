# Check Merged vs. Combined TARGET Pilot Junction files
Cindy Liang (celiang@ucsc.edu)
2026-06-11

## Background

This notebook checks junctions.bed files from two Shiba run methods for
differences. The goal is to get a sense of how close files that go into
calculating PSI values for events are to each other in the merged and
combined shiba methods. Junctions files that are different would impact
PSI calculation, so I am mainly looking to see if mismatched junctions
are present and if present, how much of the data do they represent.

Files compared in this notebook are junction count bedfiles produced
from two Shiba runs from 88 TARGET pilot samples. “Merged” refers to
shiba splice analysis run with merge_results.smk, while “combined”
refers to shiba splice analysis run with shiba in a canonical way.

Analysis outline:

Check differences in full junctions files - Check how many rows
(junction IDs) are in each bedfile. Also calculate fraction of junctions
in the merged junctions table that are not in the combined table; the
fraction of junctions in the combined table that are missed in the
merged junctions table

Check differences in junctions with at least 10 junctions in
min_samples - Check if differences change if we filter for junctions
with higher counts in at least 5 samples. Filter merged and combined
junctions tables for junctions with counts of at least 10 in at least 5
samples. Calculate fraction of junctions in the merged junctions table
that are not in the combined table; the fraction of junctions in the
combined table that are missed in the merged junctions table

## Read in directories and files

Define directories and file paths

``` r
### Directories ###

# find the root-level repo directory so we can access the other files
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

## combined shiba run results on target pilot samples ##
exploration_dir <- file.path(repo_root, "exploration")

# merged table eval dir
exploration_eval_dir <- file.path(exploration_dir, "merged-method-target-pilot-eval")

# shiba results dir
# target pilot explroation dir
target_pilot_dir <- file.path(exploration_dir, "merging-psi-tables", "target_pilot")
combined_results_dir <- file.path(target_pilot_dir, "shiba_combined")
# junctions dir
combined_junctions_dir <- file.path(combined_results_dir, "junctions")

## merged shiba run results on target pilot samples ##
# shiba results dir
merged_results_dir <- file.path(repo_root, "results", "merged_shiba", "target_pilot")

### Files ###
## combined shiba target pilot files
target_pilot_sample_sheet <- file.path(target_pilot_dir, "experiment.tsv")
# deduplicated combined junctions file - made by running bash scripts/deduplicated_shiba_junctions.sh
deduplicated_combined_junctions_file <- file.path(combined_junctions_dir, "deduplicated_junctions.bed")

## merged shiba target pilot files
merged_junctions_file <- file.path(merged_results_dir, "merged_junctions.bed")
```

Read in files

``` r
# list of target pilot sample IDs
target_pilot_samples <- readr::read_tsv(target_pilot_sample_sheet) |>
  dplyr::pull(sample)
```

    Rows: 88 Columns: 4
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): sample, bam_path, group, technology

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# minimum number of samples to filter minimum number of junctions by
min_samples <- 5

# read in combined junctions table that have been deduplicated
deduplicated_combined_junctions <- readr::read_tsv(
  deduplicated_combined_junctions_file,
  col_types = readr::cols(
    .default = "d",
    ID = "c",
    chr = "c",
    start = "d",
    end = "d"
  )
)

# read in files generated from merged shiba run
merged_junctions <- readr::read_tsv(
  merged_junctions_file,
  col_types = readr::cols(
    .default = "d",
    ID = "c",
    chr = "c",
    start = "d",
    end = "d"
  )
)
```

## Check differences in junctions bedfiles

### Check differences in full junctions files

``` r
# print number of junctions in only in combined junctions bedfile
combined_only <- dplyr::setdiff(deduplicated_combined_junctions, merged_junctions)
nrow(combined_only)
```

    [1] 53732

``` r
# check fraction of junctions that are only in the combined junctions bedfile
frac_combined_only <- nrow(combined_only) / nrow(deduplicated_combined_junctions)
frac_combined_only
```

    [1] 0.008153533

``` r
# print number of junctions only in merged junctions bedfile
merged_only <- dplyr::setdiff(merged_junctions, deduplicated_combined_junctions)
nrow(merged_only)
```

    [1] 64211

``` r
# check fraction of junctions that are only in merged junctions
frac_merged_only <- nrow(merged_only) / nrow(merged_junctions)
frac_merged_only
```

    [1] 0.009728194

There are 10479 more junctions in the merged junction file compared to
the combined junction file.

0.8153533% of the combined method junctions are missed in the merged
junctions file.

0.9728194% of the merged method junctions are not in the combined
junctions file.

### Check differences in junctions with at least 10 junctions in min_samples

``` r
filtered_merged_junctions <- merged_junctions |> dplyr::filter(
  # count number of columns per row with junction count of at least 10
  rowSums(
    dplyr::across(
      dplyr::all_of(target_pilot_samples), ~ . >= 10
    )
    # filter for rows where at least 5 columns (samples) have counts of at least 10
  ) >= 5
  )

filtered_merged_junctions |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1559044 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1559052 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1559054 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1559075 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1559100 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR1559105 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR1559133 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 | SRR1559134 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chrM | 193 | 6201 | chrM:193-6201 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 6 | 0 | 8 | 30 | 9 | 28 | 54 | 50 | 82 | 19 | 24 | 7 | 0 | 13 | 16 | 4 | 49 | 2 | 1 | 11 | 24 | 24 | 41 | 0 | 11 | 1 | 28 | 23 | 7 | 0 | 0 | 91 | 37 | 34 | 0 | 1 | 25 | 21 | 7 | 12 | 17 | 66 | 46 | 72 | 21 | 0 | 17 | 88 | 5 | 126 | 1356 | 17 | 66 | 154 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chrM | 193 | 6404 | chrM:193-6404 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 10 | 11 | 16 | 3 | 3 | 2 | 0 | 1 | 2 | 0 | 9 | 1 | 0 | 2 | 0 | 2 | 3 | 0 | 0 | 0 | 3 | 1 | 1 | 0 | 0 | 15 | 2 | 5 | 0 | 1 | 3 | 2 | 1 | 0 | 0 | 3 | 0 | 7 | 2 | 0 | 10 | 9 | 1 | 27 | 205 | 20 | 15 | 15 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| GL000224.1 | 347 | 12972 | GL000224.1:347-12972 | 0 | 6 | 0 | 0 | 0 | 0 | 3 | 1 | 2 | 0 | 11 | 0 | 4 | 5 | 0 | 1 | 8 | 41 | 4 | 1 | 4 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 148 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 51 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 11 | 26 | 0 | 0 | 0 | 0 | 56 | 0 | 4 | 0 | 0 | 0 | 0 | 7 | 0 | 0 | 0 | 10 | 1 | 0 | 6 | 0 | 4 | 0 | 0 | 0 | 0 | 0 |
| chrM | 381 | 13855 | chrM:381-13855 | 0 | 0 | 0 | 0 | 0 | 0 | 105 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 152 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1167 | 0 | 0 | 0 | 0 | 2744 | 0 | 0 | 0 | 29 | 0 | 339 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chrM | 381 | 8811 | chrM:381-8811 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 41 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 182 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 77 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 60 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 20 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| KI270744.1 | 701 | 6071 | KI270744.1:701-6071 | 0 | 21 | 3 | 1 | 0 | 106 | 0 | 1 | 5 | 3 | 2 | 0 | 3 | 5 | 2 | 13 | 0 | 1 | 0 | 9 | 0 | 0 | 9 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 35 | 0 | 3 | 0 | 0 | 0 | 22 | 41 | 4 | 2 | 0 | 3 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 5 | 7 | 1 | 5 | 0 | 0 | 16 | 3 | 1 | 0 | 1 | 0 | 0 | 13 | 1 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 5 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

``` r
filtered_combined_junctions <- deduplicated_combined_junctions |> dplyr::filter(
  # count number of columns per row with junction count of at least 10
  rowSums(
    dplyr::across(
      dplyr::all_of(target_pilot_samples), ~ . >= 10
    )
    # filter for rows where at least 5 columns (samples) have counts of at least 10
  ) >= min_samples
  )

filtered_combined_junctions |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chrM | 193 | 6201 | chrM:193-6201 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 6 | 8 | 30 | 9 | 28 | 54 | 50 | 82 | 19 | 24 | 7 | 13 | 16 | 4 | 49 | 2 | 1 | 11 | 24 | 24 | 41 | 11 | 1 | 28 | 23 | 7 | 0 | 0 | 91 | 37 | 34 | 1 | 25 | 21 | 7 | 12 | 17 | 66 | 46 | 72 | 21 | 17 | 88 | 5 | 126 | 1356 | 17 | 66 | 154 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chrM | 193 | 6404 | chrM:193-6404 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 10 | 11 | 16 | 3 | 3 | 2 | 1 | 2 | 0 | 9 | 1 | 0 | 2 | 0 | 2 | 3 | 0 | 0 | 3 | 1 | 1 | 0 | 0 | 15 | 2 | 5 | 1 | 3 | 2 | 1 | 0 | 0 | 3 | 0 | 7 | 2 | 10 | 9 | 1 | 27 | 205 | 20 | 15 | 15 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| GL000224.1 | 347 | 12972 | GL000224.1:347-12972 | 0 | 0 | 0 | 148 | 51 | 0 | 56 | 10 | 0 | 6 | 0 | 0 | 0 | 0 | 3 | 1 | 2 | 0 | 11 | 4 | 5 | 0 | 1 | 8 | 41 | 4 | 1 | 4 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 11 | 26 | 0 | 0 | 0 | 0 | 0 | 4 | 0 | 0 | 0 | 0 | 7 | 0 | 0 | 0 | 1 | 0 | 6 | 0 | 4 | 0 | 0 | 0 | 0 |
| chrM | 381 | 13855 | chrM:381-13855 | 0 | 0 | 0 | 0 | 152 | 0 | 1167 | 339 | 0 | 0 | 0 | 0 | 0 | 0 | 105 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2744 | 0 | 0 | 0 | 29 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chrM | 381 | 8811 | chrM:381-8811 | 0 | 0 | 41 | 0 | 182 | 77 | 60 | 20 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| KI270744.1 | 701 | 6071 | KI270744.1:701-6071 | 0 | 0 | 9 | 35 | 3 | 7 | 0 | 5 | 0 | 21 | 3 | 1 | 0 | 106 | 0 | 1 | 5 | 3 | 2 | 3 | 5 | 2 | 13 | 0 | 1 | 0 | 9 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 3 | 0 | 0 | 0 | 22 | 41 | 4 | 2 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 5 | 1 | 5 | 0 | 0 | 16 | 3 | 1 | 0 | 1 | 0 | 13 | 1 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Count differences between counts-filtered junction files

Count number of junctions in each filtered file

``` r
nrow(filtered_combined_junctions)
```

    [1] 322558

``` r
nrow(filtered_merged_junctions)
```

    [1] 284998

Calculate fraction of filtered junctions that are only in the combined
or merged files

``` r
# print number of junctions only in combined junctions bedfile after filtering for junctions with at least 10 counts in at least 5 samples
filtered_combined_only <- dplyr::setdiff(filtered_combined_junctions, filtered_merged_junctions)
nrow(filtered_combined_only)
```

    [1] 46792

``` r
# check fraction of filtered_combined_only junctions
nrow(filtered_combined_only) / nrow(deduplicated_combined_junctions)
```

    [1] 0.007100426

``` r
# print number of junctions only in merged junctions bedfile after filtering for junctions with at least 10 counts in at least 5 samples
filtered_merged_only <- dplyr::setdiff(filtered_merged_junctions, filtered_combined_junctions)
nrow(filtered_merged_only)
```

    [1] 9232

``` r
# check fraction of merged_chr1 junctions > 10 counts
nrow(filtered_merged_only) / nrow(merged_junctions)
```

    [1] 0.001398681

### Print examples of junctions only in merged or combined tables

Print examples of junctions only in the combined junctions file

``` r
filtered_combined_only |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| KI270721.1 | 11643 | 11644 | KI270721.1:11643-11644 | 0 | 0 | 0 | 0 | 29 | 0 | 0 | 0 | 3 | 0 | 1 | 0 | 7 | 3 | 16 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 5 | 133 | 51 | 3 | 0 | 2 | 3 | 1 | 124 | 1 | 0 | 2 | 0 | 7 | 0 | 0 | 0 | 193 | 0 | 81 | 0 | 0 | 0 | 167 | 1 | 0 | 219 | 0 | 27 | 0 | 0 | 3 | 74 | 36 | 48 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 6 | 0 | 0 | 0 | 0 | 0 |
| KI270721.1 | 11736 | 11737 | KI270721.1:11736-11737 | 1 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 20 | 2 | 0 | 0 | 0 | 0 | 0 | 9 | 0 | 0 | 0 | 0 | 5 | 0 | 0 | 0 | 40 | 0 | 24 | 0 | 0 | 0 | 27 | 0 | 0 | 20 | 0 | 3 | 0 | 2 | 0 | 13 | 35 | 37 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 8 | 0 | 0 | 0 | 0 | 0 |
| chr9 | 16061 | 16062 | chr9:16061-16062 | 73 | 8 | 31 | 43 | 8 | 29 | 17 | 106 | 12 | 62 | 30 | 23 | 15 | 19 | 12 | 5 | 0 | 1 | 2 | 3 | 2 | 0 | 2 | 2 | 1 | 3 | 0 | 5 | 2 | 2 | 25 | 1 | 0 | 50 | 24 | 7 | 9 | 3 | 23 | 9 | 4 | 7 | 2 | 3 | 11 | 2 | 13 | 3 | 0 | 37 | 6 | 24 | 4 | 0 | 2 | 0 | 0 | 3 | 5 | 0 | 4 | 0 | 7 | 17 | 0 | 2 | 176 | 2 | 1 | 5 | 2 | 2 | 25 | 24 | 7 | 3 | 29 | 1 | 4 | 5 | 0 | 13 | 4 | 6 | 1 | 2 | 1 | 5 |
| chr16 | 17750 | 17751 | chr16:17750-17751 | 3 | 0 | 5 | 26 | 35 | 15 | 28 | 17 | 15 | 20 | 16 | 15 | 6 | 2 | 3 | 0 | 1 | 0 | 1 | 2 | 2 | 1 | 2 | 2 | 2 | 0 | 0 | 0 | 4 | 4 | 2 | 0 | 0 | 0 | 10 | 1 | 2 | 4 | 0 | 0 | 4 | 2 | 1 | 1 | 2 | 0 | 0 | 2 | 0 | 12 | 0 | 8 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 1 | 1 | 0 | 0 | 4 | 0 | 0 | 0 | 2 | 1 | 2 | 6 | 6 | 5 | 8 | 1 | 1 | 0 | 1 | 1 | 4 | 1 | 2 | 3 | 14 | 2 |
| chr12 | 17868 | 17869 | chr12:17868-17869 | 19 | 2 | 11 | 9 | 24 | 19 | 24 | 18 | 27 | 20 | 14 | 12 | 3 | 1 | 10 | 0 | 6 | 8 | 6 | 1 | 0 | 5 | 1 | 4 | 8 | 2 | 0 | 0 | 17 | 8 | 10 | 8 | 1 | 3 | 45 | 5 | 14 | 22 | 12 | 2 | 12 | 12 | 9 | 5 | 5 | 3 | 4 | 1 | 1 | 5 | 4 | 9 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 2 | 2 | 10 | 11 | 32 | 2 | 1 | 5 | 0 | 1 | 4 | 13 | 15 | 54 | 71 | 4 | 0 | 3 | 13 | 3 | 1 | 8 | 7 | 6 | 31 | 3 |
| chr1 | 18369 | 18370 | chr1:18369-18370 | 7 | 10 | 9 | 4 | 7 | 5 | 7 | 7 | 12 | 10 | 2 | 10 | 2 | 1 | 2 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 2 | 2 | 7 | 1 | 0 | 1 | 2 | 1 | 2 | 1 | 6 | 11 | 8 | 1 | 7 | 1 | 0 | 4 | 5 | 0 | 2 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 1 | 0 | 3 | 0 | 1 | 4 | 2 | 0 | 2 | 0 | 1 | 4 | 0 | 0 | 3 | 1 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 3 | 4 | 1 |

Print examples of junctions only in the merged junctions file

``` r
filtered_merged_only |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1559044 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1559052 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1559054 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1559075 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1559100 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR1559105 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR1559133 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 | SRR1559134 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GL000218.1 | 42149 | 42150 | GL000218.1:42149-42150 | 2 | 16 | 21 | 88 | 10 | 12 | 5 | 0 | 0 | 0 | 18 | 5 | 0 | 5 | 0 | 10 | 0 | 0 | 0 | 0 | 7 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 29 | 0 | 0 | 0 | 3 | 1 | 0 | 0 | 0 | 1 | 0 | 25 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 36 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 76 | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 33 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 78 |
| GL000218.1 | 42323 | 42324 | GL000218.1:42323-42324 | 12 | 12 | 27 | 77 | 18 | 12 | 2 | 0 | 0 | 0 | 20 | 5 | 0 | 4 | 0 | 20 | 0 | 0 | 0 | 0 | 14 | 0 | 8 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 18 | 0 | 0 | 0 | 8 | 1 | 0 | 0 | 0 | 3 | 0 | 22 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 10 | 32 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 75 | 0 | 0 | 4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 62 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 92 |
| KI270721.1 | 52446 | 52447 | KI270721.1:52446-52447 | 0 | 3 | 2 | 29 | 35 | 81 | 1 | 4 | 3 | 0 | 2 | 6 | 0 | 4 | 1 | 0 | 2 | 0 | 1 | 2 | 1 | 99 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 64 | 4 | 185 | 5 | 0 | 1 | 0 | 2 | 2 | 1 | 1 | 10 | 9 | 0 | 0 | 0 | 0 | 0 | 0 | 7 | 5 | 3 | 2 | 11 | 39 | 19 | 6 | 26 | 3 | 0 | 0 | 82 | 0 | 6 | 67 | 2127 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 16 |
| chr16 | 60539 | 60540 | chr16:60539-60540 | 0 | 0 | 0 | 46 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 26 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 58 | 0 | 0 | 0 | 0 | 0 | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 186 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 53 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr16 | 61119 | 61120 | chr16:61119-61120 | 0 | 0 | 0 | 15 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 14 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 46 | 0 | 0 | 0 | 0 | 0 | 14 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 87 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 21 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 5 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| GL000219.1 | 77669 | 77670 | GL000219.1:77669-77670 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 80 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 538 | 0 | 0 | 0 | 0 | 267 | 0 | 0 | 730 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 65 | 60 | 0 | 0 | 0 | 0 | 0 | 76 | 0 | 0 | 0 | 35 | 0 | 52 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
