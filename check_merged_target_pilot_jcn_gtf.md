# Check Merge vs. Combined TARGET Pilot GTF and Junction files
Cindy Liang (celiang@ucsc.edu)
2026-06-08

## Read in directories and files

Define directories and file paths

``` r
### Directories ###
## combined shiba run results on target pilot samples ##
# merged table eval dir
exploration_eval_dir <- file.path("exploration", "merged-method-target-pilot-eval")

# shiba results dir
combined_results_dir <- file.path("exploration", "merging-psi-tables", "target_pilot", "shiba_combined")
# junctions dir
combined_junctions_dir <- file.path(combined_results_dir, "junctions")
# gtf dir
combined_gtf_dir <- file.path(combined_results_dir, "annotation")

## merged shiba run results on target pilot samples ##
# shiba results dir
merged_results_dir <- file.path("results", "merged_shiba", "target_pilot")

### Files ###
## combined shiba target pilot files
combined_gtf_file <- file.path(combined_gtf_dir, "assembled_annotation.gtf")
combined_junctions_file <- file.path(combined_junctions_dir, "junctions.bed")
# deduplicated combined junctions file - made by running bash scripts/deduplicated_shiba_junctions.sh
deduplicated_combined_junctions_file <- file.path(combined_junctions_dir, "deduplicated_junctions.bed")

## merged shiba target pilot files
merged_gtf_file <- file.path(merged_results_dir, "merged_gtf.gtf")
merged_junctions_file <- file.path(merged_results_dir, "merged_junctions.bed")

## output files
only_in_comb_junc <- file.path(exploration_eval_dir, "jcn_only_in_combined.bed")
only_in_merged_junc <- file.path(exploration_eval_dir, "jcn_only_in_merged.bed")
```

Create output dir if it does not exist

``` r
# Check if it exists; if not, create it
if (!dir.exists(exploration_eval_dir)) {
  dir.create(exploration_eval_dir, recursive = TRUE)
}
```

Read in files

``` r
# read in files generated from combined shiba run
combined_junctions <- readr::read_tsv(
  combined_junctions_file,
  col_types = readr::cols(
    .default = "d",
    ID = "c",
    chr = "c",
    start = "d",
    end = "d"
  )
)
```

    Warning: One or more parsing issues, call `problems()` on your data frame for details,
    e.g.:
      dat <- vroom(...)
      problems(dat)

``` r
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

Check tables for parsing errors

``` r
combined_junctions_dat <- vroom::vroom(combined_junctions_file)
```

    Rows: 6590027 Columns: 92
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (2): chr, ID
    dbl (90): start, end, SRR1559043, SRR1559044, SRR1559052, SRR1559054, SRR155...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
combined_junctions_problems <- vroom::problems(combined_junctions_dat)
```

    Warning: One or more parsing issues, call `problems()` on your data frame for details,
    e.g.:
      dat <- vroom(...)
      problems(dat)

``` r
merged_junctions_dat <- vroom::vroom(merged_junctions_file)
```

    Rows: 6600506 Columns: 92
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (2): chr, ID
    dbl (90): start, end, SRR1559043, SRR1559145, SRR1559160, SRR1559164, SRR155...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
deduplicated_comb_junctions_dat <- vroom::vroom(deduplicated_combined_junctions_file)
```

    Rows: 6590027 Columns: 92
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (2): chr, ID
    dbl (90): start, end, SRR1559043, SRR1559044, SRR1559052, SRR1559054, SRR155...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
deduplicated_comb_junctions_problems <- vroom::problems(deduplicated_combined_junctions)
```

No problems are present in merged_junctions file

## Check if files are identical

### Check if combined and merged junctions.bed files are identical

``` r
identical(combined_junctions, merged_junctions)
```

    [1] FALSE

### Check if deduplicated combined and merged junctions are identical

``` r
identical(deduplicated_combined_junctions, merged_junctions)
```

    [1] FALSE

## Check differences in junctions bedfiles

Dimensions of each bedfile

``` r
dim(combined_junctions)
```

    [1] 6590027      92

``` r
dim(deduplicated_combined_junctions)
```

    [1] 6590027      92

``` r
dim(merged_junctions)
```

    [1] 6600506      92

Take random 1000 junctions for ease of debugging

``` r
# create list of dataframes consisting of merged and combined junctions file
junctions_list <- list(merged_junctions, deduplicated_combined_junctions)

# subset junctions dfs by randomly taking 100 rows
junctions_list <- lapply(junctions_list, dplyr::filter, chr == "chr1")
```

Elements unique to each junction bedfile

``` r
# elements in deduplicated combined junctions df but not in merged
junctions_unique_to_combined <- setdiff(deduplicated_combined_junctions, merged_junctions)

# print number of junctions that are unique to combined bedfile
nrow(junctions_unique_to_combined)
```

    NULL

``` r
# elements in combined junctions df but not in merged
junctions_unique_to_merged <- setdiff(merged_junctions, deduplicated_combined_junctions)

# print number of junctions that are unique to merged bedfile
nrow(junctions_unique_to_merged)
```

    NULL

# check what class the setdiff output IDs

``` r
class(junctions_unique_to_merged)
```

    [1] "list"
