# Check Merge vs. Combined TARGET Pilot GTF and Junction files
Cindy Liang (celiang@ucsc.edu)
2026-06-10

## Read in directories and files

Define directories and file paths

``` r
### Directories ###
## combined shiba run results on target pilot samples ##
# merged table eval dir
exploration_eval_dir <- file.path("exploration", "merged-method-target-pilot-eval")

# shiba results dir
# target pilot explroation dir
target_pilot_dir <- file.path("exploration", "merging-psi-tables", "target_pilot")
combined_results_dir <- file.path(target_pilot_dir, "shiba_combined")
# junctions dir
combined_junctions_dir <- file.path(combined_results_dir, "junctions")
# gtf dir
combined_gtf_dir <- file.path(combined_results_dir, "annotation")

## merged shiba run results on target pilot samples ##
# shiba results dir
merged_results_dir <- file.path("results", "merged_shiba", "target_pilot")

### Files ###
## combined shiba target pilot files
target_pilot_sample_sheet <- file.path(target_pilot_dir, "experiment.tsv")
# combined target pilot outputs
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

No problems are present in merged_junctions file or deduplicated
junctions file

## Check if files are identical

Deduplicated and merged junctions are still not identical, so need to
drill down on what junctions are different

``` r
# ensure column order is the same before checking if data frames are identical
deduplicated_combined_junctions <- deduplicated_combined_junctions |>
  dplyr::relocate(dplyr::all_of(target_pilot_samples), .after = dplyr::last_col())

merged_junctions <- merged_junctions |>
  dplyr::relocate(dplyr::all_of(target_pilot_samples), .after = dplyr::last_col())

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

### Subsample junctions to those in chromosome 1

``` r
merged_junctions_chr1 <- merged_junctions |> dplyr::filter(chr == "chr1")
deduplicated_combined_junctions_chr1 <- deduplicated_combined_junctions |> dplyr::filter(chr == "chr1")
```

Check head of merged chromosome 1 junctions

``` r
merged_junctions_chr1 |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 10059 | 10252 | chr1:10059-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10065 | 10252 | chr1:10065-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10071 | 10252 | chr1:10071-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10077 | 10252 | chr1:10077-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10083 | 10252 | chr1:10083-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10089 | 10252 | chr1:10089-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Check head of deduplicated combined chromosome 1 junctions

``` r
deduplicated_combined_junctions_chr1 |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 10059 | 10252 | chr1:10059-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10065 | 10252 | chr1:10065-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10071 | 10252 | chr1:10071-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10077 | 10252 | chr1:10077-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10083 | 10252 | chr1:10083-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 10089 | 10252 | chr1:10089-10252 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

From quickly eyeballing chromosome 1, the column orders of samples are
not the same, but according to [setdiff
documentation](https://dplyr.tidyverse.org/reference/setops.html), the
order of the columns should not matter.

Check of chromosomes 1 of the combined and merged junctions bedfiles are
identical

``` r
identical(deduplicated_combined_junctions_chr1, merged_junctions_chr1)
```

    [1] FALSE

Many junctions also have counts below 10 across all samples. These
junctions are not expected to make it into PSI calculation based off
default Shiba filters (\> 10 junctions), but may factor into splice
event coordinate definition.

Remove junctions in chromosome 1 whose counts are \< 10 in all samples
and check for differences

``` r
merged_junctions_chr1_10 <- merged_junctions_chr1 |> dplyr::filter(
  dplyr::if_all(
    dplyr::all_of(target_pilot_samples), ~ . > 10
    )
  )

merged_junctions_chr1_10 |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 14829 | 14970 | chr1:14829-14970 | 1093 | 326 | 502 | 459 | 2985 | 1570 | 2181 | 2385 | 3673 | 1648 | 1947 | 1662 | 450 | 947 | 1135 | 268 | 216 | 124 | 123 | 287 | 284 | 228 | 273 | 448 | 212 | 242 | 156 | 158 | 440 | 443 | 302 | 181 | 207 | 560 | 689 | 219 | 520 | 513 | 1097 | 570 | 731 | 591 | 1630 | 1100 | 610 | 798 | 781 | 730 | 271 | 1119 | 187 | 526 | 128 | 58 | 54 | 98 | 420 | 972 | 345 | 636 | 448 | 455 | 513 | 815 | 285 | 412 | 668 | 205 | 113 | 313 | 776 | 774 | 769 | 972 | 554 | 788 | 705 | 84 | 133 | 181 | 121 | 243 | 203 | 369 | 31 | 90 | 238 | 99 |
| chr1 | 15038 | 15796 | chr1:15038-15796 | 266 | 114 | 111 | 148 | 960 | 477 | 816 | 963 | 1041 | 629 | 774 | 499 | 147 | 276 | 366 | 68 | 48 | 43 | 26 | 96 | 93 | 57 | 75 | 119 | 67 | 63 | 51 | 39 | 176 | 179 | 130 | 49 | 71 | 229 | 216 | 63 | 169 | 165 | 245 | 165 | 220 | 147 | 393 | 285 | 172 | 159 | 214 | 216 | 109 | 279 | 49 | 250 | 60 | 34 | 14 | 40 | 112 | 285 | 90 | 209 | 156 | 98 | 132 | 178 | 91 | 182 | 215 | 53 | 30 | 93 | 176 | 224 | 318 | 317 | 166 | 260 | 207 | 19 | 41 | 39 | 32 | 100 | 69 | 157 | 17 | 31 | 81 | 30 |
| chr1 | 17055 | 17233 | chr1:17055-17233 | 1426 | 272 | 763 | 581 | 2649 | 1442 | 2438 | 2161 | 2368 | 1417 | 1735 | 1953 | 452 | 660 | 841 | 267 | 268 | 202 | 198 | 399 | 288 | 231 | 259 | 476 | 217 | 224 | 185 | 225 | 554 | 425 | 312 | 246 | 260 | 507 | 462 | 218 | 517 | 621 | 939 | 495 | 607 | 626 | 947 | 782 | 512 | 556 | 624 | 571 | 188 | 974 | 169 | 609 | 209 | 101 | 21 | 21 | 353 | 630 | 352 | 686 | 349 | 362 | 387 | 721 | 216 | 417 | 704 | 202 | 42 | 252 | 874 | 749 | 357 | 728 | 420 | 873 | 472 | 48 | 76 | 140 | 130 | 289 | 253 | 347 | 57 | 197 | 274 | 187 |
| chr1 | 17742 | 17915 | chr1:17742-17915 | 847 | 242 | 572 | 455 | 1911 | 1137 | 1559 | 2135 | 1973 | 1328 | 1577 | 1671 | 462 | 541 | 651 | 226 | 165 | 138 | 171 | 197 | 302 | 238 | 216 | 491 | 262 | 203 | 124 | 193 | 359 | 363 | 266 | 137 | 205 | 608 | 494 | 140 | 425 | 385 | 768 | 382 | 423 | 357 | 533 | 353 | 390 | 398 | 625 | 352 | 124 | 728 | 110 | 448 | 182 | 98 | 11 | 17 | 296 | 400 | 250 | 540 | 236 | 321 | 370 | 565 | 186 | 321 | 1282 | 167 | 12 | 184 | 593 | 603 | 524 | 849 | 476 | 537 | 401 | 35 | 64 | 132 | 205 | 198 | 225 | 330 | 45 | 125 | 222 | 155 |
| chr1 | 18061 | 18268 | chr1:18061-18268 | 661 | 350 | 327 | 276 | 1373 | 628 | 918 | 1571 | 976 | 640 | 685 | 863 | 240 | 161 | 275 | 43 | 40 | 40 | 54 | 83 | 87 | 70 | 73 | 155 | 94 | 52 | 57 | 66 | 221 | 261 | 191 | 132 | 122 | 334 | 182 | 133 | 288 | 264 | 614 | 272 | 254 | 245 | 322 | 225 | 314 | 252 | 366 | 321 | 54 | 454 | 68 | 240 | 82 | 40 | 14 | 25 | 186 | 336 | 150 | 293 | 149 | 197 | 168 | 333 | 171 | 213 | 168 | 91 | 11 | 101 | 375 | 323 | 222 | 447 | 154 | 259 | 140 | 19 | 47 | 33 | 44 | 48 | 73 | 133 | 27 | 102 | 190 | 150 |
| chr1 | 185350 | 185491 | chr1:185350-185491 | 1142 | 264 | 422 | 445 | 3090 | 1663 | 2323 | 2661 | 4140 | 1848 | 2170 | 1821 | 519 | 976 | 1253 | 307 | 233 | 132 | 128 | 304 | 333 | 260 | 305 | 514 | 221 | 279 | 180 | 182 | 425 | 377 | 266 | 171 | 221 | 538 | 636 | 218 | 498 | 531 | 976 | 480 | 603 | 548 | 1265 | 865 | 499 | 753 | 627 | 722 | 291 | 1044 | 173 | 507 | 138 | 71 | 46 | 95 | 373 | 774 | 301 | 529 | 374 | 379 | 350 | 680 | 275 | 337 | 753 | 173 | 100 | 242 | 696 | 680 | 659 | 872 | 491 | 670 | 629 | 75 | 100 | 210 | 153 | 298 | 238 | 448 | 30 | 90 | 209 | 83 |

``` r
deduplicated_combined_junctions_chr1_10 <- deduplicated_combined_junctions_chr1 |> dplyr::filter(
  dplyr::if_all(
    dplyr::all_of(target_pilot_samples), ~ . > 10
    )
  )

deduplicated_combined_junctions_chr1_10 |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 14829 | 14970 | chr1:14829-14970 | 1093 | 326 | 502 | 459 | 2985 | 1570 | 2181 | 2385 | 3673 | 1648 | 1947 | 1662 | 450 | 947 | 1135 | 268 | 216 | 124 | 123 | 287 | 284 | 228 | 273 | 448 | 212 | 242 | 156 | 158 | 440 | 443 | 302 | 181 | 207 | 560 | 689 | 219 | 520 | 513 | 1097 | 570 | 731 | 591 | 1630 | 1100 | 610 | 798 | 781 | 730 | 271 | 1119 | 187 | 526 | 128 | 58 | 54 | 98 | 420 | 972 | 345 | 636 | 448 | 455 | 513 | 815 | 285 | 412 | 668 | 205 | 113 | 313 | 776 | 774 | 769 | 972 | 554 | 788 | 705 | 84 | 133 | 181 | 121 | 243 | 203 | 369 | 31 | 90 | 238 | 99 |
| chr1 | 15038 | 15796 | chr1:15038-15796 | 266 | 114 | 111 | 148 | 960 | 477 | 816 | 963 | 1041 | 629 | 774 | 499 | 147 | 276 | 366 | 68 | 48 | 43 | 26 | 96 | 93 | 57 | 75 | 119 | 67 | 63 | 51 | 39 | 176 | 179 | 130 | 49 | 71 | 229 | 216 | 63 | 169 | 165 | 245 | 165 | 220 | 147 | 393 | 285 | 172 | 159 | 214 | 216 | 109 | 279 | 49 | 250 | 60 | 34 | 14 | 40 | 112 | 285 | 90 | 209 | 156 | 98 | 132 | 178 | 91 | 182 | 215 | 53 | 30 | 93 | 176 | 224 | 318 | 317 | 166 | 260 | 207 | 19 | 41 | 39 | 32 | 100 | 69 | 157 | 17 | 31 | 81 | 30 |
| chr1 | 17055 | 17233 | chr1:17055-17233 | 1426 | 272 | 763 | 581 | 2649 | 1442 | 2438 | 2161 | 2368 | 1417 | 1735 | 1953 | 452 | 660 | 841 | 267 | 268 | 202 | 198 | 399 | 288 | 231 | 259 | 476 | 217 | 224 | 185 | 225 | 554 | 425 | 312 | 246 | 260 | 507 | 462 | 218 | 517 | 621 | 939 | 495 | 607 | 626 | 947 | 782 | 512 | 556 | 624 | 571 | 188 | 974 | 169 | 609 | 209 | 101 | 21 | 21 | 353 | 630 | 352 | 686 | 349 | 362 | 387 | 721 | 216 | 417 | 704 | 202 | 42 | 252 | 874 | 749 | 357 | 728 | 420 | 873 | 472 | 48 | 76 | 140 | 130 | 289 | 253 | 347 | 57 | 197 | 274 | 187 |
| chr1 | 17742 | 17915 | chr1:17742-17915 | 847 | 242 | 572 | 455 | 1911 | 1137 | 1559 | 2135 | 1973 | 1328 | 1577 | 1671 | 462 | 541 | 651 | 226 | 165 | 138 | 171 | 197 | 302 | 238 | 216 | 491 | 262 | 203 | 124 | 193 | 359 | 363 | 266 | 137 | 205 | 608 | 494 | 140 | 425 | 385 | 768 | 382 | 423 | 357 | 533 | 353 | 390 | 398 | 625 | 352 | 124 | 728 | 110 | 448 | 182 | 98 | 11 | 17 | 296 | 400 | 250 | 540 | 236 | 321 | 370 | 565 | 186 | 321 | 1282 | 167 | 12 | 184 | 593 | 603 | 524 | 849 | 476 | 537 | 401 | 35 | 64 | 132 | 205 | 198 | 225 | 330 | 45 | 125 | 222 | 155 |
| chr1 | 18061 | 18268 | chr1:18061-18268 | 661 | 350 | 327 | 276 | 1373 | 628 | 918 | 1571 | 976 | 640 | 685 | 863 | 240 | 161 | 275 | 43 | 40 | 40 | 54 | 83 | 87 | 70 | 73 | 155 | 94 | 52 | 57 | 66 | 221 | 261 | 191 | 132 | 122 | 334 | 182 | 133 | 288 | 264 | 614 | 272 | 254 | 245 | 322 | 225 | 314 | 252 | 366 | 321 | 54 | 454 | 68 | 240 | 82 | 40 | 14 | 25 | 186 | 336 | 150 | 293 | 149 | 197 | 168 | 333 | 171 | 213 | 168 | 91 | 11 | 101 | 375 | 323 | 222 | 447 | 154 | 259 | 140 | 19 | 47 | 33 | 44 | 48 | 73 | 133 | 27 | 102 | 190 | 150 |
| chr1 | 185350 | 185491 | chr1:185350-185491 | 1142 | 264 | 422 | 445 | 3090 | 1663 | 2323 | 2661 | 4140 | 1848 | 2170 | 1821 | 519 | 976 | 1253 | 307 | 233 | 132 | 128 | 304 | 333 | 260 | 305 | 514 | 221 | 279 | 180 | 182 | 425 | 377 | 266 | 171 | 221 | 538 | 636 | 218 | 498 | 531 | 976 | 480 | 603 | 548 | 1265 | 865 | 499 | 753 | 627 | 722 | 291 | 1044 | 173 | 507 | 138 | 71 | 46 | 95 | 373 | 774 | 301 | 529 | 374 | 379 | 350 | 680 | 275 | 337 | 753 | 173 | 100 | 242 | 696 | 680 | 659 | 872 | 491 | 670 | 629 | 75 | 100 | 210 | 153 | 298 | 238 | 448 | 30 | 90 | 209 | 83 |

The two tables look identical so far - check with setdiff

``` r
identical(deduplicated_combined_junctions_chr1, merged_junctions_chr1)
```

    [1] FALSE

``` r
identical(deduplicated_combined_junctions_chr1_10, merged_junctions_chr1_10)
```

    [1] FALSE

Check what’s different in chromosome 1 junctions \> 10 counts

``` r
combined_only_chr1 <- dplyr::setdiff(deduplicated_combined_junctions_chr1, merged_junctions_chr1)
nrow(combined_only_chr1)
```

    [1] 5173

``` r
# check fraction of combined_only_chr1 junctions that are only in chr1 of deduplicated combined junctions
nrow(combined_only_chr1) / nrow(deduplicated_combined_junctions_chr1)
```

    [1] 0.009350406

``` r
merged_only_chr1 <- dplyr::setdiff(merged_junctions_chr1, deduplicated_combined_junctions_chr1)
nrow(merged_only_chr1)
```

    [1] 6599

``` r
# check fraction of merged_only_chr1 junctions that are only in chr1 of merged combined junctions
nrow(merged_only_chr1) / nrow(merged_junctions_chr1)
```

    [1] 0.01189729

``` r
combined_only_chr1_10 <- dplyr::setdiff(deduplicated_combined_junctions_chr1_10, merged_junctions_chr1_10)
nrow(combined_only_chr1_10)
```

    [1] 61

``` r
# check fraction of combined_only_chr1 junctions > 10 counts that are only in chr1 of deduplicated combined junctions
nrow(combined_only_chr1_10) / nrow(deduplicated_combined_junctions_chr1)
```

    [1] 0.00011026

``` r
merged_only_chr1_10 <- dplyr::setdiff(merged_junctions_chr1_10, deduplicated_combined_junctions_chr1_10)
nrow(merged_only_chr1_10)
```

    [1] 0

``` r
# check fraction of combined_only_chr1 junctions > 10 counts that are only in chr1 of deduplicated combined junctions
nrow(merged_only_chr1_10) / nrow(merged_junctions_chr1)
```

    [1] 0

Overall, only 1% of all chromosome 1 junctions are only found in either
the combined or merged junctions files. If we filter for junctions with
more than 10 counts, we miss a small (0.1%) fraction of junctions in the
combined file using the merged method, but we do not have any junctions
not found in the combined method, which is great.

Potential sources of differences in the merged method that may result in
junctions or junction values not present in the combined junctions file:

- After merging separate junction.bed files, I manually convert all NAs
  to 0. I don’t know how Shiba handles junctions not found in all
  samples in the combined version

- I don’t know how shiba “decides” what junctions get the duplication
  bug and how the bug impacts how junction coordinates are defined in
  each run

Another note - because the concatenated junctions bug impacted one
junction with reads \> 10 in the combined junctions bedfile, I expect
this will result in differenes in the PSI table
