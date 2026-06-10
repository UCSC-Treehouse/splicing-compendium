# Check Merge vs. Combined TARGET Pilot Junction files (chromosome 1)
Cindy Liang (celiang@ucsc.edu)
2026-06-10

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
# gtf dir
combined_gtf_dir <- file.path(combined_results_dir, "annotation")

## merged shiba run results on target pilot samples ##
# shiba results dir
merged_results_dir <- file.path(repo_root, "results", "merged_shiba", "target_pilot")

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
# minimum number of samples to filter minimum number of junctions by
min_samples <- 5

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

## Check tables for parsing errors

### Print parsing errors in combined junctions table

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
combined_junctions_problems
```

| row | col | expected | actual | file |
|---:|---:|:---|:---|:---|
| 1277076 | 2 | a double | 6609298;6609298 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 1277076 | 3 | a double | 6609299;6609299 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 1719156 | 2 | a double | 37009809;37009809 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 1719156 | 3 | a double | 37009810;37009810 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 1907988 | 2 | a double | 24440926;24440926 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 1907988 | 3 | a double | 24440927;24440927 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 1907989 | 2 | a double | 99500871;99500871 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 1907989 | 3 | a double | 99500872;99500872 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 2114125 | 2 | a double | 43409741;43409741 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 2114125 | 3 | a double | 43409742;43409742 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 2442719 | 2 | a double | 30245611;30245611 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 2442719 | 3 | a double | 30245612;30245612 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 2442720 | 2 | a double | 89739553;89739553 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 2442720 | 3 | a double | 89739554;89739554 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 3336056 | 2 | a double | 155212126;155212126 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 3336056 | 3 | a double | 155212127;155212127 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344706 | 2 | a double | 21959232;21959232 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344706 | 3 | a double | 21959233;21959233 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344707 | 2 | a double | 241094432;241094432 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344707 | 3 | a double | 241094433;241094433 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344708 | 2 | a double | 47809157;47809157 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344708 | 3 | a double | 47809158;47809158 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344709 | 2 | a double | 61483936;61483936 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344709 | 3 | a double | 61483937;61483937 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344710 | 2 | a double | 74532819;74532819 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 4344710 | 3 | a double | 74532820;74532820 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 5837918 | 2 | a double | 72948980;72948980 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 5837918 | 3 | a double | 72948981;72948981 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 5837919 | 2 | a double | 72949102;72949102 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 5837919 | 3 | a double | 72949103;72949103 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 5837920 | 2 | a double | 99486475;99486475 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 5837920 | 3 | a double | 99486476;99486476 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 6033864 | 2 | a double | 143581328;143581328 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 6033864 | 3 | a double | 143581329;143581329 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 6271661 | 2 | a double | 35660647;35660647 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |
| 6271661 | 3 | a double | 35660648;35660648 | /home/ubuntu/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/junctions.bed |

We have junctions in the combined junctions file that are impacted by
the concatenation bug in Shiba. PSI files for events involving these
junctions are expected to be problematic in the combined PSI file.

### Check for parsing errors in the merged and deduplicated junctions tables

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

``` r
# ensure column order is the same before checking if data frames are identical
deduplicated_combined_junctions <- deduplicated_combined_junctions |>
  dplyr::relocate(dplyr::all_of(target_pilot_samples), .after = dplyr::last_col())

merged_junctions <- merged_junctions |>
  dplyr::relocate(dplyr::all_of(target_pilot_samples), .after = dplyr::last_col())

identical(deduplicated_combined_junctions, merged_junctions)
```

    [1] FALSE

Deduplicated and merged junctions are still not identical, so need to
drill down on what junctions are different

## Check differences in junctions bedfiles

Dimensions of each bedfile

``` r
dim(deduplicated_combined_junctions)
```

    [1] 6590027      92

``` r
dim(merged_junctions)
```

    [1] 6600506      92

There are 10479 more junctions in the merged junction file compared to
the combined junction file.

## Subsample junctions to those in chromosome 1

To allow for running on my laptop, subsample junctions files for those
in chromosome 1 prior to checking for differences

``` r
deduplicated_combined_junctions_chr1 <- deduplicated_combined_junctions |> dplyr::filter(chr == "chr1")
merged_junctions_chr1 <- merged_junctions |> dplyr::filter(chr == "chr1")
```

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

Many junctions also have counts below 10 across all samples. These
junctions are not expected to make it into PSI calculation based off
default Shiba filters (\> 10 junctions).

Check number of junctions in each subsetted dataframe

``` r
nrow(deduplicated_combined_junctions_chr1)
```

    [1] 553238

``` r
nrow(merged_junctions_chr1)
```

    [1] 554664

There are 1426 more junctions in the merged bedfile, compared to the
combined bedfile. The chromosome 1-subsetted junction bedfiles are not
identical.

### Filter for junctions with at least 10 junctions in min_samples

Filter for junctions in chromosome 1 whose counts are \< 10 in 5 samples
and check for differences

``` r
merged_junctions_chr1_10 <- merged_junctions_chr1 |> dplyr::filter(
  # count number of columns per row with junction count of at least 10
  rowSums(
    dplyr::across(
      dplyr::all_of(target_pilot_samples), ~ . >= 10
    )
    # filter for rows where at least 5 columns (samples) have counts of at least 10
  ) >= 5
  )

merged_junctions_chr1_10 |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 14614 | 16858 | chr1:14614-16858 | 0 | 0 | 0 | 0 | 22 | 10 | 9 | 18 | 6 | 9 | 14 | 20 | 0 | 0 | 3 | 0 | 3 | 0 | 0 | 0 | 3 | 0 | 0 | 5 | 0 | 2 | 0 | 2 | 18 | 11 | 0 | 2 | 7 | 11 | 6 | 5 | 7 | 15 | 10 | 17 | 6 | 8 | 22 | 12 | 14 | 9 | 10 | 9 | 1 | 18 | 4 | 17 | 7 | 2 | 0 | 0 | 10 | 10 | 0 | 14 | 0 | 6 | 11 | 2 | 4 | 10 | 28 | 3 | 1 | 4 | 20 | 5 | 7 | 8 | 25 | 16 | 6 | 0 | 0 | 0 | 3 | 2 | 1 | 5 | 4 | 16 | 17 | 15 |
| chr1 | 14737 | 14738 | chr1:14737-14738 | 653 | 302 | 409 | 362 | 2100 | 874 | 853 | 835 | 1482 | 361 | 899 | 914 | 200 | 542 | 359 | 72 | 117 | 66 | 28 | 92 | 114 | 113 | 121 | 176 | 78 | 83 | 59 | 57 | 92 | 117 | 75 | 35 | 30 | 112 | 135 | 42 | 88 | 129 | 236 | 136 | 154 | 140 | 740 | 365 | 170 | 115 | 181 | 102 | 49 | 243 | 40 | 71 | 30 | 20 | 16 | 39 | 151 | 214 | 107 | 189 | 114 | 128 | 157 | 172 | 129 | 152 | 141 | 47 | 36 | 80 | 114 | 175 | 157 | 193 | 113 | 177 | 115 | 2 | 9 | 62 | 21 | 46 | 58 | 48 | 5 | 1 | 17 | 9 |
| chr1 | 14764 | 24846 | chr1:14764-24846 | 54 | 88 | 17 | 31 | 129 | 64 | 51 | 65 | 108 | 29 | 62 | 61 | 8 | 25 | 20 | 0 | 4 | 0 | 0 | 6 | 5 | 6 | 3 | 0 | 7 | 7 | 5 | 5 | 8 | 10 | 16 | 6 | 5 | 11 | 26 | 8 | 11 | 18 | 31 | 17 | 42 | 15 | 112 | 48 | 20 | 22 | 24 | 10 | 14 | 26 | 13 | 12 | 0 | 0 | 0 | 7 | 21 | 25 | 13 | 17 | 14 | 15 | 18 | 26 | 11 | 16 | 0 | 12 | 6 | 11 | 23 | 23 | 27 | 34 | 22 | 16 | 18 | 4 | 13 | 4 | 1 | 3 | 5 | 0 | 0 | 1 | 1 | 3 |
| chr1 | 14829 | 14970 | chr1:14829-14970 | 1093 | 326 | 502 | 459 | 2985 | 1570 | 2181 | 2385 | 3673 | 1648 | 1947 | 1662 | 450 | 947 | 1135 | 268 | 216 | 124 | 123 | 287 | 284 | 228 | 273 | 448 | 212 | 242 | 156 | 158 | 440 | 443 | 302 | 181 | 207 | 560 | 689 | 219 | 520 | 513 | 1097 | 570 | 731 | 591 | 1630 | 1100 | 610 | 798 | 781 | 730 | 271 | 1119 | 187 | 526 | 128 | 58 | 54 | 98 | 420 | 972 | 345 | 636 | 448 | 455 | 513 | 815 | 285 | 412 | 668 | 205 | 113 | 313 | 776 | 774 | 769 | 972 | 554 | 788 | 705 | 84 | 133 | 181 | 121 | 243 | 203 | 369 | 31 | 90 | 238 | 99 |
| chr1 | 14829 | 15021 | chr1:14829-15021 | 6 | 0 | 0 | 0 | 10 | 5 | 19 | 18 | 25 | 2 | 5 | 7 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 5 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 1 | 5 | 10 | 0 | 3 | 2 | 1 | 1 | 0 | 3 | 2 | 0 | 0 | 0 | 0 | 2 | 1 | 0 | 0 | 7 | 3 | 0 | 0 | 4 | 0 | 3 | 0 | 1 | 2 | 0 | 2 | 0 | 8 | 1 | 0 | 0 | 0 | 2 | 2 | 3 | 5 | 9 | 2 | 4 | 3 | 1 | 0 | 0 | 0 | 3 | 0 | 4 | 0 | 0 | 0 | 0 |
| chr1 | 14829 | 14830 | chr1:14829-14830 | 61 | 14 | 26 | 27 | 109 | 29 | 19 | 44 | 57 | 22 | 51 | 58 | 4 | 13 | 6 | 2 | 5 | 2 | 1 | 2 | 4 | 3 | 4 | 9 | 5 | 2 | 0 | 2 | 2 | 7 | 5 | 3 | 4 | 7 | 14 | 3 | 4 | 20 | 10 | 6 | 15 | 5 | 43 | 24 | 5 | 6 | 4 | 4 | 3 | 24 | 1 | 10 | 0 | 0 | 1 | 1 | 22 | 20 | 0 | 4 | 2 | 8 | 1 | 7 | 12 | 6 | 40 | 1 | 3 | 3 | 2 | 3 | 9 | 6 | 9 | 14 | 0 | 0 | 0 | 2 | 1 | 5 | 6 | 1 | 1 | 1 | 3 | 1 |

``` r
deduplicated_combined_junctions_chr1_10 <- deduplicated_combined_junctions_chr1 |> dplyr::filter(
  # count number of columns per row with junction count of at least 10
  rowSums(
    dplyr::across(
      dplyr::all_of(target_pilot_samples), ~ . >= 10
    )
    # filter for rows where at least 5 columns (samples) have counts of at least 10
  ) >= 5
  )

deduplicated_combined_junctions_chr1_10 |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 14614 | 16858 | chr1:14614-16858 | 0 | 0 | 0 | 0 | 22 | 10 | 9 | 18 | 6 | 9 | 14 | 20 | 0 | 0 | 3 | 0 | 3 | 0 | 0 | 0 | 3 | 0 | 0 | 5 | 0 | 2 | 0 | 2 | 18 | 11 | 0 | 2 | 7 | 11 | 6 | 5 | 7 | 15 | 10 | 17 | 6 | 8 | 22 | 12 | 14 | 9 | 10 | 9 | 1 | 18 | 4 | 17 | 7 | 2 | 0 | 0 | 10 | 10 | 0 | 14 | 0 | 6 | 11 | 2 | 4 | 10 | 28 | 3 | 1 | 4 | 20 | 5 | 7 | 8 | 25 | 16 | 6 | 0 | 0 | 0 | 3 | 2 | 1 | 5 | 4 | 16 | 17 | 15 |
| chr1 | 14737 | 14738 | chr1:14737-14738 | 653 | 302 | 409 | 362 | 2100 | 874 | 853 | 835 | 1482 | 361 | 899 | 914 | 200 | 542 | 359 | 72 | 117 | 66 | 28 | 92 | 114 | 113 | 121 | 176 | 78 | 83 | 59 | 57 | 92 | 117 | 75 | 35 | 30 | 112 | 135 | 42 | 88 | 129 | 236 | 136 | 154 | 140 | 740 | 365 | 170 | 115 | 181 | 102 | 49 | 243 | 40 | 71 | 30 | 20 | 16 | 39 | 151 | 214 | 107 | 189 | 114 | 128 | 157 | 172 | 129 | 152 | 141 | 47 | 36 | 80 | 114 | 175 | 157 | 193 | 113 | 177 | 115 | 2 | 9 | 62 | 21 | 46 | 58 | 48 | 5 | 1 | 17 | 9 |
| chr1 | 14764 | 24846 | chr1:14764-24846 | 54 | 88 | 17 | 31 | 129 | 64 | 51 | 65 | 108 | 29 | 62 | 61 | 8 | 25 | 20 | 0 | 4 | 0 | 0 | 6 | 5 | 6 | 3 | 0 | 7 | 7 | 5 | 5 | 8 | 10 | 16 | 6 | 5 | 11 | 26 | 8 | 11 | 18 | 31 | 17 | 42 | 15 | 112 | 48 | 20 | 22 | 24 | 10 | 14 | 26 | 13 | 12 | 0 | 0 | 0 | 7 | 21 | 25 | 13 | 17 | 14 | 15 | 18 | 26 | 11 | 16 | 0 | 12 | 6 | 11 | 23 | 23 | 27 | 34 | 22 | 16 | 18 | 4 | 13 | 4 | 1 | 3 | 5 | 0 | 0 | 1 | 1 | 3 |
| chr1 | 14829 | 14970 | chr1:14829-14970 | 1093 | 326 | 502 | 459 | 2985 | 1570 | 2181 | 2385 | 3673 | 1648 | 1947 | 1662 | 450 | 947 | 1135 | 268 | 216 | 124 | 123 | 287 | 284 | 228 | 273 | 448 | 212 | 242 | 156 | 158 | 440 | 443 | 302 | 181 | 207 | 560 | 689 | 219 | 520 | 513 | 1097 | 570 | 731 | 591 | 1630 | 1100 | 610 | 798 | 781 | 730 | 271 | 1119 | 187 | 526 | 128 | 58 | 54 | 98 | 420 | 972 | 345 | 636 | 448 | 455 | 513 | 815 | 285 | 412 | 668 | 205 | 113 | 313 | 776 | 774 | 769 | 972 | 554 | 788 | 705 | 84 | 133 | 181 | 121 | 243 | 203 | 369 | 31 | 90 | 238 | 99 |
| chr1 | 14829 | 15021 | chr1:14829-15021 | 6 | 0 | 0 | 0 | 10 | 5 | 19 | 18 | 25 | 2 | 5 | 7 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 5 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 1 | 5 | 10 | 0 | 3 | 2 | 1 | 1 | 0 | 3 | 2 | 0 | 0 | 0 | 0 | 2 | 1 | 0 | 0 | 7 | 3 | 0 | 0 | 4 | 0 | 3 | 0 | 1 | 2 | 0 | 2 | 0 | 8 | 1 | 0 | 0 | 0 | 2 | 2 | 3 | 5 | 9 | 2 | 4 | 3 | 1 | 0 | 0 | 0 | 3 | 0 | 4 | 0 | 0 | 0 | 0 |
| chr1 | 14829 | 15796 | chr1:14829-15796 | 0 | 3 | 0 | 0 | 4 | 8 | 18 | 33 | 20 | 10 | 33 | 21 | 4 | 3 | 1 | 1 | 0 | 1 | 0 | 0 | 1 | 3 | 0 | 1 | 1 | 1 | 1 | 0 | 7 | 4 | 1 | 0 | 2 | 7 | 43 | 7 | 11 | 4 | 10 | 3 | 0 | 7 | 28 | 6 | 6 | 2 | 0 | 4 | 5 | 0 | 1 | 9 | 0 | 2 | 2 | 10 | 0 | 1 | 0 | 3 | 0 | 10 | 5 | 0 | 5 | 4 | 57 | 0 | 2 | 2 | 4 | 2 | 29 | 16 | 7 | 7 | 39 | 6 | 1 | 0 | 0 | 2 | 0 | 8 | 0 | 0 | 2 | 0 |

Check for differences between counts-filtered junction files

``` r
nrow(deduplicated_combined_junctions_chr1_10)
```

    [1] 31832

``` r
nrow(merged_junctions_chr1_10)
```

    [1] 28332

There are more junctions in the combined bedfile, compared to the merged
bedfile, after filtering for junctions with at least 10 counts in 5
samples.

### Check what’s different in chromosome 1 junctions

``` r
# print number of junctions in chromosome 1 only in combined junctions bedfile
combined_only_chr1 <- dplyr::setdiff(deduplicated_combined_junctions_chr1, merged_junctions_chr1)
nrow(combined_only_chr1)
```

    [1] 5173

``` r
# check fraction of junctions that are only in chr1 of combined junctions bedfile
frac_combined_only_chr1 <- nrow(combined_only_chr1) / nrow(deduplicated_combined_junctions_chr1)
frac_combined_only_chr1
```

    [1] 0.009350406

``` r
# print number of junctions in chromosome 1 only in merged junctions bedfile
merged_only_chr1 <- dplyr::setdiff(merged_junctions_chr1, deduplicated_combined_junctions_chr1)
nrow(merged_only_chr1)
```

    [1] 6599

``` r
# check fraction of merged_only_chr1 junctions that are only in chr1 of merged combined junctions
frac_merged_only_chr1 <- nrow(merged_only_chr1) / nrow(merged_junctions_chr1)
frac_merged_only_chr1
```

    [1] 0.01189729

Overall, only 0.9350406% of all chromosome 1 junctions are only found in
the combined junctions file and 1.1897293% of all chromosome 1 junctions
are only found in the in the merged junctions file.

``` r
# print number of junctions in chromosome 1 only in combined junctions bedfile after filtering for junctions with at least 10 counts in at least 5 samples
combined_only_chr1_10 <- dplyr::setdiff(deduplicated_combined_junctions_chr1_10, merged_junctions_chr1_10)
nrow(combined_only_chr1_10)
```

    [1] 4508

``` r
# check fraction of combined_only_chr1 junctions > 10 counts
nrow(combined_only_chr1_10) / nrow(deduplicated_combined_junctions_chr1)
```

    [1] 0.008148392

``` r
# print number of junctions in chromosome 1 only in merged junctions bedfile after filtering for junctions with at least 10 counts in at least 5 samples
merged_only_chr1_10 <- dplyr::setdiff(merged_junctions_chr1_10, deduplicated_combined_junctions_chr1_10)
nrow(merged_only_chr1_10)
```

    [1] 1008

``` r
# check fraction of merged_chr1 junctions > 10 counts
nrow(merged_only_chr1_10) / nrow(merged_junctions_chr1)
```

    [1] 0.001817316

If we filter for junctions with more than 10 counts, the fraction of
events we miss in the merged junctions file is even smaller in
chromosome 1.

### Print examples of junctions only in merged or combined tables

Print examples of junctions only in the combined junctions file

``` r
combined_only_chr1_10 |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 18369 | 18370 | chr1:18369-18370 | 7 | 10 | 9 | 4 | 7 | 5 | 7 | 7 | 12 | 10 | 2 | 10 | 2 | 1 | 2 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 2 | 2 | 7 | 1 | 0 | 1 | 2 | 1 | 2 | 1 | 6 | 11 | 8 | 1 | 7 | 1 | 0 | 4 | 5 | 0 | 2 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 1 | 0 | 3 | 0 | 1 | 4 | 2 | 0 | 2 | 0 | 1 | 4 | 0 | 0 | 3 | 1 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 3 | 4 | 1 |
| chr1 | 188847 | 188848 | chr1:188847-188848 | 10 | 9 | 24 | 4 | 8 | 10 | 9 | 28 | 11 | 12 | 3 | 19 | 2 | 0 | 9 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 2 | 0 | 0 | 6 | 5 | 6 | 2 | 2 | 1 | 5 | 13 | 2 | 4 | 3 | 3 | 1 | 11 | 7 | 8 | 1 | 0 | 9 | 0 | 0 | 0 | 1 | 0 | 0 | 2 | 1 | 0 | 9 | 3 | 4 | 0 | 4 | 0 | 5 | 6 | 2 | 0 | 0 | 0 | 6 | 1 | 2 | 1 | 2 | 2 | 0 | 1 | 0 | 0 | 1 | 5 | 6 | 2 | 2 | 7 | 2 |
| chr1 | 188889 | 188890 | chr1:188889-188890 | 1 | 2 | 18 | 1 | 17 | 7 | 14 | 13 | 11 | 21 | 5 | 16 | 2 | 3 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 2 | 0 | 0 | 1 | 2 | 0 | 0 | 0 | 0 | 2 | 3 | 0 | 1 | 0 | 2 | 0 | 5 | 3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 2 | 0 |
| chr1 | 190809 | 190810 | chr1:190809-190810 | 36 | 16 | 50 | 57 | 73 | 37 | 35 | 112 | 31 | 79 | 37 | 98 | 17 | 1 | 6 | 14 | 5 | 12 | 16 | 6 | 5 | 8 | 7 | 16 | 11 | 1 | 3 | 14 | 17 | 6 | 5 | 12 | 4 | 12 | 18 | 5 | 9 | 18 | 10 | 7 | 5 | 4 | 2 | 0 | 1 | 0 | 8 | 0 | 0 | 34 | 0 | 17 | 0 | 2 | 0 | 0 | 3 | 0 | 0 | 4 | 0 | 8 | 4 | 2 | 0 | 0 | 11 | 5 | 0 | 3 | 7 | 15 | 24 | 2 | 28 | 15 | 43 | 1 | 3 | 4 | 14 | 5 | 21 | 32 | 4 | 11 | 23 | 11 |
| chr1 | 268020 | 268021 | chr1:268020-268021 | 40 | 4 | 18 | 31 | 5 | 22 | 4 | 7 | 5 | 7 | 13 | 4 | 4 | 7 | 6 | 0 | 8 | 5 | 0 | 3 | 10 | 3 | 8 | 3 | 4 | 1 | 7 | 1 | 3 | 0 | 0 | 4 | 7 | 2 | 16 | 0 | 11 | 3 | 8 | 6 | 1 | 0 | 9 | 7 | 6 | 12 | 27 | 3 | 0 | 6 | 1 | 7 | 0 | 0 | 2 | 2 | 16 | 6 | 5 | 6 | 8 | 3 | 6 | 6 | 6 | 6 | 0 | 0 | 0 | 1 | 13 | 0 | 5 | 5 | 11 | 4 | 2 | 0 | 1 | 0 | 1 | 0 | 3 | 0 | 0 | 0 | 24 | 13 |
| chr1 | 268121 | 268122 | chr1:268121-268122 | 7 | 4 | 18 | 37 | 35 | 24 | 5 | 6 | 18 | 8 | 14 | 8 | 3 | 4 | 6 | 5 | 10 | 14 | 0 | 4 | 26 | 6 | 21 | 9 | 3 | 5 | 11 | 8 | 8 | 1 | 0 | 4 | 4 | 1 | 33 | 0 | 20 | 5 | 11 | 9 | 0 | 0 | 14 | 4 | 1 | 5 | 29 | 5 | 3 | 10 | 3 | 0 | 3 | 0 | 1 | 10 | 29 | 6 | 4 | 8 | 8 | 18 | 2 | 16 | 1 | 13 | 0 | 1 | 0 | 3 | 33 | 0 | 4 | 4 | 5 | 11 | 11 | 1 | 0 | 0 | 0 | 4 | 4 | 0 | 1 | 0 | 18 | 17 |

For comparison, examine merged chromosome 1 counts-filtered junction
table for junctions near chr1:18369-18370

``` r
merged_junctions_chr1_10 |>
  dplyr::arrange(start) |>
  dplyr::filter((18369 -  10) < start &
                  start < (18369 + 10) ) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 18362 | 18497 | chr1:18362-18497 | 22 | 3 | 9 | 9 | 23 | 10 | 20 | 44 | 29 | 15 | 16 | 38 | 3 | 5 | 6 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 16 | 10 | 0 | 6 | 25 | 0 | 0 | 13 | 0 | 25 | 12 | 16 | 14 | 0 | 16 | 17 | 0 | 18 | 14 | 6 | 25 | 5 | 6 | 0 | 0 | 0 | 0 | 8 | 25 | 0 | 15 | 0 | 19 | 5 | 11 | 13 | 0 | 0 | 5 | 0 | 0 | 19 | 0 | 16 | 14 | 7 | 0 | 7 | 0 | 3 | 0 | 2 | 0 | 0 | 0 | 0 | 2 | 6 | 2 |
| chr1 | 18366 | 18913 | chr1:18366-18913 | 184 | 95 | 102 | 90 | 317 | 110 | 308 | 459 | 228 | 179 | 305 | 162 | 51 | 26 | 33 | 16 | 5 | 1 | 8 | 8 | 9 | 6 | 4 | 7 | 6 | 3 | 2 | 4 | 23 | 36 | 38 | 19 | 23 | 66 | 32 | 16 | 46 | 20 | 47 | 23 | 10 | 38 | 113 | 30 | 45 | 13 | 42 | 29 | 6 | 53 | 3 | 65 | 18 | 9 | 0 | 0 | 23 | 43 | 11 | 41 | 6 | 11 | 10 | 13 | 12 | 11 | 15 | 14 | 4 | 11 | 31 | 32 | 39 | 60 | 36 | 23 | 17 | 5 | 16 | 7 | 34 | 27 | 19 | 29 | 1 | 15 | 39 | 31 |
| chr1 | 18366 | 24738 | chr1:18366-24738 | 177 | 174 | 152 | 130 | 449 | 341 | 375 | 542 | 354 | 500 | 401 | 479 | 108 | 179 | 325 | 86 | 93 | 53 | 62 | 199 | 153 | 141 | 102 | 175 | 122 | 108 | 50 | 113 | 131 | 117 | 121 | 82 | 66 | 205 | 175 | 70 | 257 | 258 | 329 | 231 | 173 | 154 | 115 | 99 | 109 | 181 | 232 | 123 | 48 | 204 | 68 | 156 | 105 | 51 | 16 | 7 | 188 | 194 | 111 | 194 | 99 | 91 | 127 | 187 | 79 | 111 | 172 | 20 | 3 | 46 | 139 | 74 | 138 | 230 | 40 | 228 | 105 | 15 | 19 | 44 | 47 | 59 | 66 | 187 | 13 | 39 | 81 | 43 |
| chr1 | 18366 | 29321 | chr1:18366-29321 | 3 | 11 | 0 | 4 | 24 | 6 | 14 | 12 | 24 | 10 | 8 | 9 | 0 | 3 | 12 | 0 | 0 | 0 | 0 | 13 | 0 | 3 | 0 | 0 | 0 | 0 | 3 | 0 | 11 | 0 | 5 | 6 | 0 | 8 | 11 | 9 | 6 | 24 | 10 | 10 | 5 | 8 | 0 | 4 | 12 | 8 | 11 | 34 | 0 | 8 | 2 | 8 | 2 | 3 | 2 | 1 | 10 | 2 | 4 | 0 | 15 | 8 | 0 | 26 | 0 | 21 | 22 | 2 | 0 | 4 | 17 | 11 | 29 | 56 | 17 | 33 | 32 | 0 | 5 | 0 | 0 | 1 | 3 | 2 | 2 | 0 | 0 | 0 |
| chr1 | 18369 | 18501 | chr1:18369-18501 | 32 | 10 | 2 | 19 | 36 | 18 | 52 | 42 | 32 | 15 | 18 | 46 | 9 | 9 | 8 | 0 | 0 | 1 | 2 | 0 | 1 | 0 | 1 | 2 | 0 | 0 | 2 | 0 | 2 | 4 | 3 | 2 | 4 | 4 | 0 | 3 | 7 | 3 | 10 | 7 | 7 | 3 | 3 | 1 | 6 | 1 | 14 | 6 | 0 | 6 | 2 | 4 | 0 | 0 | 0 | 0 | 2 | 0 | 1 | 11 | 4 | 2 | 3 | 4 | 0 | 3 | 3 | 1 | 0 | 1 | 5 | 7 | 4 | 10 | 1 | 5 | 0 | 1 | 3 | 0 | 0 | 2 | 0 | 3 | 0 | 1 | 3 | 7 |
| chr1 | 18369 | 18913 | chr1:18369-18913 | 40 | 50 | 77 | 43 | 163 | 56 | 174 | 256 | 117 | 76 | 191 | 94 | 24 | 18 | 26 | 0 | 2 | 0 | 6 | 4 | 8 | 4 | 4 | 0 | 10 | 12 | 0 | 0 | 9 | 6 | 15 | 8 | 9 | 36 | 9 | 19 | 19 | 5 | 18 | 8 | 3 | 18 | 45 | 16 | 25 | 8 | 15 | 7 | 4 | 33 | 2 | 19 | 5 | 3 | 0 | 3 | 0 | 17 | 5 | 10 | 0 | 8 | 4 | 4 | 2 | 4 | 0 | 2 | 1 | 2 | 13 | 11 | 12 | 18 | 39 | 9 | 0 | 3 | 5 | 9 | 15 | 5 | 11 | 16 | 3 | 4 | 12 | 17 |

This junction is completely missed, with the closest coordinate being
chr1:18369-18913.

Print examples of junctions only in the merged junctions file

``` r
merged_only_chr1_10 |>
  dplyr::arrange(start) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 268020 | 268021 | chr1:268020-268021 | 40 | 4 | 18 | 31 | 5 | 22 | 4 | 7 | 5 | 7 | 13 | 4 | 4 | 7 | 6 | 0 | 8 | 5 | 0 | 0 | 10 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 4 | 7 | 2 | 16 | 0 | 11 | 3 | 8 | 6 | 1 | 0 | 9 | 7 | 6 | 12 | 27 | 3 | 0 | 6 | 1 | 7 | 0 | 0 | 2 | 2 | 0 | 6 | 0 | 6 | 8 | 3 | 6 | 6 | 6 | 6 | 0 | 0 | 0 | 1 | 13 | 0 | 5 | 5 | 11 | 4 | 2 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 24 | 13 |
| chr1 | 268121 | 268122 | chr1:268121-268122 | 7 | 4 | 18 | 37 | 35 | 24 | 5 | 6 | 18 | 8 | 14 | 8 | 3 | 4 | 6 | 0 | 10 | 14 | 0 | 0 | 26 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 8 | 1 | 0 | 4 | 4 | 1 | 33 | 0 | 20 | 5 | 11 | 9 | 0 | 0 | 14 | 4 | 1 | 5 | 29 | 5 | 0 | 10 | 3 | 0 | 3 | 0 | 1 | 10 | 0 | 6 | 0 | 8 | 8 | 18 | 2 | 16 | 1 | 13 | 0 | 1 | 0 | 3 | 33 | 0 | 4 | 4 | 5 | 11 | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 18 | 17 |
| chr1 | 961552 | 961553 | chr1:961552-961553 | 0 | 0 | 0 | 25 | 0 | 12 | 121 | 96 | 75 | 0 | 0 | 0 | 20 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 28 | 0 | 0 | 14 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 961628 | 961629 | chr1:961628-961629 | 0 | 0 | 0 | 25 | 0 | 21 | 157 | 107 | 77 | 0 | 0 | 0 | 21 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 13 | 0 | 0 | 12 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 999613 | 999614 | chr1:999613-999614 | 0 | 78 | 80 | 43 | 210 | 341 | 157 | 132 | 168 | 274 | 589 | 317 | 100 | 0 | 0 | 215 | 292 | 140 | 185 | 282 | 322 | 565 | 674 | 489 | 580 | 102 | 134 | 123 | 62 | 29 | 0 | 0 | 0 | 13 | 0 | 0 | 0 | 0 | 30 | 265 | 110 | 65 | 1834 | 89 | 53 | 51 | 33 | 65 | 18 | 26 | 50 | 0 | 0 | 0 | 0 | 0 | 71 | 46 | 13 | 182 | 105 | 84 | 24 | 16 | 106 | 115 | 9 | 36 | 4 | 0 | 46 | 29 | 33 | 0 | 0 | 0 | 0 | 0 | 0 | 67 | 189 | 0 | 177 | 259 | 0 | 0 | 0 | 0 |
| chr1 | 999691 | 999692 | chr1:999691-999692 | 0 | 276 | 103 | 68 | 165 | 243 | 149 | 149 | 141 | 266 | 540 | 306 | 104 | 0 | 0 | 160 | 375 | 132 | 299 | 375 | 322 | 558 | 740 | 446 | 627 | 111 | 112 | 136 | 26 | 13 | 0 | 0 | 0 | 24 | 0 | 0 | 0 | 0 | 26 | 196 | 76 | 56 | 1300 | 65 | 53 | 33 | 9 | 41 | 17 | 23 | 54 | 0 | 0 | 0 | 0 | 0 | 57 | 36 | 7 | 114 | 66 | 38 | 18 | 6 | 61 | 61 | 22 | 27 | 3 | 0 | 27 | 20 | 22 | 0 | 0 | 0 | 0 | 0 | 3 | 68 | 122 | 0 | 188 | 173 | 0 | 0 | 0 | 0 |

For comparison, examine combined chromosome 1 counts-filtered junction
table for junctions near chr1:268020-268021

``` r
deduplicated_combined_junctions_chr1_10 |>
  dplyr::arrange(start) |>
  dplyr::filter((268020 -  10) < start &
                  start < (268020 + 10) ) |>
  head()
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 268020 | 268021 | chr1:268020-268021 | 40 | 4 | 18 | 31 | 5 | 22 | 4 | 7 | 5 | 7 | 13 | 4 | 4 | 7 | 6 | 0 | 8 | 5 | 0 | 3 | 10 | 3 | 8 | 3 | 4 | 1 | 7 | 1 | 3 | 0 | 0 | 4 | 7 | 2 | 16 | 0 | 11 | 3 | 8 | 6 | 1 | 0 | 9 | 7 | 6 | 12 | 27 | 3 | 0 | 6 | 1 | 7 | 0 | 0 | 2 | 2 | 16 | 6 | 5 | 6 | 8 | 3 | 6 | 6 | 6 | 6 | 0 | 0 | 0 | 1 | 13 | 0 | 5 | 5 | 11 | 4 | 2 | 0 | 1 | 0 | 1 | 0 | 3 | 0 | 0 | 0 | 24 | 13 |

Position chr1:268020-268021 is actually a junction recorded in the
combined junctions table. The difference come in the counts going into
each sample:

| sample | junction count in merged table | junction count in combined table |
|----|----|----|
| **SRR1712457** | 0 | 3 |
| **SRR1712459** | 0 | 3 |
| **SRR1712460** | 0 | 8 |
| **SRR1712461** | 0 | 3 |
| **SRR1712462** | 0 | 4 |
| … | … | … |
| **SRR1799041** | 0 | 16 |

From scanning chr1:268020-268021, the counts differences arise from
cases where the merged junctions counts have a value of 0 for a sample
and the combined file is not zero for the same sample. For all cases
except for SRR1799041, the counts that are dropped in the merged table
are \< 10. Because the values in the merged table are 0 for the counts
that differ from the combined table, I suspect these differences stem
from junctions that are not shared across samples (these appear as NA
values which we replace with 0 in the final merged bedfile). I think
this means that junctions defined differently when Shiba runs on one
vs. multiple samples in a way that I don’t understand yet. But since
this difference impacts a very small fraction of junctions in chromosome
1, maybe it is not something worth trying to fix.

Another note - because the concatenated junctions bug impacted one
junction with reads \> 10 in the combined junctions bedfile, I expect
this will result in differences in the PSI table
