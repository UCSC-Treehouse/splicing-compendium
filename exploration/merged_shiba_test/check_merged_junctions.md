# Debug merge junctions warning
Cindy Liang (celiang@ucsc.edu)
2026-05-26

## Files and directories

``` r
# find the root-level repo directory so we can access the other files
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

## define the data directories ##
exploration_dir <- file.path(repo_root, "exploration")
results_dir <- file.path(repo_root, "results", "merged_shiba", "test") # manually merged junction and gtf files and shiba results
merged_test_dir <- file.path(exploration_dir, "merged_shiba_test")
shiba_results_dir <- file.path(merged_test_dir, "shiba_results")

## data files ##
# combined run shiba results, junctions, and gtf ("control" results and files)
combined_junctions_file <- file.path(shiba_results_dir, "combined_run", "junctions", "junctions.bed") 

# junction bedfiles generated from separate shiba runs
# the way the paths are read in could probably be optimized but I don't know if that's worth doing now
# since we plan to read in paths from snakemake
SRR4376025_bedfile <- file.path(repo_root, "results/target/shiba/SRR4376025/junctions/junctions.bed")
SRR1559031_bedfile <- file.path(repo_root, "results/target/shiba/SRR1559031/junctions/junctions.bed")

# Read in junctions bed paths into a vector
junction_paths <- c(SRR4376025_bedfile, 
                   SRR1559031_bedfile)
```

## Print problems in each of the junction files

### Check separate junction bedfiles for problems

Check SRR4376025 bedfile for problems:

``` r
SRR4376025_problems <- vroom::vroom(SRR4376025_bedfile)
```

    Rows: 498794 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): chr, ID
    dbl (3): start, end, SRR4376025

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
vroom::problems(SRR4376025_problems) |> table()
```

    < table of extent 0 x 0 x 0 x 0 x 0 >

No problems in the SRR4376025 bedfile, though we may want to set the
default column type for numeric values as double (d), instead of integer
(i)

Check SRR1559031 bedfile for problems:

``` r
SRR1559031_problems <- vroom::vroom(SRR1559031_bedfile)
```

    Rows: 452467 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): chr, ID
    dbl (3): start, end, SRR1559031

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
vroom::problems(SRR1559031_problems)
```

    Warning: One or more parsing issues, call `problems()` on your data frame for details,
    e.g.:
      dat <- vroom(...)
      problems(dat)

| row | col | expected | actual | file |
|---:|---:|:---|:---|:---|
| 434393 | 2 | a double | 35660647;35660647 | /Users/cindyliang/repos/splicing-compendium/results/target/shiba/SRR1559031/junctions/junctions.bed |
| 434393 | 3 | a double | 35660648;35660648 | /Users/cindyliang/repos/splicing-compendium/results/target/shiba/SRR1559031/junctions/junctions.bed |

One source of the merging problem is in row 434393 in SRR1559031, where
the position coordinates seem to be duplicated

Check what line looks like in the bedfile:

``` r
readLines(file.path(repo_root, "results/target/shiba/SRR1559031/junctions/junctions.bed"))[434393]
```

    [1] "chr9;chr9\t35660647;35660647\t35660648;35660648\tchr9:35660647-35660648\t57"

Check if problem junction also had a “correct” entry:

``` r
SRR1559031_problems |>
  dplyr::filter(start == "35660647",
                end == "35660648")
```

| chr | start | end | ID  | SRR1559031 |
|:----|------:|----:|:----|-----------:|

No, when the junction is in incorrect bed format, there is no
corresponding duplicate ‘fixed’ version.

In the problem junction, the “chr” field also has duplicated values.
This is event ID chr9:35660647-35660648, which has “chr9;chr9” for the
chromosome

``` r
SRR1559031_problems |> dplyr::slice(434391, 434392, 434393)
```

| chr       |    start |      end | ID                     | SRR1559031 |
|:----------|---------:|---------:|:-----------------------|-----------:|
| chr9      | 99886674 | 99886675 | chr9:99886674-99886675 |         10 |
| chr9;chr9 |       NA |       NA | chr9:35660647-35660648 |         57 |
| chrM      |       26 |    11767 | chrM:26-11767          |          2 |

chr9:35660647-35660648, our problem junction, also has high read counts
( \> 10 ) so it seems like an important junction to be able to keep.
Maybe a workaround would be to manually fix the chr/start/end
coordinates of these cases by extracting them from the ID field,
assuming these cases are always duplicates…

## Check if the combined junctions.bed has this weird junction coordinate

Check for problems in combined junctions.bed

``` r
combined_jcn_problems <- vroom::vroom(combined_junctions_file) 
```

    Rows: 741591 Columns: 6
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): chr, ID
    dbl (4): start, end, SRR1559031, SRR4376025

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
vroom::problems(combined_jcn_problems) |> table()
```

    < table of extent 0 x 0 x 0 x 0 x 0 >

The problem is not in the combined junction bedfile.

Check for the specific junction in the combined junctions.bed

``` r
combined_jcn_problems |>
  dplyr::filter(ID == "chr9:35660647-35660648")
```

| chr  |    start |      end | ID                     | SRR1559031 | SRR4376025 |
|:-----|---------:|---------:|:-----------------------|-----------:|-----------:|
| chr9 | 35660647 | 35660648 | chr9:35660647-35660648 |         57 |          1 |

The event looks normal in the combined file. Importantly, the ID is the
same: chr9:35660647-35660648

## Check problem junction before fix

``` r
# read in junctions.bed files created from separate Shiba runs
merged_junctions <- purrr::map(junction_paths, \(file) {
    readr::read_tsv(file,
               # make sure columns are types that we expect
               col_types = readr::cols(.default = "d", ID = "c", chr = "c", start = "c", end = "c"))
  })|>
    # merge junctions tables from multiple samples
    # the resulting table separates junction counts from each sample by columns with the sample ID
    purrr::reduce(\(x, y) dplyr::full_join(x, y, by = c("chr", "start", "end", "ID")))
```

    Registered S3 methods overwritten by 'readr':
      method                    from 
      as.data.frame.spec_tbl_df vroom
      as_tibble.spec_tbl_df     vroom
      format.col_spec           vroom
      print.col_spec            vroom
      print.collector           vroom
      print.date_names          vroom
      print.locale              vroom
      str.col_spec              vroom

``` r
# check problem junction
merged_junctions |>
  dplyr::filter(ID == "chr9:35660647-35660648")
```

| chr | start | end | ID | SRR4376025 | SRR1559031 |
|:---|:---|:---|:---|---:|---:|
| chr9 | 35660647 | 35660648 | chr9:35660647-35660648 | 1 | NA |
| chr9;chr9 | 35660647;35660647 | 35660648;35660648 | chr9:35660647-35660648 | NA | 57 |

With how we currently merge junctions, the problem junction appears as
two separate events that are missing values in the other samples

## Proposed new junction merging method:

Use ID to extract chr, start, and end fields in events with weird
concatenated fields

``` r
# read in junctions.bed files created from separate Shiba runs
merged_junctions <- purrr::map(junction_paths, \(file) {
  # read in separate bedfiles for each sample
  withCallingHandlers(
  readr::read_tsv(
    file,
    # make sure columns are types that we expect
    col_types = readr::cols(
      .default = "d",
      ID = "c",
      chr = "c",
      start = "d",
      end = "d"
    )
  ),
  # convert warnings in reading bedfiles into a failure and stop
  # warnings here indicate there are still junctions remaining in invalid format (e.g. "12345;12345") 
  warning = \(w) stop(w)
  )
  }
  ) |>
  # merge junctions tables from multiple samples
  # the resulting table separates junction counts from each sample by columns with the sample ID
  purrr::reduce(
    \(x, y) dplyr::full_join(
      x, y,
      by = c("chr", "start", "end", "ID")
    )
  )
```

    Warning: One or more parsing issues, call `problems()` on your data frame for details,
    e.g.:
      dat <- vroom(...)
      problems(dat)

Check position ID of fixed chr/start/end entry

``` r
merged_junctions |>
  dplyr::filter(ID == "chr9:35660647-35660648")
```

| chr       |    start |      end | ID                     | SRR4376025 | SRR1559031 |
|:----------|---------:|---------:|:-----------------------|-----------:|-----------:|
| chr9      | 35660647 | 35660648 | chr9:35660647-35660648 |          1 |         NA |
| chr9;chr9 |       NA |       NA | chr9:35660647-35660648 |         NA |         57 |

In this case the junctions stay unmerged because the merging code
stopped with a warning from the invalid entry
