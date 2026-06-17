# Check Merged vs. Combined TARGET Pilot Junction files
Cindy Liang (celiang@ucsc.edu)
2026-06-17

## Background

This notebook checks junctions.bed files from two Shiba run methods for
differences. The goal is to get a sense of how close files that go into
calculating PSI values for events are to each other in the merged and
combined Shiba methods. Junctions files that are different would impact
PSI calculation, so I am mainly looking to see if mismatched junctions
are present and if present, how much of the data do they represent.

Files compared in this notebook are junction count bedfiles of 88 TARGET
pilot samples, produced from running Shiba with two different methods.

“Combined” refers to shiba splice analysis run with shiba in a canonical
way. Inputs of this method are .bam files of each sample and a reference
annotation GTF. In the combined method, unmodified Shiba scripts are run
on the input files using shiba.py, which calls on Shiba scripts that
generate intermediate files and perform the final PSI calculation. The
intermediate files that go into PSI calculation from this method are:

- GTF, which is made with
  \[bam2gtf.py\](https://github.com/Sika-Zheng-Lab/Shiba/blob/v0.8.1/src/bam2gtf.py.
  This script uses a stringtie merge command as part of that merges GTFs
  made from individual samples together with the reference GTF

- junctions.bed, which is produced by
  [bam2junc.py](https://github.com/Sika-Zheng-Lab/Shiba/blob/v0.8.1/src/bam2junc.py).
  This script counts exon-exon (regtools) and exon-intron
  (featureCounts) junctions from each bam file given as input.

  - To count exon-intron junctions, Shiba uses `create_saf_file()`
    within `bam2junc.py` to create annotations of intron/exon boundaries
    to pass to featureCounts. The exon-intron boundaries are taken from
    the `ri_event` file generated from `gtf2events.py`.

  Then, `bam2junc.py` checks if exon-exon junctions are duplicated and
  deduplicates them. There is a bug in the deduplication step, as this
  notebook will later show junctions with invalid coordinates in the
  combined junctions file. Finally, junctions files containing exon-exon
  and exon-intron junctions from all samples are merged together.

“Merged” refers to shiba splice analysis run with `merge_results.smk`.
The inputs of this method are junctions.bed and GTF files created by
Shiba run on each sample separately. `merge_results.smk` merges the GTFs
using the same stringtie merge command Shiba uses, then manually
combines the junction files after deduplicating any invalid junction
entries resulting from the `bam2junc.py` bug. The jobs for the merged
Shiba method were run in the following order, to generate the merged GTF
and junctions.bed files:

    [Mon Jun  1 03:46:24 2026]
    Finished jobid: 4 (Rule: merge_junctions)
    1 of 5 steps (20%) done
    [Mon Jun  1 04:07:03 2026]
    Finished jobid: 3 (Rule: merge_gtfs)
    2 of 5 steps (40%) done
    Select jobs to execute...
    Execute 1 jobs...
    [Mon Jun  1 04:07:03 2026]
    localrule gtf_to_events:
        input: results/merged_shiba/target_pilot/merged_gtf.gtf, references/gencode.v47.primary_assembly.annotation.gtf
        output: results/merged_shiba/target_pilot/events
        jobid: 2
        reason: Missing output files: results/merged_shiba/target_pilot/events; Input files updated by another job: results/merged_shiba/target_pilot/merged_gtf.gtf
        priority: 1
        threads: 10
        resources: tmpdir=/tmp
    [Mon Jun  1 04:43:48 2026]
    Finished jobid: 2 (Rule: gtf_to_events)

Because the junctions.bed files have exon-intron coordinates defined by
retained intron events files obtained from one sample only, some samples
will not have exon-intron counts from exon-intron boundaries that are
only present in other samples.

### Analysis outline:

Check differences in full junctions files - Check how many rows
(junction IDs) are in each bedfile. Also calculate fraction of junctions
in the merged junctions table that are not in the combined table; the
fraction of junctions in the combined table that are missed in the
merged junctions table

Check differences in junctions with at least 10 junctions in
min_samples - check if differences change if we filter for junctions
with higher counts in at least 5 samples. Filter merged and combined
junctions tables for junctions with counts of at least 10 in at least 5
samples. Calculate fraction of junctions in the merged junctions table
that are not in the combined table; the fraction of junctions in the
combined table that are missed in the merged junctions table

We discovered that the largest source of difference between the merged
and combined method lies in how exon-intron junctions are defined and
counted. Shiba uses the retained intron events coordinates generated
from gtf2events.py to define intron-exon coordinates for counting. Since
the merged method only uses the events coordinates defined from one
sample at a time, junctions counted with this method are expected to
miss all exon-intron junctions that are not shared between samples.
Coming into this analysis, I expect all junction IDs that are different
between the methods to be exon-intron junctions.

To see if this is the case, this notebook also:

- Checks if junctions that are different have a length of 1

- Checks if all junctions whose counts values are different (but IDs are
  the same) betweeen the tables are different because the value in the
  merged table is 0 and the value in the combined table is nonzero

## Define functions

``` r
# pivot junctions dataframe long and add summary column of number of samples with counts exceeding min_counts
pivot_junctions_long <- function(df) {
 df |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("SRR"),
    names_to = "sample",
    values_to = "count"
  ) |>
  dplyr::group_by(ID) |>
  # add column that counts how many samples have count values over 10 for each position
  dplyr::mutate(
    samples_over_min_count = sum(count >= params$min_counts)
  )
}
```

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
# target pilot exploration dir
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

Check that column names are the same in both dataframes

``` r
colnames(deduplicated_combined_junctions)
```

     [1] "chr"        "start"      "end"        "ID"         "SRR1559043"
     [6] "SRR1559044" "SRR1559052" "SRR1559054" "SRR1559075" "SRR1559100"
    [11] "SRR1559105" "SRR1559133" "SRR1559134" "SRR1559145" "SRR1559160"
    [16] "SRR1559164" "SRR1559177" "SRR1559183" "SRR1559184" "SRR1712453"
    [21] "SRR1712454" "SRR1712455" "SRR1712456" "SRR1712457" "SRR1712458"
    [26] "SRR1712459" "SRR1712460" "SRR1712461" "SRR1712462" "SRR1712463"
    [31] "SRR1712464" "SRR1712465" "SRR1784865" "SRR1784867" "SRR1791016"
    [36] "SRR1791028" "SRR1791108" "SRR1796863" "SRR1796867" "SRR1796893"
    [41] "SRR1796906" "SRR1796912" "SRR1796939" "SRR1796967" "SRR1796990"
    [46] "SRR1797014" "SRR1797024" "SRR1797033" "SRR1797034" "SRR1797035"
    [51] "SRR1797039" "SRR1797052" "SRR1797053" "SRR1797055" "SRR1797057"
    [56] "SRR1797087" "SRR1797107" "SRR1797111" "SRR1799022" "SRR1799025"
    [61] "SRR1799041" "SRR1799042" "SRR1799057" "SRR1799058" "SRR1799059"
    [66] "SRR1799061" "SRR1799062" "SRR1799067" "SRR1799069" "SRR1799081"
    [71] "SRR1810588" "SRR2042833" "SRR2042845" "SRR2042853" "SRR2042854"
    [76] "SRR2042856" "SRR2083154" "SRR2083162" "SRR2083171" "SRR2083176"
    [81] "SRR2083188" "SRR2239703" "SRR2239717" "SRR3162160" "SRR3162195"
    [86] "SRR3162212" "SRR3162237" "SRR3162253" "SRR4376029" "SRR4416297"
    [91] "SRR4419554" "SRR4419565"

``` r
colnames(merged_junctions)
```

     [1] "chr"        "start"      "end"        "ID"         "SRR1559043"
     [6] "SRR1559145" "SRR1559160" "SRR1559164" "SRR1559177" "SRR1559183"
    [11] "SRR1559184" "SRR1712453" "SRR1712454" "SRR1712455" "SRR1712456"
    [16] "SRR1559044" "SRR1712457" "SRR1712458" "SRR1712459" "SRR1712460"
    [21] "SRR1712461" "SRR1712462" "SRR1712463" "SRR1712464" "SRR1712465"
    [26] "SRR1784865" "SRR1559052" "SRR1784867" "SRR1791016" "SRR1791028"
    [31] "SRR1791108" "SRR1796863" "SRR1796867" "SRR1796893" "SRR1796906"
    [36] "SRR1796912" "SRR1796939" "SRR1559054" "SRR1796967" "SRR1796990"
    [41] "SRR1797014" "SRR1797024" "SRR1797033" "SRR1797034" "SRR1797035"
    [46] "SRR1797039" "SRR1797052" "SRR1797053" "SRR1559075" "SRR1797055"
    [51] "SRR1797057" "SRR1797087" "SRR1797107" "SRR1797111" "SRR1799022"
    [56] "SRR1799025" "SRR1799041" "SRR1799042" "SRR1799057" "SRR1559100"
    [61] "SRR1799058" "SRR1799059" "SRR1799061" "SRR1799062" "SRR1799067"
    [66] "SRR1799069" "SRR1799081" "SRR1810588" "SRR2042833" "SRR2042845"
    [71] "SRR1559105" "SRR2042853" "SRR2042854" "SRR2042856" "SRR2083154"
    [76] "SRR2083162" "SRR2083171" "SRR2083176" "SRR2083188" "SRR2239703"
    [81] "SRR2239717" "SRR1559133" "SRR3162160" "SRR3162195" "SRR3162212"
    [86] "SRR3162237" "SRR3162253" "SRR4376029" "SRR4416297" "SRR4419554"
    [91] "SRR4419565" "SRR1559134"

Column names are the same, save for the order of the sample junction
counts columns

## Check tables for parsing errors

### Check for parsing errors in the merged and deduplicated junctions tables

``` r
readr::problems(merged_junctions)
```

| row | col | expected | actual | file |
|----:|----:|:---------|:-------|:-----|

``` r
readr::problems(deduplicated_combined_junctions)
```

| row | col | expected | actual | file |
|----:|----:|:---------|:-------|:-----|

No problems are present in merged_junctions file or deduplicated
junctions file (readr::problems prints out nothing)

## Check differences in junctions bedfiles

### Subset junctions files

``` r
# Subset junctions files if a list of chromosomes are provided, or if one chromosome is provided in params

if(length(params$chromosome) > 1 || tolower(params$chromosome) != "all"){
  deduplicated_combined_junctions <- deduplicated_combined_junctions |>
    dplyr::filter(chr %in% params$chromosome)
}

if(length(params$chromosome) > 1 || tolower(params$chromosome) != "all"){
  merged_junctions <- merged_junctions |>
    dplyr::filter(chr %in% params$chromosome)
}

# Ensure sample columns are in the same order
merged_junctions <- merged_junctions |> dplyr::relocate(chr, start, end, ID, dplyr::all_of(target_pilot_samples))
```

### Check differences in junctions files

#### Count all differences between each data frame

``` r
# print number of junctions in only in combined junctions bedfile
combined_only <- dplyr::setdiff(deduplicated_combined_junctions, merged_junctions)
num_combined_only <- nrow(combined_only)

# check fraction of junctions that are only in the combined junctions bedfile
frac_combined_only <- num_combined_only / nrow(deduplicated_combined_junctions)

# print number of junctions only in merged junctions bedfile
merged_only <- dplyr::setdiff(merged_junctions, deduplicated_combined_junctions)
num_merged_only <- nrow(merged_only)

# check fraction of junctions that are only in merged junctions
frac_merged_only <- num_merged_only / nrow(merged_junctions)

# print results
summary_message <- c(
  paste0("Number of rows different in combined junctions: ", num_combined_only),
  paste0("Fraction of rows different in combined junctions: ", frac_combined_only),
  paste0("Number of rows different in merged junctions: ", num_merged_only),
  paste0("Fraction of rows different in merged junctions: ", frac_merged_only)
  )

summary_message
```

    [1] "Number of rows different in combined junctions: 5173"                 
    [2] "Fraction of rows different in combined junctions: 0.00935040615431333"
    [3] "Number of rows different in merged junctions: 6599"                   
    [4] "Fraction of rows different in merged junctions: 0.011897292775446"    

There are 1426 more junctions in the merged junction file compared to
the combined junction file.

0.9350406% of the combined method junctions are missed in the merged
junctions file.

1.1897293% of the merged method junctions are not in the combined
junctions file.

#### Count differences in junction position IDs

``` r
# print number of junction IDs only in combined junctions bedfile
combined_only_ids <- setdiff(deduplicated_combined_junctions$ID, merged_junctions$ID)
# print number of junctions only in the combined dataframe
num_combined_only <- length(combined_only_ids)

# check fraction of combined_only junction IDs
frac_combined_only <- num_combined_only / length(deduplicated_combined_junctions$ID)

# print number of junction IDs only in merged junctions bedfile
merged_only_ids <- setdiff(merged_junctions$ID, deduplicated_combined_junctions$ID)
# print number of unique junctions with any differences between the two dataframes
num_merged_only <- length(merged_only_ids)

# check fraction of filtered merged_only junction IDs
frac_merged_only <- num_merged_only / length(merged_junctions$ID)

# print results
summary_message <- c(
  paste0("Number of rows different in combined junctions: ", num_combined_only),
  paste0("Fraction of rows different in combined junctions: ", frac_combined_only),
  paste0("Number of rows different in merged junctions: ", num_merged_only),
  paste0("Fraction of rows different in merged junctions: ", frac_merged_only)
  )

summary_message
```

    [1] "Number of rows different in combined junctions: 974"                  
    [2] "Fraction of rows different in combined junctions: 0.00176054428654575"
    [3] "Number of rows different in merged junctions: 2400"                   
    [4] "Fraction of rows different in merged junctions: 0.00432694387953788"  

0.1760544% of the combined method junction IDs are missed in the merged
junctions file.

0.4326944% of the merged method junction IDs are not in the combined
junctions file.

### Check differences in junctions with at least 10 junctions in min_samples

Pivot data frame longer and count the number samples with junctions
above min_junctions that are present

``` r
long_combined_junctions <- pivot_junctions_long(deduplicated_combined_junctions)
long_merged_junctions <- pivot_junctions_long(merged_junctions)

# filter data frames for IDs where there are min_samples with min_count
filtered_long_combined <- long_combined_junctions |>
  dplyr::filter(samples_over_min_count >= params$min_samples)
filtered_long_merged <- long_merged_junctions |>
  dplyr::filter(samples_over_min_count >= params$min_samples)

# spot-check
filtered_long_combined |>
  head()
```

| chr  | start |   end | ID               | sample     | count | samples_over_min_count |
|:-----|------:|------:|:-----------------|:-----------|------:|-----------------------:|
| chr1 | 14614 | 16858 | chr1:14614-16858 | SRR1559043 |     0 |                     29 |
| chr1 | 14614 | 16858 | chr1:14614-16858 | SRR1559044 |     0 |                     29 |
| chr1 | 14614 | 16858 | chr1:14614-16858 | SRR1559052 |     0 |                     29 |
| chr1 | 14614 | 16858 | chr1:14614-16858 | SRR1559054 |     0 |                     29 |
| chr1 | 14614 | 16858 | chr1:14614-16858 | SRR1559075 |    22 |                     29 |
| chr1 | 14614 | 16858 | chr1:14614-16858 | SRR1559100 |    10 |                     29 |

## Count differences between counts-filtered junction files

Count number of junctions in each filtered file

``` r
length(unique(filtered_long_combined$ID))
```

    [1] 31832

``` r
length(unique(filtered_long_merged$ID))
```

    [1] 28332

Calculate fraction of filtered junctions that are different in any way
in the combined or merged files

``` r
# print number of junctions only in combined junctions bedfile after filtering for junctions with at least 10 counts in at least 5 samples
filtered_combined_only <- dplyr::setdiff(filtered_long_combined, filtered_long_merged)
# print number of unique junctions with any differences between the two dataframes
length(unique(filtered_combined_only$ID))
```

    [1] 4508

``` r
# check fraction of filtered_combined_only junctions
length(unique(filtered_combined_only$ID)) / length(unique(long_combined_junctions$ID))
```

    [1] 0.008148392

``` r
# print number of junctions only in merged junctions bedfile after filtering for junctions with at least 10 counts in at least 5 samples
filtered_merged_only <- dplyr::setdiff(filtered_long_merged, filtered_long_combined)
# print number of unique junctions with any differences between the two dataframes
length(unique(filtered_merged_only$ID))
```

    [1] 1008

``` r
# check fraction of filtered merged_only junctions
length(unique(filtered_merged_only$ID)) / length(unique(long_merged_junctions$ID))
```

    [1] 0.001817316

Check if junctions only found in combined or merged dataframes are
always 1bp in length (are they always exon-intron junctions?)

``` r
filtered_combined_only_ids <- setdiff(filtered_long_combined$ID, filtered_long_merged$ID)
long_combined_only_ids <- filtered_long_combined |>
  dplyr::filter(ID %in% filtered_combined_only_ids) |>
  # count length of junctions only found in combined table
  dplyr::mutate(junction_length = abs(end - start))

# check if the junctions only in combined df are all always 1bp in length (are they all exon-intron junctions?)
unique(long_combined_only_ids$junction_length)
```

    [1] 1

``` r
filtered_merged_only_ids <- setdiff(filtered_long_merged$ID, filtered_long_combined$ID)
long_merged_only_ids <- filtered_long_merged |>
  dplyr::filter(ID %in% filtered_merged_only_ids) |>
  # count length of junctions only found in combined table
  dplyr::mutate(junction_length = abs(end - start))

# check if the junctions only in merged df are all always 1bp in length
unique(long_combined_only_ids$junction_length)
```

    [1] 1

Junctions that are only found in either the combined or merged
dataframes are always 1bp in length, supporting our hypothesis that the
difference stems from exon-intron junctions

### Print examples of junctions only in merged or combined tables

#### Check junctions only in combined junctions file

Print examples of junctions only in the combined junctions file

``` r
filtered_combined_only |>
  dplyr::arrange(start) |>
  head()
```

| chr  | start |   end | ID               | sample     | count | samples_over_min_count |
|:-----|------:|------:|:-----------------|:-----------|------:|-----------------------:|
| chr1 | 18369 | 18370 | chr1:18369-18370 | SRR1559043 |     7 |                      5 |
| chr1 | 18369 | 18370 | chr1:18369-18370 | SRR1559044 |    10 |                      5 |
| chr1 | 18369 | 18370 | chr1:18369-18370 | SRR1559052 |     9 |                      5 |
| chr1 | 18369 | 18370 | chr1:18369-18370 | SRR1559054 |     4 |                      5 |
| chr1 | 18369 | 18370 | chr1:18369-18370 | SRR1559075 |     7 |                      5 |
| chr1 | 18369 | 18370 | chr1:18369-18370 | SRR1559100 |     5 |                      5 |

For comparison, examine merged chromosome 1 counts-filtered junction
table for junctions near chr1:18369-18370

``` r
merged_junctions |>
  dplyr::arrange(start) |>
  dplyr::filter(
    # make sure position filtered is chromosome 1 even if we analyze full bedfile
    chr == "chr1",
    (18369 -  10) < start &
      start < (18369 + 10) )
```

| chr | start | end | ID | SRR1559043 | SRR1559044 | SRR1559052 | SRR1559054 | SRR1559075 | SRR1559100 | SRR1559105 | SRR1559133 | SRR1559134 | SRR1559145 | SRR1559160 | SRR1559164 | SRR1559177 | SRR1559183 | SRR1559184 | SRR1712453 | SRR1712454 | SRR1712455 | SRR1712456 | SRR1712457 | SRR1712458 | SRR1712459 | SRR1712460 | SRR1712461 | SRR1712462 | SRR1712463 | SRR1712464 | SRR1712465 | SRR1784865 | SRR1784867 | SRR1791016 | SRR1791028 | SRR1791108 | SRR1796863 | SRR1796867 | SRR1796893 | SRR1796906 | SRR1796912 | SRR1796939 | SRR1796967 | SRR1796990 | SRR1797014 | SRR1797024 | SRR1797033 | SRR1797034 | SRR1797035 | SRR1797039 | SRR1797052 | SRR1797053 | SRR1797055 | SRR1797057 | SRR1797087 | SRR1797107 | SRR1797111 | SRR1799022 | SRR1799025 | SRR1799041 | SRR1799042 | SRR1799057 | SRR1799058 | SRR1799059 | SRR1799061 | SRR1799062 | SRR1799067 | SRR1799069 | SRR1799081 | SRR1810588 | SRR2042833 | SRR2042845 | SRR2042853 | SRR2042854 | SRR2042856 | SRR2083154 | SRR2083162 | SRR2083171 | SRR2083176 | SRR2083188 | SRR2239703 | SRR2239717 | SRR3162160 | SRR3162195 | SRR3162212 | SRR3162237 | SRR3162253 | SRR4376029 | SRR4416297 | SRR4419554 | SRR4419565 |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 18360 | 188494 | chr1:18360-188494 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18362 | 18497 | chr1:18362-18497 | 22 | 3 | 9 | 9 | 23 | 10 | 20 | 44 | 29 | 15 | 16 | 38 | 3 | 5 | 6 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 16 | 10 | 0 | 6 | 25 | 0 | 0 | 13 | 0 | 25 | 12 | 16 | 14 | 0 | 16 | 17 | 0 | 18 | 14 | 6 | 25 | 5 | 6 | 0 | 0 | 0 | 0 | 8 | 25 | 0 | 15 | 0 | 19 | 5 | 11 | 13 | 0 | 0 | 5 | 0 | 0 | 19 | 0 | 16 | 14 | 7 | 0 | 7 | 0 | 3 | 0 | 2 | 0 | 0 | 0 | 0 | 2 | 6 | 2 |
| chr1 | 18362 | 24738 | chr1:18362-24738 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18362 | 18913 | chr1:18362-18913 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18362 | 188788 | chr1:18362-188788 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18362 | 188791 | chr1:18362-188791 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18362 | 29321 | chr1:18362-29321 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
| chr1 | 18363 | 188137 | chr1:18363-188137 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
| chr1 | 18364 | 188488 | chr1:18364-188488 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 18913 | chr1:18366-18913 | 184 | 95 | 102 | 90 | 317 | 110 | 308 | 459 | 228 | 179 | 305 | 162 | 51 | 26 | 33 | 16 | 5 | 1 | 8 | 8 | 9 | 6 | 4 | 7 | 6 | 3 | 2 | 4 | 23 | 36 | 38 | 19 | 23 | 66 | 32 | 16 | 46 | 20 | 47 | 23 | 10 | 38 | 113 | 30 | 45 | 13 | 42 | 29 | 6 | 53 | 3 | 65 | 18 | 9 | 0 | 0 | 23 | 43 | 11 | 41 | 6 | 11 | 10 | 13 | 12 | 11 | 15 | 14 | 4 | 11 | 31 | 32 | 39 | 60 | 36 | 23 | 17 | 5 | 16 | 7 | 34 | 27 | 19 | 29 | 1 | 15 | 39 | 31 |
| chr1 | 18366 | 24738 | chr1:18366-24738 | 177 | 174 | 152 | 130 | 449 | 341 | 375 | 542 | 354 | 500 | 401 | 479 | 108 | 179 | 325 | 86 | 93 | 53 | 62 | 199 | 153 | 141 | 102 | 175 | 122 | 108 | 50 | 113 | 131 | 117 | 121 | 82 | 66 | 205 | 175 | 70 | 257 | 258 | 329 | 231 | 173 | 154 | 115 | 99 | 109 | 181 | 232 | 123 | 48 | 204 | 68 | 156 | 105 | 51 | 16 | 7 | 188 | 194 | 111 | 194 | 99 | 91 | 127 | 187 | 79 | 111 | 172 | 20 | 3 | 46 | 139 | 74 | 138 | 230 | 40 | 228 | 105 | 15 | 19 | 44 | 47 | 59 | 66 | 187 | 13 | 39 | 81 | 43 |
| chr1 | 18366 | 29321 | chr1:18366-29321 | 3 | 11 | 0 | 4 | 24 | 6 | 14 | 12 | 24 | 10 | 8 | 9 | 0 | 3 | 12 | 0 | 0 | 0 | 0 | 13 | 0 | 3 | 0 | 0 | 0 | 0 | 3 | 0 | 11 | 0 | 5 | 6 | 0 | 8 | 11 | 9 | 6 | 24 | 10 | 10 | 5 | 8 | 0 | 4 | 12 | 8 | 11 | 34 | 0 | 8 | 2 | 8 | 2 | 3 | 2 | 1 | 10 | 2 | 4 | 0 | 15 | 8 | 0 | 26 | 0 | 21 | 22 | 2 | 0 | 4 | 17 | 11 | 29 | 56 | 17 | 33 | 32 | 0 | 5 | 0 | 0 | 1 | 3 | 2 | 2 | 0 | 0 | 0 |
| chr1 | 18366 | 29534 | chr1:18366-29534 | 0 | 0 | 0 | 0 | 1 | 0 | 3 | 2 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 2 | 1 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 29824 | chr1:18366-29824 | 0 | 0 | 0 | 0 | 7 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 23731 | chr1:18366-23731 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 29603 | chr1:18366-29603 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 8 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 18367 | chr1:18366-18367 | 0 | 0 | 0 | 0 | 27 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 267303 | chr1:18366-267303 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 356685 | chr1:18366-356685 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 90726 | chr1:18366-90726 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 188869 | chr1:18366-188869 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 188791 | chr1:18366-188791 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18366 | 18444 | chr1:18366-18444 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18367 | 188793 | chr1:18367-188793 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18367 | 188805 | chr1:18367-188805 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18369 | 18501 | chr1:18369-18501 | 32 | 10 | 2 | 19 | 36 | 18 | 52 | 42 | 32 | 15 | 18 | 46 | 9 | 9 | 8 | 0 | 0 | 1 | 2 | 0 | 1 | 0 | 1 | 2 | 0 | 0 | 2 | 0 | 2 | 4 | 3 | 2 | 4 | 4 | 0 | 3 | 7 | 3 | 10 | 7 | 7 | 3 | 3 | 1 | 6 | 1 | 14 | 6 | 0 | 6 | 2 | 4 | 0 | 0 | 0 | 0 | 2 | 0 | 1 | 11 | 4 | 2 | 3 | 4 | 0 | 3 | 3 | 1 | 0 | 1 | 5 | 7 | 4 | 10 | 1 | 5 | 0 | 1 | 3 | 0 | 0 | 2 | 0 | 3 | 0 | 1 | 3 | 7 |
| chr1 | 18369 | 18913 | chr1:18369-18913 | 40 | 50 | 77 | 43 | 163 | 56 | 174 | 256 | 117 | 76 | 191 | 94 | 24 | 18 | 26 | 0 | 2 | 0 | 6 | 4 | 8 | 4 | 4 | 0 | 10 | 12 | 0 | 0 | 9 | 6 | 15 | 8 | 9 | 36 | 9 | 19 | 19 | 5 | 18 | 8 | 3 | 18 | 45 | 16 | 25 | 8 | 15 | 7 | 4 | 33 | 2 | 19 | 5 | 3 | 0 | 3 | 0 | 17 | 5 | 10 | 0 | 8 | 4 | 4 | 2 | 4 | 0 | 2 | 1 | 2 | 13 | 11 | 12 | 18 | 39 | 9 | 0 | 3 | 5 | 9 | 15 | 5 | 11 | 16 | 3 | 4 | 12 | 17 |
| chr1 | 18369 | 24738 | chr1:18369-24738 | 30 | 20 | 8 | 9 | 54 | 26 | 52 | 83 | 32 | 62 | 44 | 71 | 11 | 15 | 31 | 15 | 19 | 3 | 12 | 20 | 11 | 9 | 6 | 38 | 17 | 17 | 0 | 13 | 13 | 11 | 7 | 8 | 7 | 28 | 17 | 9 | 6 | 7 | 29 | 13 | 28 | 5 | 5 | 11 | 13 | 12 | 35 | 10 | 12 | 21 | 2 | 0 | 6 | 3 | 0 | 0 | 13 | 46 | 0 | 12 | 6 | 4 | 4 | 18 | 9 | 7 | 9 | 1 | 1 | 2 | 9 | 8 | 7 | 23 | 17 | 23 | 10 | 2 | 1 | 2 | 8 | 3 | 0 | 65 | 0 | 0 | 11 | 9 |
| chr1 | 18369 | 29321 | chr1:18369-29321 | 0 | 1 | 4 | 3 | 15 | 8 | 20 | 12 | 14 | 15 | 0 | 5 | 0 | 4 | 9 | 0 | 0 | 1 | 0 | 7 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 4 | 0 | 0 | 1 | 0 | 0 | 6 | 3 | 7 | 4 | 7 | 2 | 0 | 7 | 0 | 5 | 6 | 5 | 0 | 18 | 0 | 7 | 1 | 2 | 1 | 0 | 0 | 0 | 4 | 2 | 2 | 6 | 3 | 0 | 0 | 10 | 0 | 0 | 3 | 1 | 0 | 3 | 6 | 3 | 14 | 11 | 6 | 13 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 9 | 1 | 0 | 0 | 3 |
| chr1 | 18369 | 18497 | chr1:18369-18497 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18369 | 189796 | chr1:18369-189796 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18369 | 19273 | chr1:18369-19273 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18369 | 18370 | chr1:18369-18370 | 0 | 0 | 0 | 0 | 7 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18369 | 194838 | chr1:18369-194838 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18369 | 24313 | chr1:18369-24313 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18369 | 29294 | chr1:18369-29294 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18369 | 29540 | chr1:18369-29540 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| chr1 | 18377 | 188198 | chr1:18377-188198 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 |

This junction is completely missed, with the closest coordinate being
chr1:18369-18913.

Check if chr1:18369-18913 is only in the merged table and not in the
combined table

``` r
"chr1:18369-18913" %in% merged_only_ids
```

    [1] FALSE

chr1:18369-18913 is found in both the merged and combined junctions
tables.

#### Check junctions only in merged junctions file

Print examples of junctions only in the merged junctions file

``` r
filtered_merged_only |>
  dplyr::arrange(start) |>
  head()
```

| chr  |  start |    end | ID                 | sample     | count | samples_over_min_count |
|:-----|-------:|-------:|:-------------------|:-----------|------:|-----------------------:|
| chr1 | 268020 | 268021 | chr1:268020-268021 | SRR1559043 |    40 |                     14 |
| chr1 | 268020 | 268021 | chr1:268020-268021 | SRR1559044 |     4 |                     14 |
| chr1 | 268020 | 268021 | chr1:268020-268021 | SRR1559052 |    18 |                     14 |
| chr1 | 268020 | 268021 | chr1:268020-268021 | SRR1559054 |    31 |                     14 |
| chr1 | 268020 | 268021 | chr1:268020-268021 | SRR1559075 |     5 |                     14 |
| chr1 | 268020 | 268021 | chr1:268020-268021 | SRR1559100 |    22 |                     14 |

For comparison, examine combined junction table for junctions near
chr1:268020-268021

``` r
deduplicated_combined_junctions |>
  dplyr::arrange(start) |>
  dplyr::filter(
    chr == "chr1",
    (268020 -  10) < start &
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
are \< 10.

Because the values in the merged table are 0 for the counts that differ
from the combined table, I suspect these differences stem from junctions
that are not shared across samples (these appear as NA values which we
replace with 0 in the final merged bedfile). In the merged method, these
0 counts arise because exon-intron boundaries are defined from an intron
events coordinates file from gtf2event.py that was constructed from the
annotations of one sample only.

## Check if differences flagged in junction IDs that are identical stem from differences in 0 counts

``` r
# Exclude IDs that are only found in each dataframe
# The "_shared" dataframe contain rows corresponding to junctions that are shared between the merged and combined methods
# The differences within these dataframes should come from the values of the samples' junction counts for these IDs

merged_junctions_shared <- merged_junctions |>
  dplyr::filter(!ID %in% merged_only_ids)

combined_junctions_shared <- deduplicated_combined_junctions |>
  dplyr::filter(!ID %in% combined_only_ids)

# check that resulting length of dataframes excluding differing IDs are the same (they should be)
nrow(merged_junctions_shared) == nrow(combined_junctions_shared)
```

    [1] TRUE

``` r
# pull out dataframe of junction IDs with values in merged_junctions_shared not in combined_junctions_shared
merged_shared_id_differences <- dplyr::setdiff(merged_junctions_shared, combined_junctions_shared)

combined_shared_id_differences <- combined_junctions_shared |> dplyr::filter(
  ID %in% merged_shared_id_differences$ID
)
```

Examine junction counts differences in one sample only. Check if there
are any junctions that are different because the junction counts of the
merged and combined dataframes are nonzero and different

``` r
# Pull one sample ID to examine
# This should be in params but I think ideally I would want to do this for all samples
# Not sure of the best way to do that at the moment - maybe do this on the long dataframe?

sample <- "SRR1559043"

# select for only this sample in the combined and merged dataframes
merged_shared_id_differences <- merged_shared_id_differences |>
  dplyr::select(chr, start, end, ID, dplyr::all_of(sample))
combined_shared_id_differences <- combined_shared_id_differences |>
  dplyr::select(chr, start, end, ID, dplyr::all_of(sample))

# Combine the dataframes and label each row by whether the junction count for the example sample are zero-to-nonzero, the same, or different due to nonzero-to-nonzero reasons

compare_junction_values_df <- dplyr::left_join(
  merged_shared_id_differences,
  combined_shared_id_differences,
  by = c("chr", "start", "end", "ID"),
  suffix = c("_merged", "_combined")
) |>
  # label junctions by nature of difference
  dplyr::mutate(
    difference_type = dplyr::case_when(SRR1559043_merged == SRR1559043_combined ~ "same",
                                       SRR1559043_merged != 0 & SRR1559043_combined != 0 ~ "nonzero-to-nonzero",
                                       SRR1559043_merged == 0 & SRR1559043_combined != 0 ~ "zero-to-nonzero"
    )
  )

# spot-check
compare_junction_values_df |> head()
```

| chr | start | end | ID | SRR1559043_merged | SRR1559043_combined | difference_type |
|:---|---:|---:|:---|---:|---:|:---|
| chr1 | 100974902 | 100974903 | chr1:100974902-100974903 | 0 | 0 | same |
| chr1 | 100975218 | 100975219 | chr1:100975218-100975219 | 1 | 1 | same |
| chr1 | 101025308 | 101025309 | chr1:101025308-101025309 | 13 | 13 | same |
| chr1 | 101025466 | 101025467 | chr1:101025466-101025467 | 7 | 7 | same |
| chr1 | 101025639 | 101025640 | chr1:101025639-101025640 | 0 | 0 | same |
| chr1 | 101025664 | 101025665 | chr1:101025664-101025665 | 0 | 0 | same |

Summarize fraction of junctions in sample SRR1559043 in each category

``` r
compare_junction_values_df |>
  dplyr::summarise(
    total = dplyr::n(),
    n_same = sum(difference_type == "same"),
    n_zero_to_nonzero = sum(difference_type == "zero-to-nonzero"),
    n_nonzero_to_nonzero = sum(difference_type == "nonzero-to-nonzeo"),
    frac_same = n_same / total,
    frac_zero_to_nonzero = n_zero_to_nonzero / total,
    frac_nonzero_to_nonzero = n_nonzero_to_nonzero / total
  )
```

| total | n_same | n_zero_to_nonzero | n_nonzero_to_nonzero | frac_same | frac_zero_to_nonzero | frac_nonzero_to_nonzero |
|---:|---:|---:|---:|---:|---:|---:|
| 4199 | 2038 | 2161 | 0 | 0.4853537 | 0.5146463 | 0 |
