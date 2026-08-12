# Compare merged PSI table from TARGET pilot samples
Cindy Liang (celiang@ucsc.edu)
2026-08-11

## Introduction

Compendium workflow v1 (workflow with updated GTF generation and merging
rules) has been run on the TARGET pilot samples. To understand what
differences remain between PSI results generated from this workflow and
the canonical Shiba run, a portion of the `compared_merged_psi_table`
notebook is run on the PSI tables generated from the v1 workflow and the
unaltered Shiba scripts.

**v1**: Refers to PSI results generated from the most recent compendium
workflow method

**Combined**: Refers to PSI results generated from the unaltered Shiba
script.

**Combined_dedup**: PSI valiues generated from running only the Shiba
PSI.py script on the following inputs and command:

- event coordinates files generated from the Combined method

- Junctions.bed files generated from Combined, then duplicated with the
  deduplication commands used in the v1 workflow

<!-- -->

    python $CONDA_PREFIX/share/shiba-0.8.1-0/src/psi.py -m \
    10 -p \
    30 -v \
    --onlypsi \
    /private/groups/treehouse/working-projects/celiang/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/deduplicated_junctions.bed \
    /private/groups/treehouse/working-projects/celiang/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/events \
    v1_dedup_target_pilot_combined

**Question:** Does merging v1 splice tables cause us to lose out on the
trustworthiness of unannotated events to an extent that it impacts many
samples?

As part of our pipeline, we plan to merge Shiba tables created for
individual samples to obtain a final PSI table of all samples. Before
doing so, we use this notebook to test that merging PSI tables from v1
runs will not introduce a large amount of untrustworthy splice events.

As part of tracking differences between the two tables, we have labeled
NA values introduced at various steps of the data processing pipeline
with different negative values:

- -1: NA values assigned by Shiba (splice events with low coverage)

We are particularly interested in the impact merging v1 tables has on
unannotated event detection, since this is one novelty of using Shiba.

**Targeted questions this notebook is trying to answer:**

- What splice event types do we miss out on with the v1 tables method?
- Of splice events with high numbers of NA values, how many of these
  events only have a numeric PSI value in one (or a very low amount of)
  samples?
- How many NA values of different NA types are there in each type of
  splice table? How are they distributed across events? I am more
  concerned with the merge NAs (-2), as there is no good way to replace
  them in the v1 tables (they could either be 0 or 1). Another concern
  is events that are only found in the combined table, or only found in
  the v1 table, as there is no easy way to correct for those.
- Can the number of NA values be reduced in the v1 tables method if we
  filter for events with numeric PSI values in a minimum number of
  samples?
- What is the relationship between PSI values and NAs in each splice
  table type?
- What is the relationship between PSI values and NAs in each splice
  table type after filtering for events with numeric PSI values in a
  minimum number of samples? Assuming the Shiba NAs (-1) are from low
  read support, do the v1 and combined tables always call the same
  events as a Shiba NA? Merge NAs (-2) should arise from events dropped
  by Shiba due to lack of alternative transcripts for a gene, so we
  expect these NAs to be 0 or 1.
- How well correlated are the number of NA values in each event type for
  each PSI table type?
- How well correlated are the PSI values of events that are shared?

## Set up

## Directories and files

``` r
## directories ##
# find the root-level repo directory so we can access the other files
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
# define the data directories
exploration_dir <- file.path(repo_root, "exploration")
merge_exploration_dir <- file.path(exploration_dir, "merged-method-target-pilot-eval")
# these contain zipped psi files
target_pilot_workflow_results_dir <- file.path(exploration_dir, "v1_workflow_target_pilot_results", "psi")

#output dir for long tables
target_pilot_output <- file.path(exploration_dir, "v1_psi_tables")

## combined / canonical shiba run results
target_pilot_dir <- file.path(exploration_dir, "merging-psi-tables", "target_pilot")
# directory of shiba results produced from the same shiba run
combined_dir <- file.path(target_pilot_dir, "shiba_combined", "results")
# psi table directory
# these contain unzipped psi files
combined_splice_results_dir <- file.path(combined_dir, "splicing")

# psi table dir of deduplicated combined pilot results
# these are unzipped files
dedup_combined_results_dir <- file.path(exploration_dir, "dedup_target_pilot_combined_results")

# Check output dir exists; if not, create it
if (!dir.exists(target_pilot_output)) {
  dir.create(target_pilot_output)
}

## files ##
# define list of PSI results files corresponding to event types quantified by Shiba bulk analysis
v1_psi_files <- c(
  se = "PSI_SE.txt.gz",
  afe = "PSI_AFE.txt.gz",
  ale = "PSI_ALE.txt.gz",
  five = "PSI_FIVE.txt.gz",
  three = "PSI_THREE.txt.gz",
  mse = "PSI_MSE.txt.gz",
  mxe = "PSI_MXE.txt.gz",
  ri = "PSI_RI.txt.gz"
)

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
shiba_psi_paths <-file.path(combined_splice_results_dir, psi_files)
names(shiba_psi_paths) <- names(psi_files)

dedup_shiba_psi_paths <-file.path(dedup_combined_results_dir, psi_files)
names(dedup_shiba_psi_paths) <- names(psi_files)

v1_psi_paths <-file.path(target_pilot_workflow_results_dir, v1_psi_files)
names(v1_psi_paths) <- names(v1_psi_files)

## output ##
# long df output file of combined and v1 PSI tables
long_df_output <- file.path(target_pilot_output, "long_psi_df.rds")
```

Read in files

``` r
# read merged splice table from v1 workflow on the TARGET pilot samples
v1_splice_results <- v1_psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
      # convert PSI values to numeric and replace NA values from shiba representing low coverage with -1
      dplyr::mutate(across(contains("_PSI"), \(x) tidyr::replace_na(as.numeric(x), -1)))
  })

# combine PSI values of samples run together into one dataframe to compare against v1 PSI dataframe
v1_psi_table <- purrr::list_rbind(v1_splice_results, names_to = "event_type")

# read in combined splice results to compare against v1 results
combined_splice_results <- shiba_psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
      # convert PSI values to numeric and replace NA values from shiba representing low coverage with -1
      dplyr::mutate(across(contains("_PSI"), \(x) tidyr::replace_na(as.numeric(x), -1)))
  })

# combine PSI values of samples run together into one dataframe to compare against v1 PSI dataframe
combined_psi_table <- purrr::list_rbind(combined_splice_results, names_to = "event_type")

# read in deduplicated splice results
dedup_splice_results <- dedup_shiba_psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
      # convert PSI values to numeric and replace NA values from shiba representing low coverage with -1
      dplyr::mutate(across(contains("_PSI"), \(x) tidyr::replace_na(as.numeric(x), -1)))
  })

# combine PSI values of samples run together into one dataframe to compare against v1 PSI dataframe
dedup_combined_psi_table <- purrr::list_rbind(dedup_splice_results, names_to = "event_type")
```

### Obtain dimensions of each PSI table

``` r
# dimensions of v1 PSI table
dim(v1_psi_table)
```

    [1] 746294     92

``` r
# dimensions of combined PSI table
dim(combined_psi_table)
```

    [1] 746601     92

``` r
# dimensions of deduplicated combined PSI table
dim(dedup_combined_psi_table)
```

    [1] 746601     92

There are ~400 more splice events in the combined PSI table than the v1
table.

``` r
identical(v1_psi_table, dedup_combined_psi_table)
```

    [1] FALSE

``` r
identical(combined_psi_table, dedup_combined_psi_table)
```

    [1] TRUE

The combined PSI tables and PSI tables from running psi.py on the
combined run’s event files and deduplicated junctions are identical This
kind of makes sense to me because the final PSI tables don’t have
duplicated junction pos_ids; so I think the bug didn’t impact the PSI
calculation since shiba just “picked” the valid junction rows?

## What splice event types do we miss out on with the v1 tables method?

To help us prioritize what splice event types we may be the most
confident in after merging v1 Shiba PSI tables, we are interested in
seeing a breakdown of what event types are most represented in PSI
tables constructed from different methods

``` r
# create df of splice events that are shared between both combined and v1 tables
# This is an event-level dataframe (each row is one unique pos_id/event_id)
# columns are each samples' PSI value with -1 indicating low read NA; NA indicating event not shared between methods
# combined_event, v1_event, and shared_event columns label what tables the event is present in
all_events <- dplyr::full_join(
  dedup_combined_psi_table,
  v1_psi_table,
  by = c("event_type", "pos_id", "gene_id", "label"),
  # label PSI values by table they came from
  suffix = c("_combined", "_v1")) |>
  # categorize events by whether they are in the combined or v1 tables
  dplyr::mutate(
    # combined events are counted if the values are not all NA
    combined_event = ! dplyr::if_all(ends_with("_combined"), is.na),
    # v1 events are counted if the values for the event are not all NA (indicating the specific position ID is only found in the combined or v1 table)
    v1_event = ! dplyr::if_all(ends_with("_v1"), is.na),
    shared_event = combined_event & v1_event
  )
```

Print summary of counts across each method, across all event types

``` r
# number of total events shared or present only in the combined and v1 tables
# this table counts events with -1 and -2 PSI values

all_methods_summary <- all_events |> dplyr::summarise(
    # count number of events in each event type and annotation category
    shared_count = sum(shared_event),
    combined_only_count = sum(combined_event) - shared_count,
    v1_only_count = sum(v1_event) - shared_count
)

all_methods_summary
```

| shared_count | combined_only_count | v1_only_count |
|-------------:|--------------------:|--------------:|
|       744659 |                1942 |          1635 |

Although the majority of splice events (inclusive of NA values) are
shared between methods, some events remain that are only in the combined
or v1 method

Print summary of counts and percentages across each method, broken down
by event type

``` r
# event-level summary of the number of each splice event type in each splice table
all_events_summary <- all_events |> dplyr::summarise(
    .by = c(event_type, label),
    # count number of events in each event type and annotation category
    total = dplyr::n(),
    shared_count = sum(shared_event),
    combined_only_count = sum(combined_event) - shared_count,
    v1_only_count = sum(v1_event) - shared_count,
    shared_percent = shared_count / total * 100,
    combined_only_percent = combined_only_count / total * 100,
    v1_only_percent = v1_only_count / total * 100,
)

all_events_summary |> dplyr::select(
  event_type,
  label,
  total,
  shared_percent,
  combined_only_percent,
  v1_only_percent)
```

| event_type | label | total | shared_percent | combined_only_percent | v1_only_percent |
|:---|:---|---:|---:|---:|---:|
| se | annotated | 103591 | 98.70838 | 0.6458090 | 0.6458090 |
| se | unannotated | 28248 | 100.00000 | 0.0000000 | 0.0000000 |
| afe | unannotated | 64044 | 100.00000 | 0.0000000 | 0.0000000 |
| afe | annotated | 158700 | 99.94077 | 0.0296156 | 0.0296156 |
| ale | annotated | 130130 | 99.97541 | 0.0122954 | 0.0122954 |
| ale | unannotated | 35648 | 100.00000 | 0.0000000 | 0.0000000 |
| five | annotated | 32471 | 99.12537 | 0.4373133 | 0.4373133 |
| five | unannotated | 13488 | 100.00000 | 0.0000000 | 0.0000000 |
| three | annotated | 38118 | 98.64106 | 0.6794690 | 0.6794690 |
| three | unannotated | 14086 | 100.00000 | 0.0000000 | 0.0000000 |
| mse | annotated | 62589 | 99.62613 | 0.1869338 | 0.1869338 |
| mse | unannotated | 16363 | 100.00000 | 0.0000000 | 0.0000000 |
| mxe | annotated | 801 | 99.50062 | 0.2496879 | 0.2496879 |
| mxe | unannotated | 226 | 100.00000 | 0.0000000 | 0.0000000 |
| ri | unannotated | 36314 | 98.18803 | 1.3438343 | 0.4681390 |
| ri | annotated | 13419 | 96.90737 | 1.5053283 | 1.5873016 |

Check what the unshared events look like

``` r
unshared_events <- all_events |> dplyr::filter(
  shared_event == FALSE,
  event_type == "se"
) |>
  dplyr::arrange(pos_id)

head(unshared_events)
```

| event_type | pos_id | gene_id | label | SRR1559043_PSI_combined | SRR1559044_PSI_combined | SRR1559052_PSI_combined | SRR1559054_PSI_combined | SRR1559075_PSI_combined | SRR1559100_PSI_combined | SRR1559105_PSI_combined | SRR1559133_PSI_combined | SRR1559134_PSI_combined | SRR1559145_PSI_combined | SRR1559160_PSI_combined | SRR1559164_PSI_combined | SRR1559177_PSI_combined | SRR1559183_PSI_combined | SRR1559184_PSI_combined | SRR1712453_PSI_combined | SRR1712454_PSI_combined | SRR1712455_PSI_combined | SRR1712456_PSI_combined | SRR1712457_PSI_combined | SRR1712458_PSI_combined | SRR1712459_PSI_combined | SRR1712460_PSI_combined | SRR1712461_PSI_combined | SRR1712462_PSI_combined | SRR1712463_PSI_combined | SRR1712464_PSI_combined | SRR1712465_PSI_combined | SRR1784865_PSI_combined | SRR1784867_PSI_combined | SRR1791016_PSI_combined | SRR1791028_PSI_combined | SRR1791108_PSI_combined | SRR1796863_PSI_combined | SRR1796867_PSI_combined | SRR1796893_PSI_combined | SRR1796906_PSI_combined | SRR1796912_PSI_combined | SRR1796939_PSI_combined | SRR1796967_PSI_combined | SRR1796990_PSI_combined | SRR1797014_PSI_combined | SRR1797024_PSI_combined | SRR1797033_PSI_combined | SRR1797034_PSI_combined | SRR1797035_PSI_combined | SRR1797039_PSI_combined | SRR1797052_PSI_combined | SRR1797053_PSI_combined | SRR1797055_PSI_combined | SRR1797057_PSI_combined | SRR1797087_PSI_combined | SRR1797107_PSI_combined | SRR1797111_PSI_combined | SRR1799022_PSI_combined | SRR1799025_PSI_combined | SRR1799041_PSI_combined | SRR1799042_PSI_combined | SRR1799057_PSI_combined | SRR1799058_PSI_combined | SRR1799059_PSI_combined | SRR1799061_PSI_combined | SRR1799062_PSI_combined | SRR1799067_PSI_combined | SRR1799069_PSI_combined | SRR1799081_PSI_combined | SRR1810588_PSI_combined | SRR2042833_PSI_combined | SRR2042845_PSI_combined | SRR2042853_PSI_combined | SRR2042854_PSI_combined | SRR2042856_PSI_combined | SRR2083154_PSI_combined | SRR2083162_PSI_combined | SRR2083171_PSI_combined | SRR2083176_PSI_combined | SRR2083188_PSI_combined | SRR2239703_PSI_combined | SRR2239717_PSI_combined | SRR3162160_PSI_combined | SRR3162195_PSI_combined | SRR3162212_PSI_combined | SRR3162237_PSI_combined | SRR3162253_PSI_combined | SRR4376029_PSI_combined | SRR4416297_PSI_combined | SRR4419554_PSI_combined | SRR4419565_PSI_combined | SRR1559043_PSI_v1 | SRR1559044_PSI_v1 | SRR1559052_PSI_v1 | SRR1559054_PSI_v1 | SRR1559075_PSI_v1 | SRR1559100_PSI_v1 | SRR1559105_PSI_v1 | SRR1559133_PSI_v1 | SRR1559134_PSI_v1 | SRR1559145_PSI_v1 | SRR1559160_PSI_v1 | SRR1559164_PSI_v1 | SRR1559177_PSI_v1 | SRR1559183_PSI_v1 | SRR1559184_PSI_v1 | SRR1712453_PSI_v1 | SRR1712454_PSI_v1 | SRR1712455_PSI_v1 | SRR1712456_PSI_v1 | SRR1712457_PSI_v1 | SRR1712458_PSI_v1 | SRR1712459_PSI_v1 | SRR1712460_PSI_v1 | SRR1712461_PSI_v1 | SRR1712462_PSI_v1 | SRR1712463_PSI_v1 | SRR1712464_PSI_v1 | SRR1712465_PSI_v1 | SRR1784865_PSI_v1 | SRR1784867_PSI_v1 | SRR1791016_PSI_v1 | SRR1791028_PSI_v1 | SRR1791108_PSI_v1 | SRR1796863_PSI_v1 | SRR1796867_PSI_v1 | SRR1796893_PSI_v1 | SRR1796906_PSI_v1 | SRR1796912_PSI_v1 | SRR1796939_PSI_v1 | SRR1796967_PSI_v1 | SRR1796990_PSI_v1 | SRR1797014_PSI_v1 | SRR1797024_PSI_v1 | SRR1797033_PSI_v1 | SRR1797034_PSI_v1 | SRR1797035_PSI_v1 | SRR1797039_PSI_v1 | SRR1797052_PSI_v1 | SRR1797053_PSI_v1 | SRR1797055_PSI_v1 | SRR1797057_PSI_v1 | SRR1797087_PSI_v1 | SRR1797107_PSI_v1 | SRR1797111_PSI_v1 | SRR1799022_PSI_v1 | SRR1799025_PSI_v1 | SRR1799041_PSI_v1 | SRR1799042_PSI_v1 | SRR1799057_PSI_v1 | SRR1799058_PSI_v1 | SRR1799059_PSI_v1 | SRR1799061_PSI_v1 | SRR1799062_PSI_v1 | SRR1799067_PSI_v1 | SRR1799069_PSI_v1 | SRR1799081_PSI_v1 | SRR1810588_PSI_v1 | SRR2042833_PSI_v1 | SRR2042845_PSI_v1 | SRR2042853_PSI_v1 | SRR2042854_PSI_v1 | SRR2042856_PSI_v1 | SRR2083154_PSI_v1 | SRR2083162_PSI_v1 | SRR2083171_PSI_v1 | SRR2083176_PSI_v1 | SRR2083188_PSI_v1 | SRR2239703_PSI_v1 | SRR2239717_PSI_v1 | SRR3162160_PSI_v1 | SRR3162195_PSI_v1 | SRR3162212_PSI_v1 | SRR3162237_PSI_v1 | SRR3162253_PSI_v1 | SRR4376029_PSI_v1 | SRR4416297_PSI_v1 | SRR4419554_PSI_v1 | SRR4419565_PSI_v1 | combined_event | v1_event | shared_event |
|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|:---|
| se | SE@chr10@100523977-100524139@100516960-100526399 | ENSG00000255339.6 | annotated | 0.0000000 | -1.0000000 | 0.6666667 | -1.0000000 | 0.2500000 | 0.2000000 | 0.1111111 | 0.2307692 | 0.2444444 | 0.1250000 | 0.1940299 | 0.3225806 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | 0.0000000 | -1.0000000 | -1.000000 | -1.0000000 | -1.0000000 | -1.00000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | 0.0000000 | 0.000000 | 0.0149254 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.000000 | -1.0000000 | -1.0000000 | 0.0000000 | -1.0000000 | 0.0000000 | -1.000000 | -1.0000000 | 0.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1 | -1.000000 | -1.0000000 | -1.0000000 | 0.0338983 | 0.0361446 | 0.0000000 | 0.0142180 | 0.2121212 | -1.00 | -1.0 | 0.0285714 | 0.0000000 | 0.0250000 | 0.0000000 | 0.0370370 | 0.0000000 | 0.0078740 | 0.0083102 | 0.0000000 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | TRUE | FALSE | FALSE |
| se | SE@chr10@100523977-100524139@100516960-100526399 | ENSG00000075826.17 | annotated | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 0.0000000 | -1.0000000 | 0.6666667 | -1.0000000 | 0.2500000 | 0.2000000 | 0.1111111 | 0.2307692 | 0.2444444 | 0.1250000 | 0.1940299 | 0.3225806 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | 0.0000000 | -1.0000000 | -1.000000 | -1.0000000 | -1.0000000 | -1.00000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | 0.0000000 | 0.000000 | 0.0149254 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.000000 | -1.0000000 | -1.0000000 | 0.0000000 | -1.0000000 | 0.0000000 | -1.000000 | -1.0000000 | 0.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1 | -1.000000 | -1.0000000 | -1.0000000 | 0.0338983 | 0.0361446 | 0.0000000 | 0.0142180 | 0.2121212 | -1.00 | -1.0 | 0.0285714 | 0.0000000 | 0.0250000 | 0.0000000 | 0.0370370 | 0.0000000 | 0.0078740 | 0.0083102 | 0.0000000 | FALSE | TRUE | FALSE |
| se | SE@chr10@100526399-100526554@100516960-100526975 | ENSG00000255339.6 | annotated | 0.9921824 | 0.9851190 | 0.9947951 | 0.9909228 | 0.9989293 | 0.9970164 | 0.9984733 | 0.9975699 | 0.9993523 | 0.9918145 | 0.9969574 | 1.0000000 | 0.9880668 | 1.0000000 | 1.0000000 | 0.9868204 | 0.9958848 | 0.9910714 | 1.0000000 | 0.9957537 | 1.0000000 | 1.0000000 | 1.000000 | 1.0000000 | 0.9830867 | 1.00000 | 0.9927140 | 0.9915789 | 0.9974969 | 1.0000000 | 0.9870968 | 1.000000 | 0.9886364 | 0.9956757 | 0.9989610 | 1.0000000 | 0.9981735 | 0.9930556 | 1.0000000 | 0.9980658 | 1.0000000 | 0.9963931 | 1.0000000 | 0.9989270 | 0.9972106 | 1.0000000 | 1.000000 | 1.0000000 | 1.0000000 | 0.9952996 | 1.0000000 | 0.9936609 | 1.000000 | 1.0000000 | 1.0000000 | 0.9964476 | 1.0000000 | 1.0000000 | 1.0000000 | 0.9978371 | 1.0000000 | 1.000000 | 0.9993275 | 0.9992221 | 0.9971469 | 1.0000000 | 1.0000000 | 0.9980159 | 1 | 1.000000 | 0.9971510 | 0.9990352 | 1.0000000 | 0.9977439 | 0.9968701 | 0.9868938 | 1.0000000 | 1.00 | 1.0 | 0.9856115 | 0.9840000 | 0.9921260 | 0.9936306 | 0.9955056 | 0.9918699 | 0.9512195 | 0.9705882 | 0.9920000 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | TRUE | FALSE | FALSE |
| se | SE@chr10@100526399-100526554@100516960-100526975 | ENSG00000075826.17 | annotated | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 0.9921824 | 0.9851190 | 0.9947951 | 0.9909228 | 0.9989293 | 0.9970164 | 0.9984733 | 0.9975699 | 0.9993523 | 0.9918145 | 0.9969574 | 1.0000000 | 0.9880668 | 1.0000000 | 1.0000000 | 0.9868204 | 0.9958848 | 0.9910714 | 1.0000000 | 0.9957537 | 1.0000000 | 1.0000000 | 1.000000 | 1.0000000 | 0.9830867 | 1.00000 | 0.9927140 | 0.9915789 | 0.9974969 | 1.0000000 | 0.9870968 | 1.000000 | 0.9886364 | 0.9956757 | 0.9989610 | 1.0000000 | 0.9981735 | 0.9930556 | 1.0000000 | 0.9980658 | 1.0000000 | 0.9963931 | 1.0000000 | 0.9989270 | 0.9972106 | 1.0000000 | 1.000000 | 1.0000000 | 1.0000000 | 0.9952996 | 1.0000000 | 0.9936609 | 1.000000 | 1.0000000 | 1.0000000 | 0.9964476 | 1.0000000 | 1.0000000 | 1.0000000 | 0.9978371 | 1.0000000 | 1.000000 | 0.9993275 | 0.9992221 | 0.9971469 | 1.0000000 | 1.0000000 | 0.9980159 | 1 | 1.000000 | 0.9971510 | 0.9990352 | 1.0000000 | 0.9977439 | 0.9968701 | 0.9868938 | 1.0000000 | 1.00 | 1.0 | 0.9856115 | 0.9840000 | 0.9921260 | 0.9936306 | 0.9955056 | 0.9918699 | 0.9512195 | 0.9705882 | 0.9920000 | FALSE | TRUE | FALSE |
| se | SE@chr10@101793913-101793998@101792943-101797980 | ENSG00000120049.20 | annotated | 0.8383838 | 0.8217822 | 0.9000000 | 0.6796715 | 0.8206278 | 0.8333333 | 0.7910864 | 0.7319588 | 0.8063439 | 0.7768987 | 0.6767896 | 0.9116719 | 0.7859532 | 0.9033233 | 0.9201774 | 0.8095238 | 0.7843943 | 0.6490728 | 0.7413074 | 0.8010657 | 0.8200837 | 0.7783595 | 0.852618 | 0.8150209 | 0.8254675 | 0.71341 | 0.8476923 | 0.7935349 | 0.6774194 | 0.8141593 | 0.9144385 | 0.836478 | 0.8200000 | 0.8677686 | 0.9371069 | 0.8987342 | 0.9099099 | 0.8295964 | 0.8181818 | 0.8633094 | 0.8470948 | 0.8164852 | 0.5949367 | 0.8942308 | 0.8401559 | 0.8680203 | 0.886202 | 0.9273828 | 0.6842105 | 0.7121951 | 0.7916667 | 0.9731286 | 0.974359 | 0.4049587 | 0.7714286 | 1.0000000 | 0.8918919 | 0.7640449 | 0.8444444 | 0.8785714 | 0.9904762 | 0.967033 | 0.9828080 | 0.9221790 | 0.9854015 | 0.9537572 | 0.7241379 | 0.9277108 | -1 | 0.981982 | 0.9553903 | 0.9126214 | 0.9276018 | 0.9469027 | 0.8917647 | 0.9467681 | 0.8818898 | 0.95 | 0.9 | 0.7538462 | 0.8666667 | 0.8111111 | 0.6608187 | 0.4479167 | 0.8319328 | 0.9047619 | 0.9605263 | 0.7878788 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | TRUE | FALSE | FALSE |
| se | SE@chr10@101793913-101793998@101792943-101797980 | ENSG00000198408.14 | annotated | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 0.8383838 | 0.8217822 | 0.9000000 | 0.6796715 | 0.8206278 | 0.8333333 | 0.7910864 | 0.7319588 | 0.8063439 | 0.7768987 | 0.6767896 | 0.9116719 | 0.7859532 | 0.9033233 | 0.9201774 | 0.8095238 | 0.7843943 | 0.6490728 | 0.7413074 | 0.8010657 | 0.8200837 | 0.7783595 | 0.852618 | 0.8150209 | 0.8254675 | 0.71341 | 0.8476923 | 0.7935349 | 0.6774194 | 0.8141593 | 0.9144385 | 0.836478 | 0.8200000 | 0.8677686 | 0.9371069 | 0.8987342 | 0.9099099 | 0.8295964 | 0.8181818 | 0.8633094 | 0.8470948 | 0.8164852 | 0.5949367 | 0.8942308 | 0.8401559 | 0.8680203 | 0.886202 | 0.9273828 | 0.6842105 | 0.7121951 | 0.7916667 | 0.9731286 | 0.974359 | 0.4049587 | 0.7714286 | 1.0000000 | 0.8918919 | 0.7640449 | 0.8444444 | 0.8785714 | 0.9904762 | 0.967033 | 0.9828080 | 0.9221790 | 0.9854015 | 0.9537572 | 0.7241379 | 0.9277108 | -1 | 0.981982 | 0.9553903 | 0.9126214 | 0.9276018 | 0.9469027 | 0.8917647 | 0.9467681 | 0.8818898 | 0.95 | 0.9 | 0.7538462 | 0.8666667 | 0.8111111 | 0.6608187 | 0.4479167 | 0.8319328 | 0.9047619 | 0.9605263 | 0.7878788 | FALSE | TRUE | FALSE |

From looking at chr 10 iof the unshared skipped exon events, these are
actually all the same position ID; they have just been assigned a
different gene So the error likely comes from random things different in
the GTF
