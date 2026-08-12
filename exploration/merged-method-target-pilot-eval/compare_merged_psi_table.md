# Compare merged PSI table from TARGET pilot samples
Cindy Liang (celiang@ucsc.edu)
2026-08-12

## Introduction

Compendium workflow v1 (workflow with updated GTF generation and merging
rules) has been run on the TARGET pilot samples. To understand what
differences remain between PSI results generated from this workflow and
the canonical Shiba run, this notebook compares what percentage of
events remain that can only be found in the current v1 workflow and the
canonical Shiba run methods.

**v1**: Refers to PSI results generated from the most recent compendium
workflow method

**rep2_v1**: Refers to PSI results generated from the most recent
compendium workflow method, generated with the same commands and
parameters as v1 but written in a different output directory. rep2_v1
will be used to test whether identical runs of the workflow will result
in the same proportion of different events, due to randomness from GTFs
being generated from scratch in each workflow run.

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

For all PSI tables loaded in, NA values assigned by Shiba within each
method listed above are replaced with -1, to distinguish NAs from splice
events with low coverage apart from NAs from splice events that are not
matched between each method.

## Set up

Define functions

``` r
# read merged splice table from v1 workflow on the TARGET pilot samples
read_psi_tables <- function(psi_paths) {
  psi_results <- psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
      # convert PSI values to numeric and replace NA values from shiba representing low coverage with -1
      dplyr::mutate(across(contains("_PSI"), \(x) tidyr::replace_na(as.numeric(x), -1)))
  })

# combine PSI values of samples run together into one dataframe to compare against other PSI tables
psi_table <- purrr::list_rbind(psi_results, names_to = "event_type")

# return the psi table
return(psi_table)
}
```

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
rep2_target_pilot_results_dir <- file.path(exploration_dir, "v1_target_pilot_rep2", "psi")

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

# construct psi table paths of combined run shiba results
shiba_psi_paths <-file.path(combined_splice_results_dir, psi_files)
names(shiba_psi_paths) <- names(psi_files)

# psi table paths of psi.py run on deduplicated combined run junctions.bed
dedup_shiba_psi_paths <-file.path(dedup_combined_results_dir, psi_files)
names(dedup_shiba_psi_paths) <- names(psi_files)

# psi table paths of v1 workflow on target pilot samples
# these psi tables are zipped so can use the 'v1_psi_files' names
v1_psi_paths <-file.path(target_pilot_workflow_results_dir, v1_psi_files)
names(v1_psi_paths) <- names(v1_psi_files)

rep2_v1_psi_paths <- file.path(rep2_target_pilot_results_dir, v1_psi_files)
names(rep2_v1_psi_paths) <- names(v1_psi_files)
```

Read in files

``` r
# read merged splice table from v1 workflow on the TARGET pilot samples
v1_psi_table <- read_psi_tables(v1_psi_paths)
rep2_v1_psi_table <- read_psi_tables(rep2_v1_psi_paths)

# read in 'negative control' canonical shiba run tables
combined_psi_table <- read_psi_tables(shiba_psi_paths)
dedup_combined_psi_table <- read_psi_tables(dedup_shiba_psi_paths)
```

## Compare PSI tables

### Obtain dimensions of each PSI table

``` r
# dimensions of v1 PSI table
dim(v1_psi_table)
```

    [1] 746294     92

``` r
# dimensions of rep2 v1 PSI table
dim(rep2_v1_psi_table)
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

There are ~400 more splice events in the combined PSI tables than the v1
tables.

### Check if PSI values are identical

``` r
# check whether identical runs of the v1 method produce identical results
identical(v1_psi_table, rep2_v1_psi_table)
```

    [1] FALSE

``` r
# check if v1 method result is identical to the deduplicated combined method result
identical(v1_psi_table, dedup_combined_psi_table)
```

    [1] FALSE

``` r
# check if psi values from the combined method before and after deduplication of junctions.bed file are identical
identical(combined_psi_table, dedup_combined_psi_table)
```

    [1] TRUE

Repeated runs of the v1 compendium workflow do not produce identical
results.

Results from the v1 compendium workflow are not identical to results
from the canonical Shiba run.

However, PSI tables from the “combined” run and just running unaltered
`psi.py` on events files from the combined run and deduplicated combined
run `junctions.bed` file are identical. This kind of makes sense to me
because the final PSI tables don’t have duplicated junction `pos_ids`;
so I think the bug didn’t impact the PSI calculation since Shiba just
“picked” the valid junction rows to calculate PSIs on.

### Check fraction of splice events that are unshared between methods

Since the `combined` and `deduplicated_combined` results are identical,
I will make the following comparisons:

- `v1` vs. `rep2_v1`: To see what fraction of events are different when
  the workflow is run multiple times

- `v1` vs. `combined` and `rep2_v1` vs. `combined`: To see if
  differences between the latest workflow and the combined run are the
  same, despite different workflow runs.

The dataframes created in the following cells are event-level dataframes
of splice events that are shared or unshared between each method. Each
row is a unique pos_id/event_id + gene combination. Columns contain each
samples’ PSI value with -1 indicating low read NA and NA indicating
event not shared between methods.

#### v1 vs. rep2 v1 PSI table event differences

``` r
# construct table of v1 and v2 events and label by whether event is shared
v1_vs_rep2_v1_events <- dplyr::full_join(
  v1_psi_table,
  rep2_v1_psi_table,
  by = c("event_type", "pos_id", "gene_id", "label"),
  # label PSI values by table they came from
  suffix = c("_v1", "_rep2")) |>
  # categorize events by whether they are in the v1 or rep2 tables
  dplyr::mutate(
    # v1 events are counted if the values are not all NA
    v1_event = ! dplyr::if_all(ends_with("_v1"), is.na),
    # rep2 events are counted if the values for the event are not all NA
    rep2_event = ! dplyr::if_all(ends_with("_rep2"), is.na),
    shared_event = v1_event & rep2_event
  )

# event-level summary of the number of each splice event type in each splice table
v1_vs_rep2_events_summary <- v1_vs_rep2_v1_events |> dplyr::summarise(
    .by = c(event_type, label),
    # count number of events in each event type and annotation category
    total = dplyr::n(),
    shared_count = sum(shared_event),
    v1_only_count = sum(v1_event) - shared_count,
    rep2_only_count = sum(rep2_event) - shared_count,
    shared_percent = shared_count / total * 100,
    v1_only_percent = v1_only_count / total * 100,
    rep2_only_percent = rep2_only_count / total * 100,
)

v1_vs_rep2_events_summary |> 
  # no need to look at event types that are all the same
  dplyr::filter(
    shared_percent != 100
  ) |>
  dplyr::select(
  event_type,
  label,
  total,
  shared_percent,
  rep2_only_percent,
  v1_only_percent)
```

| event_type | label | total | shared_percent | rep2_only_percent | v1_only_percent |
|:---|:---|---:|---:|---:|---:|
| se | annotated | 103600 | 98.69112 | 0.6544402 | 0.6544402 |
| afe | annotated | 158723 | 99.91180 | 0.0441020 | 0.0441020 |
| ale | annotated | 130135 | 99.96773 | 0.0161371 | 0.0161371 |
| five | annotated | 32474 | 99.10698 | 0.4465111 | 0.4465111 |
| three | annotated | 38117 | 98.64627 | 0.6768633 | 0.6768633 |
| mse | annotated | 62597 | 99.60062 | 0.1996901 | 0.1996901 |
| mxe | annotated | 801 | 99.50062 | 0.2496879 | 0.2496879 |
| ri | unannotated | 36004 | 98.98622 | 0.4943895 | 0.5193867 |
| ri | annotated | 13406 | 97.24750 | 1.4098165 | 1.3426824 |

Even when the same workflow is run twice, a small fraction of
differences remain in annotated event types. The exception is retained
intron events, which we expect to be more different due to how
exon-intron junctions are constructed in this workflow.

Check what the unshared events between `v1` and `rep2` tables look like

``` r
unshared_events <- v1_vs_rep2_v1_events |> dplyr::filter(
  shared_event == FALSE,
  event_type == "se"
) |>
  dplyr::arrange(pos_id)

head(unshared_events)
```

| event_type | pos_id | gene_id | label | SRR1559043_PSI_v1 | SRR1559044_PSI_v1 | SRR1559052_PSI_v1 | SRR1559054_PSI_v1 | SRR1559075_PSI_v1 | SRR1559100_PSI_v1 | SRR1559105_PSI_v1 | SRR1559133_PSI_v1 | SRR1559134_PSI_v1 | SRR1559145_PSI_v1 | SRR1559160_PSI_v1 | SRR1559164_PSI_v1 | SRR1559177_PSI_v1 | SRR1559183_PSI_v1 | SRR1559184_PSI_v1 | SRR1712453_PSI_v1 | SRR1712454_PSI_v1 | SRR1712455_PSI_v1 | SRR1712456_PSI_v1 | SRR1712457_PSI_v1 | SRR1712458_PSI_v1 | SRR1712459_PSI_v1 | SRR1712460_PSI_v1 | SRR1712461_PSI_v1 | SRR1712462_PSI_v1 | SRR1712463_PSI_v1 | SRR1712464_PSI_v1 | SRR1712465_PSI_v1 | SRR1784865_PSI_v1 | SRR1784867_PSI_v1 | SRR1791016_PSI_v1 | SRR1791028_PSI_v1 | SRR1791108_PSI_v1 | SRR1796863_PSI_v1 | SRR1796867_PSI_v1 | SRR1796893_PSI_v1 | SRR1796906_PSI_v1 | SRR1796912_PSI_v1 | SRR1796939_PSI_v1 | SRR1796967_PSI_v1 | SRR1796990_PSI_v1 | SRR1797014_PSI_v1 | SRR1797024_PSI_v1 | SRR1797033_PSI_v1 | SRR1797034_PSI_v1 | SRR1797035_PSI_v1 | SRR1797039_PSI_v1 | SRR1797052_PSI_v1 | SRR1797053_PSI_v1 | SRR1797055_PSI_v1 | SRR1797057_PSI_v1 | SRR1797087_PSI_v1 | SRR1797107_PSI_v1 | SRR1797111_PSI_v1 | SRR1799022_PSI_v1 | SRR1799025_PSI_v1 | SRR1799041_PSI_v1 | SRR1799042_PSI_v1 | SRR1799057_PSI_v1 | SRR1799058_PSI_v1 | SRR1799059_PSI_v1 | SRR1799061_PSI_v1 | SRR1799062_PSI_v1 | SRR1799067_PSI_v1 | SRR1799069_PSI_v1 | SRR1799081_PSI_v1 | SRR1810588_PSI_v1 | SRR2042833_PSI_v1 | SRR2042845_PSI_v1 | SRR2042853_PSI_v1 | SRR2042854_PSI_v1 | SRR2042856_PSI_v1 | SRR2083154_PSI_v1 | SRR2083162_PSI_v1 | SRR2083171_PSI_v1 | SRR2083176_PSI_v1 | SRR2083188_PSI_v1 | SRR2239703_PSI_v1 | SRR2239717_PSI_v1 | SRR3162160_PSI_v1 | SRR3162195_PSI_v1 | SRR3162212_PSI_v1 | SRR3162237_PSI_v1 | SRR3162253_PSI_v1 | SRR4376029_PSI_v1 | SRR4416297_PSI_v1 | SRR4419554_PSI_v1 | SRR4419565_PSI_v1 | SRR1559043_PSI_rep2 | SRR1559044_PSI_rep2 | SRR1559052_PSI_rep2 | SRR1559054_PSI_rep2 | SRR1559075_PSI_rep2 | SRR1559100_PSI_rep2 | SRR1559105_PSI_rep2 | SRR1559133_PSI_rep2 | SRR1559134_PSI_rep2 | SRR1559145_PSI_rep2 | SRR1559160_PSI_rep2 | SRR1559164_PSI_rep2 | SRR1559177_PSI_rep2 | SRR1559183_PSI_rep2 | SRR1559184_PSI_rep2 | SRR1712453_PSI_rep2 | SRR1712454_PSI_rep2 | SRR1712455_PSI_rep2 | SRR1712456_PSI_rep2 | SRR1712457_PSI_rep2 | SRR1712458_PSI_rep2 | SRR1712459_PSI_rep2 | SRR1712460_PSI_rep2 | SRR1712461_PSI_rep2 | SRR1712462_PSI_rep2 | SRR1712463_PSI_rep2 | SRR1712464_PSI_rep2 | SRR1712465_PSI_rep2 | SRR1784865_PSI_rep2 | SRR1784867_PSI_rep2 | SRR1791016_PSI_rep2 | SRR1791028_PSI_rep2 | SRR1791108_PSI_rep2 | SRR1796863_PSI_rep2 | SRR1796867_PSI_rep2 | SRR1796893_PSI_rep2 | SRR1796906_PSI_rep2 | SRR1796912_PSI_rep2 | SRR1796939_PSI_rep2 | SRR1796967_PSI_rep2 | SRR1796990_PSI_rep2 | SRR1797014_PSI_rep2 | SRR1797024_PSI_rep2 | SRR1797033_PSI_rep2 | SRR1797034_PSI_rep2 | SRR1797035_PSI_rep2 | SRR1797039_PSI_rep2 | SRR1797052_PSI_rep2 | SRR1797053_PSI_rep2 | SRR1797055_PSI_rep2 | SRR1797057_PSI_rep2 | SRR1797087_PSI_rep2 | SRR1797107_PSI_rep2 | SRR1797111_PSI_rep2 | SRR1799022_PSI_rep2 | SRR1799025_PSI_rep2 | SRR1799041_PSI_rep2 | SRR1799042_PSI_rep2 | SRR1799057_PSI_rep2 | SRR1799058_PSI_rep2 | SRR1799059_PSI_rep2 | SRR1799061_PSI_rep2 | SRR1799062_PSI_rep2 | SRR1799067_PSI_rep2 | SRR1799069_PSI_rep2 | SRR1799081_PSI_rep2 | SRR1810588_PSI_rep2 | SRR2042833_PSI_rep2 | SRR2042845_PSI_rep2 | SRR2042853_PSI_rep2 | SRR2042854_PSI_rep2 | SRR2042856_PSI_rep2 | SRR2083154_PSI_rep2 | SRR2083162_PSI_rep2 | SRR2083171_PSI_rep2 | SRR2083176_PSI_rep2 | SRR2083188_PSI_rep2 | SRR2239703_PSI_rep2 | SRR2239717_PSI_rep2 | SRR3162160_PSI_rep2 | SRR3162195_PSI_rep2 | SRR3162212_PSI_rep2 | SRR3162237_PSI_rep2 | SRR3162253_PSI_rep2 | SRR4376029_PSI_rep2 | SRR4416297_PSI_rep2 | SRR4419554_PSI_rep2 | SRR4419565_PSI_rep2 | v1_event | rep2_event | shared_event |
|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|:---|
| se | SE@chr10@100509316-100509511@100509102-100516096 | ENSG00000255339.6 | annotated | 1 | 0.8461538 | 0.9245283 | 0.8490566 | 0.9411765 | 0.9622642 | 0.920000 | 1.0000000 | 1.0000000 | 0.8750000 | 0.8800000 | 1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | 0.8500000 | -1.0000000 | 0.8823529 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | 0.9310345 | 1.0000000 | -1.0000000 | -1.0000000 | 0.8461538 | 1.000000 | 1.0000000 | 0.8550725 | 0.6850394 | 0.9069767 | 0.7200000 | 0.8888889 | 0.500000 | 0.8888889 | 0.8333333 | -1.0000000 | -1.0000000 | -1.0000000 | -1 | -1 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1 | 0.9230769 | -1 | 0.7608696 | -1.0000000 | 0.8518519 | 0.8666667 | -1.0000000 | -1.0000000 | -1 | -1.0000000 | -1.0000000 | -1 | -1.0000000 | -1 | -1.0000000 | -1 | -1.0000000 | 0.7435897 | -1 | -1 | -1 | 0.7142857 | -1.0000000 | 0.8909091 | 0.796875 | 0.9640719 | 0.8980392 | 0.7435897 | -1 | -1 | 0.972973 | 0.8382353 | 1.0000000 | 0.9178082 | 0.8222222 | 0.8235294 | 0.8086957 | 0.8776978 | 0.6530612 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | TRUE | FALSE | FALSE |
| se | SE@chr10@100509316-100509511@100509102-100516096 | ENSG00000075826.17 | annotated | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 1 | 0.8461538 | 0.9245283 | 0.8490566 | 0.9411765 | 0.9622642 | 0.920000 | 1.0000000 | 1.0000000 | 0.8750000 | 0.8800000 | 1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | 0.8500000 | -1.0000000 | 0.8823529 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | 0.9310345 | 1.0000000 | -1.0000000 | -1.0000000 | 0.8461538 | 1.000000 | 1.0000000 | 0.8550725 | 0.6850394 | 0.9069767 | 0.7200000 | 0.8888889 | 0.500000 | 0.8888889 | 0.8333333 | -1.0000000 | -1.0000000 | -1.0000000 | -1 | -1 | -1.0000000 | -1.0000000 | -1.0000000 | -1.0000000 | -1 | 0.9230769 | -1 | 0.7608696 | -1.0000000 | 0.8518519 | 0.8666667 | -1.0000000 | -1.0000000 | -1 | -1.0000000 | -1.0000000 | -1 | -1.0000000 | -1 | -1.0000000 | -1 | -1.0000000 | 0.7435897 | -1 | -1 | -1 | 0.7142857 | -1.0000000 | 0.8909091 | 0.796875 | 0.9640719 | 0.8980392 | 0.7435897 | -1 | -1 | 0.972973 | 0.8382353 | 1.0000000 | 0.9178082 | 0.8222222 | 0.8235294 | 0.8086957 | 0.8776978 | 0.6530612 | FALSE | TRUE | FALSE |
| se | SE@chr10@101793686-101793744@101792943-101797980 | ENSG00000120049.20 | annotated | 0 | -1.0000000 | 0.1304348 | 0.0545455 | 0.0551181 | 0.0000000 | 0.044586 | 0.0000000 | 0.0333333 | 0.0000000 | 0.0000000 | 0.0117647 | 0.0000000 | 0.0588235 | 0.0526316 | 0.0000000 | 0.0000000 | 0.0040486 | 0.0000000 | 0.0000000 | 0.0038610 | 0.0000000 | 0.0043668 | 0.0148148 | 0.0200000 | 0.0000000 | 0.0000000 | 0.0050251 | 0.0000000 | 0.000000 | -1.0000000 | 0.2121212 | 0.2173913 | 0.4920635 | 0.6774194 | -1.0000000 | -1.000000 | 0.0000000 | 0.0263158 | 0.0000000 | 0.0000000 | 0.0000000 | 0 | 0 | 0.0000000 | 0.0000000 | 0.0476190 | 0.0400000 | -1 | 0.0000000 | -1 | -1.0000000 | -1.0000000 | 0.0136986 | -1.0000000 | -1.0000000 | 0.0000000 | 0 | 0.0000000 | 0.0000000 | -1 | -1.0000000 | -1 | 0.0000000 | -1 | -1.0000000 | 0.1724138 | -1 | -1 | -1 | 0.0000000 | 0.0181818 | -1.0000000 | -1.000000 | 0.4102564 | -1.0000000 | 0.1666667 | -1 | -1 | 0.030303 | 0.0400000 | 0.1500000 | 0.0333333 | 0.0000000 | 0.0476190 | -1.0000000 | -1.0000000 | 0.0344828 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | TRUE | FALSE | FALSE |
| se | SE@chr10@101793686-101793744@101792943-101797980 | ENSG00000198408.14 | annotated | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 0 | -1.0000000 | 0.1304348 | 0.0545455 | 0.0551181 | 0.0000000 | 0.044586 | 0.0000000 | 0.0333333 | 0.0000000 | 0.0000000 | 0.0117647 | 0.0000000 | 0.0588235 | 0.0526316 | 0.0000000 | 0.0000000 | 0.0040486 | 0.0000000 | 0.0000000 | 0.0038610 | 0.0000000 | 0.0043668 | 0.0148148 | 0.0200000 | 0.0000000 | 0.0000000 | 0.0050251 | 0.0000000 | 0.000000 | -1.0000000 | 0.2121212 | 0.2173913 | 0.4920635 | 0.6774194 | -1.0000000 | -1.000000 | 0.0000000 | 0.0263158 | 0.0000000 | 0.0000000 | 0.0000000 | 0 | 0 | 0.0000000 | 0.0000000 | 0.0476190 | 0.0400000 | -1 | 0.0000000 | -1 | -1.0000000 | -1.0000000 | 0.0136986 | -1.0000000 | -1.0000000 | 0.0000000 | 0 | 0.0000000 | 0.0000000 | -1 | -1.0000000 | -1 | 0.0000000 | -1 | -1.0000000 | 0.1724138 | -1 | -1 | -1 | 0.0000000 | 0.0181818 | -1.0000000 | -1.000000 | 0.4102564 | -1.0000000 | 0.1666667 | -1 | -1 | 0.030303 | 0.0400000 | 0.1500000 | 0.0333333 | 0.0000000 | 0.0476190 | -1.0000000 | -1.0000000 | 0.0344828 | FALSE | TRUE | FALSE |
| se | SE@chr10@101806045-101806143@101804019-101807730 | ENSG00000120049.20 | annotated | 1 | 1.0000000 | 0.9411765 | 0.9891008 | 0.9628483 | 0.9653680 | 1.000000 | 0.9735099 | 0.9831224 | 0.9719626 | 0.9522546 | 0.9591837 | 0.9156627 | 0.9350649 | 0.9815668 | 0.9654747 | 0.9467041 | 0.9652997 | 0.9577465 | 0.9583333 | 0.9427966 | 0.9423299 | 0.9630656 | 1.0000000 | 0.9670200 | 0.9565647 | 0.9457831 | 0.9761194 | 0.9502262 | 0.976699 | 0.9501558 | 1.0000000 | 0.9640103 | 0.9661017 | 0.9760192 | 1.0000000 | 0.974359 | 0.9746835 | 0.9756839 | 0.9647218 | 0.9641256 | 0.9656652 | 1 | 1 | 0.9812207 | 0.9613527 | 0.9348148 | 0.9622642 | -1 | 0.9560976 | 1 | 0.9650757 | 0.9402985 | 0.9416058 | 1.0000000 | 0.9864407 | 0.9130435 | 1 | 0.9694656 | 0.9675325 | 1 | 0.9517241 | 1 | 0.9243697 | 1 | 0.9724771 | 1.0000000 | 1 | -1 | 1 | 0.9621053 | 0.9786856 | 0.9759519 | 1.000000 | 0.9693141 | 1.0000000 | 0.9422383 | 1 | 1 | 0.928934 | 0.9433962 | 0.9722222 | 0.8839590 | 0.9338843 | 0.9710145 | 0.9785523 | 0.9795918 | 1.0000000 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | TRUE | FALSE | FALSE |
| se | SE@chr10@101806045-101806143@101804019-101807730 | ENSG00000198408.14 | annotated | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 1 | 1.0000000 | 0.9411765 | 0.9891008 | 0.9628483 | 0.9653680 | 1.000000 | 0.9735099 | 0.9831224 | 0.9719626 | 0.9522546 | 0.9591837 | 0.9156627 | 0.9350649 | 0.9815668 | 0.9654747 | 0.9467041 | 0.9652997 | 0.9577465 | 0.9583333 | 0.9427966 | 0.9423299 | 0.9630656 | 1.0000000 | 0.9670200 | 0.9565647 | 0.9457831 | 0.9761194 | 0.9502262 | 0.976699 | 0.9501558 | 1.0000000 | 0.9640103 | 0.9661017 | 0.9760192 | 1.0000000 | 0.974359 | 0.9746835 | 0.9756839 | 0.9647218 | 0.9641256 | 0.9656652 | 1 | 1 | 0.9812207 | 0.9613527 | 0.9348148 | 0.9622642 | -1 | 0.9560976 | 1 | 0.9650757 | 0.9402985 | 0.9416058 | 1.0000000 | 0.9864407 | 0.9130435 | 1 | 0.9694656 | 0.9675325 | 1 | 0.9517241 | 1 | 0.9243697 | 1 | 0.9724771 | 1.0000000 | 1 | -1 | 1 | 0.9621053 | 0.9786856 | 0.9759519 | 1.000000 | 0.9693141 | 1.0000000 | 0.9422383 | 1 | 1 | 0.928934 | 0.9433962 | 0.9722222 | 0.8839590 | 0.9338843 | 0.9710145 | 0.9785523 | 0.9795918 | 1.0000000 | FALSE | TRUE | FALSE |

From looking at chr 10 iof the unshared skipped exon events between `v1`
and `rep2` of the v1 workflow results, these are actually all the same
position ID; they have just been assigned a different gene. The event
differences therefore likely come from random differences in the GTF
that occur when GTFs are regenerated from scratch.

#### v1 vs. combined run PSI table event differences

``` r
# construct table of v1 and combined events and label by whether event is shared
v1_vs_combined_events <- dplyr::full_join(
  v1_psi_table,
  combined_psi_table,
  by = c("event_type", "pos_id", "gene_id", "label"),
  # label PSI values by table they came from
  suffix = c("_v1", "_combined")) |>
  # categorize events by whether they are in the v1 or combined tables
  dplyr::mutate(
    # v1 events are counted if the values are not all NA
    v1_event = ! dplyr::if_all(ends_with("_v1"), is.na),
    # combined events are counted if the values for the event are not all NA
    combined_event = ! dplyr::if_all(ends_with("_combined"), is.na),
    shared_event = v1_event & combined_event
  )

# event-level summary of the number of each splice event type in each splice table
v1_vs_combined_events_summary <- v1_vs_combined_events |> dplyr::summarise(
    .by = c(event_type, label),
    # count number of events in each event type and annotation category
    total = dplyr::n(),
    shared_count = sum(shared_event),
    v1_only_count = sum(v1_event) - shared_count,
    combined_only_count = sum(combined_event) - shared_count,
    shared_percent = shared_count / total * 100,
    v1_only_percent = v1_only_count / total * 100,
    combined_only_percent = combined_only_count / total * 100,
)

v1_vs_combined_events_summary |> 
  # no need to look at event types that are all the same
  dplyr::filter(
    shared_percent != 100
  ) |>
  dplyr::select(
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
| afe | annotated | 158700 | 99.94077 | 0.0296156 | 0.0296156 |
| ale | annotated | 130130 | 99.97541 | 0.0122954 | 0.0122954 |
| five | annotated | 32471 | 99.12537 | 0.4373133 | 0.4373133 |
| three | annotated | 38118 | 98.64106 | 0.6794690 | 0.6794690 |
| mse | annotated | 62589 | 99.62613 | 0.1869338 | 0.1869338 |
| mxe | annotated | 801 | 99.50062 | 0.2496879 | 0.2496879 |
| ri | unannotated | 36314 | 98.18803 | 1.3438343 | 0.4681390 |
| ri | annotated | 13419 | 96.90737 | 1.5053283 | 1.5873016 |

#### rep2_v1 vs. combined run event differences

``` r
# construct table of rep2_v1 and combined events and label by whether event is shared
rep2_vs_combined_events <- dplyr::full_join(
  rep2_v1_psi_table,
  combined_psi_table,
  by = c("event_type", "pos_id", "gene_id", "label"),
  # label PSI values by table they came from
  suffix = c("_rep2", "_combined")) |>
  # categorize events by whether they are in the v1 or combined tables
  dplyr::mutate(
    # rep2 v1 events are counted if the values are not all NA
    rep2_event = ! dplyr::if_all(ends_with("_rep2"), is.na),
    # combined events are counted if the values for the event are not all NA
    combined_event = ! dplyr::if_all(ends_with("_combined"), is.na),
    shared_event = rep2_event & combined_event
  )

# event-level summary of the number of each splice event type in each splice table
rep2_vs_combined_events_summary <- rep2_vs_combined_events |> dplyr::summarise(
    .by = c(event_type, label),
    # count number of events in each event type and annotation category
    total = dplyr::n(),
    shared_count = sum(shared_event),
    rep2_only_count = sum(rep2_event) - shared_count,
    combined_only_count = sum(combined_event) - shared_count,
    shared_percent = shared_count / total * 100,
    rep2_only_percent = rep2_only_count / total * 100,
    combined_only_percent = combined_only_count / total * 100,
)

rep2_vs_combined_events_summary |> 
  # no need to look at event types that are all the same
  dplyr::filter(
    shared_percent != 100
  ) |>
  dplyr::select(
  event_type,
  label,
  total,
  shared_percent,
  combined_only_percent,
  rep2_only_percent)
```

| event_type | label | total | shared_percent | combined_only_percent | rep2_only_percent |
|:---|:---|---:|---:|---:|---:|
| se | annotated | 103563 | 98.76211 | 0.6189469 | 0.6189469 |
| afe | annotated | 158754 | 99.87276 | 0.0636204 | 0.0636204 |
| ale | annotated | 130121 | 99.98924 | 0.0053796 | 0.0053796 |
| five | annotated | 32475 | 99.10085 | 0.4495766 | 0.4495766 |
| three | annotated | 38103 | 98.71926 | 0.6403695 | 0.6403695 |
| mse | annotated | 62606 | 99.57193 | 0.2140370 | 0.2140370 |
| mxe | annotated | 803 | 99.00374 | 0.4981320 | 0.4981320 |
| ri | unannotated | 36322 | 98.11960 | 1.3903419 | 0.4900611 |
| ri | annotated | 13440 | 96.66667 | 1.5922619 | 1.7410714 |

Both `v1` and `rep2_v1` methods result in events that are different from
the `combined` method’s events. Similar to differences between the `v1`
and `rep2_v1` methods, these differences are only present in annotated
event types. Although the fraction of these differences are all low,
they tend to differ when the v1 workflow is run different times
(e.g. fraction of events unanswered between `v1` and `combined`, or
between `rep2_v1` and `combined` differ slightly) As expected, in all
comparisons, the retained intron event type is the only category wehre
the percentage of unshared events exceeds 1% and unshared events are
present in both annotated and unannotated event types. We suspect this
stems from exon-intron junctions in the compendium workflow being
defined from GTFs with fewer possibel transcripts.

## Concluding remarks

- Differences in results remaining in the v1 compendium workflow are
  negligible and appear to stem from randomness in GTF generation and
  exon-intron junction definition
- Even with the differences in the retained introns, the fraction of
  unshared events is so low that I think we can keep retained introns in
  the compendium, with a caveat that they are slightly less trustworthy
  than other events types.
