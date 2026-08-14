# Compare merged PSI table from TARGET pilot samples
Cindy Liang (celiang@ucsc.edu)
2026-08-13

## Introduction

Compendium workflow v1 (workflow with updated GTF generation and merging
rules) has been run on the TARGET pilot samples. To understand what
differences remain between PSI results generated from this workflow and
the canonical Shiba run, this notebook compares what percentage of
events remain that can only be found in the current v1 workflow and the
canonical Shiba run methods.

**rep1**: Refers to PSI results generated from the most recent
compendium v1 workflow method

**rep2**: Refers to PSI results generated from the most recent
compendium v1 workflow method, generated with the same commands and
parameters as rep1 but written in a different output directory. rep2_v1
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
# construct paths to psi tables
construct_psi_paths <- function(psi_dir, psi_file_list) {
  # construct paths to individual psi files
  psi_paths <-file.path(psi_dir, psi_file_list)
  # name each path according to event type
  names(psi_paths) <- names(psi_file_list)
  
  # returned the named vector of psi paths
  return(psi_paths)
}

# read in and merge splice tables
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
workflow_psi_files <- c(
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
combined_psi_paths <- construct_psi_paths(combined_splice_results_dir, psi_files)

# psi table paths of psi.py run on deduplicated combined run junctions.bed
dedup_comb_psi_paths <- construct_psi_paths(dedup_combined_results_dir, psi_files)

# psi table paths of v1 workflow on target pilot samples
# these psi tables are zipped so can use the 'v1_psi_files' names
rep1_psi_paths <- construct_psi_paths(target_pilot_workflow_results_dir, workflow_psi_files)
rep2_psi_paths <- construct_psi_paths(rep2_target_pilot_results_dir, workflow_psi_files)
```

Read in files

``` r
# read merged splice table from v1 workflow on the TARGET pilot samples
rep1_psi_table <- read_psi_tables(rep1_psi_paths)
rep2_psi_table <- read_psi_tables(rep2_psi_paths)

# read in 'negative control' canonical shiba run tables
combined_psi_table <- read_psi_tables(combined_psi_paths)
dedup_combined_psi_table <- read_psi_tables(dedup_comb_psi_paths)
```

## Compare PSI tables

### Obtain dimensions of each PSI table

``` r
# dimensions of v1 PSI table
dim(rep1_psi_table)
```

    [1] 746294     92

``` r
# dimensions of rep2 v1 PSI table
dim(rep2_psi_table)
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

### Check fraction of splice events that are unshared between methods

Since the `combined` and `deduplicated_combined` results are identical,
I will make the following comparisons:

- `rep1` vs. `rep2`: To see what fraction of events are different when
  the workflow is run multiple times

The dataframes created in the following cells are event-level dataframes
of splice events that are shared or unshared between each method. Each
row is a unique pos_id/event_id + gene combination. Columns contain each
samples’ PSI value with -1 indicating low read NA and NA indicating
event not shared between methods.

## rep1 vs. rep2 PSI table event differences

### print fraction of events that are unshared

``` r
# construct table of rep1 and rep2 events and label by whether event is shared
rep1_vs_rep2_events <- dplyr::full_join(
  rep1_psi_table,
  rep2_psi_table,
  by = c("event_type", "pos_id", "label"),
  # label PSI values by table they came from
  suffix = c("_rep1", "_rep2")) |>
  # categorize events by whether they are in the v1 or rep2 tables
  dplyr::mutate(
    # v1 events are counted if the values are not all NA
    rep1_event = ! dplyr::if_all(ends_with("_rep1"), is.na),
    # rep2 events are counted if the values for the event are not all NA
    rep2_event = ! dplyr::if_all(ends_with("_rep2"), is.na),
    shared_event = rep1_event & rep2_event
  )

# event-level summary of the number of each splice event type in each splice table
rep1_vs_rep2_events_summary <- rep1_vs_rep2_events |> dplyr::summarise(
    .by = c(event_type, label),
    # count number of events in each event type and annotation category
    total = dplyr::n(),
    shared_count = sum(shared_event),
    rep1_only_count = sum(rep1_event) - shared_count,
    rep2_only_count = sum(rep2_event) - shared_count,
    shared_percent = shared_count / total * 100,
    rep1_only_percent = rep1_only_count / total * 100,
    rep2_only_percent = rep2_only_count / total * 100,
)

rep1_vs_rep2_events_summary |> 
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
  rep1_only_percent)
```

| event_type | label | total | shared_percent | rep2_only_percent | rep1_only_percent |
|:---|:---|---:|---:|---:|---:|
| ri | unannotated | 35933 | 99.37940 | 0.2977764 | 0.3228230 |
| ri | annotated | 13333 | 98.32746 | 0.8700218 | 0.8025201 |

### Check if PSI values are identical

``` r
# extract position ids of events that are shared between tables
rep1_rep2_shared_psi <- rep1_vs_rep2_events |>
  dplyr::filter(
    shared_event == TRUE,
    # let's not look at retained introns for now
    event_type != "ri"
    ) |>
  dplyr::pull(pos_id)

# filter tables for events in shared tables
shared_only_rep1_table <- rep1_psi_table |>
  dplyr::filter(pos_id %in% rep1_rep2_shared_psi) |>
  dplyr::select(pos_id, ends_with("_PSI"))

shared_only_rep2_table <- rep2_psi_table |>
  dplyr::filter(pos_id %in% rep1_rep2_shared_psi) |>
  dplyr::select(pos_id, ends_with("_PSI"))

all.equal(shared_only_rep1_table, shared_only_rep2_table)
```

     [1] "Component \"pos_id\": 69715 string mismatches"                   
     [2] "Component \"SRR1559043_PSI\": Mean relative difference: 1.483315"
     [3] "Component \"SRR1559044_PSI\": Mean relative difference: 1.431269"
     [4] "Component \"SRR1559052_PSI\": Mean relative difference: 1.400423"
     [5] "Component \"SRR1559054_PSI\": Mean relative difference: 1.364795"
     [6] "Component \"SRR1559075_PSI\": Mean relative difference: 1.199781"
     [7] "Component \"SRR1559100_PSI\": Mean relative difference: 1.233691"
     [8] "Component \"SRR1559105_PSI\": Mean relative difference: 1.175447"
     [9] "Component \"SRR1559133_PSI\": Mean relative difference: 1.142479"
    [10] "Component \"SRR1559134_PSI\": Mean relative difference: 1.233245"
    [11] "Component \"SRR1559145_PSI\": Mean relative difference: 1.185312"
    [12] "Component \"SRR1559160_PSI\": Mean relative difference: 1.219263"
    [13] "Component \"SRR1559164_PSI\": Mean relative difference: 1.196658"
    [14] "Component \"SRR1559177_PSI\": Mean relative difference: 1.367918"
    [15] "Component \"SRR1559183_PSI\": Mean relative difference: 1.284216"
    [16] "Component \"SRR1559184_PSI\": Mean relative difference: 1.310272"
    [17] "Component \"SRR1712453_PSI\": Mean relative difference: 1.337926"
    [18] "Component \"SRR1712454_PSI\": Mean relative difference: 1.326938"
    [19] "Component \"SRR1712455_PSI\": Mean relative difference: 1.345445"
    [20] "Component \"SRR1712456_PSI\": Mean relative difference: 1.390126"
    [21] "Component \"SRR1712457_PSI\": Mean relative difference: 1.306888"
    [22] "Component \"SRR1712458_PSI\": Mean relative difference: 1.288512"
    [23] "Component \"SRR1712459_PSI\": Mean relative difference: 1.283701"
    [24] "Component \"SRR1712460_PSI\": Mean relative difference: 1.281702"
    [25] "Component \"SRR1712461_PSI\": Mean relative difference: 1.282"   
    [26] "Component \"SRR1712462_PSI\": Mean relative difference: 1.327855"
    [27] "Component \"SRR1712463_PSI\": Mean relative difference: 1.319016"
    [28] "Component \"SRR1712464_PSI\": Mean relative difference: 1.326263"
    [29] "Component \"SRR1712465_PSI\": Mean relative difference: 1.347159"
    [30] "Component \"SRR1784865_PSI\": Mean relative difference: 1.273438"
    [31] "Component \"SRR1784867_PSI\": Mean relative difference: 1.347157"
    [32] "Component \"SRR1791016_PSI\": Mean relative difference: 1.24397" 
    [33] "Component \"SRR1791028_PSI\": Mean relative difference: 1.386342"
    [34] "Component \"SRR1791108_PSI\": Mean relative difference: 1.335521"
    [35] "Component \"SRR1796863_PSI\": Mean relative difference: 1.293248"
    [36] "Component \"SRR1796867_PSI\": Mean relative difference: 1.195812"
    [37] "Component \"SRR1796893_PSI\": Mean relative difference: 1.397051"
    [38] "Component \"SRR1796906_PSI\": Mean relative difference: 1.291464"
    [39] "Component \"SRR1796912_PSI\": Mean relative difference: 1.25922" 
    [40] "Component \"SRR1796939_PSI\": Mean relative difference: 1.297364"
    [41] "Component \"SRR1796967_PSI\": Mean relative difference: 1.266642"
    [42] "Component \"SRR1796990_PSI\": Mean relative difference: 1.330303"
    [43] "Component \"SRR1797014_PSI\": Mean relative difference: 1.317037"
    [44] "Component \"SRR1797024_PSI\": Mean relative difference: 1.41884" 
    [45] "Component \"SRR1797033_PSI\": Mean relative difference: 1.417557"
    [46] "Component \"SRR1797034_PSI\": Mean relative difference: 1.3245"  
    [47] "Component \"SRR1797035_PSI\": Mean relative difference: 1.345886"
    [48] "Component \"SRR1797039_PSI\": Mean relative difference: 1.30156" 
    [49] "Component \"SRR1797052_PSI\": Mean relative difference: 1.243718"
    [50] "Component \"SRR1797053_PSI\": Mean relative difference: 1.501428"
    [51] "Component \"SRR1797055_PSI\": Mean relative difference: 1.322663"
    [52] "Component \"SRR1797057_PSI\": Mean relative difference: 1.417556"
    [53] "Component \"SRR1797087_PSI\": Mean relative difference: 1.181426"
    [54] "Component \"SRR1797107_PSI\": Mean relative difference: 1.44048" 
    [55] "Component \"SRR1797111_PSI\": Mean relative difference: 1.351465"
    [56] "Component \"SRR1799022_PSI\": Mean relative difference: 1.354374"
    [57] "Component \"SRR1799025_PSI\": Mean relative difference: 1.524506"
    [58] "Component \"SRR1799041_PSI\": Mean relative difference: 1.367856"
    [59] "Component \"SRR1799042_PSI\": Mean relative difference: 1.414652"
    [60] "Component \"SRR1799057_PSI\": Mean relative difference: 1.434093"
    [61] "Component \"SRR1799058_PSI\": Mean relative difference: 1.254455"
    [62] "Component \"SRR1799059_PSI\": Mean relative difference: 1.385652"
    [63] "Component \"SRR1799061_PSI\": Mean relative difference: 1.311708"
    [64] "Component \"SRR1799062_PSI\": Mean relative difference: 1.356123"
    [65] "Component \"SRR1799067_PSI\": Mean relative difference: 1.273098"
    [66] "Component \"SRR1799069_PSI\": Mean relative difference: 1.481294"
    [67] "Component \"SRR1799081_PSI\": Mean relative difference: 1.337429"
    [68] "Component \"SRR1810588_PSI\": Mean relative difference: 1.378175"
    [69] "Component \"SRR2042833_PSI\": Mean relative difference: 1.448634"
    [70] "Component \"SRR2042845_PSI\": Mean relative difference: 1.564858"
    [71] "Component \"SRR2042853_PSI\": Mean relative difference: 1.348048"
    [72] "Component \"SRR2042854_PSI\": Mean relative difference: 1.272542"
    [73] "Component \"SRR2042856_PSI\": Mean relative difference: 1.214568"
    [74] "Component \"SRR2083154_PSI\": Mean relative difference: 1.183345"
    [75] "Component \"SRR2083162_PSI\": Mean relative difference: 1.198706"
    [76] "Component \"SRR2083171_PSI\": Mean relative difference: 1.17097" 
    [77] "Component \"SRR2083176_PSI\": Mean relative difference: 1.132067"
    [78] "Component \"SRR2083188_PSI\": Mean relative difference: 1.193011"
    [79] "Component \"SRR2239703_PSI\": Mean relative difference: 1.47903" 
    [80] "Component \"SRR2239717_PSI\": Mean relative difference: 1.530352"
    [81] "Component \"SRR3162160_PSI\": Mean relative difference: 1.359017"
    [82] "Component \"SRR3162195_PSI\": Mean relative difference: 1.39001" 
    [83] "Component \"SRR3162212_PSI\": Mean relative difference: 1.356906"
    [84] "Component \"SRR3162237_PSI\": Mean relative difference: 1.283128"
    [85] "Component \"SRR3162253_PSI\": Mean relative difference: 1.195055"
    [86] "Component \"SRR4376029_PSI\": Mean relative difference: 1.37848" 
    [87] "Component \"SRR4416297_PSI\": Mean relative difference: 1.360608"
    [88] "Component \"SRR4419554_PSI\": Mean relative difference: 1.335413"
    [89] "Component \"SRR4419565_PSI\": Mean relative difference: 1.3535"  

Even though position IDs of non-RI event types are the same, the PSI
values in rep 1 and rep 2 are different. based on the previous version
of this notebook where I printed out examples of cases where positions
IDs are the same but being assigned different genes, I wonder how much
of differing PSI values stem from cases where junctions in the same
position IDs are assigned different genes when Shiba calculates PSI.

One idea I have is to look at two categories of PSI values for whether
they are identical - PSI values from position IDs that are identical and
assigned to different genes, and PSI values from truly identical pos_ids
assigned to the same gene.
