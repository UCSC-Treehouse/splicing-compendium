# Compare merged PSI table from TARGET pilot samples
Cindy Liang (celiang@ucsc.edu)
2026-08-20

## Introduction

Treehouse Compendium version 1 has been run on the TARGET pilot samples.
This workflow includes generating GTFs from alignments, merging GTFs of
all samples, merging separately-produced junction bedfiles for each
sample, and running Shiba on the merged GTF and junction file as inputs
to produce event coordinate and PSI files.

To understand what differences remain between PSI results generated from
Treehouse Compendium v1 and the canonical Shiba run, this notebook
compares what percentage of events remain that can only be found in the
workflow and the canonical Shiba method.

Below are names used to refer to different workflow runs on TARGET pilot
samples:

### PSI tables produced using canonical (unaltered) Shiba scripts

**Combined**: Refers to PSI results generated from the unaltered Shiba
scripts, using `shiba.py`:

    shiba.py -p 15 exploration/merging-psi-tables/target_pilot/shiba_config.yaml

There is a bug in Shiba v0.8.1 that causes some `junctions.bed` position
IDs (coordinates) to be duplicated. To check if duplicated junctions
impact PSI calculation, I also compare the PSI tables below:

**Combined_dedup**: PSI values generated from running only the Shiba
`PSI.py` script on the following inputs and command:

- event coordinates files generated from the Combined (canonical Shiba)
  method above

- `junctions.bed` files generated from Combined, then duplicated with
  the deduplication commands used in the v1 snakemake workflow in the
  `merge_jucntions` rule

<!-- -->

    python $CONDA_PREFIX/share/shiba-0.8.1-0/src/psi.py -m \
    10 -p \
    30 -v \
    --onlypsi \
    /private/groups/treehouse/working-projects/celiang/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/junctions/deduplicated_junctions.bed \
    /private/groups/treehouse/working-projects/celiang/splicing-compendium/exploration/merging-psi-tables/target_pilot/shiba_combined/events \
    v1_dedup_target_pilot_combined

### PSI table results produced by running Treehouse Compendium v1:

Two replicates of Treehosue Compendium v1 workflow were run to check
whether repeated runs of the workflow (same command, different output
directory) produce different PSI tables.

**rep1**: Refers to PSI results generated from the most recent
compendium v1 workflow method **rep2**: Refers to PSI results generated
from the most recent compendium v1 workflow method, generated with the
same commands and parameters as rep1 but written in a different output
directory.

### Comparisons in notebook

Are combined and deduplicated results identical?

- combined vs. combined_dedup

Are there different events in Compendium vs. canonical Shiba methods?

- rep1 vs. combined

Are there differences in PSI tables if compendium workflow is run twice
(are two replicate runs different?)

- rep1 vs. rep2

Are GTFs generated from the same sample in rep1 and rep2 identical?

For all PSI tables loaded in, NA values assigned by Shiba within each
method listed above are replaced with -1, to distinguish NAs from splice
events with low coverage apart from NAs from splice events that are not
matched between each method.

## Set up

### Define functions

``` r
## Reading in PSI tables ##

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

## Manipulating PSI tables for comparison ##

# merge tables by event_type", "pos_id", and "label" only (consider event the same if position is the same regardless of gene)
join_tables <- function(table1, table2) {
  
  # create the merged table of PSI tables to compare
  joined_table <- dplyr::full_join(
  table1,
  table2,
  by = c("event_type", "pos_id", "label"),
  # label PSI values by table they came from
  suffix = c("_table1", "_table2")) |>
  # categorize events by whether they are in the v1 or rep2 tables
  dplyr::mutate(
    # Count events in one table or other if the gene in the other table is NA
    table1_only = is.na(gene_id_table2),
    table2_only = is.na(gene_id_table1),
    shared_pos_id = !(table1_only | table2_only),
    gene_difference = shared_pos_id & gene_id_table1 != gene_id_table2,
    # Label events that are completely identical if event_type/pos_id/label fields are the same and gene_ids are the same
    identical_pos_id_and_gene_id = shared_pos_id & ! gene_difference
  )
  
  return(joined_table)

}

# Summarize events shared and unshared between tables
summarize_shared_events <- function(joined_table) {
  # calculate fraction of events only present in either table and fraction events that are shared
  shared_events_summary <- joined_table |> dplyr::summarise(
    .by = c(event_type, label),
    # count number of events in each event type and annotation category
    total = dplyr::n(),
    # count event types that have share position IDs (includes events with different gene_ids)
    shared_pos_id_count = sum(shared_pos_id),
    # count events in each share category
    shared_pos_and_gene_id_count = sum(identical_pos_id_and_gene_id),
    shared_pos_diff_gene_count = sum(gene_difference),
    table1_only_count = sum(table1_only),
    table2_only_count = sum(table2_only),
    # calculate percentages of each category
    shared_pos_and_gene_id_percent = shared_pos_and_gene_id_count / total * 100,
    shared_pos_diff_gene_percent = shared_pos_diff_gene_count / total * 100,
    table1_only_percent = table1_only_count / total * 100,
    table2_only_percent = table2_only_count / total * 100,
  )
  
  # return summary table
  return(shared_events_summary)
}

# compare if PSI values are identical for a set of event IDs
# compare shared events between two psi tables
compare_shared_events <- function(
    comparison_table, 
    psi_table1, 
    psi_table2, 
    event_id_filter_expr, 
    event_type_filter_expr
    ) {
  shared_ids <- comparison_table |>
    dplyr::filter(
      {{ event_id_filter_expr }},
      {{ event_type_filter_expr }}
    ) |>
    dplyr::pull(pos_id)

  # extract pos_ids from filtered IDs to compare
  extract_shared <- function(psi_table) {
    psi_table |>
      dplyr::filter(pos_id %in% shared_ids) |>
      # filter only for fields to compare (pos_id, psi values of samples)
      dplyr::select(pos_id, ends_with("_PSI")) |>
      # sort by pos_id so order is the same
      dplyr::arrange(pos_id)
  }

  all.equal(extract_shared(psi_table1), extract_shared(psi_table2))
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

### Combined vs. combined_dedup

Are these PSI tables the same?

#### Compare dimensions

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

Dimensions between tables are identical.

#### Check whether combined and deduplicated combined tables are identical

``` r
identical(combined_psi_table, dedup_combined_psi_table)
```

    [1] TRUE

`combined` and `deduplicated_combined` results are identical (duplicated
junctions in the canonical junctions.bed file do not alter PSI
calculations). Moving forward, the combined table is compared to the
other methods.

### Treehouse Compendium workflow rep1 vs. Combined table

#### Compare dimensions

``` r
dim(rep1_psi_table)
```

    [1] 746294     92

``` r
dim(combined_psi_table)
```

    [1] 746601     92

There are ~400 more unique splice events (rows) in the combined table,
compared to the rep1 table;  
the canonical Shiba and Treehouse Compendium methods produce different
results.

#### Check fraction of splice events that are unshared between methods

First, quantify differences present in each table. To compare
differences, rep1 and combined tables are merged on “event_type”,
“pos_id”, “gene_id”, and “label” columns. (data tables are joined so
each row corresponds to a splice event, defined by the combination of
event_type, pos_id, gene_id, and label columns) Columns correspond to
PSI values for a sample and are suffixed by whether they appear in the
rep1 table, combined table, or both. -1 PSI values indicate missing PSI
values from low junction readsl. An NA indicates that a splice event is
not shared between the tables (splice event is only detected in one
method).

``` r
# merge rep1 and combined tables and label columns by whether the event is shared between tables
rep1_vs_combined_table <- join_tables(rep1_psi_table, combined_psi_table)

# summarize table
summary_with_genes_rep1_vs_combined <- summarize_shared_events(rep1_vs_combined_table)

summary_with_genes_rep1_vs_combined |> 
  # no need to look at event types that are all the same
  dplyr::filter(
    shared_pos_and_gene_id_percent != 100
  ) |>
  dplyr::select(
    event_type,
    label,
    total,
    shared_pos_and_gene_id_percent,
    shared_pos_diff_gene_percent,
    rep1_only =  table1_only_percent,
    combined_only = table2_only_percent
  )
```

<div id="tbl-rep1_vs_combined_event_comparison">

Table 1: Percentage of events only in the rep1 or combined PSI tables

<div class="cell-output-display">

| event_type | label | total | shared_pos_and_gene_id_percent | shared_pos_diff_gene_percent | rep1_only | combined_only |
|:---|:---|---:|---:|---:|---:|---:|
| se | annotated | 102922 | 99.34999 | 0.6500068 | 0.0000000 | 0.0000000 |
| afe | annotated | 158653 | 99.97038 | 0.0296244 | 0.0000000 | 0.0000000 |
| ale | annotated | 130114 | 99.98770 | 0.0122969 | 0.0000000 | 0.0000000 |
| five | annotated | 32329 | 99.56077 | 0.4392341 | 0.0000000 | 0.0000000 |
| three | annotated | 37859 | 99.31588 | 0.6841174 | 0.0000000 | 0.0000000 |
| mse | annotated | 62472 | 99.81272 | 0.1872839 | 0.0000000 | 0.0000000 |
| mxe | annotated | 799 | 99.74969 | 0.2503129 | 0.0000000 | 0.0000000 |
| ri | unannotated | 36246 | 98.37223 | 0.1876069 | 0.2814104 | 1.1587486 |
| ri | annotated | 13334 | 97.52512 | 0.6374681 | 0.9599520 | 0.8774561 |

</div>

</div>

Except for retained intron event types, all differences between the v1
compendium table and canonical Shiba table are from cases where the
position ID is the same, but gene IDs are different. These differences
are under 1% of the total events detected.

Based on how exon-intron junctions (which are used to calculate retained
intron PSI values) are handled in the compendium workflow, we expect
retained intron types to have differences between our method and the
canonical Shiba method.

#### Do non-RI position IDs that are shared between tables (including same gene ID) have the same PSI values?

PSI values are identical in events where `pos_id`, `event_type`,
`gene_id`, and `label` fields are the same.

``` r
compare_shared_events(
    rep1_vs_combined_table, 
    rep1_psi_table, 
    combined_psi_table, 
    identical_pos_id_and_gene_id == TRUE, 
    event_type != "ri"
    )
```

    [1] TRUE

#### Do RI position IDs that are shared between tables (including same gene ID) have the same PSI values?

``` r
compare_shared_events(
    rep1_vs_combined_table, 
    rep1_psi_table, 
    combined_psi_table, 
    identical_pos_id_and_gene_id == TRUE
    )
```

     [1] "Component \"SRR1559043_PSI\": Mean relative difference: 2.129928"
     [2] "Component \"SRR1559044_PSI\": Mean relative difference: 2.101152"
     [3] "Component \"SRR1559052_PSI\": Mean relative difference: 2.103862"
     [4] "Component \"SRR1559054_PSI\": Mean relative difference: 2.068997"
     [5] "Component \"SRR1559075_PSI\": Mean relative difference: 2.083076"
     [6] "Component \"SRR1559100_PSI\": Mean relative difference: 2.078521"
     [7] "Component \"SRR1559105_PSI\": Mean relative difference: 2.084496"
     [8] "Component \"SRR1559133_PSI\": Mean relative difference: 2.075777"
     [9] "Component \"SRR1559134_PSI\": Mean relative difference: 2.08605" 
    [10] "Component \"SRR1559145_PSI\": Mean relative difference: 2.066501"
    [11] "Component \"SRR1559160_PSI\": Mean relative difference: 2.082591"
    [12] "Component \"SRR1559164_PSI\": Mean relative difference: 2.059539"
    [13] "Component \"SRR1559177_PSI\": Mean relative difference: 2.053915"
    [14] "Component \"SRR1559183_PSI\": Mean relative difference: 2.069945"
    [15] "Component \"SRR1559184_PSI\": Mean relative difference: 2.070285"
    [16] "Component \"SRR1712453_PSI\": Mean relative difference: 2.084977"
    [17] "Component \"SRR1712454_PSI\": Mean relative difference: 2.089746"
    [18] "Component \"SRR1712455_PSI\": Mean relative difference: 2.090793"
    [19] "Component \"SRR1712456_PSI\": Mean relative difference: 2.085432"
    [20] "Component \"SRR1712457_PSI\": Mean relative difference: 2.091937"
    [21] "Component \"SRR1712458_PSI\": Mean relative difference: 2.098455"
    [22] "Component \"SRR1712459_PSI\": Mean relative difference: 2.096288"
    [23] "Component \"SRR1712460_PSI\": Mean relative difference: 2.084065"
    [24] "Component \"SRR1712461_PSI\": Mean relative difference: 2.09238" 
    [25] "Component \"SRR1712462_PSI\": Mean relative difference: 2.089491"
    [26] "Component \"SRR1712463_PSI\": Mean relative difference: 2.106456"
    [27] "Component \"SRR1712464_PSI\": Mean relative difference: 2.099494"
    [28] "Component \"SRR1712465_PSI\": Mean relative difference: 2.104983"
    [29] "Component \"SRR1784865_PSI\": Mean relative difference: 2.08046" 
    [30] "Component \"SRR1784867_PSI\": Mean relative difference: 2.064966"
    [31] "Component \"SRR1791016_PSI\": Mean relative difference: 2.10264" 
    [32] "Component \"SRR1791028_PSI\": Mean relative difference: 2.07308" 
    [33] "Component \"SRR1791108_PSI\": Mean relative difference: 2.098791"
    [34] "Component \"SRR1796863_PSI\": Mean relative difference: 2.155389"
    [35] "Component \"SRR1796867_PSI\": Mean relative difference: 2.197845"
    [36] "Component \"SRR1796893_PSI\": Mean relative difference: 2.115253"
    [37] "Component \"SRR1796906_PSI\": Mean relative difference: 2.125816"
    [38] "Component \"SRR1796912_PSI\": Mean relative difference: 2.11745" 
    [39] "Component \"SRR1796939_PSI\": Mean relative difference: 2.081003"
    [40] "Component \"SRR1796967_PSI\": Mean relative difference: 2.086844"
    [41] "Component \"SRR1796990_PSI\": Mean relative difference: 2.080954"
    [42] "Component \"SRR1797014_PSI\": Mean relative difference: 2.072259"
    [43] "Component \"SRR1797024_PSI\": Mean relative difference: 2.091501"
    [44] "Component \"SRR1797033_PSI\": Mean relative difference: 2.092129"
    [45] "Component \"SRR1797034_PSI\": Mean relative difference: 2.079456"
    [46] "Component \"SRR1797035_PSI\": Mean relative difference: 2.086948"
    [47] "Component \"SRR1797039_PSI\": Mean relative difference: 2.07848" 
    [48] "Component \"SRR1797052_PSI\": Mean relative difference: 2.08184" 
    [49] "Component \"SRR1797053_PSI\": Mean relative difference: 2.093787"
    [50] "Component \"SRR1797055_PSI\": Mean relative difference: 2.075565"
    [51] "Component \"SRR1797057_PSI\": Mean relative difference: 2.098048"
    [52] "Component \"SRR1797087_PSI\": Mean relative difference: 2.235139"
    [53] "Component \"SRR1797107_PSI\": Mean relative difference: 2.148125"
    [54] "Component \"SRR1797111_PSI\": Mean relative difference: 2.072445"
    [55] "Component \"SRR1799022_PSI\": Mean relative difference: 2.148893"
    [56] "Component \"SRR1799025_PSI\": Mean relative difference: 2.175363"
    [57] "Component \"SRR1799041_PSI\": Mean relative difference: 2.110227"
    [58] "Component \"SRR1799042_PSI\": Mean relative difference: 2.081482"
    [59] "Component \"SRR1799057_PSI\": Mean relative difference: 2.10009" 
    [60] "Component \"SRR1799058_PSI\": Mean relative difference: 2.090935"
    [61] "Component \"SRR1799059_PSI\": Mean relative difference: 2.097025"
    [62] "Component \"SRR1799061_PSI\": Mean relative difference: 2.103617"
    [63] "Component \"SRR1799062_PSI\": Mean relative difference: 2.090659"
    [64] "Component \"SRR1799067_PSI\": Mean relative difference: 2.088598"
    [65] "Component \"SRR1799069_PSI\": Mean relative difference: 2.112548"
    [66] "Component \"SRR1799081_PSI\": Mean relative difference: 2.104585"
    [67] "Component \"SRR1810588_PSI\": Mean relative difference: 2.228146"
    [68] "Component \"SRR2042833_PSI\": Mean relative difference: 2.090853"
    [69] "Component \"SRR2042845_PSI\": Mean relative difference: 2.049418"
    [70] "Component \"SRR2042853_PSI\": Mean relative difference: 2.075527"
    [71] "Component \"SRR2042854_PSI\": Mean relative difference: 2.076128"
    [72] "Component \"SRR2042856_PSI\": Mean relative difference: 2.083871"
    [73] "Component \"SRR2083154_PSI\": Mean relative difference: 2.21093" 
    [74] "Component \"SRR2083162_PSI\": Mean relative difference: 2.212325"
    [75] "Component \"SRR2083171_PSI\": Mean relative difference: 2.205065"
    [76] "Component \"SRR2083176_PSI\": Mean relative difference: 2.23737" 
    [77] "Component \"SRR2083188_PSI\": Mean relative difference: 2.237047"
    [78] "Component \"SRR2239703_PSI\": Mean relative difference: 2.147759"
    [79] "Component \"SRR2239717_PSI\": Mean relative difference: 2.13347" 
    [80] "Component \"SRR3162160_PSI\": Mean relative difference: 2.125868"
    [81] "Component \"SRR3162195_PSI\": Mean relative difference: 2.094416"
    [82] "Component \"SRR3162212_PSI\": Mean relative difference: 2.156903"
    [83] "Component \"SRR3162237_PSI\": Mean relative difference: 2.120659"
    [84] "Component \"SRR3162253_PSI\": Mean relative difference: 2.16378" 
    [85] "Component \"SRR4376029_PSI\": Mean relative difference: 2.076315"
    [86] "Component \"SRR4416297_PSI\": Mean relative difference: 2.113849"
    [87] "Component \"SRR4419554_PSI\": Mean relative difference: 2.128789"
    [88] "Component \"SRR4419565_PSI\": Mean relative difference: 2.080833"

Even in events where `pos_id`, `event_type`, `gene_id`, and `label`
fields are the same, retained intron events between both tables are not
identical. Based on these results, exclude RI event types from
compendium v1 release

#### Do shared position IDs with different gene IDs have the same PSI values?

``` r
compare_shared_events(
    rep1_vs_combined_table, 
    rep1_psi_table, 
    combined_psi_table, 
    shared_pos_id & gene_difference, 
    event_type != "ri"
    )
```

    [1] TRUE

When retained intron events are excluded, PSI values of tables produced
by the compendium workflow and canonical Shiba method in cases where the
position IDs are the same but gene IDs are different are identical.

### Rep1 vs. rep2

Because all differences in the PSI tables stem from `pos_id`s being
assigned different genes, I suspect the difference comes from the GTF
generation or merge steps in the workflow.

To check this, see if the same kind of differences (same pos_ids present
but assigned to different genes) are present in PSI tables generated
from identical Compendium workflow runs (rep1 and rep2).

#### Check fraction of splice events that are unshared between methods

First, quantify differences present in each table. To compare
differences, rep1 and combined tables are merged on “event_type”,
“pos_id”, “gene_id”, and “label” columns. (data tables are joined so
each row corresponds to a splice event, defined by the combination of
event_type, pos_id, gene_id, and label columns) Columns correspond to
PSI values for a sample and are suffixed by whether they appear in the
rep1 table, combined table, or both. -1 PSI values indicate missing PSI
values from low junction readsl. An NA indicates that a splice event is
not shared between the tables (splice event is only detected in one
method).

``` r
# merge rep1 and combined tables and label columns by whether the event is shared between tables
rep1_vs_rep2_table <- join_tables(rep1_psi_table, rep2_psi_table)

# summarize table
summary_rep1_vs_rep2 <- summarize_shared_events(rep1_vs_rep2_table)

summary_rep1_vs_rep2 |> 
  # no need to look at event types that are all the same
  dplyr::filter(
    shared_pos_and_gene_id_percent != 100
  ) |>
  dplyr::select(
  event_type,
  label,
  total,
  shared_pos_and_gene_id_percent,
  shared_pos_diff_gene_percent,
  rep1_only = table1_only_percent,
  rep2_only = table2_only_percent)
```

<div id="tbl-rep1_vs_rep2_event_comparison">

Table 2: Percentage of events only in the rep1 (table1) or rep2 (table2)
PSI tables

<div class="cell-output-display">

| event_type | label | total | shared_pos_and_gene_id_percent | shared_pos_diff_gene_percent | rep1_only | rep2_only |
|:---|:---|---:|---:|---:|---:|---:|
| se | annotated | 102922 | 99.34125 | 0.6587513 | 0.0000000 | 0.0000000 |
| afe | annotated | 158653 | 99.95588 | 0.0441214 | 0.0000000 | 0.0000000 |
| ale | annotated | 130114 | 99.98386 | 0.0161397 | 0.0000000 | 0.0000000 |
| five | annotated | 32329 | 99.55149 | 0.4485137 | 0.0000000 | 0.0000000 |
| three | annotated | 37859 | 99.31852 | 0.6814760 | 0.0000000 | 0.0000000 |
| mse | annotated | 62472 | 99.79991 | 0.2000896 | 0.0000000 | 0.0000000 |
| mxe | annotated | 799 | 99.74969 | 0.2503129 | 0.0000000 | 0.0000000 |
| ri | unannotated | 35933 | 99.18181 | 0.1975900 | 0.3228230 | 0.2977764 |
| ri | annotated | 13333 | 97.77994 | 0.5475137 | 0.8025201 | 0.8700218 |

</div>

</div>

Similar to @rep1_vs_combined_frac_unshared, excluding retained intron
events, all event ID differences stem from cases where the gene ID is
different for the same `pos_id`.

#### Do all position IDs that are shared between tables (including same gene ID) have the same PSI values?

``` r
compare_shared_events(
    rep1_vs_rep2_table, 
    rep1_psi_table, 
    rep2_psi_table,
    identical_pos_id_and_gene_id == TRUE
    )
```

    [1] TRUE

Including retained intron event types, PSI values are identical in
events where pos_id, event_type, gene_id, and label are the same.

#### Do shared position IDs with different gene IDs have the same PSI values?

First, check non-retained-intron event types

``` r
compare_shared_events(
    rep1_vs_combined_table, 
    rep1_psi_table, 
    combined_psi_table, 
    shared_pos_id & gene_difference, 
    event_type != "ri"
    )
```

    [1] TRUE

Next, check if retained intron events with same pos_ids but different
gene_ids have identical PSI values

``` r
compare_shared_events(
    rep1_vs_combined_table, 
    rep1_psi_table, 
    combined_psi_table, 
    shared_pos_id & gene_difference, 
    event_type == "ri"
    )
```

     [1] "Component \"SRR1559044_PSI\": Mean absolute difference: 0.1123596" 
     [2] "Component \"SRR1559052_PSI\": Mean absolute difference: 0.1448763" 
     [3] "Component \"SRR1559054_PSI\": Mean absolute difference: 0.12"      
     [4] "Component \"SRR1559075_PSI\": Mean relative difference: 2"         
     [5] "Component \"SRR1559105_PSI\": Mean absolute difference: 0.1746032" 
     [6] "Component \"SRR1559133_PSI\": Mean absolute difference: 0.1350806" 
     [7] "Component \"SRR1559134_PSI\": Mean absolute difference: 0.199115"  
     [8] "Component \"SRR1559145_PSI\": Mean relative difference: 2.386223"  
     [9] "Component \"SRR1559164_PSI\": Mean absolute difference: 0.1886848" 
    [10] "Component \"SRR1559177_PSI\": Mean absolute difference: 0.08108108"
    [11] "Component \"SRR1559183_PSI\": Mean absolute difference: 0.1349481" 
    [12] "Component \"SRR1559184_PSI\": Mean absolute difference: 0.1162791" 
    [13] "Component \"SRR1712453_PSI\": Mean absolute difference: 0.1149817" 
    [14] "Component \"SRR1712454_PSI\": Mean absolute difference: 0.182308"  
    [15] "Component \"SRR1712455_PSI\": Mean absolute difference: 0.1958143" 
    [16] "Component \"SRR1712456_PSI\": Mean absolute difference: 0.1344086" 
    [17] "Component \"SRR1712457_PSI\": Mean absolute difference: 0.295082"  
    [18] "Component \"SRR1712458_PSI\": Mean absolute difference: 0.2158696" 
    [19] "Component \"SRR1712459_PSI\": Mean absolute difference: 0.2628205" 
    [20] "Component \"SRR1712462_PSI\": Mean absolute difference: 0.3793103" 
    [21] "Component \"SRR1712464_PSI\": Mean relative difference: 2.392025"  
    [22] "Component \"SRR1712465_PSI\": Mean absolute difference: 0.1546985" 
    [23] "Component \"SRR1784865_PSI\": Mean absolute difference: 0.1284916" 
    [24] "Component \"SRR1784867_PSI\": Mean absolute difference: 0.2144397" 
    [25] "Component \"SRR1791016_PSI\": Mean absolute difference: 0.2584909" 
    [26] "Component \"SRR1791028_PSI\": Mean absolute difference: 0.2625"    
    [27] "Component \"SRR1791108_PSI\": Mean absolute difference: 0.3927348" 
    [28] "Component \"SRR1796863_PSI\": Mean absolute difference: 0.3565219" 
    [29] "Component \"SRR1796867_PSI\": Mean absolute difference: 0.4068219" 
    [30] "Component \"SRR1796893_PSI\": Mean absolute difference: 0.2241087" 
    [31] "Component \"SRR1796906_PSI\": Mean absolute difference: 0.06775068"
    [32] "Component \"SRR1796912_PSI\": Mean absolute difference: 0.09090909"
    [33] "Component \"SRR1796939_PSI\": Mean absolute difference: 0.2889564" 
    [34] "Component \"SRR1796967_PSI\": Mean absolute difference: 0.105802"  
    [35] "Component \"SRR1796990_PSI\": Mean absolute difference: 0.2141058" 
    [36] "Component \"SRR1797014_PSI\": Mean absolute difference: 0.1676768" 
    [37] "Component \"SRR1797024_PSI\": Mean relative difference: 2.524862"  
    [38] "Component \"SRR1797033_PSI\": Mean absolute difference: 0.222973"  
    [39] "Component \"SRR1797034_PSI\": Mean absolute difference: 0.08521739"
    [40] "Component \"SRR1797035_PSI\": Mean absolute difference: 0.1229508" 
    [41] "Component \"SRR1797039_PSI\": Mean absolute difference: 0.1456311" 
    [42] "Component \"SRR1797052_PSI\": Mean relative difference: 2.005747"  
    [43] "Component \"SRR1797053_PSI\": Mean absolute difference: 0.2"       
    [44] "Component \"SRR1797055_PSI\": Mean relative difference: 2.789051"  
    [45] "Component \"SRR1797057_PSI\": Mean absolute difference: 0.2984293" 
    [46] "Component \"SRR1797087_PSI\": Mean absolute difference: 0.1381302" 
    [47] "Component \"SRR1797107_PSI\": Mean absolute difference: 0.08333333"
    [48] "Component \"SRR1797111_PSI\": Mean absolute difference: 0.1666667" 
    [49] "Component \"SRR1799025_PSI\": Mean absolute difference: 0.09090909"
    [50] "Component \"SRR1799041_PSI\": Mean relative difference: 1.917407"  
    [51] "Component \"SRR1799042_PSI\": Mean absolute difference: 0.1407035" 
    [52] "Component \"SRR1799057_PSI\": Mean absolute difference: 0.06493506"
    [53] "Component \"SRR1799058_PSI\": Mean absolute difference: 0.05074627"
    [54] "Component \"SRR1799059_PSI\": Mean absolute difference: 0.02604167"
    [55] "Component \"SRR1799061_PSI\": Mean absolute difference: 0.1173134" 
    [56] "Component \"SRR1799062_PSI\": Mean absolute difference: 0.1777778" 
    [57] "Component \"SRR1799067_PSI\": Mean relative difference: 1.766059"  
    [58] "Component \"SRR1799069_PSI\": Mean relative difference: 2.092593"  
    [59] "Component \"SRR1799081_PSI\": Mean absolute difference: 0.04115226"
    [60] "Component \"SRR1810588_PSI\": Mean absolute difference: 0.06725664"
    [61] "Component \"SRR2042833_PSI\": Mean absolute difference: 0.0877193" 
    [62] "Component \"SRR2042853_PSI\": Mean relative difference: 1.915433"  
    [63] "Component \"SRR2042854_PSI\": Mean absolute difference: 0.05184981"
    [64] "Component \"SRR2042856_PSI\": Mean absolute difference: 0.2171972" 
    [65] "Component \"SRR2083154_PSI\": Mean relative difference: 2.410384"  
    [66] "Component \"SRR2083162_PSI\": Mean relative difference: 1.777648"  
    [67] "Component \"SRR2083171_PSI\": Mean absolute difference: 0.2473451" 
    [68] "Component \"SRR2083176_PSI\": Mean absolute difference: 0.09090909"
    [69] "Component \"SRR2083188_PSI\": Mean relative difference: 2.006406"  
    [70] "Component \"SRR2239703_PSI\": Mean absolute difference: 0.05555556"
    [71] "Component \"SRR2239717_PSI\": Mean absolute difference: 0.06306306"
    [72] "Component \"SRR3162160_PSI\": Mean absolute difference: 0.207986"  
    [73] "Component \"SRR3162195_PSI\": Mean relative difference: 2.552941"  
    [74] "Component \"SRR3162212_PSI\": Mean absolute difference: 0.1813538" 
    [75] "Component \"SRR3162237_PSI\": Mean absolute difference: 0.1555556" 
    [76] "Component \"SRR3162253_PSI\": Mean relative difference: 1.784615"  
    [77] "Component \"SRR4376029_PSI\": Mean relative difference: 2.311277"  
    [78] "Component \"SRR4416297_PSI\": Mean absolute difference: 0.1598837" 
    [79] "Component \"SRR4419554_PSI\": Mean relative difference: 2.228387"  
    [80] "Component \"SRR4419565_PSI\": Mean relative difference: 1.923677"  

In replicate compendium runs, retained intron events with the same
`gene_id` different `pos_id` fields have different PSI values.
