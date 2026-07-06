# Compare splice event coordinates between merged and combined Shiba runs
Cindy Liang (celiang@ucsc.edu)
2026-07-06

## Background

From analyzing junction count files produced by the merged and combined
Shiba methods in `exploration/check_merged_junctions.qmd`, only
exon-intron junctions were different.

Junctions.bed file exon-intron junctions are counted from exon-intron
boundaries defined by `EVENT_RI.txt`, produced from the GTF of all
samples’ transcripts by `gtf2event.py` in Shiba. Initial analysis of the
GTFs produced by the merged and combined methods also revealed
differences in the GTFs used to create the event coordinate files. We
suspect these differences stem from differences in the order in which
files were passed into `stringtie --merge`, and the use of
multithreading in merging.

We next wanted to answer how much GTF differences impacted similarity in
the coordinates of splice events in the event coordinate files (like
`EVENT_RI.txt`). For instance, if only one event type’s coordinates is
disproportionately impacted by GTF differences, we may decide to not
include quantification of that event type in the first release of the
splice compendium. Alternatively, if differences are present throughout
all event types, we will need to alter how the GTF is processed to make
it more similar to an unaltered Shiba run (“combined” method).

Sources of differences we want to test in this notebook are:

- Use of multithreading in the GTF merge command

- Order of GTFs that are merged in

- Merging in of reference GTF multiple times in the workflow

## Analysis outline

This notebook compares the similarity (Jaccard index) between splice
event coordinate position IDs from event coordinate files produced from
GTFs run with different StringTie conditions:

- **Negative control comparison**: Comparing event coordinates generated
  from GTFs produced with the same commands (expected to be identical):

  - Two GTFs were generated with the unaltered Shiba v 0.8.1
    `bam2gtf.py` script, with one thread, with the same commands
  - Events coordinates files for both GTFs were generated with unaltered
    `gtf2events.py` using 10 threads (the amount passed to this script
    in the splice compendium merged snakemake workflow) with the same
    commands. Position IDs from these coordinate files are compared for
    similarity with Jaccard indices in this notebook.
  - Replicate 1 files and commands used to generate them
    - **Input GTF replicate 1:**
      `python $CONDA_PREFIX/share/shiba-0.8.1-0/src/bam2gtf.py  -i exploration/merging-psi-tables/target_pilot/experiment.tsv  -r references/gencode.v47.primary_assembly.annotation.gtf  -o results/merged_shiba/target_pilot_gtf_tests/one_thread_pilot_gtf.gtf  -p 1  -v`
    - **Events files generated from input GTF 1:**
      `python $CONDA_PREFIX/share/shiba-0.8.1-0/src/gtf2event.py  -i target_pilot_gtf_tests/one_thread_pilot_gtf.gtf  -r references/gencode.v47.primary_assembly.annotation.gtf  -o shiba_pilot_gtf_events  -p 10  -v`
  - Replicate 2 files and commands used to generate them
    - **Input GTF replicate 2:**
      `python $CONDA_PREFIX/share/shiba-0.8.1-0/src/bam2gtf.py  -i exploration/merging-psi-tables/target_pilot/experiment.tsv  -r references/gencode.v47.primary_assembly.annotation.gtf  -o results/merged_shiba/target_pilot_gtf_tests_rep2/one_thread_pilot_gtf.gtf  -p 1  -v`
    - **Events files generated from input GTF 2:**
      `python $CONDA_PREFIX/share/shiba-0.8.1-0/src/gtf2event.py  -i /scratch/celiang/target_pilot_gtfs/target_pilot_gtf_tests_rep2/one_thread_pilot_gtf.gtf  -r /private/groups/treehouse/working-projects/celiang/splicing-compendium/references/gencode.v47.primary_assembly.annotation.gtf  -o /scratch/celiang/shiba_event_to_gtf_tests/shiba_pilot_gtf_events_rep2  -p 10  -v`

- **Multiple threads test:** To what extent does using multiple threads
  change events coordinates defined by GTFs made with the same commands?

  - Two GTFs were generated with the unaltered Shiba v 0.8.1
    `bam2gtf.py` script, with 15 threads, with the same commands

  - Events coordinates files for both GTFs were generated with unaltered
    `gtf2events.py` using 10 threads (the amount passed to this script
    in the splice compendium merged snakemake workflow) with the same
    commands. Position IDs from these coordinate files are compared for
    similarity with Jaccard indices in this notebook.

  - Replicate 1 files and commands used to generate them:

    - **Input GTF replicate 1:**
      `python $CONDA_PREFIX/share/shiba-0.8.1-0/src/bam2gtf.py  -i experiment_scratch.tsv  -r /private/groups/treehouse/working-projects/celiang/splicing-compendium/references/gencode.v47.primary_assembly.annotation.gtf  -o target_pilot_gtfs/15_threads_rep1/15_threads_pilot_gtf.gtf  -p 15  -v`

    - **Events files generated from input GTF 1:**
      `python $CONDA_PREFIX/share/shiba-0.8.1-0/src/gtf2event.py  -i target_pilot_gtfs/15_threads_rep1/15_threads_pilot_gtf.gtf  -r /private/groups/treehouse/working-projects/celiang/splicing-compendium/references/gencode.v47.primary_assembly.annotation.gtf  -o shiba_gtf_to_event_pilot_tests/15threads_events_rep1  -p 10  -v`

  - Replicate 2 files and commands used to generate them:

    - **Input GTF replicate 2:**
      `python $CONDA_PREFIX/share/shiba-0.8.1-0/src/bam2gtf.py  -i experiment_scratch.tsv  -r /private/groups/treehouse/working-projects/celiang/splicing-compendium/references/gencode.v47.primary_assembly.annotation.gtf  -o target_pilot_gtfs/15_threads_rep2/15_threads_pilot_gtf.gtf  -p 15  -v`

    - **Events files generated from input GTF 2:**
      `python $CONDA_PREFIX/share/shiba-0.8.1-0/src/gtf2event.py  -i target_pilot_gtfs/15_threads_rep2/15_threads_pilot_gtf.gtf  -r /private/groups/treehouse/working-projects/celiang/splicing-compendium/references/gencode.v47.primary_assembly.annotation.gtf  -o shiba_gtf_to_event_pilot_tests/15threads_events_rep2  -p 10  -v`

- **Different GTF merge order:** To what extent does merging GTFs in
  different orders (as defined by the sample order in `experiment.tsv`)
  change event coordinate definition?

  - Scrambled merge order files and commands used to generate them:
    - **Input GTF 1 (scrambled merge order):** Was generated by first
      randomly scrambling the row order of `experiment.tsv` with
      `exploration/merged-method-target-pilot-eval/scramble_experiment_file_sample_order.R`
      , then running the unaltered Shiba v 0.8.1 `bam2gtf.py` with the
      scrambled experiment.tsv file using the following command:
      `python $CONDA_PREFIX/share/shiba-0.8.1-0/src/bam2gtf.py -i scrambled_experiment.tsv -r /private/groups/treehouse/working-projects/celiang/splicing-compendium/references/gencode.v47.primary_assembly.annotation.gtf -o target_pilot_gtfs/15_threads_scrambled/15_threads_scrambled_gtf.gtf -p 15 -v`
    - **Events files generated from input GTF 1:**
      `$CONDA_PREFIX/share/shiba-0.8.1-0/src/gtf2event.py -i target_pilot_gtfs/15_threads_scrambled/15_threads_scrambled_gtf.gtf -r /private/groups/treehouse/working-projects/celiang/splicing-compendium/references/gencode.v47.primary_assembly.annotation.gtf -o shiba_gtf_to_event_pilot_tests/15_threads_scrambled -p 10 -v`
  - “Control” merge order files:
    - Input GTF replicate 1 from Multiple Threads Test
    - Events files derived from input GTF replicate 1 from Multiple
      Threads Test

- Merged vs. combined method: To what extent are event coordinates
  defined with the merged and combined methods different?

  - Merged method files and commands used to generate them:

    - **Merged input GTF:** Individual GTFs for each sample were first
      generated with `shiba-0.8.1-0/src/bam2gtf.py` with `Snakefile` by
      merging one sample’s StringTie GTF and the reference GTF. Then,
      these GTFs were merged together, along with the reference with a
      final `stringtie —merge` command with `merge_results.smk`. 4
      threads were used in the merge step.

    - **Events files generated from input merged GTF:** Events
      coordinates files were made using
      `$CONDA_PREFIX/share/shiba-0.8.1-0/src/gtf2event.py` called by
      `merge_results.smk` , with 10 threads

  - Combined method files and commands used to generate them:

    - **Combined input GTF:** Separate GTFs for 88 TARGET pilot samples
      were greated with StringTie, then merged together with the
      reference using `stringtie -–merge` using the unaltered
      `shiba-0.8.1-0/src/bam2gtf` and 15 threads.

    - **Events files generated from input combined GTF:** Events
      coordinates files were made using
      `$CONDA_PREFIX/share/shiba-0.8.1-0/src/gtf2event.py` as part of
      the unaltered `shiba.py` workflow, with 15 threads.

## Setup

### Define functions

``` r
# read in list of splice event coordinate files and bind them into one list of dataframes
read_in_events <- function(file_paths) {
event_coords <- file_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c"))
  })
}

# Calculate list of Jaccard indices for event coordinate file pairs
# Similarity scores between events files made from combined vs. merged method
calculate_jaccard_indices_of_events <- function(
  list1,
  list2
) {
  purrr::map2(
  # read in lists of event coordinate dataframes to compare and iterate the two simultaneously
  list1,
  list2,
  # take the matching two dataframes from the input lists
  \(df1, df2){
    # obtain set of position IDs 
    set1 <- df1$pos_id
    set2 <- df2$pos_id
    
    # calculate Jaccard index of position IDs
    length(intersect(set1, set2)) / length(union(set1, set2))
  }
  )
  }
```

### Read in file paths and files

``` r
# find the root-level repo directory so results files can be accessed
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

## parent directories ##
exploration_dir <- file.path(repo_root, "exploration")
method_tests_dir <- file.path(exploration_dir, "merged_shiba_target_pilot")
splice_event_results_dir <- file.path(method_tests_dir, "shiba_gtf_to_event_pilot_tests")
target_pilot_dir <- file.path(exploration_dir, "merging-psi-tables", "target_pilot")
# merged vs. combined TARGET pilot result dirs
combined_results_dir <- file.path(target_pilot_dir, "shiba_combined")
merged_results_dir <- file.path(repo_root, "results", "merged_shiba", "target_pilot")

# negative control event directories
one_thread_events_rep1 <- file.path(splice_event_results_dir, "shiba_pilot_gtf_events")
one_thread_events_rep2 <- file.path(splice_event_results_dir, "shiba_pilot_gtf_events_rep2")

# multithreading test event directories
multithread_events_rep1 <- file.path(splice_event_results_dir, "15threads_events_rep1")
multithread_events_rep2 <- file.path(splice_event_results_dir, "15threads_events_rep2")

# scrambled merge order test event directory
scrambled_merge_order_events <- file.path(splice_event_results_dir, "15_threads_scrambled")

# merged vs. combined method shiba run events directories
combined_events_dir <- file.path(combined_results_dir, "events")
merged_events_dir <- file.path(merged_results_dir, "events")

### Files ###

# list of events coordinates files to compare between the merged and combined methods
event_files <- c(
  se = "EVENT_SE.txt",
  afe = "EVENT_AFE.txt",
  ale = "EVENT_ALE.txt",
  five = "EVENT_FIVE.txt",
  three = "EVENT_THREE.txt",
  mse = "EVENT_MSE.txt",
  mxe = "EVENT_MXE.txt",
  ri = "EVENT_RI.txt"
)

## construct events coordinate file paths for event files to read in ##
# negative control
neg_ctrl_rep1_paths <- file.path(one_thread_events_rep1, event_files)
names(neg_ctrl_rep1_paths) <- names(event_files)

neg_ctrl_rep2_paths <- file.path(one_thread_events_rep2, event_files)
names(neg_ctrl_rep2_paths) <- names(event_files)

# multithreading
multithread_rep1_paths <- file.path(multithread_events_rep1, event_files)
names(multithread_rep1_paths) <- names(event_files)

multithread_rep2_paths <- file.path(multithread_events_rep2, event_files)
names(multithread_rep2_paths) <- names(event_files)

# scrambled merge order
scrambled_merge_order_paths <- file.path(scrambled_merge_order_events, event_files)
names(scrambled_merge_order_paths) <- names(event_files)

# merged vs. combined method 
merged_events_paths <-file.path(merged_events_dir, event_files)
names(merged_events_paths) <- names(event_files)

combined_events_paths <-file.path(combined_events_dir, event_files)
names(combined_events_paths) <- names(event_files)
```

Read in files

``` r
# read in each list of files into a list of dataframes

# negative control 
neg_ctrl_rep1_list <- read_in_events(neg_ctrl_rep1_paths)
neg_ctrl_rep2_list <- read_in_events(neg_ctrl_rep2_paths)

# multithreading
multithread_rep1_list <- read_in_events(multithread_rep1_paths)
multithread_rep2_list <- read_in_events(multithread_rep2_paths)

# scrambled gtf merge prder
scrambled_list <- read_in_events(scrambled_merge_order_paths)

# merged vs. combined method
merged_events_list <- read_in_events(merged_events_paths)
combined_events_list <- read_in_events(combined_events_paths)
```

## Quantify similarity of position IDs in each event coordinate file pair

Calculate Jaccard similarity index for each set of position IDs for each
splice event type identified from each Shiba run method

``` r
# negative control comparison
jaccard_indices_neg_ctrl <- calculate_jaccard_indices_of_events(
  neg_ctrl_rep1_list,
  neg_ctrl_rep2_list)

# multithreading comparison
jaccard_indices_multithreading <- calculate_jaccard_indices_of_events(
  multithread_rep1_list,
  multithread_rep2_list
)

# scrambled GTF order comparison (same threads)
jaccard_indices_scrambled <- calculate_jaccard_indices_of_events(
  multithread_rep1_list,
  scrambled_list
)

# merged vs. combined method
jaccard_indices_merged_vs_combined <- calculate_jaccard_indices_of_events(
  merged_events_list,
  combined_events_list)
```

Print results

``` r
# convert list of jaccard indices into a dataframe
jaccard_df <- data.frame(
  event_type = names(jaccard_indices_merged_vs_combined),
  neg_ctrl = unlist(jaccard_indices_neg_ctrl, use.names = FALSE),
  same_threads = unlist(jaccard_indices_multithreading, use.names = FALSE),
  scrambled = unlist(jaccard_indices_scrambled, use.names = FALSE),
  merged_vs_combined = unlist(jaccard_indices_merged_vs_combined, use.names = FALSE)
)

jaccard_df
```

<div id="tbl-jaccard">

Table 1: Jaccard indices for comparisons among methods, separated by
event type.

<div class="cell-output-display">

| event_type | neg_ctrl | same_threads | scrambled | merged_vs_combined |
|:-----------|---------:|-------------:|----------:|-------------------:|
| se         |        1 |            1 |         1 |          0.9480923 |
| afe        |        1 |            1 |         1 |          0.8393962 |
| ale        |        1 |            1 |         1 |          0.8365145 |
| five       |        1 |            1 |         1 |          0.8570423 |
| three      |        1 |            1 |         1 |          0.8871048 |
| mse        |        1 |            1 |         1 |          0.9120187 |
| mxe        |        1 |            1 |         1 |          0.8599684 |
| ri         |        1 |            1 |         1 |          0.7510809 |

</div>

</div>

## Conclusions

- Retained intron position IDs are the most different between the events
  files produced by the merged and combined methods, compared to other
  splice event types.

- Differences still exist in the other event types in the merged
  vs. combined method comparison, where similarity scores range from
  0.83-0.94 and likely stem from GTF differences from the different
  Shiba runs.

- Events coordinates derived from GTF files made without multithreading
  with the same merge order, with multithreading and the same merge
  order, and with different merge orders are the same.

- Remaining sources of differences between the methods are: Merging in
  the reference GTF multiple times due to GTFs being generated for each
  sample by running the full `bam2gtf.py` script on one sample at a
  time, or another unknown difference between the merged and combined
  methods.
