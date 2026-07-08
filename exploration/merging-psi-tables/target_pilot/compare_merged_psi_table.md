# Compare merged PSI table from TARGET pilot samples
Cindy Liang (celiang@ucsc.edu)
2026-07-08

**Question:** Does merging separate splice tables cause us to lose out
on the trustworthiness of unannotated events to an extent that it
impacts many samples?

As part of our pipeline, we plan to merge Shiba tables created for
individual samples to obtain a final PSI table of all samples. Before
doing so, we use this notebook to test that merging PSI tables from
separate runs will not introduce a large amount of untrustworthy splice
events.

As part of tracking differences between the two tables, we have labeled
NA values introduced at various steps of the data processing pipeline
with different negative values:

- -1: NA values assigned by Shiba (splice events with low coverage)

- -2: NA values from merging the separate tables together (samples have
  a NA value after this step if there are genes dropped by Shiba for
  those samples because there is only one transcript detected in the GTF
  for that gene)

- NA: NA values from merging the separate and combined tables. These NAs
  represent splice events that are only found in one analysis method
  (combined vs. separate)

We are particularly interested in the impact merging separate tables has
on unannotated event detection, since this is one novelty of using
Shiba.

**Targeted questions this notebook is trying to answer:**

- What splice event types do we miss out on with the separate tables
  method?
- Of splice events with high numbers of NA values, how many of these
  events only have a numeric PSI value in one (or a very low amount of)
  samples?
- How many NA values of different NA types are there in each type of
  splice table? How are they distributed across events? I am more
  concerned with the merge NAs (-2), as there is no good way to replace
  them in the separate tables (they could either be 0 or 1). Another
  concern is events that are only found in the combined table, or only
  found in the separate table, as there is no easy way to correct for
  those.
- Can the number of NA values be reduced in the separate tables method
  if we filter for events with numeric PSI values in a minimum number of
  samples?
- What is the relationship between PSI values and NAs in each splice
  table type?
- What is the relationship between PSI values and NAs in each splice
  table type after filtering for events with numeric PSI values in a
  minimum number of samples? Assuming the Shiba NAs (-1) are from low
  read support, do the separate and combined tables always call the same
  events as a Shiba NA? Merge NAs (-2) should arise from events dropped
  by Shiba due to lack of alternative transcripts for a gene, so we
  expect these NAs to be 0 or 1.
- How well correlated are the number of NA values in each event type for
  each PSI table type?
- How well correlated are the PSI values of events that are shared?

## Set up

### Define functions

``` r
# PSI distribution plots of events in each method
psi_distributions <- function(
    df,
    event_type_value) {
  df |>
    dplyr::filter(event_type == event_type_value) |>
    ggplot(aes(PSI, fill = event_status)) +
    geom_histogram(bins = 20) +
    facet_wrap(vars(event_type, label, event_status), scales = "free_y") +
    plot_theme +
    # make facet labels bigger
    theme(strip.text.x = element_text(size = global_size),
          # rotate x axis labels so they don't overlap
          axis.text.x = element_text(angle =45))
}

# create summary table of events detected
make_event_summary <- function (events_table) {
  events_table |>
    dplyr::summarise(
      .by = c(event_type, label),
      # count number of events in each event type and annotation category
      total = dplyr::n(),
      shared_count = sum(shared_event),
      combined_only_count = sum(combined_event) - shared_count,
      separate_only_count = sum(separate_event) - shared_count,
      shared_percent = shared_count / total * 100,
      combined_only_percent = combined_only_count / total * 100,
      separate_only_percent = separate_only_count / total * 100,
    )
}

# plot event summary frequencies as percent and raw counts bar plots
plot_event_summary <- function(
    summary_df) {
  # pivot summary df longer for plotting
  long_summary_df <- summary_df |>
    tidyr::pivot_longer(
      cols = c(
        shared_count,
        shared_percent,
        combined_only_count,
        combined_only_percent,
        separate_only_count,
        separate_only_percent
      ),
      # extract percents and counts into separate columns
      names_to = c("category", ".value"),
      names_pattern = "(.+)_(percent|count)"
    )

  # create percent stacked barplots
  percent_plot <-
    ggplot(long_summary_df, aes(fill = category, x = event_type, y = percent)) +
    geom_bar(position = "stack", stat = "identity") +
    plot_theme +
    labs(
      title = "% of splice events in combined, separate, or shared groups",
      x = "Annotation status of event",
      y = "Percentage of events"
    ) +
    facet_wrap(vars(label)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    plot_theme

  # create raw counts grouped barplot
  counts_plot <-
    ggplot(long_summary_df, aes(fill = category, x = event_type, y = count)) +
    geom_bar(position = "dodge", stat = "identity") +
    plot_theme +
    labs(
      title = "Number of splice events only in combined, separate, or shared groups",
      x = "Annotation status of event",,
      y = "Number of events"
    ) +
    facet_wrap(vars(label)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    plot_theme

  # arrange plots for printing
  percent_plot /
    counts_plot

}

# Histogram of quantified events
number_events_dist <- function(
    df,
    psi_type, title) {
    # construct histogram
    ggplot(df,  aes(x = .data[[psi_type]], fill = event_type)) +
    geom_histogram(binwidth = 1) +
    facet_wrap(vars(label)) +
    scale_fill_manual(values = cbPalette) +
    plot_theme +
    # make facet labels bigger
    theme(strip.text.x = element_text(size = global_size),
          # rotate x axis labels so they don't overlap
          axis.text.x = element_text(angle =45)) +
    ggtitle(title)
}

## Sample-level functions

# pivot event-level df long
pivot_long <- function (wide_event_df) {
  wide_event_df |>
    dplyr::select(starts_with("SRR"), event_type, pos_id, label) |>
    tidyr::pivot_longer(
      cols = matches("_PSI_(combined|separate)$"),
      names_to = c("sample", ".value"),
      names_pattern = "(.*)_PSI_(combined|separate)"
  )}

# event summary plots of frequency of events in each category
event_summary <- function(
    long_df,
    min_sample_value) {
  # filter samples by complete PSI value in min number of samples
  all_events_n_filtered <- long_df |>
    dplyr::group_by(pos_id) |>
    # count number of events are present in each shared status
    dplyr::mutate(
      combined_count = sum(method == "combined"),
      separate_count = sum(method == "separate")
    ) |>
    # filter for events with complete PSI values for minimum number of samples
    dplyr::filter(
      combined_count >= min_sample_value,
      separate_count >= min_sample_value
    ) |>
    dplyr::ungroup() |>
    # create summary table for plotting
    dplyr::summarise(
      .by = c(event_type, label),
      # count number of events in each event type and annotation category
      total = dplyr::n(),
      shared_count = sum(shared_event),
      combined_only_count = sum(combined_event) - shared_count,
      separate_only_count = sum(separate_event) - shared_count,
      shared_percent = shared_count / total * 100,
      combined_only_percent = combined_only_count / total * 100,
      separate_only_percent = separate_only_count / total * 100
    )
}

# Create summary table of number and percentage of events in each PSI-NA matchup category
na_comparison_summary <- function(long_df) {
  long_df |>
    dplyr::mutate(
      match_category =
        dplyr::case_when(
          separate == -1 & combined == -1 ~ "-1 in both",
          separate == -1 & is.na(combined) ~ "NA in combined, -1 in separate",
          separate == -2 & combined == -1 ~ "-2 in separate, -1 in combined",
          separate == -2 & combined >= 0 ~ "PSI quantified in combined, -2 in separate table",
          is.na(separate) & combined >= 0 ~ "NA in separate, quantified in combined",
          is.na(separate) & combined == -1 ~ "NA in separate, -1 in combined",
          is.na(separate) & combined == -2 ~ "NA in separate, -2 in combined",
          separate >= 0 & is.na(combined) ~ "NA in combined, quantified in separate",
          separate == -2 & is.na(combined) ~ "NA in combined, -2 in separate",
          separate >= 0 & combined >= 0 ~ "PSI quantified in both"
        )
    ) |>
    dplyr::summarise(
      .by = c(match_category),
      count = dplyr::n()
    ) |>
    dplyr::mutate(
      total = sum(count),
      percent = count / total * 100
    )
}

# plot PSI distributions of events with numeric PSI values in combined tables but are dropped in separate tables
psi_values_lost_in_tables <- function(na_comparison_df, match_queried) {
  na_comparison_df |>
    # drop NAs from putting together the separate and combined tables
    tidyr::drop_na() |>
    # assign match categories
    dplyr::mutate(
      match_category =
        dplyr::case_when(
          separate == -1 & combined == -1 ~ "both Shiba NAs",
          separate == -2 & combined == -1 ~ "separate dropped, combined shiba NA",
          separate >= 0 & combined >= 0 ~ "PSI quantified in both",
          separate == -2 & combined >= 0 ~ "PSI quantified in combined, dropped in separate table",
          is.na(separate) & combined >= 0 ~ "PSI only quantified in separate table",
          separate >= 0 & is.na(combined) ~ "PSI only quantified in combined table"
        )
    ) |>
    # plot only events in combined table with numeric PSIs that are dropped in the separate table
    dplyr::filter(match_category == match_queried,) |>
      ggplot(aes(combined, fill = event_type)) +
      geom_histogram(bins = 20) +
      facet_wrap(vars(label)) +
    scale_fill_manual(values = cbPalette) +
      plot_theme +
      # make facet labels bigger
      theme(strip.text.x = element_text(size = global_size),
            # rotate x axis labels so they don't overlap
            axis.text.x = element_text(angle =45))
}
```

## Directories and files

``` r
## directories ##
# find the root-level repo directory so we can access the other files
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
# define the data directories
exploration_dir <- file.path(repo_root, "exploration")
merge_exploration_dir <- file.path(exploration_dir, "merging-psi-tables")
target_pilot_dir <- file.path(merge_exploration_dir, "target_pilot")
#output dir for long tables
target_pilot_output <- file.path(exploration_dir, "psi_tables")
# directory of shiba results produced from the same shiba run
combined_dir <- file.path(target_pilot_dir, "shiba_combined", "results")
# psi table directory
combined_splice_results_dir <- file.path(combined_dir, "splicing")
# Directory of PSI table, merged from separate shiba runs
separate_psi_table_dir <- file.path(target_pilot_dir, "merged_results")

# Check output dir exists; if not, create it
if (!dir.exists(target_pilot_output)) {
  dir.create(target_pilot_output)
}

## files ##
# merged PSI file of two samples obtained from merge-shiba-psi-tables.qmd
separate_psi_file <- file.path(separate_psi_table_dir, "merged_psi_table.tsv")
# define list of PSI results files corresponding to event types quantified by Shiba bulk analysis
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

## output ##
# long df output file of combined and separate PSI tables
long_df_output <- file.path(target_pilot_output, "long_psi_df.rds")
```

Read in files

``` r
# read merged splice table from separate runs of Shiba on the TARGET pilot samples
separate_psi_table <- readr::read_tsv(separate_psi_file, col_types=readr::cols(.default = "c")) |>
  # replace remaining NAs from this table representing dropped genes from one sample with no detected splicing for that gene with -2
  dplyr::mutate(across(contains("_PSI"), \(x) tidyr::replace_na(as.numeric(x), -2)))

# read in combined splice results to compare against separate results
combined_splice_results <- shiba_psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
      # convert PSI values to numeric and replace NA values from shiba representing low coverage with -1
      dplyr::mutate(across(contains("_PSI"), \(x) tidyr::replace_na(as.numeric(x), -1)))
  })

# combine PSI values of samples run together into one dataframe to compare against separate PSI dataframe
combined_psi_table <- purrr::list_rbind(combined_splice_results, names_to = "event_type")
```

### Obtain dimensions of each PSI table

``` r
# dimensions of separate PSI table
dim(separate_psi_table)
```

    [1] 818412     92

``` r
# dimensions of combined PSI table
dim(combined_psi_table)
```

    [1] 746601     92

There are 71,811 more splice events in the separate PSI table than the
combined table.

## What splice event types do we miss out on with the separate tables method?

To help us prioritize what splice event types we may be the most
confident in after merging separate Shiba PSI tables, we are interested
in seeing a breakdown of what event types are most represented in PSI
tables constructed from different methods

``` r
# create df of splice events that are shared between both combined and separate tables
# This is an event-level dataframe (each row is one unique pos_id/event_id)
# columns are each samples' PSI value (with -1/-2/NA values indicating different NA categories)
# combined_event, separate_event, and shared_event columns label what tables the event is present in
all_events <- dplyr::full_join(
  combined_psi_table,
  separate_psi_table,
  by = c("event_type", "pos_id", "gene_id", "label"),
  # label PSI values by table they came from
  suffix = c("_combined", "_separate")) |>
  # categorize events by whether they are in the combined or separate tables
  dplyr::mutate(
    # combined events are counted if the values are not all NA
    combined_event = ! dplyr::if_all(ends_with("_combined"), is.na),
    # separate events are counted if the values for the event are not all NA (indicating the specific position ID is only found in the combined or separate table)
    separate_event = ! dplyr::if_all(ends_with("_separate"), is.na),
    shared_event = combined_event & separate_event
  )
```

Print summary of counts across each method, across all event types

``` r
# number of total events shared or present only in the combined and separate tables
# this table counts events with -1 and -2 PSI values

all_methods_summary <- all_events |> dplyr::summarise(
    # count number of events in each event type and annotation category
    shared_count = sum(shared_event),
    combined_only_count = sum(combined_event) - shared_count,
    separate_only_count = sum(separate_event) - shared_count
)

all_methods_summary
```

| shared_count | combined_only_count | separate_only_count |
|-------------:|--------------------:|--------------------:|
|       667022 |               79579 |              151390 |

Although the majority of splice events (inclusive of NA values) are
shared between methods, there are about twice as many splice events only
in the separate method’s table than the combined method.

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
    separate_only_count = sum(separate_event) - shared_count,
    shared_percent = shared_count / total * 100,
    combined_only_percent = combined_only_count / total * 100,
    separate_only_percent = separate_only_count / total * 100,
)

all_events_summary |> dplyr::select(
  total,
  shared_percent,
  combined_only_percent,
  separate_only_percent)
```

|  total | shared_percent | combined_only_percent | separate_only_percent |
|-------:|---------------:|----------------------:|----------------------:|
| 104872 |       97.13556 |             1.0050347 |              1.859410 |
|  32912 |       65.35914 |            20.4697375 |             14.171123 |
|  97616 |       42.63133 |            22.9767661 |             34.391903 |
| 180632 |       85.74339 |             2.0887772 |             12.167833 |
| 156353 |       81.48740 |             1.7306991 |             16.781897 |
|  62334 |       35.55042 |            21.6382712 |             42.811307 |
|  35556 |       85.59737 |             5.3268084 |              9.075824 |
|  16825 |       50.06835 |            30.0980684 |             19.833581 |
|  40691 |       89.44730 |             3.5929321 |              6.959770 |
|  17214 |       53.20669 |            28.6220518 |             18.171256 |
|  64124 |       96.40540 |             1.0183395 |              2.576258 |
|  19434 |       58.05804 |            26.1397551 |             15.802202 |
|    963 |       82.24299 |             0.7268951 |             17.030114 |
|    481 |       22.86902 |            24.1164241 |             53.014553 |
|  51635 |       50.57229 |            19.4267454 |             30.000968 |
|  16349 |       79.80916 |             0.9664200 |             19.224417 |

Out of all the event types, annotated skipped exons and annotated
multiple skipped exon types are the most shared between the two methods
(above 90%). In contrast, unannotated events of these types are only
shared at 65% and 58% respectively.

Plot event-level (one count per unique event) summary of event frequency
across methods and event types

``` r
plot_event_summary(all_events_summary)
```

<div id="fig-splice_events_breakdown">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-splice_events_breakdown-1.png"
id="fig-splice_events_breakdown" />

Figure 1

</div>

The majority of annotated event types are shared between both tables,
although there are more splice events unique to the separate table
(green) than the combined table (red). Inflated values in the separate
tables may be due to -2 NAs which arise in the Shiba processing of
separate samples when when genes without multiple transcripts are
dropped.

### How many events remain if we filter for events with PSI values in a minimum number of samples?

- Of splice events with high numbers of NA values, how many of these
  events only have a numeric PSI value in one (or a very low amount of)
  samples?

Add summary columns to dataframe

``` r
all_events <- all_events |>
  dplyr::mutate(
    ## count PSI value types in separate table ##
    # count the number of PSI values that are not one of our NA types
    n_quantified_separate = rowSums(
      # select sample PSI values that are not one of our NA types
      dplyr::pick(dplyr::ends_with("_PSI_separate")) >= 0,
      # ignore NA values (events only in separate or combined tables)
      na.rm = TRUE),
    # count the number of PSI values dropped due to insufficient reads
    n_shiba_na_separate = rowSums(
      # select sample PSI values that have -1 NA type
      dplyr::pick(dplyr::ends_with("_PSI_separate")) == -1,
      na.rm = TRUE),
    # count number of PSI values dropped due to lack of transcript diversity in single samples
    n_dropped_separate = rowSums(
      # elect sample PSI values that have -2 NA type
      dplyr::pick(dplyr::ends_with("_PSI_separate")) == -2
    ),
    ## count PSI value types in combined table ##
    # count the number of PSI values that are not one of our NA types
    n_quantified_combined = rowSums(
      # select sample PSI values that are not one of our NA types
      dplyr::pick(dplyr::ends_with("_PSI_combined")) >= 0,
      # ignore NA values (events only in separate or combined tables)
      na.rm = TRUE),
    # count the number of PSI values dropped due to insufficient reads
    n_shiba_na_combined = rowSums(
      # select sample PSI values that have -1 NA type
      dplyr::pick(dplyr::ends_with("_PSI_combined")) == -1,
      na.rm = TRUE),
    # count number of PSI values dropped due to lack of transcript diversity in single samples
    n_dropped_combined = rowSums(
      # elect sample PSI values that have -2 NA type
      dplyr::pick(dplyr::ends_with("_PSI_combined")) == -2
    )
  )
```

Filter for minimum samples with numeric PSI values in the separate table
and check n_quantified. The resulting summaries of this filtered table
will be incomplete, as it will not contain information of events in the
combined table that are present in min_samples but are missed completely
in the separate tables.

``` r
# event-level summary of the number of each splice event type in each splice table
all_events_min_samples <- all_events |>
  # filter for events where there are numeric PSI values in at least 5 samples
  dplyr::filter(n_quantified_separate >= min_samples)

# check what summary columns look like
all_events_min_samples  |>
  dplyr::select(pos_id, n_quantified_separate, n_shiba_na_separate, n_dropped_separate,
                n_quantified_combined, n_shiba_na_combined, n_dropped_combined) |>
  dplyr::arrange((n_quantified_separate)) |>
  head()
```

| pos_id | n_quantified_separate | n_shiba_na_separate | n_dropped_separate | n_quantified_combined | n_shiba_na_combined | n_dropped_combined |
|:---|---:|---:|---:|---:|---:|---:|
| SE@GL000008.2@197587-197618@194536-198476 | 5 | 23 | 60 | 5 | 83 | 0 |
| SE@GL000008.2@88636-88695@85625-129985 | 5 | 23 | 60 | 5 | 83 | 0 |
| SE@GL000195.1@143414-143629@142240-172694 | 5 | 46 | 37 | 5 | 83 | 0 |
| SE@GL000195.1@143515-143629@142240-172694 | 5 | 45 | 38 | 5 | 83 | 0 |
| SE@GL000195.1@149032-149164@142240-172694 | 5 | 45 | 38 | 5 | 83 | 0 |
| SE@GL000195.1@151955-152034@142240-172694 | 5 | 53 | 30 | 5 | 83 | 0 |

`n_quantified_separate` and `n_quantified_combined` columns have min
values of 5, so filtering appears to be correct.

Plot event distributions of the filtered table

``` r
# event-level summary of the number of each splice event type in each splice table
all_events_summary_min_samples <- all_events_min_samples |>
  # this data frame still may include -1 and -2 NA values for other splice events but the sample-level analysis later will tell us more
  dplyr::summarise(
    .by = c(event_type, label),
    # count number of events in each event type and annotation category
    total = dplyr::n(),
    shared_count = sum(shared_event),
    combined_only_count = sum(combined_event) - shared_count,
    separate_only_count = sum(separate_event) - shared_count,
    shared_percent = shared_count / total * 100,
    combined_only_percent = combined_only_count / total * 100,
    separate_only_percent = separate_only_count / total * 100,
)

plot_event_summary(all_events_summary_min_samples)
```

<div id="fig-splice_events_breakdown_5_samples">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-splice_events_breakdown_5_samples-1.png"
id="fig-splice_events_breakdown_5_samples" />

Figure 2

</div>

Observations:

- Although most events are shared after applying this filter and we
  still have unannotated events, we also still have events only found in
  the separate table. These events are harder to deal with because they
  may be real events mislabeled by Shiba, or they may be unreal events
  (we cannot easily tell)

- Additionally, we still see that most unannotated events are no longer
  present after this filtering.

Print table of splice event counts of filtered table

``` r
all_events_summary_min_samples |>
  dplyr::select(event_type, label, shared_count, separate_only_count, combined_only_count)
```

| event_type | label       | shared_count | separate_only_count | combined_only_count |
|:-----------|:------------|-------------:|--------------------:|--------------------:|
| se         | annotated   |        55020 |                 813 |                   0 |
| se         | unannotated |         2820 |                 312 |                   0 |
| afe        | annotated   |        70759 |               10601 |                   0 |
| afe        | unannotated |         2148 |                1566 |                   0 |
| ale        | annotated   |        43524 |               13828 |                   0 |
| ale        | unannotated |         1123 |                1289 |                   0 |
| five       | annotated   |        13783 |                 423 |                   0 |
| five       | unannotated |          756 |                 145 |                   0 |
| three      | annotated   |        18014 |                 472 |                   0 |
| three      | unannotated |          941 |                 140 |                   0 |
| mse        | annotated   |        38267 |                 229 |                   0 |
| mse        | unannotated |         1064 |                 123 |                   0 |
| mxe        | annotated   |          346 |                 115 |                   0 |
| mxe        | unannotated |            8 |                   9 |                   0 |
| ri         | annotated   |         9585 |                1675 |                   0 |
| ri         | unannotated |         6692 |                1139 |                   0 |

Observations:

- Once the “splice event with numeric value in at least 5 samples”
  filter is applied, there are no longer splice events only in the
  combined method’s PSI table.

- The combined method is expected capture more events found in fewer
  samples, due to calculating splice events usage defined from a more
  complete set of events.

## How many of each splice event types are quantified?

This section plots the distributions of splice events quantified in each
method.

### Histogram of quantified events

``` r
number_events_dist(all_events, "n_quantified_separate", "Number of quantified PSI values in separate table") /

number_events_dist(all_events, "n_quantified_combined", "Number of quantified PSI values in combined table") + plot_layout(guides = "collect")
```

<div id="fig-num_quantified_psi">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-num_quantified_psi-1.png"
id="fig-num_quantified_psi" />

Figure 3

</div>

Check skipped exon quantified PSI histograms

``` r
all_events |>
  dplyr::filter(event_type == "se") |>
  number_events_dist("n_quantified_separate", "Number of quantified PSI values in separate table") /

all_events |>
  dplyr::filter(event_type == "se") |>
  number_events_dist("n_quantified_combined", "Number of quantified PSI values in combined table") + plot_layout(guides = "collect")
```

<div id="fig-num_quantified_psi_se">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-num_quantified_psi_se-1.png"
id="fig-num_quantified_psi_se" />

Figure 4

</div>

As expected, we have more observations of quantified unannotated PSI
values in the combined table (higher bars at the end of the X axis) in
both sets of histograms. More critically, the unannotated observations
in the combined table that we lose in the separate tables seem to be
quantified in most of the TARGET pilot samples.

### Histogram of events dropped by Shiba due to insufficient read counts

The following histograms show the counts of PSI values that are DROPPED
in each method.

``` r
number_events_dist(all_events, "n_shiba_na_separate", "Number of PSIs dropped due to insufficient read counts in separate table") /

number_events_dist(all_events, "n_shiba_na_combined", "Number of PSIs dropped due to insufficient read counts in combined table") + plot_layout(guides = "collect")
```

<div id="fig-num_low_read_na">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-num_low_read_na-1.png"
id="fig-num_low_read_na" />

Figure 5

</div>

There appear to be more events dropped across more samples in the
combined results compared to the separate results. However, I think this
is because the combined method captures more “rare” (e.g. present in
only one sample) splice events. This would increase the overall pool of
possible loci that could be sequenced so sparsely that they do not
actually make it into PSI calculation.

### Histogram of events dropped due to genes being thrown out in separate method

I do not plot this for the combined table because they are all 0 (we
only define what is dropped due to insufficient transcript
representation based off what transcripts are missing when individual
separate tables are merged).

``` r
number_events_dist(all_events, "n_dropped_separate", "# PSIs dropped due to gene being thrown out in separate table")
```

    Warning: Removed 79579 rows containing non-finite outside the scale range
    (`stat_bin()`).

<div id="fig-num_low_transcript_na">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-num_low_transcript_na-1.png"
id="fig-num_low_transcript_na" />

Figure 6

</div>

This category of missing values are inflated in unannotated events in
the separate table.

## How many of each splice event types are quantified?

## Examine sample-level splice events

### Make long df of events

Pivot longer for plotting sample-level splice event info

Check if result exists in cache

``` r
use_cached_long <- file.exists(long_df_output) && params$use_cache
```

``` r
# Pivot all events df longer
# Each row corresponds to a splice event found in one sample, so there are multiple rows of splice events with the same sample ID
if (use_cached_long) {
  long_all_events <- readr::read_rds(long_df_output)
} else {
long_all_events <- pivot_long(all_events)
}

# print column names
colnames(long_all_events)
```

    [1] "event_type" "pos_id"     "label"      "sample"     "combined"  
    [6] "separate"  

### Write output

Write long table into output

``` r
# long df
saveRDS(long_all_events, file = long_df_output)
```

### Examine relationships between of PSI and NA values of combined and separate tables

#### Are PSI values of events that are shared between each method identical?

``` r
long_psi_in_both_df <- long_all_events |> dplyr::filter(combined > 0 & separate > 0)
identical(long_psi_in_both_df$combined, long_psi_in_both_df$combined)
```

    [1] TRUE

Yes, PSI values of shared events in both methods are identical.

#### How well correlated are the PSI values and NAs in each splice table type?

Plot the different types of matchup categories

``` r
# plot basic scatterplot to see the different PSI value / NA matchup categories in the data
ggplot(long_all_events, aes(x = combined, y = separate)) +
  geom_point()
```

    Warning: Removed 20325272 rows containing missing values or values outside the scale
    range (`geom_point()`).

![](compare_merged_psi_table_files/figure-commonmark/psi_na_matchup_categories-1.png)

From this plot, I see four different categories of splice event matches
between the separate and combined tables:

- Upper dot: We have some splice events that are dropped by Shiba due to
  no alternative splicing in the separate tables (Y = -2) but are also
  dropped in the combined tables because there were too few reads anyway
  (X = -1).
- Lower dot: We have some splice events that had too few reads for a
  numeric PSI in both tables (x = -1, y = -1)
- Diagonal line: We also see events whose PSI values are numeric and
  have a linear relationship in both combined and separate methods
  (these are likley the events that were identical and found in both)
- The most concerning case is the horizontal line at Y=-2, corresponding
  to splice events dropped in the separate tables due to lack of
  transcripts for a gene but was captured in the combined splice table.
  In this section, we can see that the PSI values calculated in the
  combined table span the full PSI value range from 0-1, so there is no
  easy way to replace the -2 NA values in the separate table.

Summarize how many of each of these match categories are present in the
pilot tables

``` r
important_categories <- c("PSI quantified in both", "NA in separate, quantified in combined", "NA in combined, quantified in separate", "PSI quantified in combined, -2 in separate table")

na_comparison_summary(long_all_events) |>
  # arrange most important categories to the front
  dplyr::arrange(match(match_category, important_categories))
```

| match_category | count | total | percent |
|:---|---:|---:|---:|
| PSI quantified in both | 12670697 | 79023208 | 16.034147 |
| NA in separate, quantified in combined | 3307421 | 79023208 | 4.185379 |
| NA in combined, quantified in separate | 1546935 | 79023208 | 1.957571 |
| PSI quantified in combined, -2 in separate table | 7205168 | 79023208 | 9.117787 |
| -2 in separate, -1 in combined | 27666041 | 79023208 | 35.010020 |
| -1 in both | 11156030 | 79023208 | 14.117410 |
| NA in separate, -1 in combined | 3695531 | 79023208 | 4.676514 |
| NA in combined, -2 in separate | 10942532 | 79023208 | 13.847238 |
| NA in combined, -1 in separate | 832853 | 79023208 | 1.053935 |

We care most about “PSI quantified in both”, “NA in separate, quantified
in combined”, “NA in combined, quantified in separate”, and “PSI
quantified in combined, -2 in separate table” categories. \* Events we
expect are real and are captured in both methods: PSI quantified in both
\* Events we expect are real and are missed or miscategorized in the
separate method: NA in separate, quantified in combined \* Events that
may be real and miscategorized (but are difficult to confirm): NA in
combined, quantified in separate \* Events we expect are real but are
missed in the separate method due to how events are defined in the GTF:
PSI quantified in combined, -2 in separate table

Summarize how many of each of these match categories are present in the
pilot tables, filtered for events with numeric PSI in min_samples

### Make filtered df long

``` r
long_min_sample_all_events <- pivot_long(all_events_min_samples)

# create summary table of match categories
na_comparison_summary(long_min_sample_all_events) |>
  # arrange most important categories to the front
  dplyr::arrange(match(match_category, important_categories))
```

| match_category | count | total | percent |
|:---|---:|---:|---:|
| PSI quantified in both | 12363779 | 26200152 | 47.189722 |
| NA in combined, quantified in separate | 1421478 | 26200152 | 5.425457 |
| PSI quantified in combined, -2 in separate table | 884696 | 26200152 | 3.376683 |
| -1 in both | 6284307 | 26200152 | 23.985765 |
| -2 in separate, -1 in combined | 3774018 | 26200152 | 14.404565 |
| NA in combined, -2 in separate | 909789 | 26200152 | 3.472457 |
| NA in combined, -1 in separate | 562085 | 26200152 | 2.145350 |

In both tables, the majority of NAs are from low counts. There is a
higher percentage of numeric PSIs that we lose out on in the separate
table, even with the min_samples filter.

#### Plot PSI distribution of PSI values lost in each method

This plot shows the distribution of the “true” values of -2 (“gene
dropped”) NAs in the combined table.

##### Quantified PSI values dropped in separate tables due to lack of transcript diversity

``` r
# plot PSI distributions of events in combined dataframe
psi_values_lost_in_tables(long_all_events, "PSI quantified in combined, dropped in separate table")
```

![](compare_merged_psi_table_files/figure-commonmark/psi_dist_of_combined_only_events-1.png)

The PSI values of the “gene dropped” NAs span the full range of PSI
values, so we miss out entirely on a wide range of values using the
separate method.

Following this analysis, we concluded that the Shiba splice compendium
workflow needed to be modified to allow for the “gene dropped” NA values
to be quantified.
