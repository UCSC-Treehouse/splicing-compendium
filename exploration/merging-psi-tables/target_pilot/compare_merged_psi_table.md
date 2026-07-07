# Compare merged PSI table from TARGET pilot samples
Cindy Liang (celiang@ucsc.edu)
2026-07-07

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
on unannotated event detection, since this is why we chose to use Shiba
over other tools.

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
      facet_wrap(vars(event_type, label)) +
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

Print summary of counts and percentages across each method, broken down
by event type

``` r
# event-level summary of the number of each splice event type in each splice table
all_events_summary <- make_event_summary(all_events)
all_events_summary
```

| event_type | label | total | shared_count | combined_only_count | separate_only_count | shared_percent | combined_only_percent | separate_only_percent |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| se | annotated | 104872 | 101868 | 1054 | 1950 | 97.13556 | 1.0050347 | 1.859410 |
| se | unannotated | 32912 | 21511 | 6737 | 4664 | 65.35914 | 20.4697375 | 14.171123 |
| afe | unannotated | 97616 | 41615 | 22429 | 33572 | 42.63133 | 22.9767661 | 34.391903 |
| afe | annotated | 180632 | 154880 | 3773 | 21979 | 85.74339 | 2.0887772 | 12.167833 |
| ale | annotated | 156353 | 127408 | 2706 | 26239 | 81.48740 | 1.7306991 | 16.781897 |
| ale | unannotated | 62334 | 22160 | 13488 | 26686 | 35.55042 | 21.6382712 | 42.811307 |
| five | annotated | 35556 | 30435 | 1894 | 3227 | 85.59737 | 5.3268084 | 9.075824 |
| five | unannotated | 16825 | 8424 | 5064 | 3337 | 50.06835 | 30.0980684 | 19.833581 |
| three | annotated | 40691 | 36397 | 1462 | 2832 | 89.44730 | 3.5929321 | 6.959770 |
| three | unannotated | 17214 | 9159 | 4927 | 3128 | 53.20669 | 28.6220518 | 18.171256 |
| mse | annotated | 64124 | 61819 | 653 | 1652 | 96.40540 | 1.0183395 | 2.576258 |
| mse | unannotated | 19434 | 11283 | 5080 | 3071 | 58.05804 | 26.1397551 | 15.802202 |
| mxe | annotated | 963 | 792 | 7 | 164 | 82.24299 | 0.7268951 | 17.030114 |
| mxe | unannotated | 481 | 110 | 116 | 255 | 22.86902 | 24.1164241 | 53.014553 |
| ri | unannotated | 51635 | 26113 | 10031 | 15491 | 50.57229 | 19.4267454 | 30.000968 |
| ri | annotated | 16349 | 13048 | 158 | 3143 | 79.80916 | 0.9664200 | 19.224417 |

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

### What specific event types are impacted (missed) by the separate tables method?

These following plots answer the question for me better than the above
summary plots.

``` r
combined_only_summary <- all_events_summary |>
  dplyr::select(event_type, label, combined_only_count, total, combined_only_percent)

long_summary_df <- combined_only_summary |>
    tidyr::pivot_longer(
      cols = c(
        combined_only_count,
        combined_only_percent
      ),
      # extract percents and counts into separate columns
      names_to = c("category", ".value"),
      names_pattern = "(.+)_(percent|count)"
    )

# create raw counts grouped barplot
counts_plot <- ggplot(long_summary_df, aes(fill = event_type, x = event_type, y = count)) +
    geom_bar(stat = "identity") +
      scale_fill_manual(values = cbPalette) +
    plot_theme +
    labs(
      title = "# annotated and unannotated splice events only in combined table",
      x = "Event Type",
      y = "Number of events"
    ) +
    facet_wrap(vars(label)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme(strip.text = element_text(size = 15))

# create percents grouped barplot
percents_plot <- ggplot(long_summary_df, aes(fill = event_type, x = event_type, y = percent)) +
    geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values = cbPalette) +
    plot_theme +
    labs(
      title = "% annotated and unannotated splice events only in combined table",
      x = "Annotation status of event",
      y = "Percent of all events"
    ) +
    facet_wrap(vars(label)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme(strip.text = element_text(size = 15))


percents_plot /
  counts_plot + plot_layout(guides = "collect")
```

<div id="fig-combined_only_splice_events">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-combined_only_splice_events-1.png"
id="fig-combined_only_splice_events" />

Figure 2

</div>

I think of these as “percent of each event category missed in the
separate tables”. In the worst cases, 30% of unannotated events might be
only in the combined table (FIVE events). Losing 30% of unannotated real
events seems rough.

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

# check what summary looks like
all_events |>
  dplyr::select(pos_id, n_quantified_separate, n_shiba_na_separate, n_dropped_separate,
                n_quantified_combined, n_shiba_na_combined, n_dropped_combined) |>
  head()
```

| pos_id | n_quantified_separate | n_shiba_na_separate | n_dropped_separate | n_quantified_combined | n_shiba_na_combined | n_dropped_combined |
|:---|---:|---:|---:|---:|---:|---:|
| SE@GL000008.2@129985-130583@85625-155430 | 2 | 24 | 62 | 2 | 86 | 0 |
| SE@GL000008.2@135134-135173@85625-155430 | 1 | 19 | 68 | 1 | 87 | 0 |
| SE@GL000008.2@154869-154964@154715-156667 | 2 | 30 | 56 | 2 | 86 | 0 |
| SE@GL000008.2@155430-155531@135173-156721 | 3 | 24 | 61 | 3 | 85 | 0 |
| SE@GL000008.2@156667-156758@154715-157528 | 1 | 25 | 62 | 1 | 87 | 0 |
| SE@GL000008.2@156667-156761@154715-157528 | 1 | 28 | 59 | 1 | 87 | 0 |

Filter for minimum samples with numeric PSI values in the separate table
and check n_quantified. The resulting summaries of this filtered table
will be incomplete, as it will not contain information of events in the
combined table that are present in min_samples but are missed completely
in the separate tables.

``` r
# event-level summary of the number of each splice event type in each splice table
all_events_min_samples <- all_events |>
  # filter for events where there are numeric PSI values in at least 5 samples
  dplyr::filter(n_quantified_separate >= min_samples | n_quantified_combined >= min_samples)

# check what summary columns look like
all_events_min_samples  |>
  dplyr::select(pos_id, n_quantified_separate, n_shiba_na_separate, n_dropped_separate,
                n_quantified_combined, n_shiba_na_combined, n_dropped_combined) |>
  dplyr::arrange((n_quantified_separate)) |>
  head()
```

| pos_id | n_quantified_separate | n_shiba_na_separate | n_dropped_separate | n_quantified_combined | n_shiba_na_combined | n_dropped_combined |
|:---|---:|---:|---:|---:|---:|---:|
| SE@GL000195.1@143417-149164@142240-151955 | 0 | 1 | 87 | 11 | 77 | 0 |
| SE@GL000195.1@143417-149164@142240-154339 | 0 | 1 | 87 | 11 | 77 | 0 |
| SE@GL000195.1@143417-149164@142240-172694 | 0 | 1 | 87 | 16 | 72 | 0 |
| SE@GL000195.1@143417-149164@142240-157171 | 0 | 1 | 87 | 11 | 77 | 0 |
| SE@GL000195.1@143417-149164@142240-149745 | 0 | 1 | 87 | 13 | 75 | 0 |
| SE@GL000195.1@143417-149164@142240-164904 | 0 | 1 | 87 | 11 | 77 | 0 |

#### Plot event distributions of the filtered separate table

``` r
# event-level summary of the number of each splice event type in min_sample filtered separate table
events_summary_min_separate_samples <- all_events_min_samples |>
  # this data frame still may include -1 and -2 NA values for other splice events but the sample-level analysis later will tell us more
  # filter for events where min_samples are in separate table only
  dplyr::filter(n_quantified_separate >= min_samples) |> 
  make_event_summary()

plot_event_summary(events_summary_min_separate_samples )
```

<div id="fig-splice_events_breakdown_5_samples_separate">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-splice_events_breakdown_5_samples_separate-1.png"
id="fig-splice_events_breakdown_5_samples_separate" />

Figure 3

</div>

Print table of splice event counts of filtered combined table

``` r
events_summary_min_separate_samples |>
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

- Although most events are shared after applying this filter and we
  still have unannotated events, we also still have events only found in
  the separate table. These events are harder to deal with because they
  may be real events mislabeled by Shiba, or they may be unreal events
  (we cannot easily tell)

#### Plot event distributions of the filtered combined table

``` r
# event-level summary of the number of each splice event type in min_sample filtered combined table
events_summary_min_combined_samples <- all_events_min_samples |>
  # this data frame still may include -1 and -2 NA values for other splice events but the sample-level analysis later will tell us more
  # filter for events where min_samples are in separate table only
  dplyr::filter(n_quantified_combined >= min_samples) |> 
  make_event_summary()

plot_event_summary(events_summary_min_combined_samples)
```

<div id="fig-splice_events_breakdown_5_samples_combined">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-splice_events_breakdown_5_samples_combined-1.png"
id="fig-splice_events_breakdown_5_samples_combined" />

Figure 4

</div>

Print table of splice event counts of filtered combined table

``` r
events_summary_min_combined_samples |>
  dplyr::select(event_type, label, shared_count, separate_only_count, combined_only_count)
```

| event_type | label       | shared_count | separate_only_count | combined_only_count |
|:-----------|:------------|-------------:|--------------------:|--------------------:|
| se         | annotated   |        56860 |                   0 |                 824 |
| se         | unannotated |        18704 |                   0 |                5643 |
| afe        | unannotated |        31078 |                   0 |               17502 |
| afe        | annotated   |        76302 |                   0 |                1784 |
| ale        | annotated   |        47278 |                   0 |                 919 |
| ale        | unannotated |        17182 |                   0 |               11613 |
| five       | annotated   |        17802 |                   0 |                1039 |
| five       | unannotated |         6805 |                   0 |                3481 |
| three      | annotated   |        21386 |                   0 |                 893 |
| three      | unannotated |         7719 |                   0 |                3475 |
| mse        | annotated   |        40865 |                   0 |                 626 |
| mse        | unannotated |        10606 |                   0 |                4864 |
| mxe        | annotated   |          360 |                   0 |                   6 |
| mxe        | unannotated |           77 |                   0 |                 105 |
| ri         | unannotated |        23495 |                   0 |                9636 |
| ri         | annotated   |        10010 |                   0 |                 130 |

Observations:

- We can see that the filtered combined table results in more shared
  events that are unannotated, compared to the filtered separate table
  (by maybe 5x, looking at the Y axis scale).

- There are also a high number of combined_only events that are missed
  in the separate tables, even after filtering for combined events that
  are present in min_samples. There are ~5,000 skipped exon events only
  in the combined table that are in min_samples that are missed in the
  separate tables. That’s more than the total number of skipped exon
  events remaining in the min_samples filtered separate table.

## How many of each splice event types are quantified?

Histogram of quantified events

``` r
number_events_dist(all_events, "n_quantified_separate", "Number of quantified PSI values in separate table") /

number_events_dist(all_events, "n_quantified_combined", "Number of quantified PSI values in combined table") + plot_layout(guides = "collect")
```

<div id="fig-num_quantified_psi">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-num_quantified_psi-1.png"
id="fig-num_quantified_psi" />

Figure 5

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

Figure 6

</div>

As expected, we have more observations of quantified unannotated PSI
values in the combined table (higher bars at the end of the X axis).
More critically, the unannotated observations in the combined table that
we lose in the separate tables seem to be quantified in most of the
TARGET pilot samples.

Histogram of events dropped by Shiba due to insufficient read counts

``` r
number_events_dist(all_events, "n_shiba_na_separate", "Number of PSIs dropped due to insufficient read counts in separate table") /

number_events_dist(all_events, "n_shiba_na_combined", "Number of PSIs dropped due to insufficient read counts in combined table") + plot_layout(guides = "collect")
```

<div id="fig-num_low_read_na">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-num_low_read_na-1.png"
id="fig-num_low_read_na" />

Figure 7

</div>

These histograms are now counting events that are lost (number of -1
values). In the histogram of PSI values dropped due to low read counts
in the combined table, most events that are dropped are dropped in
almost all the samples. We see fewer events dropped due to low read
counts in general in the separate tables, but this might be because the
events that would have been dropped are missing from the separate are
also missing due to the GTF (-2 NAs).

Check PSI values dropped from insufficient read histograms for skipped
exons

``` r
all_events |>
  dplyr::filter(event_type == "se") |>
  number_events_dist("n_shiba_na_separate", "Number of PSIs dropped due to insufficient read counts in separate table") /

all_events |>
  dplyr::filter(event_type == "se") |>
  number_events_dist("n_shiba_na_combined", "Number of PSIs dropped due to insufficient read counts in combined table") + plot_layout(guides = "collect")
```

<div id="fig-num_low_read_na_se">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-num_low_read_na_se-1.png"
id="fig-num_low_read_na_se" />

Figure 8

</div>

The same trend regarding unannotated events is seen here as in
<a href="#fig-num_low_read_na" class="quarto-xref">Figure 7</a>. We also
see a higher number of annotated events dropped in all samples in the
combined tables. I suspect the reason for these two observations is
because the combined tables have a larger dictionary of possible splice
events due to the merged GTF. So all these splice events not found in
the separate tables have the potential to have -1 values.

Check histograms of events dropped in separate tables (-2 values). These
histograms plot the number of splice events dropped in splice tables
made form the separate Shiba runs method. These events should correspond
to splice events that only have one transcript annotated in the sample’s
GTF and so are dropped by Shiba.

``` r
number_events_dist(all_events, "n_dropped_separate", "# PSIs dropped due to insufficient transcripts in separate table")
```

    Warning: Removed 79579 rows containing non-finite outside the scale range
    (`stat_bin()`).

``` r
# I do not plot this for the combined table because they are all 0 (we only define what is dropped in the separate tables based off what is missing when individual separate tables are merged)
```

<div id="fig-num_separate_table_na">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-num_separate_table_na-1.png"
id="fig-num_separate_table_na" />

Figure 9

</div>

Histogram of events missing in separate tables for skipped exon events

``` r
all_events |>
  dplyr::filter(event_type == "se") |>
  number_events_dist("n_dropped_separate", "# PSIs dropped in separate table")
```

    Warning: Removed 7791 rows containing non-finite outside the scale range
    (`stat_bin()`).

<div id="fig-num_separate_table_na_se">

<img
src="compare_merged_psi_table_files/figure-commonmark/fig-num_separate_table_na_se-1.png"
id="fig-num_separate_table_na_se" />

Figure 10

</div>

Missing values are present in both annotated and unannotated events. The
annotated events that are missed in the separate method are mainly in
only a few ( n \< 5 ) samples. In contrast, almost all unannotated
events that are dropped are present in over 75/88 samples. In both
<a href="#fig-num_separate_table_na" class="quarto-xref">Figure 9</a>
and <a href="#fig-num_separate_table_na_se"
class="quarto-xref">Figure 10</a>, the events that are dropped are
absent in most samples.

Because of the extent of missing splice events in the separate method,
the conclusion from this analysis is that we need to revise the method
to allow for splice events that have only one transcript and events that
are only found in one sample to be quantified in the final workflow.
