# Compare merged PSI table vs. Shiba PSI tables
Cindy Liang (celiang@ucsc.edu)
2026-03-25

As part of our pipeline, we plan to merge Shiba tables created for
individual samples to obtain a final PSI table of all samples. Before
doing so, we use this notebook to test that merging PSI tables and
junction bedfiles from separate runs will not introduce a large amount
of untrustworthy splice events or junction counts.

## Set up

## Directories and files

``` r
# define samples list
samples <- c("SRR601500", "SRR604528")

# define the data directories
# shiba results dir
shiba_dir <- file.path("shiba_results")

# directory of shiba results produced from the same shiba run
combined_dir <- file.path(shiba_dir, "combined_run")
# psi table directory
combined_splice_results_dir <- file.path(combined_dir, "splicing")
# junctions bed file directory
combined_junctions_dir <- file.path(combined_dir, "junctions")
# junctions bed file
combined_junctions_file <- file.path(combined_junctions_dir, "junctions.bed")

# directory of shiba results produced from separate shiba runs
separate_dir <- file.path(shiba_dir, "separate_runs")
# make sample paths to the junctions.bed files for separate splice runs
junction_paths <- file.path(
  separate_dir,
  samples,
  "junctions",
  "junctions.bed"
)
# name the separate junction file paths
names(junction_paths) <- names(samples)

# Directory of PSI table, merged from separate shiba runs
separate_psi_table_dir <- file.path("merged_results")
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
```

Define function to read in junctions.bed files

``` r
# Make function for reading junctions from multiple samples run separately
read_junctions <- function(junction_paths) {
  # read in files
  purrr::map(junction_paths, \(file) {
    read.table(file, 
               header = TRUE, 
               sep="\t",
               stringsAsFactors=FALSE, 
               quote="")
  }) |>
    # merge junctions tables from multiple samples 
    # the resulting table separates junction counts from each sample by columns with the sample ID
    purrr::reduce(\(x, y) dplyr::full_join(x, y, by = c("ID", "start", "end", "chr")))
}
```

Read in files

``` r
## read in splice data from separate shiba run

# read merged splice table from separate runs
separate_psi_table <- readr::read_tsv(separate_psi_file)
```

    Rows: 352893 Columns: 6
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): event_type, pos_id, gene_id, label
    dbl (2): SRR601500_PSI, SRR604528_PSI

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read and merge junctions bedfile
separate_junctions <- read_junctions(junction_paths)

## read in splice data from combined shiba run

# read junctions file from combined run
combined_junctions <- read.table(combined_junctions_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

# read in combined splice results to compare against separate results
combined_splice_results_paths <- shiba_psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
      dplyr::mutate(across(contains("_PSI"), as.numeric))
  })

# combine psi values of samples run together into one dataframe to compare against separate PSI dataframe
combined_psi_table <- purrr::list_rbind(combined_splice_results_paths, names_to = "event_type")
```

## Check differences between separate and combined junction files

``` r
# print head of one chromosome's junctions for separate bed table
separate_junctions |>
  dplyr::filter(chr == "chr1") |>
  dplyr::arrange(start) |>
  head() 
```

| chr  | start     | end       | ID                       | SRR601500 | SRR604528 |
|:-----|:----------|:----------|:-------------------------|----------:|----------:|
| chr1 | 100007156 | 100011365 | chr1:100007156-100011365 |         9 |         7 |
| chr1 | 100007600 | 100011365 | chr1:100007600-100011365 |        NA |         2 |
| chr1 | 100011533 | 100015302 | chr1:100011533-100015302 |         9 |         7 |
| chr1 | 100015420 | 100017682 | chr1:100015420-100017682 |         7 |         6 |
| chr1 | 100017815 | 100022386 | chr1:100017815-100022386 |        12 |         4 |
| chr1 | 100017815 | 100024645 | chr1:100017815-100024645 |         2 |        NA |

Check junction counts of combined junctions.bed file

``` r
# print head of one chromosome's junctions for merged bed table
combined_junctions |>
  dplyr::filter(chr == "chr1") |>
  dplyr::arrange(start) |>
  head()
```

| chr  | start     | end       | ID                       | SRR601500 | SRR604528 |
|:-----|:----------|:----------|:-------------------------|----------:|----------:|
| chr1 | 100007156 | 100011365 | chr1:100007156-100011365 |         9 |         7 |
| chr1 | 100007600 | 100011365 | chr1:100007600-100011365 |         0 |         2 |
| chr1 | 100011533 | 100015302 | chr1:100011533-100015302 |         9 |         7 |
| chr1 | 100015420 | 100017682 | chr1:100015420-100017682 |         7 |         6 |
| chr1 | 100017815 | 100022386 | chr1:100017815-100022386 |        12 |         4 |
| chr1 | 100017815 | 100024645 | chr1:100017815-100024645 |         2 |         0 |

Generally, when there is a NA value for a junction on the separate
junctions bedfile, the combined junctions bedfile will have a value of 0
for that sample and junction.

### Obtain dimensions of the separate and combined junction count tables to compare values

``` r
# dimensions of separate PSI table (separate)
dim(separate_junctions)
```

    [1] 383376      6

``` r
# dimensions of combined PSI table
dim(combined_junctions)
```

    [1] 383322      6

There are only ~50 more junctions in the separate junctions table than
the combined table. One possible explanation is that there are extra
junctions represent positions that are very close together (in other
words, a fewer number of junctions would have been called if there were
better read support from all possible junction regions like in the
combined bedfile).

### Check if all NAs in one sample turn into 0 in the combined table

First, filter out junctions that are too short to be incorporated into
PSI calculation

Shiba uses the following default junction length thresholds to determine
what junctions are used for PSI calculation
([source](https://sika-zheng-lab.github.io/Shiba/quickstart/diff_splicing_bulk/?h=junction+length#1-prepare-inputs)):

    minimum_anchor_length:
      6 
    minimum_intron_length:
      70 
    maximum_intron_length:
      500000 
    strand:
      XS 

Filter for junctions that are greater than 1 bp for comparison. Here, I
am operating on the assumption that these junctions would get filtered
out anyway once PSI values are calculated.

``` r
# filter separate junctions table for junctions > 6bp and convert NAs to 0 for comparison
filtered_separate_junctions <- separate_junctions |>
  # calculate junction length
  dplyr::mutate(
    # calculate junction length
    junc_len = abs(as.numeric(start) - as.numeric(end))
    ) |>
  # replace NA junction counts with 0
  tidyr::replace_na(list(SRR604528 = 0, SRR601500 = 0)) |>
  # filter out junctions < 1bp
  dplyr::filter(
    junc_len > 1
  ) 
```

    Warning: There were 2 warnings in `dplyr::mutate()`.
    The first warning was:
    ℹ In argument: `junc_len = abs(as.numeric(start) - as.numeric(end))`.
    Caused by warning:
    ! NAs introduced by coercion
    ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.

``` r
# filter combined junctions table for junctions > 1bp
filtered_combined_junctions <- combined_junctions |>
  # calculate junction length
  dplyr::mutate(
    # calculate junction length
    junc_len = abs(as.numeric(start) - as.numeric(end))
  ) |>
  # filter out junctions < 1bp
  dplyr::filter(
    junc_len > 1
  ) 
```

    Warning: There were 2 warnings in `dplyr::mutate()`.
    The first warning was:
    ℹ In argument: `junc_len = abs(as.numeric(start) - as.numeric(end))`.
    Caused by warning:
    ! NAs introduced by coercion
    ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.

Check above warning - what are the NAs in the dataframes?

``` r
sum(is.na(filtered_combined_junctions))
```

    [1] 0

``` r
sum(is.na(filtered_separate_junctions))
```

    [1] 0

Despite the warning about NAs being introduced by coercion, no NA values
are present in either filtered dataframe.

Obtain list of junction IDs only in the separate or combined junctions
tables

``` r
combined_jcn_only <- setdiff(filtered_combined_junctions$ID, filtered_separate_junctions$ID)
separate_jcn_only <- setdiff(filtered_separate_junctions$ID, filtered_combined_junctions$ID)
```

After filtering for junctions longer than 6bp and replacing NAs with 0,
are there still junction IDs only in the combined or separate table?

``` r
# obtain fraction of events only found in the combined or separate junctions table
# Filter junction tables for junctions that are > 6bp
compare_junctions <- dplyr::full_join(
  filtered_separate_junctions,
  filtered_combined_junctions,
  by = c("ID", "chr", "start", "end"),
  suffix = c("_separate", "_combined")
  ) |>
  # label junction IDs by whether they are only in the combined or separate junctions table
  dplyr::mutate(
    combined_only = ID %in% combined_jcn_only,
    separate_only = ID %in% separate_jcn_only,
    shared = !ID %in% combined_jcn_only & !ID %in% separate_jcn_only
    )

# Print summary table of number of values in each category
compare_junctions |>
  dplyr::summarise(
    total = dplyr::n(),
    combined_only = sum(combined_only),
    separate_only = sum(separate_only),
    shared = sum(shared),
    frac_combined_only = combined_only / total,
    frac_separate_only = separate_only / total,
    frac_shared = shared / total
  )
```

| total | combined_only | separate_only | shared | frac_combined_only | frac_separate_only | frac_shared |
|---:|---:|---:|---:|---:|---:|---:|
| 350239 | 0 | 0 | 350239 | 0 | 0 | 1 |

Once we filter for junctions \> 1 bp , all IDs are shared in both
junction tables.

### Check similarity of counts in junctions that are matched

Of junctions with IDs that are the same between the tables, how many
have the same counts once we replace NAs with 0?

``` r
# identical returns false but all returns true
# perhaps because the classes of the values are not the same
all.equal(compare_junctions$SRR601500_combined, compare_junctions$SRR601500_separate)
```

    [1] TRUE

``` r
all(compare_junctions$SRR604528_combined == compare_junctions$SRR604528_separate)
```

    [1] TRUE

In conclusion, for this test case of 2 samples, once we filter out
junctions smaller than 1bp in length and replace NAs with 0, the counts
for both combined and separate tables become identical

## Examine splice events that are consistently present in both combined and separate splice tables

To help us prioritize what splice event types we may be the most
confident in after merging separate Shiba PSI tables, we are interested
in seeing a breakdown of what event types are most represented in both
PSI tables

``` r
# obtain df of splice events that are shared between both combined and separate tables
shared_events <- dplyr::full_join(
  combined_psi_table, 
  separate_psi_table,
  by = c("event_type", "pos_id", "gene_id", "label"),
  suffix = c("_combined", "_separate")) # label PSI values by table they came from

head(shared_events)
```

| event_type | pos_id | gene_id | label | SRR601500_PSI_combined | SRR604528_PSI_combined | SRR601500_PSI_separate | SRR604528_PSI_separate |
|:---|:---|:---|:---|---:|---:|---:|---:|
| se | SE@GL000008.2@156667-156758@154715-157528 | ENSG00000296775.1 | annotated | NA | NA | NA | NA |
| se | SE@GL000008.2@83860-84014@83545-85567 | ENSG00000296732.1 | annotated | NA | NA | NA | NA |
| se | SE@GL000008.2@83860-84014@83545-85457 | ENSG00000296732.1 | annotated | NA | NA | NA | NA |
| se | SE@GL000009.2@124033-124133@122622-124755 | ENSG00000306721.1 | annotated | NA | NA | NA | NA |
| se | SE@GL000009.2@124164-124272@122622-124755 | ENSG00000306721.1 | annotated | NA | NA | NA | NA |
| se | SE@GL000009.2@22153-22202@22061-28179 | ENSG00000297619.1 | annotated | NA | NA | NA | NA |

``` r
# recording a note that we may want to pivot this table longer so a "method" column tells us whether the values are from the combined or separate table
```

Print summary of event types in the shared events table

``` r
shared_events |>
  dplyr::summarise(.by = c(event_type, label),
                   count = dplyr::n()) |>
  dplyr::mutate(percent = count / sum(count) * 100)
```

| event_type | label       |  count |    percent |
|:-----------|:------------|-------:|-----------:|
| se         | annotated   |  66379 | 18.7667159 |
| se         | unannotated |    988 |  0.2793280 |
| afe        | annotated   | 100422 | 28.3913759 |
| afe        | unannotated |   3041 |  0.8597536 |
| ale        | annotated   |  78337 | 22.1474897 |
| ale        | unannotated |   1156 |  0.3268251 |
| five       | annotated   |  15221 |  4.3032914 |
| five       | unannotated |    406 |  0.1147846 |
| three      | annotated   |  19710 |  5.5724246 |
| three      | unannotated |    572 |  0.1617162 |
| mse        | annotated   |  46169 | 13.0529310 |
| mse        | unannotated |    428 |  0.1210045 |
| mxe        | annotated   |    586 |  0.1656743 |
| mxe        | unannotated |     17 |  0.0048063 |
| ri         | annotated   |  14108 |  3.9886233 |
| ri         | unannotated |   6166 |  1.7432557 |

Of the events that are the same between the combined and separate PSI
tables, annotated SE, AFE, ALE, and MSE events are the most abundant.

## Examine splice events only found in combined or separate splice tables

### Obtain dimensions of each PSI table

``` r
# dimensions of separate PSI table 
dim(separate_psi_table)
```

    [1] 352893      6

``` r
# dimensions of combined PSI table
dim(combined_psi_table)
```

    [1] 351307      6

There are 1,586 more splice events in the separte PSI table than the
combined table. However, it is still a small fraction of the total
number of events.

### Examine splice events only present in each PSI table

``` r
# produce dataframe of elements in combined splice results reference dataframe
combined_only_list <- setdiff(combined_psi_table$pos_id, separate_psi_table$pos_id)

# produce dataframe of elements in merged splice results dataframe but not in the reference df
separate_only_list <- setdiff(separate_psi_table$pos_id, combined_psi_table$pos_id)
```

### Obtain the number of splice events only found in each table

``` r
# events only present in separate tables
# Number of splice events only in separate PSI table 
length(separate_only_list)
```

    [1] 1839

``` r
# Number of splice events only in combined PSI table 
length(combined_only_list)
```

    [1] 657

### Examine percentage of each splice event type in each PSI table

Print the number of annotated vs. unannotated events across all event
types in each table

Fraction of annotated vs. unannotated events in the combined table

``` r
# filter for pos_ids unique to the shiba dataframe
psi_summary <- shared_events |>
  dplyr::mutate(
    combined_only = pos_id %in% combined_only_list,
    separate_only = pos_id %in% separate_only_list
  ) |>
# print summary of how many unique events with values that are novel vs. unannotated
  dplyr::summarise(.by = c(label, event_type),
                   combined_only = sum(combined_only),
                   separate_only = sum(separate_only),
                   # dplyr::n() gives the size of the group (annotated events or unannotated events)
                   total = dplyr::n(),
                   combined_frac = combined_only / total,
                   separate_frac = separate_only / total)

psi_summary            
```

| label | event_type | combined_only | separate_only | total | combined_frac | separate_frac |
|:---|:---|---:|---:|---:|---:|---:|
| annotated | se | 7 | 4 | 66379 | 0.0001055 | 0.0000603 |
| unannotated | se | 25 | 5 | 988 | 0.0253036 | 0.0050607 |
| annotated | afe | 103 | 680 | 100422 | 0.0010257 | 0.0067714 |
| unannotated | afe | 78 | 59 | 3041 | 0.0256495 | 0.0194015 |
| annotated | ale | 211 | 834 | 78337 | 0.0026935 | 0.0106463 |
| unannotated | ale | 54 | 123 | 1156 | 0.0467128 | 0.1064014 |
| annotated | five | 47 | 17 | 15221 | 0.0030878 | 0.0011169 |
| unannotated | five | 22 | 4 | 406 | 0.0541872 | 0.0098522 |
| annotated | three | 28 | 11 | 19710 | 0.0014206 | 0.0005581 |
| unannotated | three | 21 | 6 | 572 | 0.0367133 | 0.0104895 |
| annotated | mse | 0 | 1 | 46169 | 0.0000000 | 0.0000217 |
| unannotated | mse | 6 | 3 | 428 | 0.0140187 | 0.0070093 |
| annotated | mxe | 0 | 8 | 586 | 0.0000000 | 0.0136519 |
| unannotated | mxe | 0 | 2 | 17 | 0.0000000 | 0.1176471 |
| annotated | ri | 1 | 13 | 14108 | 0.0000709 | 0.0009215 |
| unannotated | ri | 54 | 69 | 6166 | 0.0087577 | 0.0111904 |

Together, the unmatched splice events in the separate PSI table make up
a small percent of the events.

## Inspect individual position IDs of unique events

We are most concerned if unannotated events have altered positions in
the separate table and are unmatched in the combined table (is Shiba
doing something strange to assign reads into position coordinates?)

Check for one event type (SE)

For skipped exon events, the first coordinate spans the exon of the
inclusion event and the second coordinate spans the intron of the
exclusion event (intron_c of [this
diagram](https://sika-zheng-lab.github.io/Shiba/output/shiba/#psi_setxt)).

``` r
# filter for unannotated SE events only in the separate table
se_only_in_separate <- separate_psi_table |>
  dplyr::filter(
    pos_id %in% separate_only_list,
    event_type == "se",
    label == "unannotated")

se_only_in_separate |>
  # print in markdown format so full values will be printed out instead of truncated
  knitr::kable(format = "markdown")
```

| event_type | pos_id | gene_id | label | SRR601500_PSI | SRR604528_PSI |
|:---|:---|:---|:---|---:|---:|
| se | SE@chr10@112448882-112449019@112447467-112460524 | ENSG00000151532.15 | unannotated | 0.200000 | NA |
| se | SE@chr12@98745524-98750748@98735634-98751355 | ENSG00000185046.21 | unannotated | 0.184466 | NA |
| se | SE@chr3@44862358-44862527@44862116-44864197 | ENSG00000169964.8 | unannotated | 0.559322 | NA |
| se | SE@chr1@155562153-155562349@155521618-155562779 | ENSG00000116539.14 | unannotated | NA | NA |
| se | SE@chr8@144510368-144510507@144510015-144510598 | ENSG00000167700.9 | unannotated | NA | 0.7380952 |

I don’t care as much about the events with NA values in both samples
(from eyeballing the junction counts tables, it seems like these
junctions would just be assigned 0), so I will check cases where there
is a PSI value in at least one sample.

## Check splice events where the PSI is low in one sample and NA in the other

### IGV view of alignments at chr10@112448882-112449019@112447467-112460524

![](images/chr10_112447467-112460524.png)

Above: IGV screenshot of this unannotated skipped exon event
(SE@chr10@112448882-112449019@112447467-112460524) sample alignments and
Refseq track. Red highlight indicates the second position coordinate
(chr10:112447467-112460524). An alternate explanation for this event if
it is not a real skipped exon is that it is just the last exon in a
ZDHHC6 isoform

Look for closest match to this event in the combined PSI table

``` r
combined_psi_table |>
  dplyr::filter(gene_id == "ENSG00000151532.15") |>
  knitr::kable(format = "markdown")
```

| event_type | pos_id | gene_id | label | SRR601500_PSI | SRR604528_PSI |
|:---|:---|:---|:---|---:|---:|
| se | SE@chr10@112448382-112448489@112447467-112448882 | ENSG00000151532.15 | annotated | NA | NA |
| se | SE@chr10@112533530-112533550@112527164-112538246 | ENSG00000151532.15 | annotated | 0.7551020 | NA |
| afe | AFE@chr10@112447467-112460524@112449019-112460524 | ENSG00000151532.15 | unannotated | 0.9230769 | 1 |
| ale | ALE@chr10@112668998-112736723;112668288-112668937;112538330-112668218;112527164-112538246;112464657-112527087@112464657-112484969 | ENSG00000151532.15 | annotated | 1.0000000 | NA |
| ale | ALE@chr10@112668998-112815290@112668998-112736723 | ENSG00000151532.15 | annotated | NA | NA |
| ale | ALE@chr10@112668998-112815290;112668288-112668937;112538330-112668218;112527164-112538246;112464657-112527087@112464657-112484969 | ENSG00000151532.15 | annotated | 1.0000000 | NA |
| ale | ALE@chr10@112668998-112815290;112668288-112668937;112538330-112668218;112533550-112538246;112527164-112533530;112464657-112527087@112464657-112484969 | ENSG00000151532.15 | annotated | 1.0000000 | NA |

Check junction counts for ENSG00000151532.15 events in the combined
junctions table close to the position start

``` r
combined_junctions |>
  dplyr::filter(chr == "chr10",
                dplyr::between(as.double(start), 112447467 - 20, 112447467 + 20)) |>
  knitr::kable(format = "markdown")
```

    Warning: There was 1 warning in `dplyr::filter()`.
    ℹ In argument: `dplyr::between(...)`.
    Caused by warning in `dplyr::between()`:
    ! NAs introduced by coercion

| chr   | start     | end       | ID                        | SRR601500 | SRR604528 |
|:------|:----------|:----------|:--------------------------|----------:|----------:|
| chr10 | 112447467 | 112448382 | chr10:112447467-112448382 |         2 |         0 |
| chr10 | 112447467 | 112448882 | chr10:112447467-112448882 |         5 |         7 |
| chr10 | 112447467 | 112460524 | chr10:112447467-112460524 |        12 |        19 |

Check junction counts for this region on the separate junctions table

``` r
separate_junctions |>
  dplyr::filter(chr == "chr10",
                dplyr::between(as.double(start), 112447467 - 20, 112447467 + 20)) |>
  knitr::kable(format = "markdown")
```

    Warning: There was 1 warning in `dplyr::filter()`.
    ℹ In argument: `dplyr::between(...)`.
    Caused by warning in `dplyr::between()`:
    ! NAs introduced by coercion

| chr   | start     | end       | ID                        | SRR601500 | SRR604528 |
|:------|:----------|:----------|:--------------------------|----------:|----------:|
| chr10 | 112447467 | 112448382 | chr10:112447467-112448382 |         2 |        NA |
| chr10 | 112447467 | 112448882 | chr10:112447467-112448882 |         5 |         7 |
| chr10 | 112447467 | 112460524 | chr10:112447467-112460524 |        12 |        19 |

Check alignments for this locus on IGV for the closest events in the
combined table

![](images/chr10_112447467-112448882.png)

Above: IGV screenshot of this skipped exon event
(SE@chr10@112448382-112448489@112447467-112448882) sample alignments and
Refseq track. Red highlight indicates the second position coordinate,
which corresponds to the intron of the exclusion isoform
(chr10:112447467-112448882).

![](images/VTI1A.png)

Above: IGV screenshot of the second kipped exon event detected for VTI1A
(SE@chr10@112533530-112533550@112527164-112538246) sample alignments and
Refseq track. Red highlight indicates the second position coordinate.

It seems like this event is not detected in the separate table, so it is
not too concerning that this event would be dropped in the separate
table (due to PSI values/junctions not being detected at that locus for
all samples).

### Check SE@chr12@98745524-98750748@98735634-98751355

![](images/chr12_98735634_98751355.png)

Above: IGV views of alignments at the highlighted locus
chr12:98735634-98751355 (flanking exons of ANKS1B event from ENSG ID)

![](images/chr12_98745524_98750748.png)

Above: IGV views of alignments at the highlighted locus
chr12:98745524-98750748 (skipped exon of ANKS1B event)

This event looks like a miscategorized splice event as well (if
anything, it looks like an unannotated retained intron in ANKS1B), so it
is less of a concern that this event would be filtered out if we select
for splice events with numeric PSI values in some minimum number of
samples.

Check junction counts for positions nearby in the combined junctions
table

``` r
combined_junctions |>
  dplyr::filter(chr == "chr12",
                dplyr::between(as.double(start), 98735634 - 20, 98735634 + 20)) |>
  knitr::kable(format = "markdown")
```

    Warning: There was 1 warning in `dplyr::filter()`.
    ℹ In argument: `dplyr::between(as.double(start), 98735634 - 20, 98735634 +
      20)`.
    Caused by warning in `dplyr::between()`:
    ! NAs introduced by coercion

| chr   | start    | end      | ID                      | SRR601500 | SRR604528 |
|:------|:---------|:---------|:------------------------|----------:|----------:|
| chr12 | 98735634 | 98751355 | chr12:98735634-98751355 |        42 |        18 |

Check closest junction counts in the separate junctions table

``` r
separate_junctions |>
  dplyr::filter(chr == "chr12",
                dplyr::between(as.double(start), 98735634 - 20, 98735634 + 20)) |>
  knitr::kable(format = "markdown")
```

    Warning: There was 1 warning in `dplyr::filter()`.
    ℹ In argument: `dplyr::between(as.double(start), 98735634 - 20, 98735634 +
      20)`.
    Caused by warning in `dplyr::between()`:
    ! NAs introduced by coercion

| chr   | start    | end      | ID                      | SRR601500 | SRR604528 |
|:------|:---------|:---------|:------------------------|----------:|----------:|
| chr12 | 98735634 | 98751355 | chr12:98735634-98751355 |        42 |        18 |

## Check event where the PSI is high in one sample but NA in another

### IGV view of alignments at chr8@144510368-144510507@144510015-144510598

![](images/chr8_144510015_144510598.png)

Above: IGV screenshot of this skipped exon event
(SE@chr8@144510368-144510507@144510015-144510598) sample alignments and
Refseq track. Red highlight indicates the second position coordinate
(chr10:144510015-144510598). The first coordinate lands squarely on the
second exon in the first, third, and fourth SLC33A2 isoform models, so
I’m unsure why this was labeled unannotated.

Look for closest match to this event in the combined PSI table

``` r
combined_psi_table |>
  dplyr::filter(gene_id == "ENSG00000167700.9") |>
  knitr::kable(format = "markdown")
```

| event_type | pos_id | gene_id | label | SRR601500_PSI | SRR604528_PSI |
|:---|:---|:---|:---|---:|---:|
| se | SE@chr8@144510360-144510507@144510015-144510598 | ENSG00000167700.9 | annotated | 0.9333333 | 0.7904762 |
| three | THREE@chr8@144510015-144510360@144510015-144510598 | ENSG00000167700.9 | annotated | 0.9142857 | 0.7027027 |
| ri | RI@chr8@144510722-144510794 | ENSG00000167700.9 | annotated | 0.1532847 | 0.1234568 |

The closest event in the combined table is
SE@chr8@144510360-144510507@144510015-144510598, which is an annotated
skipped exon event.

We can actually see that the coordinates are very close to the
unannotated event in the separate PSI table
(SE@chr8144510368-144510507@144510015-144510598). The event has the same
second position (junctions of the flanking exons), but the first
position (junctions of the exon that is skipped) is just off by 8bp,
causing it to be called as unannotated

Check junction counts for ENSG00000167700.9 events in the combined
junctions table close to the position start

``` r
combined_junctions |>
  dplyr::filter(chr == "chr8",
                dplyr::between(as.double(start), 144510015 - 20, 144510015 + 20)) |>
  knitr::kable(format = "markdown")
```

    Warning: There was 1 warning in `dplyr::filter()`.
    ℹ In argument: `dplyr::between(...)`.
    Caused by warning in `dplyr::between()`:
    ! NAs introduced by coercion

| chr  | start     | end       | ID                       | SRR601500 | SRR604528 |
|:-----|:----------|:----------|:-------------------------|----------:|----------:|
| chr8 | 144510014 | 144510397 | chr8:144510014-144510397 |         2 |         1 |
| chr8 | 144510015 | 144510360 | chr8:144510015-144510360 |        64 |        26 |
| chr8 | 144510015 | 144510368 | chr8:144510015-144510368 |        12 |         5 |
| chr8 | 144510015 | 144510425 | chr8:144510015-144510425 |         9 |         0 |
| chr8 | 144510015 | 144510598 | chr8:144510015-144510598 |         6 |        11 |

Check junction counts for this region on the separate junctions file

``` r
separate_junctions |>
  dplyr::filter(chr == "chr8",
                dplyr::between(as.double(start), 144510015 - 20, 144510015 + 20)) |>
  knitr::kable(format = "markdown")
```

    Warning: There was 1 warning in `dplyr::filter()`.
    ℹ In argument: `dplyr::between(...)`.
    Caused by warning in `dplyr::between()`:
    ! NAs introduced by coercion

| chr  | start     | end       | ID                       | SRR601500 | SRR604528 |
|:-----|:----------|:----------|:-------------------------|----------:|----------:|
| chr8 | 144510014 | 144510397 | chr8:144510014-144510397 |         2 |         1 |
| chr8 | 144510015 | 144510360 | chr8:144510015-144510360 |        64 |        26 |
| chr8 | 144510015 | 144510368 | chr8:144510015-144510368 |        12 |         5 |
| chr8 | 144510015 | 144510425 | chr8:144510015-144510425 |         9 |        NA |
| chr8 | 144510015 | 144510598 | chr8:144510015-144510598 |         6 |        11 |

The junction counts for events near the start position of this event are
the same for both the combined and separate junction tables
(specifically chr8:144510015-144510360 and chr8:144510015-144510368). It
seems that when Shiba is run one sample at a time, it will call
chr8:144510015-144510368 as its own splice event instead of assigning it
to the chr8:144510015-144510360 event.

The chr8:144510015-144510368 event is only found in the separate table,
which might be a good sign since it is not a real event and would be
dropped from our analysis if we filtered only for events with numeric
PSI values in all (or some minimum number of) samples.

## Check IGV alignments at SE@chr3@44862358-44862527@44862116-44864197

|     |                                             |
|-----|---------------------------------------------|
|     | SE@chr3@44862358-44862527@44862116-44864197 |

![](images/chr3_44862358_44862527.png)

Above: IGV views of alignments at the highlighted locus
chr3:44862358-44862527 (skipped exon of TMEM42 event)

![](images/chr3_44862116_44864197.png)

Above: IGV views of alignments at the highlighted locus
chr3:44862116-44864197 (flanking exon of TMEM42 event)

This actually looks like a real unannotated exon skipping event that
occurs in both samples, so I’m surprised that the PSI value was NA for
one of the samples.

Look for closest match in combined PSI table - surprisingly, there are
no events for this gene

``` r
combined_psi_table |>
  dplyr::filter(gene_id == "ENSG00000169964.8") |>
  knitr::kable(format = "markdown")
```

| event_type | pos_id | gene_id | label | SRR601500_PSI | SRR604528_PSI |
|:-----------|:-------|:--------|:------|--------------:|--------------:|

Check junction counts in the separate and combined junctions tables.
This time, I am looking for junctions matching the start and/or end
position of the unannotated exon

``` r
combined_junctions |>
  dplyr::filter(chr == "chr3",
                dplyr::between(as.double(start), 44862527 - 20, 44862527 + 20)) |>
  knitr::kable(format = "markdown")
```

    Warning: There was 1 warning in `dplyr::filter()`.
    ℹ In argument: `dplyr::between(as.double(start), 44862527 - 20, 44862527 +
      20)`.
    Caused by warning in `dplyr::between()`:
    ! NAs introduced by coercion

| chr  | start    | end      | ID                     | SRR601500 | SRR604528 |
|:-----|:---------|:---------|:-----------------------|----------:|----------:|
| chr3 | 44862527 | 44864197 | chr3:44862527-44864197 |        26 |        19 |
| chr3 | 44862527 | 44864222 | chr3:44862527-44864222 |         1 |         1 |
| chr3 | 44862527 | 44873329 | chr3:44862527-44873329 |         0 |         1 |

Check the junctions in the separate junctions table

``` r
separate_junctions |>
  dplyr::filter(chr == "chr3",
                dplyr::between(as.double(start), 44862527 - 20, 44862527 + 20)) |>
  knitr::kable(format = "markdown")
```

    Warning: There was 1 warning in `dplyr::filter()`.
    ℹ In argument: `dplyr::between(as.double(start), 44862527 - 20, 44862527 +
      20)`.
    Caused by warning in `dplyr::between()`:
    ! NAs introduced by coercion

| chr  | start    | end      | ID                     | SRR601500 | SRR604528 |
|:-----|:---------|:---------|:-----------------------|----------:|----------:|
| chr3 | 44862527 | 44864197 | chr3:44862527-44864197 |        26 |        19 |
| chr3 | 44862527 | 44864222 | chr3:44862527-44864222 |         1 |         1 |
| chr3 | 44862527 | 44873329 | chr3:44862527-44873329 |        NA |         1 |

Considering that there are high junction counts for this exon, I am
surprised that there are no splice events for this gene in the combined
PSI table. But since this was the only skipped exon event that looked
real and missed in the combined table, we should continue with running
Shiba separately on our files.
