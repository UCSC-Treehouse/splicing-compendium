# Compare merged PSI table vs. Shiba PSI tables
Cindy Liang (celiang@ucsc.edu)
2026-03-12

## Set up

## Directories and files

``` r
# define samples list
samples <- c("SRR601500", "SRR604528")

# define repo directory to access merged psi table dir
base_dir = here::here()

# define the data directories
# psi table comparison dir
merging_dir <- file.path(base_dir, "exploration/merging-psi-tables")
# splice results dir
shiba_dir <- shiba_dir <- file.path(merging_dir, "shiba_results")

# separate runs dir
separate_dir <- file.path(shiba_dir, "separate_runs")

# shiba splice results directory of two samples run together
shiba_combined_dir <- file.path(shiba_dir, "combined_run")
shiba_splice_results_dir <- file.path(shiba_combined_dir, "splicing")
shiba_junctions_dir <- file.path(shiba_combined_dir, "junctions")

# merged PSI table of two samples run separately with shiba
merged_psi_table_dir <- file.path(merging_dir, "merged_results")

# Read in file paths

# merged PSI file of two samples obtained from merge-shiba-psi-tables.qmd
merged_psi_file <- file.path(merged_psi_table_dir, "merged_psi_table.tsv")

# shiba junctions file
shiba_junctions_file <- file.path(shiba_junctions_dir, "junctions.bed")

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
shiba_psi_paths <-file.path(shiba_splice_results_dir, psi_files)
names(shiba_psi_paths) <- names(psi_files)

# make sample paths to the junctions.bed files for separate splice runs
junction_paths <- file.path(
  
  separate_dir,
  samples,
  "junctions"
)
```

Define function to read in junctions.bed files

``` r
# Make function for reading event types for a single sample
read_junctions <- function(sample_id, junction_paths) {
  # construct file path for each sample
  file_paths <- file.path(junction_paths, "junctions.bed")
  # name each PSI table file path by event type
  names(file_paths) <- names(sample_id)
  
  # read in files
  purrr::map(file_paths, \(file) {
    read.table(file, 
               header = TRUE, 
               sep="\t",
               stringsAsFactors=FALSE, 
               quote="")
  }) |>
    # merge individual event type PSI tables to get one PSI table per sample
    purrr::reduce(\(x, y) dplyr::full_join(x, y, by = c("ID", "start", "end", "chr")))
}
```

Read in files

``` r
# read in manually combined splice table
merged_psi_table <- readr::read_tsv(merged_psi_file)
```

    Rows: 352893 Columns: 6
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): event_type, pos_id, gene_id, label
    dbl (2): SRR601500_PSI, SRR604528_PSI

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read in manually merged junctions table
merged_junctions <- read_junctions(samples, junction_paths)

shiba_junctions <- read.table(shiba_junctions_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

# read in shiba splice results to compare against manually combined results
shiba_splice_results_paths <- shiba_psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
      dplyr::mutate(across(contains("_PSI"), as.numeric))
  })

# combine psi values of samples run together into one dataframe to compare against manually combined dataframe
shiba_psi_table <- purrr::list_rbind(shiba_splice_results_paths, names_to = "event_type")
```

## Spot check differences between merged and shiba junction files

``` r
# print head of one chromosome's junctions for merged bed table
merged_junctions |>
  dplyr::filter(chr == "chr1") |>
  head() 
```

| chr  | start | end    | ID                | SRR601500 | SRR604528 |
|:-----|:------|:-------|:------------------|----------:|----------:|
| chr1 | 14599 | 185016 | chr1:14599-185016 |         1 |        NA |
| chr1 | 14614 | 16858  | chr1:14614-16858  |         2 |         1 |
| chr1 | 14697 | 185125 | chr1:14697-185125 |         1 |        NA |
| chr1 | 14764 | 24846  | chr1:14764-24846  |        12 |        13 |
| chr1 | 14788 | 185230 | chr1:14788-185230 |         1 |        NA |
| chr1 | 14829 | 14930  | chr1:14829-14930  |         1 |         3 |

Check junction counts of shiba junctions.bed file

``` r
# print head of one chromosome's junctions for merged bed table
shiba_junctions |>
  dplyr::filter(chr == "chr1") |>
  head()
```

| chr  | start | end    | ID                | SRR601500 | SRR604528 |
|:-----|:------|:-------|:------------------|----------:|----------:|
| chr1 | 11671 | 12010  | chr1:11671-12010  |         0 |         1 |
| chr1 | 14599 | 185016 | chr1:14599-185016 |         1 |         0 |
| chr1 | 14614 | 16858  | chr1:14614-16858  |         2 |         1 |
| chr1 | 14695 | 185175 | chr1:14695-185175 |         0 |         3 |
| chr1 | 14697 | 185125 | chr1:14697-185125 |         1 |         0 |
| chr1 | 14720 | 185148 | chr1:14720-185148 |         0 |         1 |

It is reassuring to see that junctions that are NA in samples in the
merged table appear as 0 in the shiba table (e.g. chr1:14599-185016)

## Obtain dimensions of the merged and shiba tables to compare values

``` r
paste0("Merged PSI table dimensions ", 
       dim(merged_psi_table)[1], 
       "x", 
       dim(merged_psi_table)[2])
```

    [1] "Merged PSI table dimensions 352893x6"

``` r
paste0("Shiba PSI table dimensions ", 
       dim(shiba_psi_table)[1],
       "x",
       dim(shiba_psi_table)[2])
```

    [1] "Shiba PSI table dimensions 351307x6"

## Obtain splice events unique to each table

``` r
# produce dataframe of elements in combined splice results reference dataframe but not in merged psi df
unique_shiba <- setdiff(shiba_psi_table$pos_id, merged_psi_table$pos_id)

# produce dataframe of elements in merged splice results dataframe but not in the reference df
unique_merged <- setdiff(merged_psi_table$pos_id, shiba_psi_table$pos_id)
```

## Obtain dimensions of the merged and shiba tables to compare values

``` r
paste0("Number of total unique merged events: ", 
       length(unique_merged))
```

    [1] "Number of total unique merged events: 1839"

``` r
paste0("Number of total unique Shiba events: ", 
       length(unique_shiba))
```

    [1] "Number of total unique Shiba events: 657"

## Examine percentage of each splice event type in each PSI table

Print the number of annotated vs. unannotated events in each event type
only in the merged table

``` r
# filter for pos_ids unique to the shiba dataframe
unique_shiba_psi_table <- shiba_psi_table |>
  dplyr::filter(pos_id %in% unique_shiba)

# print summary of how many unique events with values in both samples are novel vs. unannotated
unique_shiba_psi_table |> 
  dplyr::summarise(.by = c(event_type), 
                   count_annotated = sum(label == "annotated"),
                   frac_annotated = count_annotated / dplyr::n(),
                   count_unannotated = sum(label == "unannotated"),
                   frac_unannotated = count_unannotated / dplyr::n()
                     )
```

| event_type | count_annotated | frac_annotated | count_unannotated | frac_unannotated |
|:-----------|----------------:|---------------:|------------------:|-----------------:|
| se         |               7 |      0.2187500 |                25 |        0.7812500 |
| afe        |             103 |      0.5690608 |                78 |        0.4309392 |
| ale        |             211 |      0.7962264 |                54 |        0.2037736 |
| five       |              47 |      0.6811594 |                22 |        0.3188406 |
| three      |              28 |      0.5714286 |                21 |        0.4285714 |
| mse        |               0 |      0.0000000 |                 6 |        1.0000000 |
| ri         |               1 |      0.0181818 |                54 |        0.9818182 |

Print the number of annotated vs. unannotated events in each event type
only in the merged table

``` r
# filter for the pos_ids in the merged df that are unique
unique_merged_psi_table <- merged_psi_table |>
  dplyr::filter(pos_id %in% unique_merged)

# print summary of how many unique events with values in both samples are novel vs unannotated
unique_merged_psi_table |> 
  dplyr::summarise(.by = c(event_type), 
                   count_annotated = sum(label == "annotated"),
                   frac_annotated = count_annotated / dplyr::n(),
                   count_unannotated = sum(label == "unannotated"),
                   frac_unannotated = count_unannotated / dplyr::n())
```

| event_type | count_annotated | frac_annotated | count_unannotated | frac_unannotated |
|:-----------|----------------:|---------------:|------------------:|-----------------:|
| se         |               4 |      0.4444444 |                 5 |        0.5555556 |
| afe        |             680 |      0.9201624 |                59 |        0.0798376 |
| ale        |             834 |      0.8714734 |               123 |        0.1285266 |
| five       |              17 |      0.8095238 |                 4 |        0.1904762 |
| three      |              11 |      0.6470588 |                 6 |        0.3529412 |
| mse        |               1 |      0.2500000 |                 3 |        0.7500000 |
| mxe        |               8 |      0.8000000 |                 2 |        0.2000000 |
| ri         |              13 |      0.1585366 |                69 |        0.8414634 |

## Inspect individual position IDs of unique events

We are most concerned if unannotated events have altered positions in
the merged table and are unmatched in the shiba table (is Shiba doing
something strange to assign reads into position coordinates?)

Check for one event type (SE)

``` r
# filter for unannotated SE events only in the merged table
se_only_in_merged <- unique_merged_psi_table |>
  dplyr::filter(event_type == "se",
                label == "unannotated")

se_only_in_merged |>
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
junctions would just be assigned 0 anyway), so I will check cases where
there is a PSI value in at least one sample.

## Check event where the PSI is low in one sample and NA in the other

### IGV view of alignments at chr10@112448882-112449019@112447467-112460524

![](images/chr10_112447467-112460524.png)

Above: IGV screenshot of this unannotated skipped exon event
(SE@chr10@112448882-112449019@112447467-112460524) sample alignments and
Refseq track. Red highlight indicates the second position coordinate
(chr10:112447467-112460524). An alternate explanation for this event if
it is not a real skipped exon is that it is just the last exon in a
ZDHHC6 isoform

Look for closest match to this event in the shiba PSI table

``` r
shiba_psi_table |>
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

Check junction counts for ENSG00000151532.15 events in the shiba table
close to the position start

``` r
shiba_junctions |>
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

Check junction counts for this region on the merged junctions file

``` r
merged_junctions |>
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
shiba table

![](images/chr10_112447467-112448882.png)

Above: IGV screenshot of this skipped exon event
(SE@chr10@112448382-112448489@112447467-112448882) sample alignments and
Refseq track. Red highlight indicates the second position coordinate
(chr10:112447467-112448882).

![](images/Screenshot%202026-03-12%20at%203.53.35%20PM.png)

Above: IGV screenshot of the second kipped exon event detected for VTI1A
(SE@chr10@112533530-112533550@112527164-112538246) sample alignments and
Refseq track. Red highlight indicates the second position coordinate.

It seems like this event is not detected in the merged table, so it is
not too concerning that this event would be dropped in the merged table
(due to PSI values/junctions not being detected at that locus for all
samples).

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

Check junction counts for positions nearby in the shiba junctions table

``` r
shiba_junctions |>
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

Check closest junction counts in merged table

``` r
merged_junctions |>
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

Look for closest match to this event in the shiba PSI table

``` r
shiba_psi_table |>
  dplyr::filter(gene_id == "ENSG00000167700.9") |>
  knitr::kable(format = "markdown")
```

| event_type | pos_id | gene_id | label | SRR601500_PSI | SRR604528_PSI |
|:---|:---|:---|:---|---:|---:|
| se | SE@chr8@144510360-144510507@144510015-144510598 | ENSG00000167700.9 | annotated | 0.9333333 | 0.7904762 |
| three | THREE@chr8@144510015-144510360@144510015-144510598 | ENSG00000167700.9 | annotated | 0.9142857 | 0.7027027 |
| ri | RI@chr8@144510722-144510794 | ENSG00000167700.9 | annotated | 0.1532847 | 0.1234568 |

The closest event in the shiba table is
SE@chr8@144510360-144510507@144510015-144510598, which is an annotated
skipped exon event.

We can actually see that the coordinates are very close to the
unannotated event in the merged table
(SE@chr8144510368-144510507@144510015-144510598). The event has the same
second position (junctions of the flanking exons), but the first
position (junctions of the exon that is skipped) is just off by 8bp,
causing it to be called as unannotated

Check junction counts for ENSG00000167700.9 events in the shiba table
close to the position start

``` r
shiba_junctions |>
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

Check junction counts for this region on the merged junctions file

``` r
merged_junctions |>
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
the same for both the shiba and merged tables (specifically
chr8:144510015-144510360 and chr8:144510015-144510368). It seems that
when Shiba is run one sample at a time, it will call
chr8:144510015-144510368 as its own splice event instead of assigning it
to the chr8:144510015-144510360 event.

The chr8:144510015-144510368 event is only found in the merged table,
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

Look for closest match in shiba PSI table - surprisingly, there are no
events for this gene

``` r
shiba_psi_table |>
  dplyr::filter(gene_id == "ENSG00000169964.8") |>
  knitr::kable(format = "markdown")
```

| event_type | pos_id | gene_id | label | SRR601500_PSI | SRR604528_PSI |
|:-----------|:-------|:--------|:------|--------------:|--------------:|

Check junction counts in merged and shiba junctions tables. This time, I
am looking for junctions matching the start and/or end position of the
unannotated exon

``` r
shiba_junctions |>
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

Check merged

``` r
merged_junctions |>
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
surprised that there are no splice events for this gene in the shiba PSI
table.
