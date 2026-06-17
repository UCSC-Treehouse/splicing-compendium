# Check Merged vs. Combined TARGET Pilot GTFs
Cindy Liang (celiang@ucsc.edu)
2026-06-16

## Introduction

## Background

This notebook checks junctions.bed files from two Shiba run methods for
differences. The goal is to get a sense of how close files that go into
calculating PSI values for events are to each other in the merged and
combined Shiba methods. Junctions files that are different would impact
PSI calculation, so I am mainly looking to see if mismatched junctions
are present and if present, how much of the data do they represent.

Files compared in this notebook are GTFs containing transcript
annotations of 88 TARGET pilot samples, produced from running Shiba with
two different methods.

“Combined” refers to shiba splice analysis run with shiba in a canonical
way. Inputs of this method are .bam files of each sample and a reference
annotation GTF. In the combined method, unmodified Shiba scripts are run
on the input files using shiba.py, which calls on Shiba scripts that
generate intermediate files and perform the final PSI calculation. The
relevant intermediate files that go into PSI calculation from this
method are:

- GTF, which is made with
  \[bam2gtf.py\](https://github.com/Sika-Zheng-Lab/Shiba/blob/v0.8.1/src/bam2gtf.py.
  This script uses a stringtie merge command as part of that merges GTFs
  made from individual samples together with the reference GTF

“Merged” refers to shiba splice analysis run with `merge_results.smk`.
The inputs of this method are junctions.bed and GTF files created by
Shiba run on each sample separately. `merge_results.smk` merges the GTFs
using the same stringtie merge command Shiba uses:

The snakemake GTF merge command:
`stringtie --merge -p {threads} -G {input.reference_gtf} -o {output.merged_gtf} $manifest`
Where `$manifest` refers to the list of GTFs to merge

The canonical/combined Shiba GTF merge command, for reference (from
[bam2gtf.py](https://github.com/Sika-Zheng-Lab/Shiba/blob/v0.8.1/src/bam2gtf.py)):

            "stringtie",
            "--merge",
            "-p", str(num_processors),
            "-G", reference_gtf,
            "-o", assembled_gtf
        ] + gtf_list

### Analysis outline:

- Check if GTFs are identical

## Setup

``` r
# import libraries
library(rtracklayer)
```

    Loading required package: GenomicRanges

    Loading required package: stats4

    Loading required package: BiocGenerics


    Attaching package: 'BiocGenerics'

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
        table, tapply, union, unique, unsplit, which.max, which.min

    Loading required package: S4Vectors


    Attaching package: 'S4Vectors'

    The following object is masked from 'package:utils':

        findMatches

    The following objects are masked from 'package:base':

        expand.grid, I, unname

    Loading required package: IRanges

    Loading required package: GenomeInfoDb

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
# gtf dir
combined_gtf_dir <- file.path(combined_results_dir, "annotation")

## merged shiba run results on target pilot samples ##
# shiba results dir
merged_results_dir <- file.path(repo_root, "results", "merged_shiba", "target_pilot")

### Files ###
# combined GTF
combined_gtf_file <- file.path(combined_gtf_dir, "assembled_annotation.gtf")

## merged shiba target pilot files
merged_gtf_file <- file.path(merged_results_dir, "merged_gtf.gtf")
```

Read in files

``` r
combined_gtf <- import(combined_gtf_file)
merged_gtf <- import(merged_gtf_file)

# quick and dirty check if they are identical
# this is a strict check and will fail if row order is different
identical(combined_gtf, merged_gtf)
```

    [1] FALSE

## Sort GTFs for comparison

The GTFs are not identical - first check if row/column orders are
different, or if there is another easy explanation why they are
different by looking at the head

``` r
head(combined_gtf)
```

    GRanges object with 6 ranges and 9 metadata columns:
            seqnames        ranges strand |    source       type     score
               <Rle>     <IRanges>  <Rle> |  <factor>   <factor> <numeric>
      [1] GL000008.2      157-5657      * | StringTie transcript      1000
      [2] GL000008.2      157-5657      * | StringTie exon            1000
      [3] GL000009.2 200967-201409      * | StringTie transcript      1000
      [4] GL000009.2 200967-201409      * | StringTie exon            1000
      [5] GL000009.2 174663-175059      * | StringTie transcript      1000
      [6] GL000009.2 174663-175059      * | StringTie exon            1000
              phase     gene_id transcript_id exon_number   gene_name ref_gene_id
          <integer> <character>   <character> <character> <character> <character>
      [1]      <NA>     MSTRG.1     MSTRG.1.1        <NA>        <NA>        <NA>
      [2]      <NA>     MSTRG.1     MSTRG.1.1           1        <NA>        <NA>
      [3]      <NA>     MSTRG.2     MSTRG.2.1        <NA>        <NA>        <NA>
      [4]      <NA>     MSTRG.2     MSTRG.2.1           1        <NA>        <NA>
      [5]      <NA>     MSTRG.3     MSTRG.3.1        <NA>        <NA>        <NA>
      [6]      <NA>     MSTRG.3     MSTRG.3.1           1        <NA>        <NA>
      -------
      seqinfo: 83 sequences from an unspecified genome; no seqlengths

``` r
head(merged_gtf)
```

    GRanges object with 6 ranges and 9 metadata columns:
            seqnames      ranges strand |    source       type     score     phase
               <Rle>   <IRanges>  <Rle> |  <factor>   <factor> <numeric> <integer>
      [1] GL000008.2    157-5657      * | StringTie transcript      1000      <NA>
      [2] GL000008.2    157-5657      * | StringTie exon            1000      <NA>
      [3] GL000009.2  5121-16309      + | StringTie transcript      1000      <NA>
      [4] GL000009.2   5121-5294      + | StringTie exon            1000      <NA>
      [5] GL000009.2 12446-12476      + | StringTie exon            1000      <NA>
      [6] GL000009.2 15157-16309      + | StringTie exon            1000      <NA>
              gene_id     transcript_id exon_number       gene_name
          <character>       <character> <character>     <character>
      [1]     MSTRG.1         MSTRG.1.1        <NA>            <NA>
      [2]     MSTRG.1         MSTRG.1.1           1            <NA>
      [3]     MSTRG.2 ENST00000749351.1        <NA> ENSG00000297619
      [4]     MSTRG.2 ENST00000749351.1           1 ENSG00000297619
      [5]     MSTRG.2 ENST00000749351.1           2 ENSG00000297619
      [6]     MSTRG.2 ENST00000749351.1           3 ENSG00000297619
                ref_gene_id
                <character>
      [1]              <NA>
      [2]              <NA>
      [3] ENSG00000297619.1
      [4] ENSG00000297619.1
      [5] ENSG00000297619.1
      [6] ENSG00000297619.1
      -------
      seqinfo: 83 sequences from an unspecified genome; no seqlengths

It looks like the row order is not the same in each GTF. Attempt to sort
the GTFs so that the rows appear in the same order, prior to comparison

``` r
sorted_combined_gtf <- sort(combined_gtf)
sorted_merged_gtf <- sort(merged_gtf)

identical(sorted_combined_gtf, sorted_merged_gtf)
```

    [1] FALSE

GTFs are still not identical, so print head to see if there is a clear
reason

``` r
head(sorted_merged_gtf)
```

    GRanges object with 6 ranges and 9 metadata columns:
          seqnames      ranges strand |    source       type     score     phase
             <Rle>   <IRanges>  <Rle> |  <factor>   <factor> <numeric> <integer>
      [1]     chr1 11121-11211      + | StringTie exon            1000      <NA>
      [2]     chr1 11121-14764      + | StringTie transcript      1000      <NA>
      [3]     chr1 11125-11211      + | StringTie exon            1000      <NA>
      [4]     chr1 11125-14764      + | StringTie transcript      1000      <NA>
      [5]     chr1 11410-11671      + | StringTie exon            1000      <NA>
      [6]     chr1 11410-14764      + | StringTie transcript      1000      <NA>
              gene_id     transcript_id exon_number   gene_name       ref_gene_id
          <character>       <character> <character> <character>       <character>
      [1]   MSTRG.391 ENST00000832824.1           1    DDX11L16 ENSG00000290825.2
      [2]   MSTRG.391 ENST00000832824.1        <NA>    DDX11L16 ENSG00000290825.2
      [3]   MSTRG.391 ENST00000832825.1           1    DDX11L16 ENSG00000290825.2
      [4]   MSTRG.391 ENST00000832825.1        <NA>    DDX11L16 ENSG00000290825.2
      [5]   MSTRG.391 ENST00000832826.1           1    DDX11L16 ENSG00000290825.2
      [6]   MSTRG.391 ENST00000832826.1        <NA>    DDX11L16 ENSG00000290825.2
      -------
      seqinfo: 83 sequences from an unspecified genome; no seqlengths

``` r
head(sorted_combined_gtf)
```

    GRanges object with 6 ranges and 9 metadata columns:
          seqnames      ranges strand |    source       type     score     phase
             <Rle>   <IRanges>  <Rle> |  <factor>   <factor> <numeric> <integer>
      [1]     chr1 11121-11211      + | StringTie exon            1000      <NA>
      [2]     chr1 11121-14764      + | StringTie transcript      1000      <NA>
      [3]     chr1 11125-11211      + | StringTie exon            1000      <NA>
      [4]     chr1 11125-14764      + | StringTie transcript      1000      <NA>
      [5]     chr1 11410-11671      + | StringTie exon            1000      <NA>
      [6]     chr1 11410-14764      + | StringTie transcript      1000      <NA>
              gene_id     transcript_id exon_number   gene_name       ref_gene_id
          <character>       <character> <character> <character>       <character>
      [1]   MSTRG.324 ENST00000832824.1           1    DDX11L16 ENSG00000290825.2
      [2]   MSTRG.324 ENST00000832824.1        <NA>    DDX11L16 ENSG00000290825.2
      [3]   MSTRG.324 ENST00000832825.1           1    DDX11L16 ENSG00000290825.2
      [4]   MSTRG.324 ENST00000832825.1        <NA>    DDX11L16 ENSG00000290825.2
      [5]   MSTRG.324 ENST00000832826.1           1    DDX11L16 ENSG00000290825.2
      [6]   MSTRG.324 ENST00000832826.1        <NA>    DDX11L16 ENSG00000290825.2
      -------
      seqinfo: 83 sequences from an unspecified genome; no seqlengths

## Check differences in junctions bedfiles

``` r
# convert sorted gtfs into dataframes
combined_df <- as.data.frame(sorted_combined_gtf)
merged_df <- as.data.frame(sorted_merged_gtf)

merged_gtf_only <- dplyr::setdiff(combined_df, merged_df)
merged_gtf_only |> head()
```

| seqnames | start | end | width | strand | source | type | score | phase | gene_id | transcript_id | exon_number | gene_name | ref_gene_id |
|:---|---:|---:|---:|:---|:---|:---|---:|---:|:---|:---|:---|:---|:---|
| chr1 | 11121 | 11211 | 91 | \+ | StringTie | exon | 1000 | NA | MSTRG.324 | ENST00000832824.1 | 1 | DDX11L16 | ENSG00000290825.2 |
| chr1 | 11121 | 14764 | 3644 | \+ | StringTie | transcript | 1000 | NA | MSTRG.324 | ENST00000832824.1 | NA | DDX11L16 | ENSG00000290825.2 |
| chr1 | 11125 | 11211 | 87 | \+ | StringTie | exon | 1000 | NA | MSTRG.324 | ENST00000832825.1 | 1 | DDX11L16 | ENSG00000290825.2 |
| chr1 | 11125 | 14764 | 3640 | \+ | StringTie | transcript | 1000 | NA | MSTRG.324 | ENST00000832825.1 | NA | DDX11L16 | ENSG00000290825.2 |
| chr1 | 11410 | 11671 | 262 | \+ | StringTie | exon | 1000 | NA | MSTRG.324 | ENST00000832826.1 | 1 | DDX11L16 | ENSG00000290825.2 |
| chr1 | 11410 | 14764 | 3355 | \+ | StringTie | transcript | 1000 | NA | MSTRG.324 | ENST00000832826.1 | NA | DDX11L16 | ENSG00000290825.2 |

Check for closest match to chr1 start position 11121 in combined gtf

``` r
combined_df |>
  dplyr::filter(
    seqnames == "chr1",
    start == "11121"
  )
```

| seqnames | start | end | width | strand | source | type | score | phase | gene_id | transcript_id | exon_number | gene_name | ref_gene_id |
|:---|---:|---:|---:|:---|:---|:---|---:|---:|:---|:---|:---|:---|:---|
| chr1 | 11121 | 11211 | 91 | \+ | StringTie | exon | 1000 | NA | MSTRG.324 | ENST00000832824.1 | 1 | DDX11L16 | ENSG00000290825.2 |
| chr1 | 11121 | 14764 | 3644 | \+ | StringTie | transcript | 1000 | NA | MSTRG.324 | ENST00000832824.1 | NA | DDX11L16 | ENSG00000290825.2 |

Both combined GTF entries for chr1:11121-any match the merged GTF
entries, so I’m unsure why it was flagged as different

Check if column names are the same between each dataframe

``` r
dplyr::setdiff(names(merged_df), names(combined_df))
```

    character(0)

They are the same (setdiff returns 0)

Check if classes in each column are the same between data frames (they
are)

``` r
sapply(merged_df, class) == sapply(combined_df, class)
```

         seqnames         start           end         width        strand 
             TRUE          TRUE          TRUE          TRUE          TRUE 
           source          type         score         phase       gene_id 
             TRUE          TRUE          TRUE          TRUE          TRUE 
    transcript_id   exon_number     gene_name   ref_gene_id 
             TRUE          TRUE          TRUE          TRUE 

Check for differences column by column

``` r
print("seqnames")
```

    [1] "seqnames"

``` r
dplyr::setdiff(merged_df$seqnames, combined_df$seqnames)
```

    character(0)

``` r
print("start")
```

    [1] "start"

``` r
dplyr::setdiff(merged_df$start, combined_df$start) |> head()
```

    [1] 1020076 1041296 1450494 1471795 1503248 1615868

``` r
dplyr::setdiff(combined_df$start, merged_df$start) |> head()
```

    [1] 629067 631814 632131 633879 633983 634200

``` r
print("end")
```

    [1] "end"

``` r
dplyr::setdiff(merged_df$end, combined_df$end) |> head()
```

    [1]  784049  908512 1280216 1509736 1450516 1626366

``` r
dplyr::setdiff(combined_df$end, merged_df$end) |> head()
```

    [1] 630538 634095 632077 784058 857835 908516

``` r
print("width")
```

    [1] "width"

``` r
dplyr::setdiff(merged_df$width, combined_df$width) |> head()
```

    [1]  50872 111760  75735  78835  76715  86499

``` r
dplyr::setdiff(combined_df$width, merged_df$width) |> head()
```

    [1]  53645  78802  86118  83625 108480 106801

``` r
print("strand")
```

    [1] "strand"

``` r
dplyr::setdiff(merged_df$strand, combined_df$strand) |> head()
```

    character(0)

``` r
dplyr::setdiff(combined_df$strand, merged_df$strand) |> head()
```

    character(0)

``` r
print("source")
```

    [1] "source"

``` r
dplyr::setdiff(merged_df$source, combined_df$source) |> head()
```

    character(0)

``` r
dplyr::setdiff(combined_df$source, merged_df$source) |> head()
```

    character(0)

``` r
print("type")
```

    [1] "type"

``` r
dplyr::setdiff(merged_df$type, combined_df$type) |> head()
```

    character(0)

``` r
dplyr::setdiff(combined_df$type, merged_df$type) |> head()
```

    character(0)

``` r
print("score")
```

    [1] "score"

``` r
dplyr::setdiff(merged_df$score, combined_df$score) |> head()
```

    numeric(0)

``` r
dplyr::setdiff(combined_df$score, merged_df$score) |> head()
```

    numeric(0)

``` r
print("phase")
```

    [1] "phase"

``` r
dplyr::setdiff(merged_df$phase, combined_df$phase) |> head()
```

    integer(0)

``` r
dplyr::setdiff(combined_df$phase, merged_df$phase) |> head()
```

    integer(0)

``` r
print("gene_id")
```

    [1] "gene_id"

``` r
dplyr::setdiff(merged_df$gene_id, combined_df$gene_id) |> head()
```

    [1] "ENSG00000251823.2" "ENSG00000283591.1" "MSTRG.70214"      
    [4] "MSTRG.70218"       "MSTRG.70220"       "MSTRG.70222"      

``` r
dplyr::setdiff(combined_df$gene_id, merged_df$gene_id) |> head()
```

    [1] "ENSG00000296088.1"  "ENSG00000223823.1"  "ENSG00000306476.1" 
    [4] "ENSG00000169885.10" "ENSG00000287396.2"  "ENSG00000302855.1" 

``` r
print("transcript_id")
```

    [1] "transcript_id"

``` r
dplyr::setdiff(merged_df$transcript_id, combined_df$transcript_id) |> head()
```

    [1] "MSTRG.399.3" "MSTRG.364.1" "MSTRG.364.2" "MSTRG.363.4" "MSTRG.363.3"
    [6] "MSTRG.371.1"

``` r
dplyr::setdiff(combined_df$transcript_id, merged_df$transcript_id) |> head()
```

    [1] "MSTRG.331.3" "MSTRG.348.2" "MSTRG.349.1" "MSTRG.349.2" "MSTRG.348.3"
    [6] "MSTRG.348.4"

``` r
print("exon_number")
```

    [1] "exon_number"

``` r
dplyr::setdiff(merged_df$exon_number, combined_df$exon_number) |> head()
```

    character(0)

``` r
dplyr::setdiff(combined_df$exon_number, merged_df$exon_number) |> head()
```

    character(0)

``` r
print("gene_name")
```

    [1] "gene_name"

``` r
dplyr::setdiff(merged_df$gene_name, combined_df$gene_name) |> head()
```

    [1] "ELOA3CP"         "ENSG00000285723"

``` r
dplyr::setdiff(combined_df$gene_name, merged_df$gene_name) |> head()
```

    [1] "TUBB8P7"         "ENSG00000291265" "ENSG00000293284"

``` r
print("ref_gene_id")
```

    [1] "ref_gene_id"

``` r
dplyr::setdiff(merged_df$ref_gene_id, combined_df$ref_gene_id ) |> head()
```

    [1] "ENSG00000275553.5" "ENSG00000285723.2"

``` r
dplyr::setdiff(combined_df$ref_gene_id, merged_df$ref_gene_id ) |> head()
```

    [1] "ENSG00000289718.2" "ENSG00000261812.7" "ENSG00000291265.2"
    [4] "ENSG00000293284.1"

Some of the differences seem like they come from stringtie-assigned
transcripts that don’t otherwise have transcript IDs

Investigate ELOA3CP - Does it really only appear in the merged dataframe
and not the combined?

``` r
merged_df |>
  dplyr::filter(gene_name == "ELOA3CP")
```

| seqnames | start | end | width | strand | source | type | score | phase | gene_id | transcript_id | exon_number | gene_name | ref_gene_id |
|:---|---:|---:|---:|:---|:---|:---|---:|---:|:---|:---|:---|:---|:---|
| chr18 | 46968695 | 46969912 | 1218 | \- | StringTie | transcript | 1000 | NA | MSTRG.36039 | ENST00000620522.2 | NA | ELOA3CP | ENSG00000275553.5 |
| chr18 | 46968695 | 46969912 | 1218 | \- | StringTie | exon | 1000 | NA | MSTRG.36039 | ENST00000620522.2 | 1 | ELOA3CP | ENSG00000275553.5 |

``` r
combined_df |>
  dplyr::filter(gene_name == "ELOA3CP")
```

| seqnames | start | end | width | strand | source | type | score | phase | gene_id | transcript_id | exon_number | gene_name | ref_gene_id |
|:---|---:|---:|---:|:---|:---|:---|---:|---:|:---|:---|:---|:---|:---|

Yes, this transcript is only in the merged dataframe for some reason
