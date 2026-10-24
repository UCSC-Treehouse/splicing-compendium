# Explore Treehouse PolyA Compendium
Cindy Liang (celiang@ucsc.edu)
2026-10-23

## Introduction

We are interested in adding more datasets to the Treehouse Splice
Compendium prior to submitting the manuscript to a journal. Currently,
the priorities are to include samples in the Treehouse polyA compendia,
with healthy comparators from dGTEx.

This notebook reads in the clinical metadata from the v25 Treehouse
PolyA Compendium and explores how many samples in each disease type are
represented. Additionally, it reads in sample metadata from the first
release of the developmental GTEx (dGTEx) project to explore how many
RNA-seq samples of what normal tissue types are represented.

## Setup

Libraries

``` r
library(ggplot2)
```

Files and directories

``` r
## Directories ##

repo_root <- rprojroot::find_root(rprojroot::is_git_root)
metadata_dir <- file.path(repo_root, "metadata")
## Files ##

# dgtex metadata with library prep and tissue type info
dgtex_metadata_file <- file.path(metadata_dir, "dGTEx_Analysis_2026-01-30_v1_Sample_Attributes_Dataset_DS_open.xlsx")
```

Read in files

``` r
## Read in files ##
dgtex_metadata_df <- readxl::read_xlsx(dgtex_metadata_file)
```

Explore comparator tissue types in dGTEx v1

RNA libraries appear to be stranded “Illumina Tru-Seq Strand Specific
Large Insert RNA Sequencing (75M pairs) v1 - DRAGEN”

Currently stranded libraries fail with Shiba:
https://github.com/Sika-Zheng-Lab/Shiba/issues/100

``` r
dgtex_rna_summary <- dgtex_metadata_df |>
  dplyr::filter(
    ANALYTE_TYPE == "RNA:Total RNA"
  ) |>
  # select for tissue type column and sample ID column only for now
  dplyr::select(
    SAMPLE_ID,
    SMTS
  ) |>
  dplyr::summarise(
    .by = SMTS,
    n = dplyr::n()
  ) |>
    dplyr::arrange(n)

dgtex_rna_summary
```

| SMTS            |   n |
|:----------------|----:|
| Thymus          |   1 |
| Pituitary       |   2 |
| Kidney          |   3 |
| Brain           |   4 |
| Uterus          |   4 |
| Ovary           |   4 |
| Pancreas        |   5 |
| Adrenal Gland   |   6 |
| Blood Vessel    |   7 |
| Heart           |   7 |
| Spleen          |   8 |
| Nerve           |  10 |
| Testis          |  10 |
| Liver           |  11 |
| Stomach         |  14 |
| Thyroid         |  15 |
| Esophagus       |  18 |
| Lymph Node      |  18 |
| Skin            |  22 |
| Small Intestine |  24 |
| Colon           |  25 |
| Lung            |  36 |
| Adipose Tissue  |  40 |
| Muscle          |  48 |
