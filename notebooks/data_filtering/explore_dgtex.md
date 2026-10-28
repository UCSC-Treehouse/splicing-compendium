# Explore Treehouse PolyA Compendium
Cindy Liang (celiang@ucsc.edu)
2026-10-28

## Introduction

We are interested in adding more datasets to the Treehouse Splice
Compendium prior to submitting the manuscript to a journal. Currently,
the priorities are to include samples in the Treehouse polyA compendia,
with healthy comparators from dGTEx.

This notebook reads sample metadata from the first release of the
developmental GTEx (dGTEx) project to explore how many RNA-seq samples
of what normal tissue types are represented. The dGTEx metadata does not
have file size information, so I will be unable to incorporate these
files into the glacier backup costs until I obtain dbGaP access to the
files.

RNA libraries from dGTEx are stranded from the metadata: “Illumina
Tru-Seq Strand Specific Large Insert RNA Sequencing (75M pairs) v1 -
DRAGEN”

## Setup

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

## Summarize comparator tissue types in dGTEx v1

``` r
dgtex_rna_summary <- dgtex_metadata_df |>
  dplyr::filter(
    ANALYTE_TYPE == "RNA:Total RNA"
  ) |>
  # select for tissue type column and sample ID column only for now
  dplyr::select(
    SAMPLE_ID,
    Tissue_type = SMTS
  ) |>
  dplyr::summarise(
    .by = Tissue_type,
    n_samples = dplyr::n()
  ) |>
    dplyr::arrange(n_samples) |>
  dplyr::mutate(
    total_dGTEx = sum(n_samples)
  )

dgtex_rna_summary
```

| Tissue_type     | n_samples | total_dGTEx |
|:----------------|----------:|------------:|
| Thymus          |         1 |         342 |
| Pituitary       |         2 |         342 |
| Kidney          |         3 |         342 |
| Brain           |         4 |         342 |
| Uterus          |         4 |         342 |
| Ovary           |         4 |         342 |
| Pancreas        |         5 |         342 |
| Adrenal Gland   |         6 |         342 |
| Blood Vessel    |         7 |         342 |
| Heart           |         7 |         342 |
| Spleen          |         8 |         342 |
| Nerve           |        10 |         342 |
| Testis          |        10 |         342 |
| Liver           |        11 |         342 |
| Stomach         |        14 |         342 |
| Thyroid         |        15 |         342 |
| Esophagus       |        18 |         342 |
| Lymph Node      |        18 |         342 |
| Skin            |        22 |         342 |
| Small Intestine |        24 |         342 |
| Colon           |        25 |         342 |
| Lung            |        36 |         342 |
| Adipose Tissue  |        40 |         342 |
| Muscle          |        48 |         342 |

The current version of dGTEx contains 342 samples, mostly represented by
muscle, lung, and adipose tissue.
