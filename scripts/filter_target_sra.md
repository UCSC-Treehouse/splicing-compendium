# filter_target_sra
Cindy Liang (celiang@ucsc.edu)
2025-07-08

## Introduction

The goals of this notebook are to filter the following SRA tables for
files to download:

- TARGET SRA table, for pediatric cancer bulk RNA-seq to download. These
  sequences will be used as a query group for performing splicing

- GTEX SRA table and metadata matched ages to subject IDs, for samples
  from the youngest donors

Additionally, obtain an estimate of how much space files will take up.
There is 2.64Tb of free space on the huge OpenStack instance.

Metadata used for this analysis are:

- `SraRunTable-TARGET.csv`, metadata file of accession IDs for dbGaP
  samples with the TARGET dataset. Created by downloading the SRA run
  table on dbGaP, with the “phs000218” filter (parent phs ID for TARGET
  dataset)

## Running Code

Read in directory and file paths

``` r
# set seed for random subsampling of TARGET dataset
set.seed(1)

# find the project directory
base_dir = here::here()

# find the root-level repo directory
repo_root = rprojroot::find_root(rprojroot::is_git_root)
 
# define metadata directory
metadata_dir <- file.path(repo_root, "data/metadata")
gtex_target_metadata_dir <- file.path(metadata_dir, "filter_target_gtex")

# define paths to input files
target_sra_file <- file.path(gtex_target_metadata_dir, "SraRunTable-TARGET.csv")

# define path to output file
target_accession_path <- file.path(gtex_target_metadata_dir, "target_accessions.tsv")

# read in files as variables
target_sra <- readr::read_csv(target_sra_file,  col_types = readr::cols(Bytes = "i", .default = "c"))
```

Filter accession metadata for TARGET RNA-seq datasets

``` r
# filter accession table for accession IDs of TARGET samples with fastqs
target_filtered_sra <- target_sra |>
  # filter only for accession IDs with sequence files associated with them
  dplyr::filter(
    # Filter for TARGET samples
    biospecimen_repository == "NCI_TARGET",
    # exclude osteosarcomas due to their complex transcriptome
    study_name != "TARGET: Osteosarcoma (OS)",
    # also exclude the pilot phase 1 ALL samples because there are far more ALL samples in pilot phase 2
    study_name != "TARGET: Acute Lymphoblastic Leukemia (ALL) Pilot Phase 1",
    # exclude cell line RNA-seq
    study_name != "TARGET: Cancer Model Systems (MDLS): Cell Lines and Xenografts (including PPTP)",
    # filter for samples with fastqs
    grepl("fastq", `DATASTORE filetype`)
    ) |>
 # select for only relevant columns
  dplyr::select(
    Run, biospecimen_repository, `DATASTORE filetype`, body_site, Bytes, study_name, `DATASTORE provider`, analyte_type, `Assay Type`, `SRA Study`, BioProject, BioSample
  ) 

# count how many rows there are
paste0("There are ", nrow(target_filtered_sra), " accessions matching the filters.")
```

    [1] "There are 1421 accessions matching the filters."

``` r
# estimate amount of space files will take up
file_space <- sum(target_filtered_sra$Bytes)/1e12
paste0("About ", file_space, " terabytes of space will be taken up by file downloads.")
```

    [1] "About 13.599346454039 terabytes of space will be taken up by file downloads."

Summarize distribution of analyte type for TARGET samples (RNA, DNA, or
both)

``` r
target_filtered_sra |> dplyr::count(analyte_type)
```

| analyte_type |    n |
|:-------------|-----:|
| RNA          | 1421 |

Summarize distribution of assay type for TARGET samples

``` r
target_filtered_sra |> dplyr::count(`Assay Type`)
```

| Assay Type |    n |
|:-----------|-----:|
| RNA-Seq    | 1421 |

Conclusions: All TARGET samples appear to be RNA-seq. Another dataset
with matched DNA & RNA sequencing will be needed to examine SNV status
of pediatric cancers.

Summarize distribution of tumor types across TARGET samples:

``` r
target_filtered_sra |> dplyr::count(study_name)
```

| study_name                                                   |   n |
|:-------------------------------------------------------------|----:|
| TARGET: Acute Lymphoblastic Leukemia (ALL) Expansion Phase 2 | 579 |
| TARGET: Acute Myeloid Leukemia (AML)                         | 450 |
| TARGET: Kidney, Clear Cell Sarcoma of the Kidney (CCSK)      |  13 |
| TARGET: Kidney, Rhabdoid Tumor (RT)                          |  81 |
| TARGET: Kidney, Wilms Tumor (WT)                             | 137 |
| TARGET: Neuroblastoma (NBL)                                  | 161 |

Conclusions: The predominant cancer type in the TARGET dataset consist
of leukemias (ALL and AML).

Subsample TARGET dataset to obtain a max of 20 samples per tumor type
(aiming to take up ~ half the space of the huge OpenStack instance,
~1.3TB

``` r
subsampled_target <- target_filtered_sra |>
  dplyr::slice_sample(n = 20, by = study_name)

# count how many rows there are
paste0("There are ", nrow(subsampled_target), " accessions matching the filters.")
```

    [1] "There are 113 accessions matching the filters."

``` r
# estimate amount of space files will take up
subsampled_file_space <- (sum(subsampled_target$Bytes)/1e12)
paste0("About ", subsampled_file_space, " terabytes of space will be taken up by file downloads.")
```

    [1] "About 1.226703850466 terabytes of space will be taken up by file downloads."

Summarize distribution of tumor types across TARGET samples, to confirm
subsampling:

``` r
subsampled_target |> dplyr::count(study_name)
```

| study_name                                                   |   n |
|:-------------------------------------------------------------|----:|
| TARGET: Acute Lymphoblastic Leukemia (ALL) Expansion Phase 2 |  20 |
| TARGET: Acute Myeloid Leukemia (AML)                         |  20 |
| TARGET: Kidney, Clear Cell Sarcoma of the Kidney (CCSK)      |  13 |
| TARGET: Kidney, Rhabdoid Tumor (RT)                          |  20 |
| TARGET: Kidney, Wilms Tumor (WT)                             |  20 |
| TARGET: Neuroblastoma (NBL)                                  |  20 |

Note: Clear Cell Sarcoma samples are the least represented in this
subsampled dataset, since it has the smallest N.

Export filtered and subsampled TARGET accession file

``` r
target_accessions <- target_filtered_sra |>
  dplyr::select(Run, `SRA Study`, BioProject, BioSample, study_name)

readr::write_tsv(target_accessions, file = file.path(target_accession_path))
```
