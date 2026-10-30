# Obtain TARGET osteosarcoma SRA accessions
Cindy Liang (celiang@ucsc.edu)
2026-10-29

## Introduction

This notebook produces a file with the TARGET osteosarcoma accession
IDs, to be used to download the data. It will also print a summary of
how much space files will take up. Does not exclude TARGET pilot samples
from filtering, because the pilot did not include osteosarcomas.

### Inputs

- `SraRunTable-TARGET.csv`, metadata file of accession IDs for dbGaP
  samples with the TARGET dataset. Created by downloading the SRA run
  table on dbGaP, with the “phs000218” filter (parent phs ID for TARGET
  dataset)
- `metadata/pilot_shiba_run/target_accessions.tsv` accessions file of
  TARGET sequences downloaded and analyzed from the pilot Shiba run, to
  exclude from this filtering

We also use the [TARGET Molecular Characterization
Platforms](https://www.cancer.gov/ccg/research/genome-sequencing/target/using-target-data/technology#all-phase-1)
page, which describes what library prep approach was used for some
samples.

## Setup

Read in directory and file paths

``` r
# find the project directory (splicing-compendium)
# working directory is in scripts/
base_dir <- here::here()

# define metadata directory
metadata_dir <- file.path(base_dir, "metadata")
pilot_dir <- file.path(metadata_dir, "pilot_shiba_run")
gtex_target_metadata_dir <- file.path(metadata_dir, "filter_target_gtex")

## define paths to input files
# TARGET accessions metadata from dbGaP
target_sra_file <- file.path(gtex_target_metadata_dir, "SraRunTable-TARGET.csv")

# define path to output files
target_accession_path <-  file.path(gtex_target_metadata_dir, "target_osteo_accessions.tsv")

# read in files
# SRA metadata of TARGET samples
target_sra <- readr::read_csv(target_sra_file, col_types = readr::cols(Bytes = "d", .default = "c"))
```

## Filter TARGET dataset for download

``` r
target_filtered_sra <- target_sra |>
  # filter for TARGET accessions
  dplyr::filter(
    # Filter for TARGET samples
    biospecimen_repository == "NCI_TARGET",
    study_name == "TARGET: Osteosarcoma (OS)",
    # filter for samples with fastqs
    grepl("fastq", `DATASTORE filetype`),
    # filter for paired end reads
    LibraryLayout == "PAIRED",
    # Filter for RNA-seq
    analyte_type == "RNA",
    # don't exclude target pilot because no osteos there
  ) |>
  # transform Bytes column to numeric for estimating space
  dplyr::mutate(Bytes = as.numeric(Bytes)) |>
  # select for only relevant columns
  dplyr::select(
    Run, biospecimen_repository, `DATASTORE filetype`, body_site, Bytes, study_name, `DATASTORE provider`, analyte_type, `Assay Type`, `SRA Study`, BioProject, BioSample, LibraryLayout, LibrarySelection, `Center Name`
  )

# count how many rows there are to count number of expected samples
paste0("There are ", nrow(target_filtered_sra), " accessions matching the filters.")
```

    [1] "There are 179 accessions matching the filters."

``` r
# estimate amount of space files will take up
file_space <- sum(target_filtered_sra$Bytes) / 1e12
paste0("About ", file_space, " terabytes of space will be taken up by file downloads.")
```

    [1] "About 0.755886850559 terabytes of space will be taken up by file downloads."

Summarize samples by sequencing center of origin. All osteosarcomas
originate from the same sequencing center.

``` r
target_filtered_sra |> dplyr::count(study_name, `Center Name`)
```

| study_name                | Center Name |   n |
|:--------------------------|:------------|----:|
| TARGET: Osteosarcoma (OS) | NCI-MELTZER | 179 |

Export filtered TARGET accession file

``` r
readr::write_tsv(target_filtered_sra, file = file.path(target_accession_path))
```

Print session info

``` r
sessionInfo()
```

    R version 4.4.3 (2025-02-28)
    Platform: aarch64-apple-darwin20
    Running under: macOS 27.0.1

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: America/Los_Angeles
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    loaded via a namespace (and not attached):
     [1] crayon_1.5.3      vctrs_0.7.3       cli_3.6.6         knitr_1.51       
     [5] rlang_1.3.0       xfun_0.60         otel_0.2.0        generics_0.1.4   
     [9] jsonlite_2.0.0    bit_4.6.0         glue_1.8.1        rprojroot_2.1.1  
    [13] htmltools_0.5.9   hms_1.1.4         rmarkdown_2.31    evaluate_1.0.5   
    [17] tibble_3.3.1      tzdb_0.5.0        fastmap_1.2.0     yaml_2.3.12      
    [21] lifecycle_1.0.5   compiler_4.4.3    dplyr_1.2.1       pkgconfig_2.0.3  
    [25] here_1.0.2        rstudioapi_0.18.0 digest_0.6.39     R6_2.6.1         
    [29] tidyselect_1.2.1  readr_2.2.0       parallel_4.4.3    vroom_1.7.1      
    [33] pillar_1.11.1     magrittr_2.0.5    withr_3.0.3       tools_4.4.3      
    [37] bit64_4.8.2      
