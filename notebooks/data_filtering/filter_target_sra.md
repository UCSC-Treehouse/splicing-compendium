# Filtering TARGET samples from SRA
Cindy Liang (celiang@ucsc.edu)
2026-04-09

## Introduction

The goals of this notebook are to filter the following SRA tables for
files to download:

- TARGET SRA table, for pediatric cancer bulk RNA-seq to download. These
  sequences will be used as a query group for performing differential
  splicing

Additionally, this notebook obtains an estimate of how much space files
will take up. There is 2.64Tb of free space on the huge OpenStack
instance. We have access to 5 more instances, which equate to 13.2TB of
space

Metadata used for this analysis are:

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

From Maryke’s work, we also learned that a subset of the TARGET ALL
phase 2 expansion is ribo-D. These samples are described in [Liu et
al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5535770/#S11). An
excel sheet containing IDs for samples that were used in this study
(which in the methods are described to have undergone ribo-depletion)
has been downloaded from the Liu et al 2017 supplemental data.

## Setup

## Read in input and output paths

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
# TARGET accessions analyzed in pilot Shiba run to exclude from final accessions list
pilot_target_file <- file.path(pilot_dir, "target_accessions.tsv")

## ALL expansion phase 2 metadata
# excel sheet that includes sample IDs of 265 ribo-D RNA-seq samples from the TARGET ALL phase 2 expansion
# sample IDs are under "RNA_barcode_D" column in the first tab
ribod_all_phase_two_samples <- file.path(gtex_target_metadata_dir, "41588_2017_BFng3909_MOESM2_ESM.xlsx")

# define path to output files
target_accession_path <-  file.path(gtex_target_metadata_dir, "target_accessions.tsv")

# read in files
target_sra <- readr::read_csv(target_sra_file, col_types = readr::cols(Bytes = "d", .default = "c"))
all_phase_two_ribod <- readxl::read_excel(ribod_all_phase_two_samples, sheet = "Table S1 cohort")
pilot_target <- readr::read_tsv(pilot_target_file, col_types = readr::cols(.default = "c"))
```

## Filter TARGET dataset for download

### Check TARGET documentation for any tissue types recorded to be ribo D

Based off the information on the [TARGET Molecular Characterization
Platforms](https://www.cancer.gov/ccg/research/genome-sequencing/target/using-target-data/technology#all-phase-1)
page, check whether there are cancer types that are recorded to be
ribo-deplete instead of polyA selected

| TARGET cancer type | Documented library selection method | Sequencing center |
|----|----|----|
| ALL phase 1 | PolyA | BCC |
| ALL phase 2 | PolyA | BCC |
| ALL phase 3 | Unspecified | St Jude |
| AML | PolyA (for ssRNA-seq) | BCC |
| AML induction failure | PolyA (for ssRNA-seq) | BCC |
| Kidney - Wilms tumor | PolyA | BCC |
| Kidney - clear cell sarcoma | **Ribo-deplete** | NCI Center for Cancer Research |
| Kidney - Rhabdoid tumor | PolyA | BCC |
| Neuroblastoma | polyA | BCC |

### Make lists of specific sample IDs that we know are ribo D

Associate sample IDs from [Liu et
al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5535770/#S11) that
were recorded to have undergone ribo-D selection with accession IDs

``` r
# Associate sample IDs that are recorded as ribo-deplete from the Liu et al 2017 study with run accessions from the TARGET dbGaP metadata
all_phase_two_ribo_d_to_remove <- all_phase_two_ribod |>
  dplyr::left_join(target_sra, by = c("RNA_barcode_D" = "biospecimen_repository_sample_id"))

# define list of rna barcode IDs to exclude
all_phase_two_exclude_list <- all_phase_two_ribo_d_to_remove$RNA_barcode_D
```

### Filter TARGET data

``` r
target_filtered_sra <- target_sra |>
  # filter for TARGET accessions
  dplyr::filter(
    # Filter for TARGET samples
    biospecimen_repository == "NCI_TARGET",
    # exclude cell line RNA-seq
    study_name != "TARGET: Cancer Model Systems (MDLS): Cell Lines and Xenografts (including PPTP)",
    # exclude osteosarcomas because they have a complex transcriptome and we have no bone GTEx comparator
    study_name != "TARGET: Osteosarcoma (OS)",
    # filter for samples with fastqs
    grepl("fastq", `DATASTORE filetype`),
    # filter for paired end reads
    LibraryLayout == "PAIRED",
    # Filter for RNA-seq
    analyte_type == "RNA",
    # exclude ALL phase 2 samples that have been documented to be ribo-D
    !biospecimen_repository_sample_id %in% all_phase_two_exclude_list,
    # exclude all clear cell sarcoma samples as they are ribo-D according to TARGET documentation
    study_name != "TARGET: Kidney\\, Clear Cell Sarcoma of the Kidney (CCSK)",
    # exclude accessions downloaded from the pilot run
    !Run %in% pilot_target$Run
  ) |>
  # transform Bytes column to numeric for estimating space
  dplyr::mutate(Bytes = as.numeric(Bytes)) |>
  # select for only relevant columns
  dplyr::select(
    Run, biospecimen_repository, `DATASTORE filetype`, body_site, Bytes, study_name, `DATASTORE provider`, analyte_type, `Assay Type`, `SRA Study`, BioProject, BioSample, LibraryLayout, LibrarySelection, `Center Name`
  )

# count how many rows there are
paste0("There are ", nrow(target_filtered_sra), " accessions matching the filters.")
```

    [1] "There are 1076 accessions matching the filters."

``` r
# estimate amount of space files will take up
file_space <- sum(target_filtered_sra$Bytes) / 1e12
paste0("About ", file_space, " terabytes of space will be taken up by file downloads.")
```

    [1] "About 10.825550981914 terabytes of space will be taken up by file downloads."

Summarize distribution of tumor types across TARGET samples:

``` r
target_filtered_sra |> dplyr::count(study_name, `Center Name`)
```

| study_name | Center Name | n |
|:---|:---|---:|
| TARGET: Acute Lymphoblastic Leukemia (ALL) Expansion Phase 2 | BCCAGSC | 304 |
| TARGET: Acute Lymphoblastic Leukemia (ALL) Pilot Phase 1 | STJUDE | 3 |
| TARGET: Acute Myeloid Leukemia (AML) | BCCAGSC | 371 |
| TARGET: Acute Myeloid Leukemia (AML) | HAIB | 64 |
| TARGET: Kidney, Rhabdoid Tumor (RT) | BCCAGSC | 66 |
| TARGET: Kidney, Wilms Tumor (WT) | BCCAGSC | 122 |
| TARGET: Neuroblastoma (NBL) | NCI-KHAN | 146 |

## Divide TARGET accessions into groups of 1.2TB samples

For partitioning data onto OpenStack instances, we need to divide the
TARGET accessions based on how much space they will take up. So that
there will be enough space for holding raw fastqs and alignments, we
will aim to download a little under half the available space on the huge
OpenStack instance (~1.2TB)

``` r
# define target sum of data
max_tb <- 1.2

target_filtered_sra <- target_filtered_sra |>
  dplyr::mutate(
    # calculate cumulative sum of Bytes column
    cumulative_tb = cumsum(Bytes) / 1e12,
    # assign batch number of data by dividing cumulative tb by max tb of interest and round to nearest integer
    batch_id = ceiling(cumulative_tb / 1.2)
  ) 
```

Summarize number of accessions in each batch

``` r
target_filtered_sra |>
  dplyr::summarise(
    .by = batch_id,
    n = dplyr::n()
  )
```

| batch_id |   n |
|---------:|----:|
|        1 | 123 |
|        2 | 107 |
|        3 | 118 |
|        4 | 115 |
|        5 | 108 |
|        6 | 109 |
|        7 | 123 |
|        8 | 108 |
|        9 | 161 |
|       10 |   4 |

Export filtered TARGET accession file

``` r
readr::write_tsv(target_filtered_sra, file = file.path(target_accession_path))
```

Print session info

``` r
sessionInfo()
```

    R version 4.4.2 (2024-10-31)
    Platform: x86_64-conda-linux-gnu
    Running under: Ubuntu 24.04 LTS

    Matrix products: default
    BLAS/LAPACK: /home/ubuntu/splicing-compendium/.pixi/envs/default/lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
     [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
     [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

    time zone: Etc/UTC
    tzcode source: system (glibc)

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    loaded via a namespace (and not attached):
     [1] crayon_1.5.3     vctrs_0.7.2      cli_3.6.5        knitr_1.51      
     [5] rlang_1.1.7      xfun_0.57        generics_0.1.4   renv_1.1.7      
     [9] jsonlite_2.0.0   glue_1.8.0       bit_4.6.0        rprojroot_2.1.1 
    [13] htmltools_0.5.9  readxl_1.4.5     hms_1.1.4        rmarkdown_2.31  
    [17] cellranger_1.1.0 evaluate_1.0.5   tibble_3.3.1     tzdb_0.5.0      
    [21] fastmap_1.2.0    yaml_2.3.12      lifecycle_1.0.5  compiler_4.4.2  
    [25] dplyr_1.2.0      pkgconfig_2.0.3  here_1.0.2       digest_0.6.39   
    [29] R6_2.6.1         tidyselect_1.2.1 readr_2.2.0      parallel_4.4.2  
    [33] vroom_1.7.0      pillar_1.11.1    magrittr_2.0.4   withr_3.0.2     
    [37] tools_4.4.2      bit64_4.6.0-1   
