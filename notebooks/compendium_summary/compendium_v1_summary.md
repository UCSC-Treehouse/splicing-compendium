# Summary of samples on splice compendium v1
Cindy Liang (celiang@ucsc.edu)
2026-07-27

## Introduction

Splice compendium v1 contains data from TARGET and GTEx datasets. This
notebook summarizes the distribution of tissue types, demographic
information, and sequencing center of origin information of samples
analyzed in the compendium.

Metadata used for this analysis are:

- [sample_sheet.tsv](https://github.com/UCSC-Treehouse/splicing-compendium/blob/main/config/sample_sheet.tsv):
  Accessions list of samples processed in splice compendium v1. Columns
  include accession ID (“sample”) and dataset (“group”).
- [gtex_accessions.tsv](https://github.com/UCSC-Treehouse/splicing-compendium/blob/main/metadata/filter_target_gtex/gtex_accessions.tsv):
  metadata file of GTEX sample accession IDs, along with age, sex of
  donor and tissue type of the sample.
- [target_accessions_clinical.tsv](https://github.com/UCSC-Treehouse/splicing-compendium/blob/celiang/new-clinical-metadata/metadata/filter_target_gtex/target_accessions_clinical.tsv):
  metadata file of TARGET sample accession IDs, along with age, sex of
  donor, cancer type of the sample, and sequencing center of origin of
  the sample.

## Setup

## Read in input and output paths

Read in directory and file paths

``` r
# find the root-level repo directory so we can access the other files
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# find the project directory (splicing-compendium)
# working directory is in scripts
base_dir <- here::here()

## Directories ##
# define metadata directory
metadata_dir <- file.path(base_dir, "metadata")
gtex_target_metadata_dir <- file.path(metadata_dir, "filter_target_gtex")

# config directory of compendium
config_dir <- file.path(repo_root, "config")

### Files ##
# gtex metadata with age info
gtex_metadata_file <- file.path(gtex_target_metadata_dir, "gtex_accessions.tsv")
# target metadata with age info
target_metadata_file <- file.path(gtex_target_metadata_dir, "target_accessions_clinical.tsv")
# sample sheet of compendium files 
compendium_sample_file <- file.path(config_dir, "sample_sheet.tsv")

## Read in files ##
gtex_metadata <- readr::read_tsv(
  gtex_metadata_file, 
  col_types = readr::cols(Bytes = "d", 
                          cumulative_tb = "d",
                          .default = "c")
)

target_metadata <- readr::read_tsv(
  target_metadata_file,
  col_types = readr::cols(
    age_at_earliest_diagnosis_in_years.diagnoses.xena_derived = "d",
    age_at_diagnosis_days = "d",
    .default = "c")
)

sample_df <- readr::read_tsv(
  compendium_sample_file,
  col_types = readr::cols(.default = "c")
)
```

## Filter TARGET and GTEx samples for what’s in the compendium

Not all samples filtered for downloading are processed as part of the
compendium due to `fasterq-dump` download issues or the file size being
too large to align in a timely manner.

``` r
# obtain list of accessions in compendium for each group

target_accessions_list <- sample_df |>
  dplyr::filter(group == "target") |>
  dplyr::pull(sample)

gtex_accessions_list <- sample_df |>
  dplyr::filter(group == "gtex") |>
  dplyr::pull(sample)

# filter gtex and target metadata for those corresponding to samples in compendium

target_compendium_df <- target_metadata |>
  dplyr::filter(Run %in% target_accessions_list) |>
  # the xena browser metadata's "name.project" has the cancer type in a prettier format 
  # (no "TARGET:" prefix) 
  # but since some ALL phase 2 samples were miscategorized as  ALL phase 3 in the xena metadata, 
  # we will use thee study_name column from dbGaP metadata
  dplyr::mutate(
    # add dataset column to facet plots by
    dataset = "TARGET",
    # remove TARGET prefix from study_name values
    study_name = stringr::str_remove(study_name, "TARGET: ")) |>
  # clean up TARGET metadata by selecting clinically relevant and batch-relevant information
  dplyr::select(
    # fields shared with gtex metadata
    dataset,
    Run, # accession ID, same format as GTEx,
    BioSample, # sample-specific accession ID from dbgap
    center_name = `Center Name`, # sequencing center of origin
    age_in_years = age_at_earliest_diagnosis_in_years.diagnoses.xena_derived, # in continuous floating point numbers
    sex = gender.demographic, # in "male" / "female" values
    # fields specific to target metadata
    target_patient_id = `_PATIENT`, # ID of patient sample came from - some samples came from the same patient
    # in cases where patient ID is missing, the ID can be derived from the biospecimen sample ID (first three fields separated by hyphens)
    target_biospecimen_sample_id = biospecimen_repository_sample_id, # specific replicate ID for a sample
    # additional demographic info only in target metadata
    target_race = race.demographic,
    target_ethnicity = ethnicity.demographic,
    # clinical metadata specific to target samples
    target_OS = OS,
    target_OS_time = OS.time,
    target_disease_type = disease_type,
    target_project_id = project_id.project,
    target_project_name = name.project,
    target_tumor_class = classification_of_tumor.diagnoses,
    target_primary_diagnosis = primary_diagnosis.diagnoses,
    target_sample_type = sample_type.samples,
    target_tissue_type = tissue_type.samples,
    target_histological_type = histological_type,
    target_body_site = body_site,
    target_study_name = study_name
  )

gtex_compendium_df <- gtex_metadata |>
  dplyr::filter(Run %in% gtex_accessions_list) |>
  dplyr::mutate(
    # add dataset column to facet plots by
    dataset = "GTEx"
  ) |>
  # clean up gtex metadata by selecting clinically relevant and batch-relevant information
  dplyr::select(
    # fields shared with target metadata
    dataset,
    Run, # accession ID
    BioSample, # sample-specific accession ID from dbgap
    age_in_years = AGE, # in 10-year bins (characters)
    sex = SEX, # coded as 1 or 2
    center_name = `Center Name`,
    # gtex-specific metadata fields
    gtex_batch_id = batch_id,
    gtex_version = version, # gtex version sample was added
    gtex_subject_id = SUBJID, # ID of the individual the sample came from - some samples come from the same subject
    gtex_body_site = body_site, # tissue of origin of sample
  )
```

Merge GTEx and TARGET metadata tables for analysis and unify format of
common column values. GTEx sex information is coded in values of 1 and
2, while TARGET values are in ‘male’ and ‘female’ format. According to
`GTEx_Analysis_v10_Annotations_SampleAttributesDD.xlsx` on the [GTEx
portal metadata download
page](https://gtexportal.org/home/downloads/adult-gtex/metadata), 1 =
Male and 2 = Female.

``` r
# ensure shared column values are in a consistent format before merging
cleaned_target_df <- target_compendium_df |>
   dplyr::mutate(
     # bin TARGET sample age in years into 10 year blocks to make data closer to GTEx representation
     age_in_years = cut(
       age_in_years,
       breaks = seq(0, 100, by = 10),
       left = TRUE,
       right = FALSE,
       labels = c(
         "0-9",
         "10-19",
         "20-29",
         "30-39",
         "40-49",
         "50-59",
         "60-69",
         "70-79",
         "80-89",
         "90-99"
       )
       )
    )

cleaned_gtex_df <- gtex_compendium_df |>
  dplyr::mutate(
    # relabel sex column values according to what they represent
    sex =
      dplyr::case_when(
        sex == 1 ~ "male",
        sex == 2 ~ "female"
  )
  )

merged_compendium_df <- dplyr::full_join(
  cleaned_target_df,
  cleaned_gtex_df,
  by = c(
    "dataset",
    "Run",
    "BioSample",
    "age_in_years",
    "sex",
    "center_name")
)

# print column names in merged df tom spot-check
colnames(merged_compendium_df)
```

     [1] "dataset"                      "Run"                         
     [3] "BioSample"                    "center_name"                 
     [5] "age_in_years"                 "sex"                         
     [7] "target_patient_id"            "target_biospecimen_sample_id"
     [9] "target_race"                  "target_ethnicity"            
    [11] "target_OS"                    "target_OS_time"              
    [13] "target_disease_type"          "target_project_id"           
    [15] "target_project_name"          "target_tumor_class"          
    [17] "target_primary_diagnosis"     "target_sample_type"          
    [19] "target_tissue_type"           "target_histological_type"    
    [21] "target_body_site"             "target_study_name"           
    [23] "gtex_batch_id"                "gtex_version"                
    [25] "gtex_subject_id"              "gtex_body_site"              

## Summaries of sample composition of GTEx and TARGET accessions

### Table of number of samples in each GTEx tissue type

Print a table summarizing number of samples in each tissue type in
compendium

``` r
# summarize number of samples in each GTEx tissue type
merged_compendium_df |>
  dplyr::filter(dataset == "GTEx") |>
  dplyr::group_by(gtex_body_site) |>
  dplyr::summarise(
    n = dplyr::n()
    )
```

| gtex_body_site                      |   n |
|:------------------------------------|----:|
| Cells - EBV-transformed lymphocytes | 144 |
| Kidney - Cortex                     |  29 |
| Muscle - Skeletal                   | 455 |
| Whole Blood                         | 449 |

Our curated GTEx compendium dataset is primarily composed of whole blood
(for leukemia comparisons) and skeletal muscle (for soft tissue
comparison) samples.

### Table of number of samples in each TARGET cancer type

``` r
# summarize number of samples in each GTEx tissue type
merged_compendium_df |>
  dplyr::filter(dataset == "TARGET") |>
  dplyr::group_by(target_study_name) |>
  dplyr::summarise(
    n = dplyr::n()
    )
```

| target_study_name                                    |   n |
|:-----------------------------------------------------|----:|
| Acute Lymphoblastic Leukemia (ALL) Expansion Phase 2 | 313 |
| Acute Lymphoblastic Leukemia (ALL) Pilot Phase 1     |   3 |
| Acute Myeloid Leukemia (AML)                         | 448 |
| Kidney, Clear Cell Sarcoma of the Kidney (CCSK)      |  13 |
| Kidney, Rhabdoid Tumor (RT)                          |  81 |
| Kidney, Wilms Tumor (WT)                             | 137 |
| Neuroblastoma (NBL)                                  | 161 |

## Print session info

``` r
sessionInfo()
```

    R version 4.4.3 (2025-02-28)
    Platform: aarch64-apple-darwin20
    Running under: macOS 26.5.2

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: America/Los_Angeles
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] patchwork_1.3.2 ggplot2_4.0.3  

    loaded via a namespace (and not attached):
     [1] bit_4.6.0          gtable_0.3.6       jsonlite_2.0.0     crayon_1.5.3      
     [5] dplyr_1.2.1        compiler_4.4.3     tidyselect_1.2.1   stringr_1.6.0     
     [9] parallel_4.4.3     scales_1.4.0       yaml_2.3.12        fastmap_1.2.0     
    [13] here_1.0.2         readr_2.2.0        R6_2.6.1           generics_0.1.4    
    [17] knitr_1.51         tibble_3.3.1       rprojroot_2.1.1    pillar_1.11.1     
    [21] RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.3.0        stringi_1.8.7     
    [25] xfun_0.60          S7_0.2.2           bit64_4.8.2        otel_0.2.0        
    [29] cli_3.6.6          withr_3.0.3        magrittr_2.0.5     digest_0.6.39     
    [33] grid_4.4.3         vroom_1.7.1        rstudioapi_0.18.0  hms_1.1.4         
    [37] lifecycle_1.0.5    vctrs_0.7.3        evaluate_1.0.5     glue_1.8.1        
    [41] farver_2.1.2       rmarkdown_2.31     tools_4.4.3        pkgconfig_2.0.3   
    [45] htmltools_0.5.9   
