# Filtering TARGET samples from SRA
Cindy Liang (celiang@ucsc.edu)
2026-10-29

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
target_osteo_accession_path <-  file.path(gtex_target_metadata_dir, "osteo_target_accessions.tsv")

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

### Filter TARGET data for 2026_09 release

``` r
target_filtered_sra <- target_sra |>
  # filter for TARGET accessions
  dplyr::filter(
    # Filter for TARGET samples
    biospecimen_repository == "NCI_TARGET",
    # exclude cell line RNA-seq
    study_name != "TARGET: Cancer Model Systems (MDLS): Cell Lines and Xenografts (including PPTP)",
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

    [1] "There are 1255 accessions matching the filters."

``` r
# estimate amount of space files will take up
file_space <- sum(target_filtered_sra$Bytes) / 1e12
paste0("About ", file_space, " terabytes of space will be taken up by file downloads.")
```

    [1] "About 11.581437832473 terabytes of space will be taken up by file downloads."

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
| TARGET: Osteosarcoma (OS) | NCI-MELTZER | 179 |

``` r
target_filtered_sra |>
  dplyr::filter(
    study_name == "TARGET: Osteosarcoma (OS)"
  ) |>
    dplyr::mutate(
      total_samples = dplyr::n(),
      total_gb = sum(Bytes) / 1000000000
    )
```

| Run | biospecimen_repository | DATASTORE filetype | body_site | Bytes | study_name | DATASTORE provider | analyte_type | Assay Type | SRA Study | BioProject | BioSample | LibraryLayout | LibrarySelection | Center Name | total_samples | total_gb |
|:---|:---|:---|:---|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|---:|---:|
| SRR1811024 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 4583405840 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736332 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811025 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 753017638 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736330 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811026 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3390290511 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736324 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811027 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3736989330 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736334 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811028 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 3959655879 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736331 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811029 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 3440214042 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736334 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811030 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 5613002965 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736324 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811031 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 4907203193 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736330 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811032 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 5254065694 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736332 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811033 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3473746940 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736356 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811034 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 4079953246 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736354 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811035 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 4326296875 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736350 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811036 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3833027895 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736345 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811037 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 4506697909 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736342 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811038 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 3908805392 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736356 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811039 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 4476627786 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736354 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811040 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 4553338517 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736350 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811041 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3950056354 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736345 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811042 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 4963443621 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736342 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811043 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 4332308845 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736338 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811044 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3962840058 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736338 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811045 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 2658850416 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736331 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811046 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3103896123 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736330 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811047 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 2544248089 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736384 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811048 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2750175990 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736377 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811049 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 2350699851 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736372 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811050 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 2680170694 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736384 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811051 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2652085837 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736377 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811052 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2520313068 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736372 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811053 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 3338638024 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736362 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811054 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 3997509228 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736375 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811055 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 3954718717 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736375 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811056 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 3257845042 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736379 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811057 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3291382031 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736379 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811058 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 5530066819 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736364 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811059 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 5515084069 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736364 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811060 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 1796727046 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736369 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811061 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 1799220610 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736368 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811062 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 1674315764 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736367 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811063 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 1855263840 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736366 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811064 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 1764207511 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736347 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811065 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 1821141792 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736340 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811066 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 1928142511 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736325 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811067 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 1103282218 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736388 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811068 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 1421801045 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260991 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811069 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 908183148 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261067 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811070 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 1036665787 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261045 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811071 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 1099302915 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260982 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811072 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 990808908 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261012 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811073 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 1146508567 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260982 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811074 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 2255258365 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736388 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811075 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 2696849142 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260991 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811076 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 1618676359 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261067 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811077 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2022570216 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261045 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811078 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 2240753490 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260982 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811079 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 2197428214 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261012 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811080 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 2346621909 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260982 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811081 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 3032567747 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736347 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811082 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 2946435942 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736366 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811083 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2892641666 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736368 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811084 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 2772728421 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736367 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811085 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 2931867327 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736370 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811086 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 2881531059 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736369 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811087 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 2214500772 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261003 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811088 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 2244376571 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261003 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811089 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2002924522 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261059 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811090 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 2071830630 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261059 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811091 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 4027182975 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261023 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811092 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3671563447 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261017 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811093 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 3063753556 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261032 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811094 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2949180675 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261020 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811095 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3276170668 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260977 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811096 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 6111078115 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443617 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811097 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3750673911 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261003 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811098 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 3265157802 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261003 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811099 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 4750037138 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736328 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811100 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3488735672 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261050 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811101 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 6510563380 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261064 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811102 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 7127623126 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260996 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811103 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 4581207650 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261029 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811104 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3352905305 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260977 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811105 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 2738705781 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261020 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811106 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 2759722349 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261032 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811107 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3047237948 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261017 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811108 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3594450622 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261023 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811109 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 3461829528 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736328 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811110 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2002027372 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261050 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811111 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 5259898240 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260996 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811112 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 4263166624 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261029 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811113 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 5782907687 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443617 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811114 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 3412748448 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260996 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811115 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 3307073310 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261050 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811116 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 3493016049 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261029 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811117 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3718138401 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736328 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811118 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 6044958275 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261064 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811119 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 4247233972 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736400 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811120 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 5113249626 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736398 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811121 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3604749307 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736396 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811122 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3588199313 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736394 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811123 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3458537309 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736381 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811124 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 3752049884 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736360 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811125 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 3620483906 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736352 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811126 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 3970999585 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736352 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811127 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 3987856185 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736360 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811128 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 3367090051 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736381 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811129 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 3882526348 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736394 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811130 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 3668417174 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736396 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811131 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 5002906305 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736398 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811132 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 4168340755 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736400 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811133 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 1240132368 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774016 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811134 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 7050396753 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774001 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811135 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 8582240950 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773995 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811136 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 8699542804 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774010 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811137 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 9951600155 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773980 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811138 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 87247262 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773989 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811139 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 8456891778 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773983 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811140 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 5606667824 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773956 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811141 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 5227202909 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773944 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811142 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 2827747495 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261067 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811143 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 3746446594 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736388 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811144 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 4236035801 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261045 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811145 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 4088883134 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00261012 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811146 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 4289560115 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00260991 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811147 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 2086151873 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736347 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811148 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 13967621787 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774013 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811149 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 9169281740 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773974 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811150 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 7472810097 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773962 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811152 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 6141281631 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02047610 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811153 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 7821560494 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773998 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811154 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 8379854543 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774022 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811155 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 8255957841 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774007 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811156 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 9056852577 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773977 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811157 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 7260911822 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774004 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811158 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 8212156035 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773965 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811159 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 11242101733 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773950 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811160 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 3520261174 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773983 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811161 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 3906873828 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773980 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811162 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 4167596697 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773977 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811163 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3401865964 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773974 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811164 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3308982383 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773983 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811165 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 3266836003 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773980 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811166 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 3682118398 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773974 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811167 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 4443891321 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773977 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811168 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 2896369408 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774016 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811169 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 3994807143 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774016 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811170 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 2854087105 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773989 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811171 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 3381366662 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773989 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811172 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 579508146 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443657 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811173 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 493736973 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443621 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811174 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 512847629 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443601 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811177 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 2352279758 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443637 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811178 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 2931690974 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443623 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811179 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 1585940569 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773989 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811180 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2299011175 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443637 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811181 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 2885695605 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443623 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811182 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 1594341721 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773989 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811183 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 23222110426 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00774016 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811184 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 1271261573 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773989 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811185 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 4543803897 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00773944 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811186 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 9476729099 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02047609 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811187 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 2742199573 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736362 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811188 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 11135139164 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443613 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811189 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 6916278228 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443623 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811190 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 6496542999 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443609 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811191 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 10009272612 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443657 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811192 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 8885994990 | TARGET: Osteosarcoma (OS) | ncbi,gs,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443621 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811193 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 6390556803 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443637 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811194 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 4809371170 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443622 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811195 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 11054874244 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443601 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811196 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 6190642885 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443633 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811197 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 10152148400 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443597 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811198 | NCI_TARGET | fastq,run.zq,sra | Primary Solid Tumor | 7567620828 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443608 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811199 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 6617522774 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN02443642 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811200 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 3772104547 | TARGET: Osteosarcoma (OS) | ncbi,s3,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736358 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811201 | NCI_TARGET | sra,run.zq,fastq | Primary Solid Tumor | 4160544214 | TARGET: Osteosarcoma (OS) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736336 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811202 | NCI_TARGET | fastq,sra,run.zq | Primary Solid Tumor | 4350422532 | TARGET: Osteosarcoma (OS) | s3,ncbi,gs | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736336 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811203 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 4137377155 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736358 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811204 | NCI_TARGET | sra,fastq,run.zq | Primary Solid Tumor | 2424562734 | TARGET: Osteosarcoma (OS) | gs,s3,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736340 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |
| SRR1811205 | NCI_TARGET | run.zq,sra,fastq | Primary Solid Tumor | 2232362645 | TARGET: Osteosarcoma (OS) | s3,gs,ncbi | RNA | RNA-Seq | SRP012003 | PRJNA89527 | SAMN00736325 | PAIRED | cDNA | NCI-MELTZER | 179 | 755.8869 |

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
|        9 | 189 |
|       10 | 155 |

### Filter TARGET data for 2027_XX release

Filter TARGET dataset for osteosarcomas only, which we plan to download
and process for the next splice compendium release.

``` r
osteo_filtered_sra <- target_sra |>
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
paste0("There are ", nrow(osteo_filtered_sra), " accessions matching the filters.")
```

    [1] "There are 179 accessions matching the filters."

``` r
# estimate amount of space files will take up
file_space <- sum(osteo_filtered_sra$Bytes) / 1e12
paste0("About ", file_space, " terabytes of space will be taken up by file downloads.")
```

    [1] "About 0.755886850559 terabytes of space will be taken up by file downloads."

Summarize samples by sequencing center of origin. All osteosarcomas
originate from the same sequencing center.

``` r
osteo_filtered_sra |> dplyr::count(study_name, `Center Name`)
```

| study_name                | Center Name |   n |
|:--------------------------|:------------|----:|
| TARGET: Osteosarcoma (OS) | NCI-MELTZER | 179 |

## Export filtered TARGET accession files

``` r
readr::write_tsv(target_filtered_sra, file = file.path(target_accession_path))
readr::write_tsv(osteo_filtered_sra, file = file.path(target_osteo_accession_path ))
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
    [13] htmltools_0.5.9   readxl_1.5.0      hms_1.1.4         rmarkdown_2.31   
    [17] cellranger_1.1.0  evaluate_1.0.5    tibble_3.3.1      tzdb_0.5.0       
    [21] fastmap_1.2.0     yaml_2.3.12       lifecycle_1.0.5   compiler_4.4.3   
    [25] dplyr_1.2.1       pkgconfig_2.0.3   here_1.0.2        rstudioapi_0.18.0
    [29] digest_0.6.39     R6_2.6.1          tidyselect_1.2.1  readr_2.2.0      
    [33] parallel_4.4.3    vroom_1.7.1       pillar_1.11.1     magrittr_2.0.5   
    [37] withr_3.0.3       tools_4.4.3       bit64_4.8.2      
