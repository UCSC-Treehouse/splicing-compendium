# Explore Treehouse PolyA Compendium
Cindy Liang (celiang@ucsc.edu)
2026-10-22

## Introduction

We are interested in adding more datasets to the Treehouse Splice
Compendium prior to submitting the manuscript to a journal. Currently,
the priorities are to include samples in the Treehouse PolyA compendia,
with healthy comparators from dGTEx.

This notebook reads in the clinical metadata from the v25 Treehouse
PolyA Compendium and filters for samples to download to process for
Splice Compendium v1.2. The input files to this notebook are:

**polya_clinical_file**: [Treehouse PolyA Compendium v25 clinical
metadata](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294351),
which contains clinical (disease, age at dx, sample source) information
associated with a unique `th_dataset_id`. This file has more samples
than recorded in `treehouse_files_on_aws`, because Treehouse does not
keep all of the TARGET files processed in the PolyA compendium on AWS.

**treehouse_files_on_aws**: An internal log of sequence files in
Treehouse storage. Fields consist of date, time, bytes, and path.
Treehouse IDs for associating files with those in the PolyA compendium
are within the values for `path`.

**TH_dataset_ID_by_collaborator_THIDs.tsv**: Internal Treehouse key of
what collaborators are associated with specific Treehouse dataset IDs

**41698_2025_852_MOESM2_ESM.xlsx**: Supplementary data table 1 from
[Vasquez et
al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11928651/#sec29)
listing Treehouse sample IDs of 33 samples used in study.

**JAMA_TH_IDs.tsv**: space-separated file of 144 sample IDs used in
[Vaske et
al. 2019](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2753519),
pasted from eTable 3b in the paper’s
[supplement](https://cdn.jamanetwork.com/ama/content_public/journal/jamanetworkopen/938206/zoi190533supp1_prod.pdf?Expires=1794613443&Signature=WLyr7rwB8MKSlrxifyxwoJ75UMi1ByrZUUi2X~k3JhloUTzZvJTDuBzkgS8mHOvixSNDp6kPS5ygqHaNk-EcwPNC6G9mG0VXY5ORDMVdrSlCGy2WujAvboi0lCObb3GEtNBXy0iuBOKPERC3FFQX7LI4xHPxm060Gzu3ColP91exQVd1w10XJYGA6cnn5GhqMCzpC~YBCPtalUlWqaGGmKaYzaUmDCWeTDN5W0~X~VsX8yZSHS8rbj18EciADgJdPawWGDI3ml0WGrgcgbD8V2tkG33t1DcE7gJwHJLogU9PfHYNvye8TOppM8lMbcZXSt40Npr7S~XfDFvv5XA4GQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

JAMA Treehouse IDs are obtained from JAMA_TH_IDs.tsv document through:

    awk 'NR > 1 {print $1}' JAMA_TH_IDs.tsv > JAMA_TH_ID_column.txt

**Jama_REDCap_Donor_Sample_20260730.tsv**: REDCap record of JAMA samples
in Treehouse’s internal database, with MEND QC information.

This notebook produces output files of accession IDs of Treehouse PolyA
Compendium v25 files, separated by dataset of origin. Accessions from
these files will be used as input for downstream downloading scripts.

Additionally, this notebook produces summary tables breaking down
individual contributors of origin for samples in the Treehouse PolyA
Compendium, and samples from the Vaske 2019 and Vasquez 2025 papers.

## Setup

Files and directories

``` r
## Directories ##
repo_root <- rprojroot::find_root(rprojroot::is_git_root)
metadata_dir <- file.path(repo_root, "metadata")

## Files ##
# Treehouse polyA compendium with disease type of samples and treehouse ID
polya_clinical_file <- file.path(metadata_dir, "GSE294351_clinical_Treehouse-Tumor-Compendium-25.01-PolyA_20250131v1.tsv")
# record of files on AWS, with file type information and file size.
# to be associated with polya_clinical_file by treehouse ID
treehouse_aws_summary_file <- file.path(metadata_dir, "archived_files_on_archive-treehouse-ucsc-edu_primary.txt")
# table associating Treehouse IDs with individual partner contributor source
thid_partner_file <- file.path(metadata_dir, "TH_dataset_ID_by_collaborator_THIDs.tsv")
jama_sample_id_file <- file.path(metadata_dir, "JAMA_TH_ID_column.txt")
vasquez_sample_id_file <- file.path(metadata_dir, "41698_2025_852_MOESM2_ESM.xlsx")

jama_redcap_file <- file.path(metadata_dir, "jama_REDCap_Donor_Sample_20260730.tsv")

# output accessions 
# list of SRA accessions fof polyA compendium files
sra_accessions <- file.path(metadata_dir, "sra_polya_compendium_accessions.txt")
ega_accessions <- file.path(metadata_dir, "ega_polya_compendium_accessions.txt")
cbtn_accessions <- file.path(metadata_dir, "cbtn_polya_compendium_accessions.txt")
# these files should only have treehouse IDs rather than "accession" IDs
treehouse_collaborator_ids <-file.path(metadata_dir, "collaborator_polya_compendium_ids.txt")
publication_only_th_collaborator_ids <- file.path(metadata_dir, "pub_only_collaborator_ids.txt")
jama_samples_by_collaborator <- file.path(metadata_dir, "jama_samples_by_collaborator_count.txt")
```

Read in files

``` r
## Read in files ##
# age_at_dx have "unknown" when age is unknown so cannot set coltype as double
polya_metadata_df <- readr::read_tsv(polya_clinical_file,
                                     col_types = readr::cols(
                                       .default = "c"
                                     ))

# internal Treehouse AWS log of files in storage, alongside file size
# this file is separated by spaces
treehouse_aws_df <- read.table(
  treehouse_aws_summary_file,
  sep = "", 
  header = FALSE, 
  col.names = c(
    "date",
    "time",
    "bytes",
    "path"),
  colClasses = c(
    bytes = "double"
  )
  )

# file associating TH dataset ID to collaboration partner
th_dataset_id_partner_df <- readr::read_tsv(thid_partner_file)
```

    Rows: 65 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): TH_ID, partner

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# vector of treehouse sample IDs in JAMA / Vaske 2019 publiation
jama_sample_vector <- readLines(
  jama_sample_id_file
) |>
  unique()

# vector of treehouse sample IDs in CCK2 / Vasquez 2025 publiation
vasquez_sample_vector <- readxl::read_xlsx(vasquez_sample_id_file) |>
  dplyr::pull(`Patient ID`)

# list of JAMA samples on redcap 
jama_redcap_data <- readr::read_tsv(
  jama_redcap_file,
  col_types = readr::cols_only(
    th_donor_id = "c", 
    site_id = "c",
    th_sampleid = "c",
    qc = "n",
    umend = "n",
    # polyA samples are denoted with a value of "10" here
    inv10_ribod  = "n")
  )
```

Note: JAMA sample IDs are sample IDs, while CKCC2 IDs are patient IDs,
which can be derived from the sample ID.

``` r
jama_sample_vector |> head()
```

    [1] "TH01_0053_S01" "TH01_0054_S01" "TH01_0055_S01" "TH01_0061_S01"
    [5] "TH01_0062_S01" "TH01_0063_S01"

``` r
vasquez_sample_vector |> head()
```

    [1] "TH34_1149" "TH34_1150" "TH34_1162" "TH34_1163" "TH34_1179" "TH34_1238"

Out of all JAMA samples, check how many have no QC data (assume they are
not in polyA compendium for this reason)

samples pass MEND QC if they have \>= 10M MEND reads : “The RNA-Seq data
must contain at least 10 million mapped exonic, non-duplicate (MEND)
reads.” [(Beale
2025)](https://www.nature.com/articles/s41597-025-05376-z)

``` r
jama_redcap_data |>
  dplyr::mutate(
    pass_mend =
    dplyr::case_when(
      umend > 10000000 ~ TRUE,
      umend < 10000000 ~ FALSE
    )
    )|>
  dplyr::summarise(.by = c(qc, pass_mend, inv10_ribod),
                   n = dplyr::n())
```

|  qc | pass_mend | inv10_ribod |   n |
|----:|:----------|------------:|----:|
|   5 | TRUE      |          10 |  88 |
|   5 | TRUE      |           9 |  10 |
|   5 | TRUE      |           1 |  15 |
|   6 | NA        |          10 |  31 |

88/144 samples pass MEND QC and are polyA, so I should expect to see 88
JAMA samples in the polyA compendium.

Extract MEND-QC-pass samples in a vector

``` r
jamna_mend_qc_pass_ids <-  jama_redcap_data |>
  dplyr::mutate(
    pass_mend =
      dplyr::case_when(
        umend > 10000000 ~ TRUE,
        umend < 10000000 ~ FALSE
      )
  )|>
  dplyr::filter(
    pass_mend == TRUE
  ) |>
  dplyr::pull(th_sampleid)
```

## Clean and merge Treehouse sequence metadata

Split “path” column in AWS file into filetype and Treehouse ID fields

Filetypes in Treehouse storage are: .sra, .sra.vdbcache,
sequence.txt.gz, markdown.txt, .fastq, .fq.gz, and .bam. .bam files
processed by Treehouse should contained unmapped reads.

``` r
cleaned_aws_df <- treehouse_aws_df  |>
  dplyr::mutate(
    filename = basename(path),
    suffix = dplyr::case_when(
      stringr::str_detect(filename, "fastq") ~ "fastq",
      stringr::str_detect(filename, "fq.gz") ~ "fastq",
      stringr::str_detect(filename, "bam") ~ "bam",
      stringr::str_detect(filename, "sra") ~ "sra",
      stringr::str_detect(filename, "swp") ~ "swp",
      stringr::str_detect(filename, "txt") ~ "txt"
    ),
    th_dataset_id = stringr::str_split_i(path, "/", 3)
  )
```

Merge treehouse polyA compendium metadata with aws file information to
associate file size with compendium samples. This information will be
used to inform download batching strategy and estimate AWS file transfer
costs.

``` r
merged_compendium_aws_df <- polya_metadata_df |> 
  dplyr::left_join(cleaned_aws_df,
                   by = "th_dataset_id") |>
  # select relevant columns for downstream filtering and summary of samples
  dplyr::select(
    th_dataset_id,
    study_id,
    age_at_dx,
    disease,
    organism,
    sex,
    study_accession,
    source_name,
    bytes,
    suffix
  ) |>
  # turn bytes column back into double for summing later on
  dplyr::mutate(
    bytes = as.double(bytes),
    # derive treehouse patient id from sample id by splitting on the second underscore
    th_patient_id = stringr::str_split_i(th_dataset_id, "(?<=^[^_]{1,100}_[^_]{0,100})_", 1),
    # label IDs based off if they come from JAMA or CKCC2 paper
    publication_origin = dplyr::case_when(
      th_patient_id %in% vasquez_sample_vector ~ "ckcc2_vasquez_2025",
      th_dataset_id %in% jama_sample_vector ~ "jama_vaske_2019"
    )
  )
```

Filter samples prior to summarizing. For now, I will exclude samples
from TARGET and TCGA, since these are datasets not unique to Treehouse.
TARGET is in Splice Compendium v1, and TCGA is primarily composed of
adult samples. We could consider downloading more of TARGET
(e.g. osteosarcomas) that were not downloaded previously for Splice
Compendium v1.

At this stage, I will not filter for pedaya, as some samples have
inconsistently annotated pedaya and/or age_at_dx (some samples with
pedaya age_at_dx are labeled “unknown” pedaya and some samples have
“unknown” age_at_dx but I assume are pediatric based off the study area
of collaborators where samples are from

``` r
filtered_compendium_aws_df <- merged_compendium_aws_df |>
  dplyr::filter(
    # exclude TARGET and TCGA for now
    ! source_name %in% c("TARGET", "TCGA"),
    # sanity check - some aws files are dummies with 0 bytes
    bytes > 0
  )
```

## Summarize filtered data

937 files pass the above filters. The below summary tables will provide
information on the sex/disease type breakdowns of the filtered samples,
how many files originate from what data source, and the file size of
each file categorized by data source.

sex/disease breakdown of filtered polyA compendium files

``` r
filtered_compendium_aws_df |>
  dplyr::summarise(
    .by = c(disease, sex),
    n_samples = dplyr::n()
  ) |>
  dplyr::arrange(desc(n_samples)) |>
  # there are 101 rows, so print head for ease of viewing
  head()
```

<div id="tbl-polya_v25_top_samples">

Table 1: Top 6 sex/disease categories with the most files

<div class="cell-output-display">

| disease         | sex    | n_samples |
|:----------------|:-------|----------:|
| osteosarcoma    | male   |       147 |
| medulloblastoma | male   |       136 |
| osteosarcoma    | female |       129 |
| medulloblastoma | female |        69 |
| Ewing sarcoma   | male   |        59 |
| Ewing sarcoma   | female |        51 |

</div>

</div>

Much of the compendium is composed of unusual sarcomas (n = 1 or n \<
5). The majority of the compendium after these filters consist of
medulloblastoma, osteosarcoma, glioma, and Ewing’s sarcoma.

Print summary of files by data source and total bytes

There are some cases where the same Treehouse sample ID appears multiple
times, but the bytes are different. These are cases where multiple fastq
pairs are associated with the same sample.  
The strategy for processing these files is to concatenate the fastqs
into a single read1 / read2 pair, prior to trimming and alignment.

``` r
filtered_compendium_aws_df |>
  dplyr::summarise(
    .by = source_name,
    n_datasets = length(unique(th_dataset_id)),
    total_bytes = sum(bytes)
  )
```

<div id="tbl-polya_v25_source_filesize_summary">

Table 2: Summary of polyA v25 compendium files, by data source and avg
bytes per file

<div class="cell-output-display">

| source_name           | n_datasets |  total_bytes |
|:----------------------|-----------:|-------------:|
| EGA                   |         69 | 1.831246e+12 |
| unavailable           |        193 | 2.109279e+12 |
| Sequence Read Archive |        684 | 5.337788e+12 |
| CBTN                  |          9 | 1.188477e+11 |

</div>

</div>

~1.8TB of data should be accessible on SRA. Some of these files are
listed as .sra in glacier, so the estimated filesize may be smaller than
the zipped fastqs. This group represents data that I can try downloading
now, since we have sra-toolkit in our pixi environment.

“unavailable” indicates the data originate from individual
collaborators, such as collaborators. These files should be downloadable
on glacier, but the sequence files may need to be deposited in a
repository prior to publication of their splice values.

## Summarize sequence files originating from individual collaborators

Merge filtered polyA compendium dataframe with information mapping
treehouse ID to collaborator

``` r
# merge filtered compendium df with file info from aws with supplement labeling individual collaborator source of Treehouse IDs

collaborator_source_compendium_df <- filtered_compendium_aws_df |>
  dplyr::left_join(
    th_dataset_id_partner_df,
    by = c("study_id" = "TH_ID"))
```

Count number of files (multiple per patient) and samples from individual
collaborators

``` r
# summarize data by partner they are from
collaborator_source_compendium_df |>
  # filter for samples from collaborators only
  dplyr::filter(source_name == "unavailable") |>
  dplyr::summarise(
    .by = c(publication_origin, `partner`),
    n_files = dplyr::n(),
    n_samples = length(unique(th_dataset_id))
  )
```

<div id="tbl-polya_v25_individual_collaborator_files">

Table 3: Number of files in the polyA v25 compendium from each
individual collaborator

<div class="cell-output-display">

| publication_origin | partner | n_files | n_samples |
|:---|:---|---:|---:|
| ckcc2_vasquez_2025 | Samples associated with Stanford Clinical Registry | 136 | 36 |
| jama_vaske_2019 | Alejandro at Stanford (incl pdx) | 6 | 3 |
| NA | Alejandro at UCSF - July 2017 on, except for new samples from TH03 cases | 240 | 120 |
| NA | Nationwide Children’s | 40 | 10 |
| NA | Dr. Aru Narendran at Alberta Children’s Hospital (part of POETIC) | 16 | 5 |
| NA | University of Arkansas Medical Sciences (Analiz Rodriguez) | 8 | 1 |
| NA | Sameer Agnihotri University of Pittsburgh | 2 | 1 |
| NA | Samples associated with Stanford Clinical Registry | 26 | 13 |
| NA | Alejandro at Stanford (incl pdx) | 8 | 4 |

</div>

</div>

In total, we have sequences from 193 samples originating from individual
collaborators. 36/34 expected samples are from the CKCC2 study, while
only 3/113 expected samples are from the JAMA study.

Count number of individual patients whose data come from individual
collaborators

``` r
# summarize number of unique patients whose samples are from collaborators
collaborator_source_compendium_df |>  
  # filter for samples from collaborators only
  dplyr::filter(source_name == "unavailable") |>
  dplyr::summarise(
    .by = c(publication_origin, `partner`),
    # note: some samples come from the same patient
    n_patients = length(unique(th_patient_id))
  )
```

<div id="tbl-polya_v25_individual_collaborator_patients">

Table 4: Number of unique patients in the polyA v25 compendium from each
individual collaborator

<div class="cell-output-display">

| publication_origin | partner | n_patients |
|:---|:---|---:|
| ckcc2_vasquez_2025 | Samples associated with Stanford Clinical Registry | 31 |
| jama_vaske_2019 | Alejandro at Stanford (incl pdx) | 3 |
| NA | Alejandro at UCSF - July 2017 on, except for new samples from TH03 cases | 97 |
| NA | Nationwide Children’s | 9 |
| NA | Dr. Aru Narendran at Alberta Children’s Hospital (part of POETIC) | 4 |
| NA | University of Arkansas Medical Sciences (Analiz Rodriguez) | 1 |
| NA | Sameer Agnihotri University of Pittsburgh | 1 |
| NA | Samples associated with Stanford Clinical Registry | 12 |
| NA | Alejandro at Stanford (incl pdx) | 2 |

</div>

</div>

Create a table of collaborator origin of JAMA samples not in polyA
compendium. These are likely collaborators we need to reach out to to
request for data to be deposited.

``` r
jama_samples_by_collaborator_df <- dplyr::tibble(jama_sample_vector) |>
  dplyr::mutate(
    TH_ID = stringr::str_split_i(jama_sample_vector, "_", 1)
  ) |>
  dplyr::left_join(
    th_dataset_id_partner_df
  ) |>
  dplyr::summarise(
    .by = c(partner, TH_ID),
    n_samples = dplyr::n()
  )
```

    Joining with `by = join_by(TH_ID)`

``` r
jama_samples_by_collaborator_df
```

<div id="tbl-jama_by_collaborator_source">

Table 5: Number of JAMA samples broken down by collaborator

<div class="cell-output-display">

| partner                          | TH_ID | n_samples |
|:---------------------------------|:------|----------:|
| CHOC and CHOC/Nant               | TH01  |        23 |
| PNOC UCSF (Mueller)              | TH02  |        24 |
| Alejandro at Stanford (incl pdx) | TH03  |        66 |
| BC-POG                           | TH06  |        31 |

</div>

</div>

Most of these samples are not in the PolyA compendium for a mixture of
reasons

CHOC are all ribo-deplete and BC-POC fastqs cannot be accessed

Double check from REDCap data which of these samples are ribo-D
(potentially a third table) – an additional column to this summary as
well as another column on how many samples pass MEND QC

## Output a list of file accessions per data source

My assumption is that I will need a separate download script for each
different data source (or at least go through a different system to
download the data).

``` r
# prep filtered df for exporting file accessions
for_export_compendium_df <- collaborator_source_compendium_df |>
  dplyr::select(
    th_dataset_id,
    study_id,
    th_patient_id,
    publication_origin,
    partner,
    study_accession,
    source_name
  )

# export accessions lists to output
sra_df <- for_export_compendium_df |>
  dplyr::filter(
    source_name == "Sequence Read Archive"
  ) |>
  readr::write_tsv(sra_accessions)

for_export_compendium_df |>
  dplyr::filter(
    source_name == "EGA"
  ) |>
  readr::write_tsv(ega_accessions)

for_export_compendium_df |>
  dplyr::filter(
    source_name == "CBTN"
  ) |>
  readr::write_tsv(cbtn_accessions)

for_export_compendium_df |>
  dplyr::filter(
    source_name == "unavailable"
  ) |>
  dplyr::arrange(publication_origin) |>
  readr::write_tsv(treehouse_collaborator_ids)

for_export_compendium_df |>
  dplyr::filter(
    source_name == "unavailable",
    !is.na(publication_origin)
  ) |>
  dplyr::arrange(publication_origin) |>
  readr::write_tsv(publication_only_th_collaborator_ids)

jama_samples_by_collaborator_df |>
  readr::write_tsv(jama_samples_by_collaborator)
```

Check more study accessions from piolya metadata for any other accession
outside of project ID for SRA downloads
