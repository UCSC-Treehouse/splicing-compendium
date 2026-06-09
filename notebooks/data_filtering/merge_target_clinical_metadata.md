# Add clinical metadata to TARGET compendium v1 sample sheet
Cindy Liang (celiang@ucsc.edu)
2026-06-09

## Introduction

This notebook adds clinical metadata to the TARGET accessions sheet.
Clinical metadata is obtained from the UCSC Xena browser:
https://xenabrowser.net/datapages/?hub=https://gdc.xenahubs.net:443

We plan to include the following metadata fields in the accessions
sheet:

- age at dx

- race/ethnicity

- sex

- survival

- Some columns with a single number or value per column (values that can
  be computed on easily)

## Setup

### Read in directories and files

``` r
## Define paths to directories ##

# find the project directory
base_dir <- here::here()

# define directory paths
# accession metadata dir
metadata_dir <- file.path(base_dir, "metadata")
# clinical TARGET metadata directory
clinical_target_dir <- file.path(metadata_dir, "clinical", "target")
# filtered accessions directory
gtex_target_metadata_dir <- file.path(metadata_dir, "filter_target_gtex")
# sample sheet dir of samples that were downloaded and processed
config_dir <- file.path(base_dir, "config")

## Define paths to input files ##

# full TARGET dataset SRA accessions and metadata
target_sra_path <- file.path(gtex_target_metadata_dir, "SraRunTable-TARGET.csv")
# final list of accessions that were downloaded and processed, including pilot accessions
compendium_samples_path <- file.path(config_dir, "sample_sheet.tsv")

# TARGET clinical metadata files
# Define list of TARGET clinical metadata files
target_clinical_files <- c(
  all_p1_clinical = "TARGET-ALL-P1.clinical.tsv",
  all_p2_clinical = "TARGET-ALL-P2.clinical.tsv",
  aml_clinical = "TARGET-AML.clinical.tsv",
  nbl_clinical = "TARGET-NBL.clinical.tsv",
  rt_clinical = "TARGET-RT.clinical.tsv",
  wt_clinical = "TARGET-WT.clinical.tsv"
)

target_survival_files <- c(
  all_p1_survival = "TARGET-ALL-P1.survival.tsv",
  all_p2_survival = "TARGET-ALL-P2.survival.tsv",
  aml_survival = "TARGET-AML.survival.tsv",
  nbl_survival = "TARGET-NBL.survival.tsv",
  rt_survival = "TARGET-RT.survival.tsv",
  wt_survival = "TARGET-WT.survival.tsv"
)

# construct TARGET clinical metadata paths
target_clinical_metadata_paths <-file.path(clinical_target_dir, target_clinical_files)
names(target_clinical_metadata_paths) <- names(target_clinical_files)
target_survival_metadata_paths <- file.path(clinical_target_dir, target_survival_files)
names(target_survival_files) <- names(target_survival_files)

## Output ##
# define path to output file
target_accessions_with_clinical_metadata <- file.path(gtex_target_metadata_dir, "target_accessions_clinical.tsv")
```

Read in files

``` r
## Sample sheets ##

# SRA accessions tables - target_sra "Run" column is to be associated with compendium_samples "sample" column
# target SRA table
target_sra <- readr::read_csv(target_sra_path, col_types = readr::cols(.default = "c"))
# sample sheet of what was downloaded and processed
compendium_samples <- readr::read_tsv(compendium_samples_path, col_types = readr::cols( .default = "c"))

## Clinical metadata ##

# read in metadata paths
target_clinical_metadata <- target_clinical_metadata_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c"))
  }) |> 
  purrr::list_rbind()

target_survival_metadata <- target_survival_metadata_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default="c"))
  }) |>
  purrr::list_rbind()
```

The “sample” column in the survival and clinical metadata tables have
entries that look like: TARGET-10-PAMDRM-09A

The “sample” column in the metadata table correspond to the
“biospecimen_repository_sample_id” column in the target accessions
table, which looks like: TARGET-10-PAMDRM-09A-01D

“biospecimen_repository_sample_id” needs to be cleaned so the last
substring separated by “-” is removed

Combine metadata with accessions sheet and filter target SRA metadata
for those corresponding to samples that were downloaded and processed

``` r
# combine target clinical and survival metadata
target_metadata <- dplyr::full_join(
  target_clinical_metadata,
  target_survival_metadata,
  by = "sample",
  relationship = "one-to-one"
)

# make a list of target accessions that were downloaded and processed
target_downloaded <- compendium_samples |>
  dplyr::filter(group == "target") |>
  dplyr::pull(sample)

# filter target accessions for accessions that are in the sample sheet
target_sra_metadata <- target_sra |>
  dplyr::filter(Run %in% target_downloaded) |>
  # select for relevant columns
  # dplyr::select(
  #   Run, biospecimen_repository, `DATASTORE filetype`, body_site,
  #   study_name, analyte_type, `Assay Type`, `SRA Study`,
  #   BioProject, BioSample, LibraryLayout, LibrarySelection, `Center Name`,
  #   biospecimen_repository_sample_id
  # ) |>
  # derive target sample ID from sample ID 
  dplyr::mutate(
    # remove the hyphen followed by characters that are not hyphens that are at the end of the string
    target_sample_id = sub("-[^-]+$", "", biospecimen_repository_sample_id)
  )

# Check target metadata contents
colnames(target_sra_metadata)
```

      [1] "Run"                              "analyte_type"                    
      [3] "Assay Type"                       "BioProject"                      
      [5] "BioSample"                        "biospecimen_repository"          
      [7] "biospecimen_repository_sample_id" "Bytes"                           
      [9] "Center Name"                      "Consent code"                    
     [11] "Consent"                          "DATASTORE filetype"              
     [13] "DATASTORE provider"               "DATASTORE region"                
     [15] "Experiment"                       "dbGaP accession"                 
     [17] "Instrument"                       "Library Name"                    
     [19] "LibraryLayout"                    "LibrarySelection"                
     [21] "LibrarySource"                    "Organism"                        
     [23] "Platform"                         "ReleaseDate"                     
     [25] "Sample Name"                      "SRA Study"                       
     [27] "study_design"                     "study_name"                      
     [29] "submitted_subject_id"             "Bases"                           
     [31] "AvgSpotLen"                       "histological_type"               
     [33] "body_site"                        "sex"                             
     [35] "AssemblyName"                     "create_date"                     
     [37] "version"                          "Is_Tumor"                        
     [39] "molecular_data_type"              "alignment_software (exp)"        
     [41] "gap_parent_phs"                   "AvgReadLength (run)"             
     [43] "coverage (run)"                   "lsid (exp)"                      
     [45] "lsid (run)"                       "run_barcode (run)"               
     [47] "project (exp)"                    "project (run)"                   
     [49] "work_request (exp)"               "work_request (run)"              
     [51] "run_name (run)"                   "read_group_platform_unit (run)"  
     [53] "molecular_indexing_scheme (run)"  "flowcell_barcode (run)"          
     [55] "read_group_id (run)"              "gssr_id (exp)"                   
     [57] "gssr_id (run)"                    "material_type (exp)"             
     [59] "root_sample_id (exp)"             "root_sample_id (run)"            
     [61] "sample_id (exp)"                  "Sample_ID (run)"                 
     [63] "analysis_type (exp)"              "analysis_type (run)"             
     [65] "library_type (exp)"               "library_type (run)"              
     [67] "lane (run)"                       "instrument_name (run)"           
     [69] "sample_type (exp)"                "research_project (exp)"          
     [71] "research_project (run)"           "RUN (run)"                       
     [73] "product_order (exp)"              "product_order (run)"             
     [75] "data_type (run)"                  "data_type (exp)"                 
     [77] "product_part_number (exp)"        "product_part_number (run)"       
     [79] "Sequencing_date (run)"            "subject_is_affected"             
     [81] "study_disease"                    "CompleteGenomics_sample_ID (run)"
     [83] "sample_barcode (exp)"             "sample_barcode (run)"            
     [85] "Library_Construction_batch (exp)" "target_set (exp)"                
     [87] "is_technical_control"             "Assembly (run)"                  
     [89] "instrument_model (run)"           "primary_disease (exp)"           
     [91] "RUNS (run)"                       "Sequencing_Dates (run)"          
     [93] "BI_GSSR_sample_ID (exp)"          "BI_GSSR_sample_ID (run)"         
     [95] "BI_GSSR_sample_LSID (exp)"        "BI_GSSR_sample_LSID (run)"       
     [97] "BI_project_name (exp)"            "BI_project_name (run)"           
     [99] "BI_run_barcode (run)"             "BI_run_name (run)"               
    [101] "BI_target_set (exp)"              "BI_work_request_ID (exp)"        
    [103] "BI_work_request_ID (run)"         "FLAG (exp)"                      
    [105] "secondary_accessions (run)"       "secondary_accessions (exp)"      
    [107] "aggregation_project (exp)"        "aggregation_project (run)"       
    [109] "library (exp)"                    "library (run)"                   
    [111] "molecular_idx_scheme (run)"       "work_request_or_pdo (exp)"       
    [113] "work_request_or_pdo (run)"        "bait_set (run)"                  
    [115] "Illumina_HiSeq_1000 (run)"        "rg_platform_unit (run)"          
    [117] "dangling_references (run)"        "missing_file (run)"              
    [119] "rg_platform_unit_lib (run)"       "TRUNCATED_DATA (exp)"            
    [121] "target_sample_id"                

There are more fields in the metadata than we planned to include, so
some discussion on what may or may not be useful as a group will be
helpful. For now, I focus on excluding fields that are duplicates of
other columns, or fields I don’t think are as useful for clinical
analysis (like vital_status, program project name).

## Merge TARGET accessions with clinical metadata

Row 2421 of the clinical metadata corresponds to one AML patient with
two RNA-seq samples in the TARGET dbGaP accessions metadata, but only
one entry in the clinical metadata. Some clinical metadata columns have
a list of values for this sample, which I think indicates it contains
information from different samples.

``` r
# print row of target_metadata that has multiple rows which match to target-sra_metadata
target_metadata[2421,] # sample column value is "TARGET-20-PARTXH-09A"
```

| sample | id | disease_type | case_id | submitter_id | primary_site | cause_of_death.demographic | race.demographic | gender.demographic | ethnicity.demographic | vital_status.demographic | age_at_index.demographic | days_to_birth.demographic | age_is_obfuscated.demographic | days_to_death.demographic | primary_site.project | project_id.project | disease_type.project | name.project | name.program.project | entity_submitter_id.annotations | notes.annotations | submitter_id.annotations | classification.annotations | entity_id.annotations | created_datetime.annotations | annotation_id.annotations | entity_type.annotations | updated_datetime.annotations | case_id.annotations | state.annotations | category.annotations | status.annotations | case_submitter_id.annotations | tissue_or_organ_of_origin.diagnoses | age_at_diagnosis.diagnoses | morphology.diagnoses | classification_of_tumor.diagnoses | icd_10_code.diagnoses | days_to_diagnosis.diagnoses | primary_diagnosis.diagnoses | year_of_diagnosis.diagnoses | diagnosis_is_primary_disease.diagnoses | site_of_resection_or_biopsy.diagnoses | age_at_earliest_diagnosis.diagnoses.xena_derived | age_at_earliest_diagnosis_in_years.diagnoses.xena_derived | protocol_identifier.treatments.diagnoses | updated_datetime.treatments.diagnoses | treatment_id.treatments.diagnoses | submitter_id.treatments.diagnoses | state.treatments.diagnoses | treatment_or_therapy.treatments.diagnoses | created_datetime.treatments.diagnoses | sample_type_id.samples | tumor_descriptor.samples | sample_id.samples | sample_type.samples | tumor_code.samples | preservation_method.samples | freezing_method.samples | tumor_code_id.samples | oct_embedded.samples | specimen_type.samples | tissue_type.samples | last_known_disease_status.diagnoses | days_to_last_follow_up.diagnoses | tumor_grade.diagnoses | progression_or_recurrence.diagnoses | sites_of_involvement.diagnoses | timepoint_category.treatments.diagnoses | treatment_type.treatments.diagnoses | course_number.treatments.diagnoses | reason_treatment_ended.treatments.diagnoses | treatment_outcome.treatments.diagnoses | therapeutic_agents.treatments.diagnoses | annotations.samples | inss_stage.diagnoses | cog_neuroblastoma_risk_group.diagnoses | mitosis_karyorrhexis_index.diagnoses | inpc_grade.diagnoses | pathology_detail_id.pathology_details.diagnoses | updated_datetime.pathology_details.diagnoses | submitter_id.pathology_details.diagnoses | state.pathology_details.diagnoses | created_datetime.pathology_details.diagnoses | necrosis_percent.pathology_details.diagnoses | percent_tumor_nuclei.pathology_details.diagnoses | pediatric_kidney_staging.diagnoses | wilms_tumor_histologic_subtype.diagnoses | OS.time | OS | \_PATIENT |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| TARGET-20-PARTXH-09A | 51af485c-dfcb-5e00-a22a-0d7d539892ad | Myeloid Leukemias | 51af485c-dfcb-5e00-a22a-0d7d539892ad | TARGET-20-PARTXH | Hematopoietic and reticuloendothelial systems | NA | white | male | hispanic or latino | Alive | 9.0 | -3386.0 | False | NA | \[‘Unknown’, ‘Hematopoietic and reticuloendothelial systems’\] | TARGET-AML | \[‘Not Applicable’, ‘Myeloid Leukemias’\] | Acute Myeloid Leukemia | TARGET | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | Bone marrow | 3386.0 | 9861/3 | primary | C92.0 | 0.0 | Acute myeloid leukemia, NOS | 2008 | True | Not Reported | 3386.0 | 9.276712328767124 | \[’‘,’‘,’‘, ’AAML0531’, ’’\] | \[‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’\] | \[‘3cd0ce34-fb69-4c9e-a1d3-fd009f016383’, ‘bc0cf9b8-33e9-4f20-a89d-90bff6bdedb0’, ‘bd7050b4-cab7-4e40-a690-2330751321f9’, ‘cf1b0cf4-0c2a-4191-8e98-62c2a6ef6ae5’, ‘f89f75cf-244f-49be-b3cb-830d798fe351’\] | \[‘TARGET-20-PARTXH_treatment3’, ‘TARGET-20-PARTXH_treatment2’, ‘TARGET-20-PARTXH_treatment5’, ‘TARGET-20-PARTXH_treatment’, ‘TARGET-20-PARTXH_treatment4’\] | \[‘released’, ‘released’, ‘released’, ‘released’, ‘released’\] | \[‘yes’, ‘no’, ‘yes’, ‘yes’, ‘yes’\] | \[‘2023-07-20T19:06:10.425901-05:00’, ‘2023-07-20T19:06:10.425901-05:00’, ‘2023-07-20T19:06:10.425901-05:00’, ‘2023-07-20T19:06:10.425901-05:00’, ‘2023-07-20T19:06:10.425901-05:00’\] | 09 | Primary | aacb5151-785c-57ec-9129-419091dc2820 | Primary Blood Derived Cancer - Bone Marrow | Acute myeloid leukemia (AML) | Unknown | None | 20 | None | Bone Marrow NOS | Tumor | NA | NA | NA | NA | NA | \[‘End of Treatment Course’, ‘First Complete Response’, ’‘,’‘, ’End of Treatment Course’\] | \[’‘, ’Stem Cell Transplantation, NOS’, ‘Pharmaceutical Therapy, NOS’, ’‘,’’\] | \[‘1.0’, ’‘,’‘,’‘, ’2.0’\] | \[‘Course of Therapy Completed’, ’‘,’‘,’‘, ’Course of Therapy Completed’\] | \[‘Complete Response’, ’‘,’‘,’‘, ’Complete Response’\] | \[’‘,’‘, ’Gemtuzumab Ozogamicin’, ’‘,’’\] | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

``` r
# print rows corresponding to sample TARGET-20-PARTXH-09A in the target_sra_metadata sheet to check if multiple rows match to the clinical metadata for explainable reasons 
target_sra_metadata |> dplyr::filter(
  target_sample_id == "TARGET-20-PARTXH-09A"
)
```

| Run | analyte_type | Assay Type | BioProject | BioSample | biospecimen_repository | biospecimen_repository_sample_id | Bytes | Center Name | Consent code | Consent | DATASTORE filetype | DATASTORE provider | DATASTORE region | Experiment | dbGaP accession | Instrument | Library Name | LibraryLayout | LibrarySelection | LibrarySource | Organism | Platform | ReleaseDate | Sample Name | SRA Study | study_design | study_name | submitted_subject_id | Bases | AvgSpotLen | histological_type | body_site | sex | AssemblyName | create_date | version | Is_Tumor | molecular_data_type | alignment_software (exp) | gap_parent_phs | AvgReadLength (run) | coverage (run) | lsid (exp) | lsid (run) | run_barcode (run) | project (exp) | project (run) | work_request (exp) | work_request (run) | run_name (run) | read_group_platform_unit (run) | molecular_indexing_scheme (run) | flowcell_barcode (run) | read_group_id (run) | gssr_id (exp) | gssr_id (run) | material_type (exp) | root_sample_id (exp) | root_sample_id (run) | sample_id (exp) | Sample_ID (run) | analysis_type (exp) | analysis_type (run) | library_type (exp) | library_type (run) | lane (run) | instrument_name (run) | sample_type (exp) | research_project (exp) | research_project (run) | RUN (run) | product_order (exp) | product_order (run) | data_type (run) | data_type (exp) | product_part_number (exp) | product_part_number (run) | Sequencing_date (run) | subject_is_affected | study_disease | CompleteGenomics_sample_ID (run) | sample_barcode (exp) | sample_barcode (run) | Library_Construction_batch (exp) | target_set (exp) | is_technical_control | Assembly (run) | instrument_model (run) | primary_disease (exp) | RUNS (run) | Sequencing_Dates (run) | BI_GSSR_sample_ID (exp) | BI_GSSR_sample_ID (run) | BI_GSSR_sample_LSID (exp) | BI_GSSR_sample_LSID (run) | BI_project_name (exp) | BI_project_name (run) | BI_run_barcode (run) | BI_run_name (run) | BI_target_set (exp) | BI_work_request_ID (exp) | BI_work_request_ID (run) | FLAG (exp) | secondary_accessions (run) | secondary_accessions (exp) | aggregation_project (exp) | aggregation_project (run) | library (exp) | library (run) | molecular_idx_scheme (run) | work_request_or_pdo (exp) | work_request_or_pdo (run) | bait_set (run) | Illumina_HiSeq_1000 (run) | rg_platform_unit (run) | dangling_references (run) | missing_file (run) | rg_platform_unit_lib (run) | TRUNCATED_DATA (exp) | target_sample_id |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| SRR2083137 | RNA | RNA-Seq | PRJNA89525 | SAMN01778140 | NCI_TARGET | TARGET-20-PARTXH-09A-03R | 11784423639 | BCCAGSC | 1 | DS-PEDCR | sra,run.zq,fastq | ncbi,gs,s3 | ncbi.dbgap,gs.us-east1,s3.us-east-1 | SRX1077694 | phs000465 | Illumina HiSeq 2500 | A12663 | PAIRED | cDNA | TRANSCRIPTOMIC | Homo sapiens | ILLUMINA | 2015-06-30T00:00:00Z | TARGET-20-PARTXH-09A-03R | SRP012000 | Tumor vs. Matched-Normal | TARGET: Acute Myeloid Leukemia (AML) | TARGET-20-PARTXH | 24494158650 | 150 | AML | Primary Blood Derived Cancer - Bone Marrow | male | NA | 2015-06-30T16:14:00Z | 1 | Yes | NA | NA | phs000218 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 147590 | NA | NA | NA | NA | NA | NA | 2014-12-01T22:44:00Z | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | TARGET-20-PARTXH-09A |
| SRR2083139 | RNA | RNA-Seq | PRJNA89525 | SAMN01778022 | NCI_TARGET | TARGET-20-PARTXH-09A-02R | 9201211556 | BCCAGSC | 1 | DS-PEDCR | fastq,sra,run.zq | s3,ncbi,gs | s3.us-east-1,ncbi.dbgap,gs.us-east1 | SRX1077619 | phs000465 | Illumina HiSeq 2500 | A12616 | PAIRED | cDNA | TRANSCRIPTOMIC | Homo sapiens | ILLUMINA | 2015-06-30T00:00:00Z | TARGET-20-PARTXH-09A-02R | SRP012000 | Tumor vs. Matched-Normal | TARGET: Acute Myeloid Leukemia (AML) | TARGET-20-PARTXH | 18870628050 | 150 | AML | Primary Blood Derived Cancer - Bone Marrow | male | NA | 2015-06-30T15:56:00Z | 1 | Yes | NA | NA | phs000218 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 147863 | NA | NA | NA | NA | NA | NA | 2014-12-08T21:16:29Z | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | TARGET-20-PARTXH-09A |

Check - are there sample-specific entries in the clinical metadata?

``` r
target_metadata |> dplyr::filter(sample == "TARGET-20-PARTXH-09A") |>
  dplyr::select(primary_site.project,
                disease_type.project,
                protocol_identifier.treatments.diagnoses,
                updated_datetime.treatments.diagnoses,
                treatment_id.treatments.diagnoses,
                submitter_id.treatments.diagnoses,
                state.treatments.diagnoses,
                treatment_or_therapy.treatments.diagnoses,
                created_datetime.treatments.diagnoses
                )
```

| primary_site.project | disease_type.project | protocol_identifier.treatments.diagnoses | updated_datetime.treatments.diagnoses | treatment_id.treatments.diagnoses | submitter_id.treatments.diagnoses | state.treatments.diagnoses | treatment_or_therapy.treatments.diagnoses | created_datetime.treatments.diagnoses |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| \[‘Unknown’, ‘Hematopoietic and reticuloendothelial systems’\] | \[‘Not Applicable’, ‘Myeloid Leukemias’\] | \[’‘,’‘,’‘, ’AAML0531’, ’’\] | \[‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’\] | \[‘3cd0ce34-fb69-4c9e-a1d3-fd009f016383’, ‘bc0cf9b8-33e9-4f20-a89d-90bff6bdedb0’, ‘bd7050b4-cab7-4e40-a690-2330751321f9’, ‘cf1b0cf4-0c2a-4191-8e98-62c2a6ef6ae5’, ‘f89f75cf-244f-49be-b3cb-830d798fe351’\] | \[‘TARGET-20-PARTXH_treatment3’, ‘TARGET-20-PARTXH_treatment2’, ‘TARGET-20-PARTXH_treatment5’, ‘TARGET-20-PARTXH_treatment’, ‘TARGET-20-PARTXH_treatment4’\] | \[‘released’, ‘released’, ‘released’, ‘released’, ‘released’\] | \[‘yes’, ‘no’, ‘yes’, ‘yes’, ‘yes’\] | \[‘2023-07-20T19:06:10.425901-05:00’, ‘2023-07-20T19:06:10.425901-05:00’, ‘2023-07-20T19:06:10.425901-05:00’, ‘2023-07-20T19:06:10.425901-05:00’, ‘2023-07-20T19:06:10.425901-05:00’\] |

There are, but there is no way to associate which items in the
list-columns correspond to which samples. Additionally, some columns
like treatment_or_therapy.treatments.diagnoses have five entries,
indicating they correspond to samples from the same patient not in the
compendium. In this case, it makes sense to me to preserve both samples
from the same patient and drop metadata columns that have these list
values. So samples from the same patient will just have the same
clinical metadata.

### Merge TARGET accessions with clinical metadata without assuming one-to-one relationship

``` r
# merge target accessions with clinical metadata 
target_clinical_accessions <- dplyr::left_join(
  target_sra_metadata,
  target_metadata, 
  by = dplyr::join_by(
    target_sample_id == sample)
    )

# print columns
colnames(target_clinical_accessions)
```

      [1] "Run"                                                      
      [2] "analyte_type"                                             
      [3] "Assay Type"                                               
      [4] "BioProject"                                               
      [5] "BioSample"                                                
      [6] "biospecimen_repository"                                   
      [7] "biospecimen_repository_sample_id"                         
      [8] "Bytes"                                                    
      [9] "Center Name"                                              
     [10] "Consent code"                                             
     [11] "Consent"                                                  
     [12] "DATASTORE filetype"                                       
     [13] "DATASTORE provider"                                       
     [14] "DATASTORE region"                                         
     [15] "Experiment"                                               
     [16] "dbGaP accession"                                          
     [17] "Instrument"                                               
     [18] "Library Name"                                             
     [19] "LibraryLayout"                                            
     [20] "LibrarySelection"                                         
     [21] "LibrarySource"                                            
     [22] "Organism"                                                 
     [23] "Platform"                                                 
     [24] "ReleaseDate"                                              
     [25] "Sample Name"                                              
     [26] "SRA Study"                                                
     [27] "study_design"                                             
     [28] "study_name"                                               
     [29] "submitted_subject_id"                                     
     [30] "Bases"                                                    
     [31] "AvgSpotLen"                                               
     [32] "histological_type"                                        
     [33] "body_site"                                                
     [34] "sex"                                                      
     [35] "AssemblyName"                                             
     [36] "create_date"                                              
     [37] "version"                                                  
     [38] "Is_Tumor"                                                 
     [39] "molecular_data_type"                                      
     [40] "alignment_software (exp)"                                 
     [41] "gap_parent_phs"                                           
     [42] "AvgReadLength (run)"                                      
     [43] "coverage (run)"                                           
     [44] "lsid (exp)"                                               
     [45] "lsid (run)"                                               
     [46] "run_barcode (run)"                                        
     [47] "project (exp)"                                            
     [48] "project (run)"                                            
     [49] "work_request (exp)"                                       
     [50] "work_request (run)"                                       
     [51] "run_name (run)"                                           
     [52] "read_group_platform_unit (run)"                           
     [53] "molecular_indexing_scheme (run)"                          
     [54] "flowcell_barcode (run)"                                   
     [55] "read_group_id (run)"                                      
     [56] "gssr_id (exp)"                                            
     [57] "gssr_id (run)"                                            
     [58] "material_type (exp)"                                      
     [59] "root_sample_id (exp)"                                     
     [60] "root_sample_id (run)"                                     
     [61] "sample_id (exp)"                                          
     [62] "Sample_ID (run)"                                          
     [63] "analysis_type (exp)"                                      
     [64] "analysis_type (run)"                                      
     [65] "library_type (exp)"                                       
     [66] "library_type (run)"                                       
     [67] "lane (run)"                                               
     [68] "instrument_name (run)"                                    
     [69] "sample_type (exp)"                                        
     [70] "research_project (exp)"                                   
     [71] "research_project (run)"                                   
     [72] "RUN (run)"                                                
     [73] "product_order (exp)"                                      
     [74] "product_order (run)"                                      
     [75] "data_type (run)"                                          
     [76] "data_type (exp)"                                          
     [77] "product_part_number (exp)"                                
     [78] "product_part_number (run)"                                
     [79] "Sequencing_date (run)"                                    
     [80] "subject_is_affected"                                      
     [81] "study_disease"                                            
     [82] "CompleteGenomics_sample_ID (run)"                         
     [83] "sample_barcode (exp)"                                     
     [84] "sample_barcode (run)"                                     
     [85] "Library_Construction_batch (exp)"                         
     [86] "target_set (exp)"                                         
     [87] "is_technical_control"                                     
     [88] "Assembly (run)"                                           
     [89] "instrument_model (run)"                                   
     [90] "primary_disease (exp)"                                    
     [91] "RUNS (run)"                                               
     [92] "Sequencing_Dates (run)"                                   
     [93] "BI_GSSR_sample_ID (exp)"                                  
     [94] "BI_GSSR_sample_ID (run)"                                  
     [95] "BI_GSSR_sample_LSID (exp)"                                
     [96] "BI_GSSR_sample_LSID (run)"                                
     [97] "BI_project_name (exp)"                                    
     [98] "BI_project_name (run)"                                    
     [99] "BI_run_barcode (run)"                                     
    [100] "BI_run_name (run)"                                        
    [101] "BI_target_set (exp)"                                      
    [102] "BI_work_request_ID (exp)"                                 
    [103] "BI_work_request_ID (run)"                                 
    [104] "FLAG (exp)"                                               
    [105] "secondary_accessions (run)"                               
    [106] "secondary_accessions (exp)"                               
    [107] "aggregation_project (exp)"                                
    [108] "aggregation_project (run)"                                
    [109] "library (exp)"                                            
    [110] "library (run)"                                            
    [111] "molecular_idx_scheme (run)"                               
    [112] "work_request_or_pdo (exp)"                                
    [113] "work_request_or_pdo (run)"                                
    [114] "bait_set (run)"                                           
    [115] "Illumina_HiSeq_1000 (run)"                                
    [116] "rg_platform_unit (run)"                                   
    [117] "dangling_references (run)"                                
    [118] "missing_file (run)"                                       
    [119] "rg_platform_unit_lib (run)"                               
    [120] "TRUNCATED_DATA (exp)"                                     
    [121] "target_sample_id"                                         
    [122] "id"                                                       
    [123] "disease_type"                                             
    [124] "case_id"                                                  
    [125] "submitter_id"                                             
    [126] "primary_site"                                             
    [127] "cause_of_death.demographic"                               
    [128] "race.demographic"                                         
    [129] "gender.demographic"                                       
    [130] "ethnicity.demographic"                                    
    [131] "vital_status.demographic"                                 
    [132] "age_at_index.demographic"                                 
    [133] "days_to_birth.demographic"                                
    [134] "age_is_obfuscated.demographic"                            
    [135] "days_to_death.demographic"                                
    [136] "primary_site.project"                                     
    [137] "project_id.project"                                       
    [138] "disease_type.project"                                     
    [139] "name.project"                                             
    [140] "name.program.project"                                     
    [141] "entity_submitter_id.annotations"                          
    [142] "notes.annotations"                                        
    [143] "submitter_id.annotations"                                 
    [144] "classification.annotations"                               
    [145] "entity_id.annotations"                                    
    [146] "created_datetime.annotations"                             
    [147] "annotation_id.annotations"                                
    [148] "entity_type.annotations"                                  
    [149] "updated_datetime.annotations"                             
    [150] "case_id.annotations"                                      
    [151] "state.annotations"                                        
    [152] "category.annotations"                                     
    [153] "status.annotations"                                       
    [154] "case_submitter_id.annotations"                            
    [155] "tissue_or_organ_of_origin.diagnoses"                      
    [156] "age_at_diagnosis.diagnoses"                               
    [157] "morphology.diagnoses"                                     
    [158] "classification_of_tumor.diagnoses"                        
    [159] "icd_10_code.diagnoses"                                    
    [160] "days_to_diagnosis.diagnoses"                              
    [161] "primary_diagnosis.diagnoses"                              
    [162] "year_of_diagnosis.diagnoses"                              
    [163] "diagnosis_is_primary_disease.diagnoses"                   
    [164] "site_of_resection_or_biopsy.diagnoses"                    
    [165] "age_at_earliest_diagnosis.diagnoses.xena_derived"         
    [166] "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"
    [167] "protocol_identifier.treatments.diagnoses"                 
    [168] "updated_datetime.treatments.diagnoses"                    
    [169] "treatment_id.treatments.diagnoses"                        
    [170] "submitter_id.treatments.diagnoses"                        
    [171] "state.treatments.diagnoses"                               
    [172] "treatment_or_therapy.treatments.diagnoses"                
    [173] "created_datetime.treatments.diagnoses"                    
    [174] "sample_type_id.samples"                                   
    [175] "tumor_descriptor.samples"                                 
    [176] "sample_id.samples"                                        
    [177] "sample_type.samples"                                      
    [178] "tumor_code.samples"                                       
    [179] "preservation_method.samples"                              
    [180] "freezing_method.samples"                                  
    [181] "tumor_code_id.samples"                                    
    [182] "oct_embedded.samples"                                     
    [183] "specimen_type.samples"                                    
    [184] "tissue_type.samples"                                      
    [185] "last_known_disease_status.diagnoses"                      
    [186] "days_to_last_follow_up.diagnoses"                         
    [187] "tumor_grade.diagnoses"                                    
    [188] "progression_or_recurrence.diagnoses"                      
    [189] "sites_of_involvement.diagnoses"                           
    [190] "timepoint_category.treatments.diagnoses"                  
    [191] "treatment_type.treatments.diagnoses"                      
    [192] "course_number.treatments.diagnoses"                       
    [193] "reason_treatment_ended.treatments.diagnoses"              
    [194] "treatment_outcome.treatments.diagnoses"                   
    [195] "therapeutic_agents.treatments.diagnoses"                  
    [196] "annotations.samples"                                      
    [197] "inss_stage.diagnoses"                                     
    [198] "cog_neuroblastoma_risk_group.diagnoses"                   
    [199] "mitosis_karyorrhexis_index.diagnoses"                     
    [200] "inpc_grade.diagnoses"                                     
    [201] "pathology_detail_id.pathology_details.diagnoses"          
    [202] "updated_datetime.pathology_details.diagnoses"             
    [203] "submitter_id.pathology_details.diagnoses"                 
    [204] "state.pathology_details.diagnoses"                        
    [205] "created_datetime.pathology_details.diagnoses"             
    [206] "necrosis_percent.pathology_details.diagnoses"             
    [207] "percent_tumor_nuclei.pathology_details.diagnoses"         
    [208] "pediatric_kidney_staging.diagnoses"                       
    [209] "wilms_tumor_histologic_subtype.diagnoses"                 
    [210] "OS.time"                                                  
    [211] "OS"                                                       
    [212] "_PATIENT"                                                 

## Filter merged TARGET clinical metadata

Below is a breakdown of columns that are excluded from the final
clinical metadata:

### Columns excluded from the metadata because they contain redundant values found in other columns

The following columns are excluded because they are redundant with the
submitter_id column

``` r
target_clinical_accessions |>
  # showing ALL phase 1 because columns ending in "annotations" only have value for this cancer type
  dplyr::filter(study_name == "TARGET: Acute Lymphoblastic Leukemia (ALL) Pilot Phase 1") |>
  dplyr::select(submitter_id, 
                entity_submitter_id.annotations, 
                case_submitter_id.annotations,
                # submitter_id is redundant because it concatenates classification.annotations and case_submitter_id.annotations
                submitter_id.annotations,
                classification.annotations) |>
  head()
```

| submitter_id | entity_submitter_id.annotations | case_submitter_id.annotations | submitter_id.annotations | classification.annotations |
|:---|:---|:---|:---|:---|
| TARGET-10-PANNGL | TARGET-10-PANNGL | TARGET-10-PANNGL | TARGET-10-PANNGL_Notification | Notification |
| TARGET-10-PANSFD | TARGET-10-PANSFD | TARGET-10-PANSFD | TARGET-10-PANSFD_Notification | Notification |
| TARGET-10-PANEHF | NA | NA | NA | NA |

The following columns are excluded because they are redundant with the
disease_type column

``` r
target_clinical_accessions |>
  dplyr::select(disease_type, 
                disease_type.project) |>
  head()
```

| disease_type       | disease_type.project |
|:-------------------|:---------------------|
| Lymphoid Leukemias | Lymphoid Leukemias   |
| Lymphoid Leukemias | Lymphoid Leukemias   |
| Lymphoid Leukemias | Lymphoid Leukemias   |
| Lymphoid Leukemias | Lymphoid Leukemias   |
| Lymphoid Leukemias | Lymphoid Leukemias   |
| Lymphoid Leukemias | Lymphoid Leukemias   |

The following columns are excluded because they are redundant with the
primary_site column

``` r
target_clinical_accessions |>
  dplyr::select(primary_site, 
                primary_site.project) |>
  head()
```

| primary_site | primary_site.project |
|:---|:---|
| Hematopoietic and reticuloendothelial systems | Hematopoietic and reticuloendothelial systems |
| Hematopoietic and reticuloendothelial systems | Hematopoietic and reticuloendothelial systems |
| Hematopoietic and reticuloendothelial systems | Hematopoietic and reticuloendothelial systems |
| Hematopoietic and reticuloendothelial systems | Hematopoietic and reticuloendothelial systems |
| Hematopoietic and reticuloendothelial systems | Hematopoietic and reticuloendothelial systems |
| Hematopoietic and reticuloendothelial systems | Hematopoietic and reticuloendothelial systems |

The following columns are excluded because they are redundant with the
age_at_diagnosis.diagnoses column

``` r
target_clinical_accessions |>
  dplyr::select(age_at_diagnosis.diagnoses,
                age_at_earliest_diagnosis.diagnoses.xena_derived) |>
  head()
```

| age_at_diagnosis.diagnoses | age_at_earliest_diagnosis.diagnoses.xena_derived |
|:---|:---|
| 4948.0 | 4948.0 |
| 1468.0 | 1468.0 |
| 5421.0 | 5421.0 |
| 5486.0 | 5486.0 |
| 2374.0 | 2374.0 |
| 551.0 | 551.0 |

The following columns are excluded because they are redundant with the
id column

``` r
target_clinical_metadata |>
  dplyr::select(id, case_id, entity_id.annotations, case_id.annotations) |>
  head()
```

| id | case_id | entity_id.annotations | case_id.annotations |
|:---|:---|:---|:---|
| 53891f08-a221-4ca8-a502-cf990bdb6020 | 53891f08-a221-4ca8-a502-cf990bdb6020 | 53891f08-a221-4ca8-a502-cf990bdb6020 | 53891f08-a221-4ca8-a502-cf990bdb6020 |
| 8706070d-5735-40f2-b77c-ce440fe3ef29 | 8706070d-5735-40f2-b77c-ce440fe3ef29 | 8706070d-5735-40f2-b77c-ce440fe3ef29 | 8706070d-5735-40f2-b77c-ce440fe3ef29 |
| 97d3cab7-80be-45fc-95bf-db16118dd95a | 97d3cab7-80be-45fc-95bf-db16118dd95a | 97d3cab7-80be-45fc-95bf-db16118dd95a | 97d3cab7-80be-45fc-95bf-db16118dd95a |
| 46cc9635-1213-4add-a704-33f25b957db2 | 46cc9635-1213-4add-a704-33f25b957db2 | 46cc9635-1213-4add-a704-33f25b957db2 | 46cc9635-1213-4add-a704-33f25b957db2 |
| 9b28de34-5805-428d-a4b3-b38fd7c66417 | 9b28de34-5805-428d-a4b3-b38fd7c66417 | 9b28de34-5805-428d-a4b3-b38fd7c66417 | 9b28de34-5805-428d-a4b3-b38fd7c66417 |
| 127b223b-5a5f-4d84-a80c-5985c4e3b146 | 127b223b-5a5f-4d84-a80c-5985c4e3b146 | 127b223b-5a5f-4d84-a80c-5985c4e3b146 | 127b223b-5a5f-4d84-a80c-5985c4e3b146 |

age_at_earliest_diagnosis_in_years.xena_derived is selected over
age_at_index.demographic because it has non-rounded age in years

``` r
target_clinical_accessions |>
  dplyr::select(age_at_index.demographic, age_at_earliest_diagnosis_in_years.diagnoses.xena_derived) |>
  head()
```

| age_at_index.demographic | age_at_earliest_diagnosis_in_years.diagnoses.xena_derived |
|:---|:---|
| 13.0 | 13.556164383561644 |
| 4.0 | 4.021917808219178 |
| 14.0 | 14.852054794520548 |
| 15.0 | 15.03013698630137 |
| 6.0 | 6.504109589041096 |
| 1.0 | 1.5095890410958903 |

There are multiple age columns, so I select age_at_diagnosis.diagnoses

``` r
target_clinical_accessions |>
  dplyr::select(age_at_earliest_diagnosis.diagnoses.xena_derived, age_at_diagnosis.diagnoses, age_at_index.demographic)
```

| age_at_earliest_diagnosis.diagnoses.xena_derived | age_at_diagnosis.diagnoses | age_at_index.demographic |
|:---|:---|:---|
| 4948.0 | 4948.0 | 13.0 |
| 1468.0 | 1468.0 | 4.0 |
| 5421.0 | 5421.0 | 14.0 |
| 5486.0 | 5486.0 | 15.0 |
| 2374.0 | 2374.0 | 6.0 |
| 551.0 | 551.0 | 1.0 |
| 1020.0 | 1020.0 | 2.0 |
| 539.0 | 539.0 | 1.0 |
| 4782.0 | 4782.0 | 13.0 |
| 795.0 | 795.0 | 2.0 |
| 3979.0 | 3979.0 | 10.0 |
| 1986.0 | 1986.0 | 5.0 |
| 2342.0 | 2342.0 | 6.0 |
| 2591.0 | 2591.0 | 7.0 |
| 7209.0 | 7209.0 | 19.0 |
| 1374.0 | 1374.0 | 3.0 |
| 2079.0 | 2079.0 | 5.0 |
| 8232.0 | 8232.0 | 22.0 |
| 5814.0 | 5814.0 | 15.0 |
| 5216.0 | 5216.0 | 14.0 |
| 5570.0 | 5570.0 | 15.0 |
| 5580.0 | 5580.0 | 15.0 |
| 4694.0 | 4694.0 | 12.0 |
| 1987.0 | 1987.0 | 5.0 |
| 5715.0 | 5715.0 | 15.0 |
| 5368.0 | 5368.0 | 14.0 |
| 3417.0 | 3417.0 | 9.0 |
| 5650.0 | 5650.0 | 15.0 |
| 4074.0 | 4074.0 | 11.0 |
| 821.0 | 821.0 | 2.0 |
| 5587.0 | 5587.0 | 15.0 |
| 4950.0 | 4950.0 | 13.0 |
| 6436.0 | 6436.0 | 17.0 |
| 4593.0 | 4593.0 | 12.0 |
| 5433.0 | 5433.0 | 14.0 |
| 6315.0 | 6315.0 | 17.0 |
| 3273.0 | 3273.0 | 8.0 |
| 4243.0 | 4243.0 | 11.0 |
| 6082.0 | 6082.0 | 16.0 |
| 3399.0 | 3399.0 | 9.0 |
| 4103.0 | 4103.0 | 11.0 |
| 4906.0 | 4906.0 | 13.0 |
| 4698.0 | 4698.0 | 12.0 |
| 1430.0 | 1430.0 | 3.0 |
| 4050.0 | 4050.0 | 11.0 |
| 4665.0 | 4665.0 | 12.0 |
| 4936.0 | 4936.0 | 13.0 |
| 4846.0 | 4846.0 | 13.0 |
| 866.0 | 866.0 | 2.0 |
| 932.0 | 932.0 | 2.0 |
| 4459.0 | 4459.0 | 12.0 |
| 705.0 | 705.0 | 1.0 |
| 6631.0 | 6631.0 | 18.0 |
| 113.0 | 113.0 | 0.0 |
| 5108.0 | 5108.0 | 13.0 |
| 6273.0 | 6273.0 | 17.0 |
| 2145.0 | 2145.0 | 5.0 |
| 2995.0 | 2995.0 | 8.0 |
| 6424.0 | 6424.0 | 17.0 |
| 2904.0 | 2904.0 | 7.0 |
| 3673.0 | 3673.0 | 10.0 |
| 567.0 | 567.0 | 1.0 |
| 5786.0 | 5786.0 | 15.0 |
| 349.0 | 349.0 | 0.0 |
| 3853.0 | 3853.0 | 10.0 |
| 3791.0 | 3791.0 | 10.0 |
| 816.0 | 816.0 | 2.0 |
| 5368.0 | 5368.0 | 14.0 |
| 6006.0 | 6006.0 | 16.0 |
| 4155.0 | 4155.0 | 11.0 |
| 6879.0 | 6879.0 | 18.0 |
| 3386.0 | 3386.0 | 9.0 |
| 4807.0 | 4807.0 | 13.0 |
| 3386.0 | 3386.0 | 9.0 |
| 277.0 | 277.0 | 0.0 |
| 2608.0 | 2608.0 | 7.0 |
| 5389.0 | 5389.0 | 14.0 |
| 10.0 | 10.0 | 0.0 |
| 4887.0 | 4887.0 | 13.0 |
| 5611.0 | 5611.0 | 15.0 |
| 1636.0 | 1636.0 | 4.0 |
| 4353.0 | 4353.0 | 11.0 |
| 5157.0 | 5157.0 | 14.0 |
| 5215.0 | 5215.0 | 14.0 |
| 5506.0 | 5506.0 | 15.0 |
| 1358.0 | 1358.0 | 3.0 |
| 5109.0 | 5109.0 | 13.0 |
| 5298.0 | 5298.0 | 14.0 |
| 2808.0 | 2808.0 | 7.0 |
| 6436.0 | 6436.0 | 17.0 |
| 2964.0 | 2964.0 | 8.0 |
| 4902.0 | 4902.0 | 13.0 |
| 6641.0 | 6641.0 | 18.0 |
| 2606.0 | 2606.0 | 7.0 |
| 3417.0 | 3417.0 | 9.0 |
| 462.0 | 462.0 | 1.0 |
| 4465.0 | 4465.0 | 12.0 |
| 1661.0 | 1661.0 | 4.0 |
| 6021.0 | 6021.0 | 16.0 |
| 4593.0 | 4593.0 | 12.0 |
| 5797.0 | 5797.0 | 15.0 |
| 4312.0 | 4312.0 | 11.0 |
| 1935.0 | 1935.0 | 5.0 |
| 1394.0 | 1394.0 | 3.0 |
| 939.0 | 939.0 | 2.0 |
| 3891.0 | 3891.0 | 10.0 |
| 5749.0 | 5749.0 | 15.0 |
| 4152.0 | 4152.0 | 11.0 |
| 5392.0 | 5392.0 | 14.0 |
| 5973.0 | 5973.0 | 16.0 |
| 2516.0 | 2516.0 | 6.0 |
| 4757.0 | 4757.0 | 13.0 |
| 986.0 | 986.0 | 2.0 |
| 6335.0 | 6335.0 | 17.0 |
| 4209.0 | 4209.0 | 11.0 |
| 2948.0 | 2948.0 | 8.0 |
| 4861.0 | 4861.0 | 13.0 |
| 557.0 | 557.0 | 1.0 |
| 6657.0 | 6657.0 | 18.0 |
| 2277.0 | 2277.0 | 6.0 |
| 1373.0 | 1373.0 | 3.0 |
| 3354.0 | 3354.0 | 9.0 |
| 4837.0 | 4837.0 | 13.0 |
| 4183.0 | 4183.0 | 11.0 |
| 2558.0 | 2558.0 | 7.0 |
| 6157.0 | 6157.0 | 16.0 |
| 5787.0 | 5787.0 | 15.0 |
| 3832.0 | 3832.0 | 10.0 |
| 6241.0 | 6241.0 | 17.0 |
| 1234.0 | 1234.0 | 3.0 |
| 7442.0 | 7442.0 | 20.0 |
| 6006.0 | 6006.0 | 16.0 |
| 2330.0 | 2330.0 | 6.0 |
| 5917.0 | 5917.0 | 16.0 |
| 6219.0 | 6219.0 | 17.0 |
| 549.0 | 549.0 | 1.0 |
| 6366.0 | 6366.0 | 17.0 |
| 2275.0 | 2275.0 | 6.0 |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| 5305.0 | 5305.0 | 14.0 |
| 590.0 | 590.0 | 1.0 |
| 4950.0 | 4950.0 | 13.0 |
| 703.0 | 703.0 | 1.0 |
| 5122.0 | 5122.0 | 14.0 |
| 358.0 | 358.0 | 0.0 |
| 1695.0 | 1695.0 | 4.0 |
| 3262.0 | 3262.0 | 8.0 |
| 4401.0 | 4401.0 | 12.0 |
| 900.0 | 900.0 | 2.0 |
| 645.0 | 645.0 | 1.0 |
| 6315.0 | 6315.0 | 17.0 |
| 3830.0 | 3830.0 | 10.0 |
| 4723.0 | 4723.0 | 12.0 |
| 2692.0 | 2692.0 | 7.0 |
| 1454.0 | 1454.0 | 3.0 |
| 1430.0 | 1430.0 | 3.0 |
| 6446.0 | 6446.0 | 17.0 |
| 6102.0 | 6102.0 | 16.0 |
| 3399.0 | 3399.0 | 9.0 |
| 491.0 | 491.0 | 1.0 |
| 4855.0 | 4855.0 | 13.0 |
| 2855.0 | 2855.0 | 7.0 |
| 359.0 | 359.0 | 0.0 |
| 2321.0 | 2321.0 | 6.0 |
| 2318.0 | 2318.0 | 6.0 |
| 2455.0 | 2455.0 | 6.0 |
| 6922.0 | 6922.0 | 18.0 |
| 733.0 | 733.0 | 2.0 |
| 367.0 | 367.0 | 1.0 |
| 3310.0 | 3310.0 | 9.0 |
| 901.0 | 901.0 | 2.0 |
| 3774.0 | 3774.0 | 10.0 |
| 2201.0 | 2201.0 | 6.0 |
| 5830.0 | 5830.0 | 15.0 |
| 5916.0 | 5916.0 | 16.0 |
| 5251.0 | 5251.0 | 14.0 |
| 227.0 | 227.0 | 0.0 |
| 283.0 | 283.0 | 0.0 |
| 283.0 | 283.0 | 0.0 |
| 3381.0 | 3381.0 | 9.0 |
| 1597.0 | 1597.0 | 4.0 |
| 6979.0 | 6979.0 | 19.0 |
| 5305.0 | 5305.0 | 14.0 |
| 2758.0 | 2758.0 | 7.0 |
| 3806.0 | 3806.0 | 10.0 |
| 2695.0 | 2695.0 | 7.0 |
| 446.0 | 446.0 | 1.0 |
| 1262.0 | 1262.0 | 3.0 |
| 339.0 | 339.0 | 0.0 |
| 6000.0 | 6000.0 | 16.0 |
| 1259.0 | 1259.0 | 3.0 |
| 5452.0 | 5452.0 | 14.0 |
| 3399.0 | 3399.0 | 9.0 |
| 292.0 | 292.0 | 0.0 |
| 3273.0 | 3273.0 | 8.0 |
| 1189.0 | 1189.0 | 3.0 |
| 664.0 | 664.0 | 1.0 |
| 748.0 | 748.0 | 2.0 |
| 6000.0 | 6000.0 | 16.0 |
| 2321.0 | 2321.0 | 6.0 |
| 6573.0 | 6573.0 | 17.0 |
| 6470.0 | 6470.0 | 17.0 |
| 2692.0 | 2692.0 | 7.0 |
| 807.0 | 807.0 | 2.0 |
| 1713.0 | 1713.0 | 4.0 |
| 1182.0 | 1182.0 | 3.0 |
| 1474.0 | 1474.0 | 4.0 |
| 2904.0 | 2904.0 | 7.0 |
| 2476.0 | 2476.0 | 6.0 |
| 1342.0 | 1342.0 | 3.0 |
| 1124.0 | 1124.0 | 3.0 |
| 1958.0 | 1958.0 | 5.0 |
| 1476.0 | 1476.0 | 4.0 |
| 2201.0 | 2201.0 | 6.0 |
| 1650.0 | 1650.0 | 4.0 |
| 824.0 | 824.0 | 2.0 |
| 2070.0 | 2070.0 | 5.0 |
| 5082.0 | 5082.0 | 13.0 |
| 405.0 | 405.0 | 1.0 |
| 1754.0 | 1754.0 | 4.0 |
| 970.0 | 970.0 | 2.0 |
| 1348.0 | 1348.0 | 3.0 |
| 160.0 | 160.0 | 0.0 |
| 1710.0 | 1710.0 | 4.0 |
| 2412.0 | 2412.0 | 6.0 |
| 1599.0 | 1599.0 | 4.0 |
| 160.0 | 160.0 | 0.0 |
| 222.0 | 222.0 | 0.0 |
| 1466.0 | 1466.0 | 4.0 |
| 800.0 | 800.0 | 2.0 |
| 1124.0 | 1124.0 | 3.0 |
| 1509.0 | 1509.0 | 4.0 |
| 581.0 | 581.0 | 1.0 |
| 2271.0 | 2271.0 | 6.0 |
| 314.0 | 314.0 | 0.0 |
| 2247.0 | 2247.0 | 6.0 |
| 707.0 | 707.0 | 1.0 |
| 1340.0 | 1340.0 | 3.0 |
| 2307.0 | 2307.0 | 6.0 |
| 2904.0 | 2904.0 | 7.0 |
| 851.0 | 851.0 | 2.0 |
| 1638.0 | 1638.0 | 4.0 |
| 2842.0 | 2842.0 | 7.0 |
| 768.0 | 768.0 | 2.0 |
| 1030.0 | 1030.0 | 2.0 |
| 979.0 | 979.0 | 2.0 |
| 3037.0 | 3037.0 | 8.0 |
| 1365.0 | 1365.0 | 3.0 |
| 670.0 | 670.0 | 1.0 |
| 1757.0 | 1757.0 | 4.0 |
| 707.0 | 707.0 | 1.0 |
| 2168.0 | 2168.0 | 5.0 |
| 325.0 | 325.0 | 0.0 |
| 4779.0 | 4779.0 | 13.0 |
| 314.0 | 314.0 | 0.0 |
| 1968.0 | 1968.0 | 5.0 |
| 1112.0 | 1112.0 | 3.0 |
| 2879.0 | 2879.0 | 7.0 |
| 1895.0 | 1895.0 | 5.0 |
| 1950.0 | 1950.0 | 5.0 |
| 2464.0 | 2464.0 | 6.0 |
| 1821.0 | 1821.0 | 4.0 |
| 2488.0 | 2488.0 | 6.0 |
| 2819.0 | 2819.0 | 7.0 |
| 281.0 | 281.0 | 0.0 |
| 961.0 | 961.0 | 2.0 |
| 1389.0 | 1389.0 | 3.0 |
| 1895.0 | 1895.0 | 5.0 |
| 521.0 | 521.0 | 1.0 |
| 2326.0 | 2326.0 | 6.0 |
| 1914.0 | 1914.0 | 5.0 |
| 1789.0 | 1789.0 | 4.0 |
| 5698.0 | 5698.0 | 15.0 |
| 797.0 | 797.0 | 2.0 |
| 547.0 | 547.0 | 1.0 |
| 490.0 | 490.0 | 1.0 |
| 1953.0 | 1953.0 | 5.0 |
| 851.0 | 851.0 | 2.0 |
| 2750.0 | 2750.0 | 7.0 |
| 1785.0 | 1785.0 | 4.0 |
| 1324.0 | 1324.0 | 3.0 |
| 2104.0 | 2104.0 | 5.0 |
| 3327.0 | 3327.0 | 9.0 |
| 1275.0 | 1275.0 | 3.0 |
| 706.0 | 706.0 | 1.0 |
| 1245.0 | 1245.0 | 3.0 |
| 1336.0 | 1336.0 | 3.0 |
| 1121.0 | 1121.0 | 3.0 |
| 1069.0 | 1069.0 | 2.0 |
| 1007.0 | 1007.0 | 2.0 |
| 1789.0 | 1789.0 | 4.0 |
| 2244.0 | 2244.0 | 6.0 |
| 1508.0 | 1508.0 | 4.0 |
| 829.0 | 829.0 | 2.0 |
| 2360.0 | 2360.0 | 6.0 |
| 691.0 | 691.0 | 1.0 |
| 1921.0 | 1921.0 | 5.0 |
| 1785.0 | 1785.0 | 4.0 |
| 1002.0 | 1002.0 | 2.0 |
| 2809.0 | 2809.0 | 7.0 |
| 4941.0 | 4941.0 | 13.0 |
| 1310.0 | 1310.0 | 3.0 |
| 787.0 | 787.0 | 2.0 |
| 699.0 | 699.0 | 1.0 |
| 5147.0 | 5147.0 | 14.0 |
| 1024.0 | 1024.0 | 2.0 |
| 1865.0 | 1865.0 | 5.0 |
| 1538.0 | 1538.0 | 4.0 |
| 972.0 | 972.0 | 2.0 |
| 3250.0 | 3250.0 | 8.0 |
| 6555.0 | 6555.0 | 17.0 |
| 506.0 | 506.0 | 1.0 |
| 2030.0 | 2030.0 | 5.0 |
| 6375.0 | 6375.0 | 17.0 |
| 3571.0 | 3571.0 | 9.0 |
| 4086.0 | 4086.0 | 11.0 |
| 6555.0 | 6555.0 | 17.0 |
| 4173.0 | 4173.0 | 11.0 |
| 1449.0 | 1449.0 | 3.0 |
| 4948.0 | 4948.0 | 13.0 |
| 3658.0 | 3658.0 | 10.0 |
| 1089.0 | 1089.0 | 2.0 |
| 10946.0 | 10946.0 | 29.0 |
| 2414.0 | 2414.0 | 6.0 |
| 5109.0 | 5109.0 | 13.0 |
| 4719.0 | 4719.0 | 12.0 |
| 2636.0 | 2636.0 | 7.0 |
| 1449.0 | 1449.0 | 3.0 |
| 2496.0 | 2496.0 | 6.0 |
| 2127.0 | 2127.0 | 5.0 |
| 1060.0 | 1060.0 | 2.0 |
| 5405.0 | 5405.0 | 14.0 |
| 4456.0 | 4456.0 | 12.0 |
| 6491.0 | 6491.0 | 17.0 |
| 3636.0 | 3636.0 | 9.0 |
| 889.0 | 889.0 | 2.0 |
| 5405.0 | 5405.0 | 14.0 |
| 3901.0 | 3901.0 | 10.0 |
| 417.0 | 417.0 | 1.0 |
| 2164.0 | 2164.0 | 5.0 |
| 4793.0 | 4793.0 | 13.0 |
| 6273.0 | 6273.0 | 17.0 |
| 5346.0 | 5346.0 | 14.0 |
| 4793.0 | 4793.0 | 13.0 |
| 5181.0 | 5181.0 | 14.0 |
| 1875.0 | 1875.0 | 5.0 |
| 1913.0 | 1913.0 | 5.0 |
| 173.0 | 173 | 0 |
| 361.0 | 361 | 0 |
| 740.0 | 740 | 2 |
| 354.0 | 354 | 0 |
| 569.0 | 569 | 1 |
| 399.0 | 399 | 1 |
| 1466.0 | 1466 | 4 |
| 352.0 | 352 | 0 |
| 400.0 | 400 | 1 |
| 220.0 | 220 | 0 |
| 286.0 | 286 | 0 |
| 54.0 | 54 | 0 |
| 290.0 | 290 | 0 |
| 332.0 | 332 | 0 |
| 146.0 | 146 | 0 |
| 272.0 | 272 | 0 |
| 5516.0 | 5516.0 | 15.0 |
| 4173.0 | 4173.0 | 11.0 |
| 2916.0 | 2916.0 | 7.0 |
| 2005.0 | 2005.0 | 5.0 |
| 1986.0 | 1986.0 | 5.0 |
| 1786.0 | 1786.0 | 4.0 |
| 4948.0 | 4948.0 | 13.0 |
| 1449.0 | 1449.0 | 3.0 |
| 2329.0 | 2329.0 | 6.0 |
| 1224.0 | 1224.0 | 3.0 |
| 2342.0 | 2342.0 | 6.0 |
| 2591.0 | 2591.0 | 7.0 |
| 3051.0 | 3051.0 | 8.0 |
| 1043.0 | 1043.0 | 2.0 |
| 2788.0 | 2788.0 | 7.0 |
| 1189.0 | 1189.0 | 3.0 |
| 2356.0 | 2356.0 | 6.0 |
| 829.0 | 829.0 | 2.0 |
| 5231.0 | 5231.0 | 14.0 |
| 1054.0 | 1054.0 | 2.0 |
| 2041.0 | 2041.0 | 5.0 |
| 1098.0 | 1098.0 | 3.0 |
| 699.0 | 699.0 | 1.0 |
| 2511.0 | 2511.0 | 6.0 |
| 812.0 | 812.0 | 2.0 |
| 3519.0 | 3519.0 | 9.0 |
| 1129.0 | 1129.0 | 3.0 |
| 1335.0 | 1335.0 | 3.0 |
| 4029.0 | 4029.0 | 11.0 |
| 2543.0 | 2543.0 | 6.0 |
| 5021.0 | 5021.0 | 13.0 |
| 1124.0 | 1124.0 | 3.0 |
| 785.0 | 785.0 | 2.0 |
| 5421.0 | 5421.0 | 14.0 |
| 1223.0 | 1223.0 | 3.0 |
| 1020.0 | 1020.0 | 2.0 |
| 5273.0 | 5273.0 | 14.0 |
| 6270.0 | 6270.0 | 17.0 |
| 6713.0 | 6713.0 | 18.0 |
| 4756.0 | 4756.0 | 13.0 |
| 837.0 | 837.0 | 2.0 |
| 2775.0 | 2775.0 | 7.0 |
| 5906.0 | 5906.0 | 16.0 |
| 2788.0 | 2788.0 | 7.0 |
| 1786.0 | 1786.0 | 4.0 |
| 2526.0 | 2526.0 | 6.0 |
| 5491.0 | 5491.0 | 15.0 |
| 886.0 | 886.0 | 2.0 |
| 1711.0 | 1711.0 | 4.0 |
| 2591.0 | 2591.0 | 7.0 |
| 651.0 | 651.0 | 1.0 |
| 4845.0 | 4845.0 | 13.0 |
| 1479.0 | 1479.0 | 4.0 |
| 3658.0 | 3658.0 | 10.0 |
| 6103.0 | 6103.0 | 16.0 |
| 396.0 | 396.0 | 1.0 |
| 2823.0 | 2823.0 | 7.0 |
| 1833.0 | 1833.0 | 5.0 |
| 978.0 | 978.0 | 2.0 |
| 6297.0 | 6297.0 | 17.0 |
| 3519.0 | 3519.0 | 9.0 |
| 1468.0 | 1468.0 | 4.0 |
| 1358.0 | 1358.0 | 3.0 |
| 1217.0 | 1217.0 | 3.0 |
| 4018.0 | 4018.0 | 11.0 |
| 5203.0 | 5203.0 | 14.0 |
| 539.0 | 539.0 | 1.0 |
| 4655.0 | 4655.0 | 12.0 |
| 7519.0 | 7519.0 | 20.0 |
| 4798.0 | 4798.0 | 13.0 |
| 1129.0 | 1129.0 | 3.0 |
| 679.0 | 679.0 | 1.0 |
| 2146.0 | 2146.0 | 5.0 |
| 1355.0 | 1355.0 | 3.0 |
| 5098.0 | 5098.0 | 13.0 |
| 6685.0 | 6685.0 | 18.0 |
| 4756.0 | 4756.0 | 13.0 |
| 2219.0 | 2219.0 | 6.0 |
| 853.0 | 853.0 | 2.0 |
| 1033.0 | 1033.0 | 2.0 |
| 4810.0 | 4810.0 | 13.0 |
| 2414.0 | 2414.0 | 6.0 |
| 648.0 | 648.0 | 1.0 |
| 2005.0 | 2005.0 | 5.0 |
| 2329.0 | 2329.0 | 6.0 |
| 903.0 | 903.0 | 2.0 |
| 3051.0 | 3051.0 | 8.0 |
| 789.0 | 789.0 | 2.0 |
| 883.0 | 883.0 | 2.0 |
| 4892.0 | 4892.0 | 13.0 |
| 5147.0 | 5147.0 | 14.0 |
| 785.0 | 785.0 | 2.0 |
| 1986.0 | 1986.0 | 5.0 |
| 4093.0 | 4093.0 | 11.0 |
| 3976.0 | 3976.0 | 10.0 |
| 1133.0 | 1133.0 | 3.0 |
| 2606.0 | 2606.0 | 7.0 |
| 4616.0 | 4616.0 | 12.0 |
| 3625.0 | 3625.0 | 9.0 |
| 759.0 | 759.0 | 2.0 |
| 845.0 | 845.0 | 2.0 |
| 1043.0 | 1043.0 | 2.0 |
| 3000.0 | 3000.0 | 8.0 |
| 2677.0 | 2677.0 | 7.0 |
| 5516.0 | 5516.0 | 15.0 |
| 3215.0 | 3215.0 | 8.0 |
| 5652.0 | 5652.0 | 15.0 |
| 1302.0 | 1302.0 | 3.0 |
| 940.0 | 940.0 | 2.0 |
| 2543.0 | 2543.0 | 6.0 |
| 673.0 | 673.0 | 1.0 |
| 2342.0 | 2342.0 | 6.0 |
| 1464.0 | 1464.0 | 4.0 |
| 3250.0 | 3250.0 | 8.0 |
| 1077.0 | 1077.0 | 2.0 |
| 679.0 | 679.0 | 1.0 |
| 2374.0 | 2374.0 | 6.0 |
| 4616.0 | 4616.0 | 12.0 |
| 2636.0 | 2636.0 | 7.0 |
| 852.0 | 852.0 | 2.0 |
| 2013.0 | 2013.0 | 5.0 |
| 1024.0 | 1024.0 | 2.0 |
| 6842.0 | 6842.0 | 18.0 |
| 431.0 | 431.0 | 1.0 |
| 585.0 | 585.0 | 1.0 |
| 5213.0 | 5213.0 | 14.0 |
| 1131.0 | 1131.0 | 3.0 |
| 6599.0 | 6599.0 | 18.0 |
| 1129.0 | 1129.0 | 3.0 |
| 3679.0 | 3679.0 | 10.0 |
| 4087.0 | 4087.0 | 11.0 |
| 837.0 | 837.0 | 2.0 |
| 3187.0 | 3187.0 | 8.0 |
| 853.0 | 853.0 | 2.0 |
| 2511.0 | 2511.0 | 6.0 |
| 815.0 | 815.0 | 2.0 |
| 6535.0 | 6535.0 | 17.0 |
| 2191.0 | 2191.0 | 5.0 |
| 549.0 | 549.0 | 1.0 |
| 163.0 | 163.0 | 0.0 |
| 5452.0 | 5452.0 | 14.0 |
| 4723.0 | 4723.0 | 12.0 |
| 3298.0 | 3298.0 | 9.0 |
| 1248.0 | 1248.0 | 3.0 |
| 1330.0 | 1330.0 | 3.0 |
| 3102.0 | 3102.0 | 8.0 |
| 1159.0 | 1159.0 | 3.0 |
| 4151.0 | 4151.0 | 11.0 |
| 5830.0 | 5830.0 | 15.0 |
| 5082.0 | 5082.0 | 13.0 |
| 658.0 | 658.0 | 1.0 |
| 3972.0 | 3972.0 | 10.0 |
| 3687.0 | 3687.0 | 10.0 |
| 5325.0 | 5325.0 | 14.0 |
| 1248.0 | 1248.0 | 3.0 |
| 282.0 | 282.0 | 0.0 |
| 677.0 | 677.0 | 1.0 |
| 3262.0 | 3262.0 | 8.0 |
| 4559.0 | 4559.0 | 12.0 |
| 1597.0 | 1597.0 | 4.0 |
| 5885.0 | 5885.0 | 16.0 |
| 5420.0 | 5420.0 | 14.0 |
| 900.0 | 900.0 | 2.0 |
| 5587.0 | 5587.0 | 15.0 |
| 664.0 | 664.0 | 1.0 |
| 1022.0 | 1022.0 | 2.0 |
| 768.0 | 768.0 | 2.0 |
| 890.0 | 890.0 | 2.0 |
| 6303.0 | 6303.0 | 17.0 |
| 6660.0 | 6660.0 | 18.0 |
| 428.0 | 428.0 | 1.0 |
| 3889.0 | 3889.0 | 10.0 |
| 2455.0 | 2455.0 | 6.0 |
| 5600.0 | 5600.0 | 15.0 |
| 2855.0 | 2855.0 | 7.0 |
| 5251.0 | 5251.0 | 14.0 |
| 5814.0 | 5814.0 | 15.0 |
| 6922.0 | 6922.0 | 18.0 |
| 1159.0 | 1159.0 | 3.0 |
| 137.0 | 137.0 | 0.0 |
| 5356.0 | 5356.0 | 14.0 |
| 3102.0 | 3102.0 | 8.0 |
| 6660.0 | 6660.0 | 18.0 |
| 282.0 | 282.0 | 0.0 |
| 3748.0 | 3748.0 | 10.0 |
| 5082.0 | 5082.0 | 13.0 |
| 757.0 | 757.0 | 2.0 |
| 1482.0 | 1482.0 | 4.0 |
| 677.0 | 677.0 | 1.0 |
| 963.0 | 963.0 | 2.0 |
| 4458.0 | 4458.0 | 12.0 |
| 161.0 | 161.0 | 0.0 |
| 4103.0 | 4103.0 | 11.0 |
| 447.0 | 447.0 | 1.0 |
| 3806.0 | 3806.0 | 10.0 |
| 2748.0 | 2748.0 | 7.0 |
| 3890.0 | 3890.0 | 10.0 |
| 4545.0 | 4545.0 | 12.0 |
| 6590.0 | 6590.0 | 18.0 |
| 5325.0 | 5325.0 | 14.0 |
| 534.0 | 534.0 | 1.0 |
| 887.0 | 887.0 | 2.0 |
| 1685.0 | 1685.0 | 4.0 |
| 2265.0 | 2265.0 | 6.0 |
| 2702.0 | 2702.0 | 7.0 |
| 1868.0 | 1868.0 | 5.0 |
| 1236.0 | 1236.0 | 3.0 |
| 1476.0 | 1476.0 | 4.0 |
| 4397.0 | 4397.0 | 12.0 |
| 1548.0 | 1548.0 | 4.0 |
| 401.0 | 401.0 | 1.0 |
| 1246.0 | 1246.0 | 3.0 |
| 1912.0 | 1912.0 | 5.0 |
| 1287.0 | 1287.0 | 3.0 |
| 1532.0 | 1532.0 | 4.0 |
| 2238.0 | 2238.0 | 6.0 |
| 3703.0 | 3703.0 | 10.0 |
| 1326.0 | 1326.0 | 3.0 |
| 5557.0 | 5557.0 | 15.0 |
| 2122.0 | 2122.0 | 5.0 |
| 1262.0 | 1262.0 | 3.0 |
| 2201.0 | 2201.0 | 6.0 |
| 156.0 | 156.0 | 0.0 |
| 1340.0 | 1340.0 | 3.0 |
| 2423.0 | 2423.0 | 6.0 |
| 1449.0 | 1449.0 | 3.0 |
| 2681.0 | 2681.0 | 7.0 |
| 1681.0 | 1681.0 | 4.0 |
| 1309.0 | 1309.0 | 3.0 |
| 2003.0 | 2003.0 | 5.0 |
| 551.0 | 551.0 | 1.0 |
| 1020.0 | 1020.0 | 2.0 |
| 1464.0 | 1464.0 | 4.0 |
| 4151.0 | 4151.0 | 11.0 |
| 4562.0 | 4562.0 | 12.0 |
| 1329.0 | 1329.0 | 3.0 |
| 5486.0 | 5486.0 | 15.0 |
| 2677.0 | 2677.0 | 7.0 |
| 2005.0 | 2005.0 | 5.0 |
| 1865.0 | 1865.0 | 5.0 |
| 596.0 | 596.0 | 1.0 |
| 795.0 | 795.0 | 2.0 |
| 5109.0 | 5109.0 | 13.0 |
| 5479.0 | 5479.0 | 15.0 |
| 648.0 | 648.0 | 1.0 |
| 5906.0 | 5906.0 | 16.0 |
| 4151.0 | 4151.0 | 11.0 |
| 4601.0 | 4601.0 | 12.0 |
| 2127.0 | 2127.0 | 5.0 |
| 6419.0 | 6419.0 | 17.0 |
| 5103.0 | 5103.0 | 13.0 |
| 1538.0 | 1538.0 | 4.0 |
| 5479.0 | 5479.0 | 15.0 |
| 4016.0 | 4016.0 | 10.0 |
| 5486.0 | 5486.0 | 15.0 |
| 3979.0 | 3979.0 | 10.0 |
| 4782.0 | 4782.0 | 13.0 |
| 3076.0 | 3076.0 | 8.0 |
| 3440.0 | 3440.0 | 9.0 |
| 6089.0 | 6089.0 | 16.0 |
| 3836.0 | 3836.0 | 10.0 |
| 31.0 | 31.0 | 0.0 |
| 430.0 | 430.0 | 1.0 |
| 6008.0 | 6008.0 | 16.0 |
| 3901.0 | 3901.0 | 10.0 |
| 417.0 | 417.0 | 1.0 |
| 6089.0 | 6089.0 | 16.0 |
| 1875.0 | 1875.0 | 5.0 |
| 450.0 | 450.0 | 1.0 |
| 6491.0 | 6491.0 | 17.0 |
| 3076.0 | 3076.0 | 8.0 |
| 430.0 | 430.0 | 1.0 |
| 6273.0 | 6273.0 | 17.0 |
| 5405.0 | 5405.0 | 14.0 |
| 6273.0 | 6273.0 | 17.0 |
| 3711.0 | 3711.0 | 10.0 |
| 3636.0 | 3636.0 | 9.0 |
| 450.0 | 450.0 | 1.0 |
| 2145.0 | 2145.0 | 5.0 |
| 8581.0 | 8581.0 | 23.0 |
| 430.0 | 430.0 | 1.0 |
| 6008.0 | 6008.0 | 16.0 |
| 3440.0 | 3440.0 | 9.0 |
| 3836.0 | 3836.0 | 10.0 |
| 4859.0 | 4859.0 | 13.0 |
| 2164.0 | 2164.0 | 5.0 |
| 8581.0 | 8581.0 | 23.0 |
| 6273.0 | 6273.0 | 17.0 |
| 318.0 | 318.0 | 0.0 |
| 5338.0 | 5338.0 | 14.0 |
| 3711.0 | 3711.0 | 10.0 |
| 5346.0 | 5346.0 | 14.0 |
| 1913.0 | 1913.0 | 5.0 |
| 318.0 | 318.0 | 0.0 |
| 4793.0 | 4793.0 | 13.0 |
| 1799.0 | 1799.0 | 4.0 |
| 3711.0 | 3711.0 | 10.0 |
| 5181.0 | 5181.0 | 14.0 |
| 318.0 | 318.0 | 0.0 |
| 3440.0 | 3440.0 | 9.0 |
| 6089.0 | 6089.0 | 16.0 |
| 889.0 | 889.0 | 2.0 |
| 3440.0 | 3440.0 | 9.0 |
| 4456.0 | 4456.0 | 12.0 |
| 8581.0 | 8581.0 | 23.0 |
| 1875.0 | 1875.0 | 5.0 |
| 4456.0 | 4456.0 | 12.0 |
| 889.0 | 889.0 | 2.0 |
| 2145.0 | 2145.0 | 5.0 |
| 1799.0 | 1799.0 | 4.0 |
| 430.0 | 430.0 | 1.0 |
| 31.0 | 31.0 | 0.0 |
| 4456.0 | 4456.0 | 12.0 |
| 2164.0 | 2164.0 | 5.0 |
| 318.0 | 318.0 | 0.0 |
| 6491.0 | 6491.0 | 17.0 |
| 450.0 | 450.0 | 1.0 |
| 6008.0 | 6008.0 | 16.0 |
| 6089.0 | 6089.0 | 16.0 |
| 5346.0 | 5346.0 | 14.0 |
| 3971.0 | 3971.0 | 10.0 |
| 3901.0 | 3901.0 | 10.0 |
| 343.0 | 343 | 0 |
| 205.0 | 205 | 0 |
| 526.0 | 526 | 1 |
| 5386.0 | 5386 | 14 |
| 361.0 | 361 | 0 |
| 580.0 | 580 | 1 |
| 214.0 | 214 | 0 |
| 464.0 | 464 | 1 |
| 141.0 | 141 | 0 |
| 306.0 | 306 | 0 |
| 332.0 | 332 | 0 |
| 312.0 | 312 | 0 |
| 147.0 | 147 | 0 |
| 54.0 | 54 | 0 |
| 64.0 | 64 | 0 |
| 286.0 | 286 | 0 |
| 102.0 | 102 | 0 |
| 138.0 | 138 | 0 |
| 592.0 | 592 | 1 |
| 948.0 | 948 | 2 |
| 214.0 | 214 | 0 |
| 63.0 | 63 | 0 |
| 239.0 | 239 | 0 |
| 1038.0 | 1038 | 2 |
| 761.0 | 761 | 2 |
| 281.0 | 281 | 0 |
| 1106.0 | 1106.0 | 3.0 |
| 503.0 | 503.0 | 1.0 |
| 4164.0 | 4164.0 | 11.0 |
| 569.0 | 569.0 | 1.0 |
| 5236.0 | 5236.0 | 14.0 |
| 3491.0 | 3491.0 | 9.0 |
| 1826.0 | 1826.0 | 4.0 |
| 948.0 | 948.0 | 2.0 |
| 949.0 | 949.0 | 2.0 |
| 425.0 | 425.0 | 1.0 |
| 383.0 | 383.0 | 1.0 |
| 4958.0 | 4958.0 | 13.0 |
| 682.0 | 682.0 | 1.0 |
| 2057.0 | 2057.0 | 5.0 |
| 3627.0 | 3627.0 | 9.0 |
| 4024.0 | 4024.0 | 11.0 |
| 1398.0 | 1398.0 | 3.0 |
| 3621.0 | 3621.0 | 9.0 |
| 1884.0 | 1884.0 | 5.0 |
| 456.0 | 456.0 | 1.0 |
| 2758.0 | 2758.0 | 7.0 |
| 1370.0 | 1370.0 | 3.0 |
| 383.0 | 383.0 | 1.0 |
| 1641.0 | 1641.0 | 4.0 |
| 1677.0 | 1677.0 | 4.0 |
| 762.0 | 762.0 | 2.0 |
| 1646.0 | 1646.0 | 4.0 |
| 5247.0 | 5247.0 | 14.0 |
| 1577.0 | 1577.0 | 4.0 |
| 1182.0 | 1182.0 | 3.0 |
| 851.0 | 851.0 | 2.0 |
| 2199.0 | 2199.0 | 6.0 |
| 1187.0 | 1187.0 | 3.0 |
| 1652.0 | 1652.0 | 4.0 |
| 928.0 | 928.0 | 2.0 |
| 956.0 | 956.0 | 2.0 |
| 1176.0 | 1176.0 | 3.0 |
| 365.0 | 365.0 | 0.0 |
| 1873.0 | 1873.0 | 5.0 |
| 719.0 | 719.0 | 1.0 |
| 912.0 | 912.0 | 2.0 |
| 5684.0 | 5684.0 | 15.0 |
| 3777.0 | 3777.0 | 10.0 |
| 1887.0 | 1887.0 | 5.0 |
| 3813.0 | 3813.0 | 10.0 |
| 6604.0 | 6604.0 | 18.0 |
| 1786.0 | 1786.0 | 4.0 |
| 1542.0 | 1542.0 | 4.0 |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| 804.0 | 804 | 2.0 |
| 1544.0 | 1544 | 4.0 |
| 1466.0 | 1466 | 4.0 |
| 1279.0 | 1279 | 3.0 |
| 823.0 | 823 | 2.0 |
| 1214.0 | 1214 | 3.0 |
| 1810.0 | 1810 | 4.0 |
| 954.0 | 954 | 2.0 |
| 699.0 | 699 | 1.0 |
| 746.0 | 746 | 2.0 |
| 1661.0 | 1661 | 4.0 |
| 1554.0 | 1554 | 4.0 |
| 1030.0 | 1030 | 2.0 |
| 614.0 | 614 | 1.0 |
| 1064.0 | 1064 | 2.0 |
| 1486.0 | 1486 | 4.0 |
| 993.0 | 993 | 2.0 |
| 1074.0 | 1074 | 2.0 |
| 710.0 | 710 | 1.0 |
| 1123.0 | 1123 | 3.0 |
| 1034.0 | 1034 | 2.0 |
| 659.0 | 659 | 1.0 |
| 728.0 | 728 | 1.0 |
| 990.0 | 990 | 2.0 |
| 1315.0 | 1315 | 3.0 |
| 679.0 | 679 | 1.0 |
| 1367.0 | 1367 | 3.0 |
| 1059.0 | 1059 | 2.0 |
| 602.0 | 602 | 1.0 |
| 1266.0 | 1266 | 3.0 |
| 1010.0 | 1010 | 2.0 |
| 1280.0 | 1280 | 3.0 |
| 1757.0 | 1757 | 4.0 |
| 1713.0 | 1713 | 4.0 |
| 1710.0 | 1710 | 4.0 |
| 993.0 | 993 | 2.0 |
| 806.0 | 806 | 2.0 |
| 898.0 | 898 | 2.0 |
| 666.0 | 666 | 1.0 |
| 1117.0 | 1117 | 3.0 |
| 515.0 | 515 | 1.0 |
| 1595.0 | 1595 | 4.0 |
| 3425.0 | 3425 | 9.0 |
| 1010.0 | 1010 | 2.0 |
| 1442.0 | 1442 | 3.0 |
| 1730.0 | 1730 | 4.0 |
| 1973.0 | 1973 | 5.0 |
| 4781.0 | 4781 | 13.0 |
| 2329.0 | 2329 | 6.0 |
| 2234.0 | 2234 | 6.0 |
| 2937.0 | 2937 | 8.0 |
| 1860.0 | 1860 | 5.0 |
| 2008.0 | 2008 | 5.0 |
| 2547.0 | 2547 | 6.0 |
| 1462.0 | 1462 | 4.0 |
| 754.0 | 754 | 2.0 |
| 1955.0 | 1955 | 5.0 |
| 1961.0 | 1961 | 5.0 |
| 1366.0 | 1366 | 3.0 |
| 550.0 | 550 | 1.0 |
| 1180.0 | 1180 | 3.0 |
| 885.0 | 885 | 2.0 |
| 3625.0 | 3625 | 9.0 |
| 2048.0 | 2048 | 5.0 |
| 1555.0 | 1555 | 4.0 |
| 2060.0 | 2060 | 5.0 |
| 2241.0 | 2241 | 6.0 |
| 1616.0 | 1616 | 4.0 |
| 2390.0 | 2390 | 6.0 |
| 2004.0 | 2004 | 5.0 |
| 6021.0 | 6021 | 16.0 |
| 2279.0 | 2279 | 6.0 |
| 2273.0 | 2273 | 6.0 |
| 1091.0 | 1091 | 2.0 |
| 1037.0 | 1037 | 2.0 |
| 1259.0 | 1259 | 3.0 |
| 1765.0 | 1765 | 4.0 |
| 4802.0 | 4802 | 13.0 |
| 1099.0 | 1099 | 3.0 |
| 753.0 | 753 | 2.0 |
| 1328.0 | 1328 | 3.0 |
| 2492.0 | 2492 | 6.0 |
| 1900.0 | 1900 | 5.0 |
| 3035.0 | 3035 | 8.0 |
| 1008.0 | 1008 | 2.0 |
| 1418.0 | 1418 | 3.0 |
| 631.0 | 631 | 1.0 |
| 727.0 | 727 | 1.0 |
| 2477.0 | 2477 | 6.0 |
| 834.0 | 834 | 2.0 |
| 1438.0 | 1438 | 3.0 |
| 583.0 | 583 | 1.0 |
| 645.0 | 645 | 1.0 |
| 1112.0 | 1112 | 3.0 |
| 1435.0 | 1435 | 3.0 |
| 688.0 | 688 | 1.0 |
| 1275.0 | 1275 | 3.0 |
| 1434.0 | 1434 | 3.0 |
| 752.0 | 752 | 2.0 |
| 1040.0 | 1040 | 2.0 |
| 648.0 | 648 | 1.0 |
| 1037.0 | 1037 | 2.0 |
| 2109.0 | 2109 | 5.0 |
| 1037.0 | 1037 | 2.0 |
| 5734.0 | 5734 | 15.0 |
| 710.0 | 710 | 1.0 |
| 564.0 | 564 | 1.0 |
| 641.0 | 641 | 1.0 |
| 2300.0 | 2300 | 6.0 |
| 21.0 | 21 | 0.0 |
| 72.0 | 72 | 0.0 |
| 58.0 | 58 | 0.0 |
| 11.0 | 11 | 0.0 |
| 60.0 | 60 | 0.0 |
| 56.0 | 56 | 0.0 |
| 41.0 | 41 | 0.0 |
| 3.0 | 3 | 0.0 |
| 18.0 | 18 | 0.0 |
| 60.0 | 60 | 0.0 |
| 51.0 | 51 | 0.0 |
| 74.0 | 74 | 0.0 |
| 25.0 | 25 | 0.0 |
| 8.0 | 8 | 0.0 |
| 62.0 | 62 | 0.0 |
| 7.0 | 7 | 0.0 |
| 4.0 | 4 | 0.0 |
| 31.0 | 31 | 0.0 |
| 52.0 | 52 | 0.0 |
| 59.0 | 59 | 0.0 |
| 7.0 | 7 | 0.0 |
| 1330.0 | 1330 | 3.0 |
| 1330.0 | 1330 | 3.0 |
| 406.0 | 406 | 1.0 |
| 35.0 | 35 | 0.0 |
| 35.0 | 35 | 0.0 |
| 341.0 | 341 | 0.0 |
| 341.0 | 341 | 0.0 |
| 242.0 | 242 | 0.0 |
| 318.0 | 318 | 0.0 |
| 318.0 | 318 | 0.0 |
| 1514.0 | 1514 | 4.0 |
| 3446.0 | 3446 | 9.0 |
| 1720.0 | 1720 | 4.0 |
| 958.0 | 958 | 2.0 |
| 837.0 | 837 | 2.0 |
| 328.0 | 328 | 0.0 |
| 328.0 | 328 | 0.0 |
| 361.0 | 361 | 0.0 |
| 361.0 | 361 | 0.0 |
| 1155.0 | 1155 | 3.0 |
| 635.0 | 635 | 1.0 |
| 1531.0 | 1531 | 4.0 |
| 1080.0 | 1080 | 2.0 |
| 696.0 | 696 | 1.0 |
| 1070.0 | 1070 | 2.0 |
| 696.0 | 696 | 1.0 |
| 5345.0 | 5345 | 14.0 |
| 2537.0 | 2537 | 6.0 |
| 1278.0 | 1278 | 3.0 |
| 1041.0 | 1041 | 2.0 |
| 1659.0 | 1659 | 4.0 |
| 408.0 | 408.0 | 1.0 |
| 4607.0 | 4607.0 | 12.0 |
| 526.0 | 526 | 1 |
| 6524.0 | 6524.0 | 17.0 |
| 5346.0 | 5346.0 | 14.0 |
| 530.0 | 530.0 | 1.0 |
| 821.0 | 821.0 | 2.0 |
| 4243.0 | 4243.0 | 11.0 |
| 3089.0 | 3089.0 | 8.0 |
| 5140.0 | 5140.0 | 14.0 |
| 5727.0 | 5727.0 | 15.0 |
| 1262.0 | 1262.0 | 3.0 |
| 6118.0 | 6118.0 | 16.0 |
| 4312.0 | 4312.0 | 11.0 |
| 5853.0 | 5853.0 | 16.0 |
| 3438.0 | 3438.0 | 9.0 |
| 6416.0 | 6416.0 | 17.0 |
| 4425.0 | 4425.0 | 12.0 |
| 618.0 | 618.0 | 1.0 |
| 590.0 | 590.0 | 1.0 |
| 5934.0 | 5934.0 | 16.0 |
| 5926.0 | 5926.0 | 16.0 |
| 1430.0 | 1430.0 | 3.0 |
| 6397.0 | 6397.0 | 17.0 |
| 261.0 | 261.0 | 0.0 |
| 3921.0 | 3921.0 | 10.0 |
| 1852.0 | 1852.0 | 5.0 |
| 4074.0 | 4074.0 | 11.0 |
| 5544.0 | 5544.0 | 15.0 |
| 5987.0 | 5987.0 | 16.0 |
| 3543.0 | 3543.0 | 9.0 |
| 3454.0 | 3454.0 | 9.0 |
| 1802.0 | 1802.0 | 4.0 |
| 2792.0 | 2792.0 | 7.0 |
| 1678.0 | 1678.0 | 4.0 |
| 4231.0 | 4231.0 | 11.0 |
| 3298.0 | 3298.0 | 9.0 |
| 766.0 | 766.0 | 2.0 |
| 5157.0 | 5157.0 | 14.0 |
| 4932.0 | 4932.0 | 13.0 |
| 4701.0 | 4701.0 | 12.0 |
| 3163.0 | 3163.0 | 8.0 |
| 147.0 | 147.0 | 0.0 |
| 690.0 | 690.0 | 1.0 |
| 3235.0 | 3235.0 | 8.0 |
| 637.0 | 637.0 | 1.0 |
| 5348.0 | 5348.0 | 14.0 |
| 5177.0 | 5177.0 | 14.0 |
| 969.0 | 969.0 | 2.0 |
| 2059.0 | 2059.0 | 5.0 |
| 4852.0 | 4852.0 | 13.0 |
| 4357.0 | 4357.0 | 11.0 |
| 2154.0 | 2154.0 | 5.0 |
| 4401.0 | 4401.0 | 12.0 |
| 3972.0 | 3972.0 | 10.0 |
| 4983.0 | 4983.0 | 13.0 |
| 1049.0 | 1049.0 | 2.0 |
| 1695.0 | 1695.0 | 4.0 |
| 227.0 | 227.0 | 0.0 |
| 2686.0 | 2686.0 | 7.0 |
| 4545.0 | 4545.0 | 12.0 |
| 5129.0 | 5129.0 | 14.0 |
| 5433.0 | 5433.0 | 14.0 |
| 6631.0 | 6631.0 | 18.0 |
| 594.0 | 594.0 | 1.0 |
| 6460.0 | 6460.0 | 17.0 |
| 4421.0 | 4421.0 | 12.0 |
| 5405.0 | 5405.0 | 14.0 |
| 137.0 | 137.0 | 0.0 |
| 3539.0 | 3539.0 | 9.0 |
| 5822.0 | 5822.0 | 15.0 |
| 3200.0 | 3200.0 | 8.0 |
| 8231.0 | 8231.0 | 22.0 |
| 6090.0 | 6090.0 | 16.0 |
| 5218.0 | 5218.0 | 14.0 |
| 1192.0 | 1192.0 | 3.0 |
| 1102.0 | 1102.0 | 3.0 |
| 1287.0 | 1287.0 | 3.0 |
| 1287.0 | 1287.0 | 3.0 |
| 1789.0 | 1789.0 | 4.0 |
| 3703.0 | 3703.0 | 10.0 |
| 1546.0 | 1546.0 | 4.0 |
| 847.0 | 847.0 | 2.0 |
| 1762.0 | 1762.0 | 4.0 |
| 2728.0 | 2728.0 | 7.0 |
| 276.0 | 276.0 | 0.0 |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| NA | NA | NA |
| 832.0 | 832 | 2 |
| 950.0 | 950 | 2 |
| 1085.0 | 1085 | 2 |
| 1141.0 | 1141 | 3 |
| 1085.0 | 1085 | 2 |
| 541.0 | 541 | 1 |
| 541.0 | 541 | 1 |
| 1141.0 | 1141 | 3 |
| 4698.0 | 4698.0 | 12.0 |
| 1661.0 | 1661.0 | 4.0 |
| 4861.0 | 4861.0 | 13.0 |
| 5108.0 | 5108.0 | 13.0 |
| 2558.0 | 2558.0 | 7.0 |
| 6366.0 | 6366.0 | 17.0 |
| 5797.0 | 5797.0 | 15.0 |
| 4837.0 | 4837.0 | 13.0 |
| 302.0 | 302.0 | 0.0 |
| 4902.0 | 4902.0 | 13.0 |
| 7599.0 | 7599.0 | 20.0 |
| 6219.0 | 6219.0 | 17.0 |
| 4312.0 | 4312.0 | 11.0 |
| 3235.0 | 3235.0 | 8.0 |
| 4425.0 | 4425.0 | 12.0 |
| 5408.0 | 5408.0 | 14.0 |
| 6460.0 | 6460.0 | 17.0 |
| 3102.0 | 3102.0 | 8.0 |
| 3748.0 | 3748.0 | 10.0 |
| 1852.0 | 1852.0 | 5.0 |
| 2948.0 | 2948.0 | 8.0 |
| 4231.0 | 4231.0 | 11.0 |
| 5822.0 | 5822.0 | 15.0 |
| 4421.0 | 4421.0 | 12.0 |
| 5472.0 | 5472.0 | 14.0 |
| 4701.0 | 4701.0 | 12.0 |
| 6273.0 | 6273.0 | 17.0 |
| 2758.0 | 2758.0 | 7.0 |
| 338.0 | 338.0 | 0.0 |
| 4780.0 | 4780.0 | 13.0 |
| 963.0 | 963.0 | 2.0 |
| 3163.0 | 3163.0 | 8.0 |
| 1967.0 | 1967.0 | 5.0 |
| 3454.0 | 3454.0 | 9.0 |
| 3381.0 | 3381.0 | 9.0 |
| 3687.0 | 3687.0 | 10.0 |
| 6986.0 | 6986.0 | 19.0 |
| 5916.0 | 5916.0 | 16.0 |
| 3670.0 | 3670.0 | 10.0 |
| 4151.0 | 4151.0 | 11.0 |
| 5727.0 | 5727.0 | 15.0 |
| 5386.0 | 5386.0 | 14.0 |
| 3682.0 | 3682.0 | 10.0 |
| 1189.0 | 1189.0 | 3.0 |
| 3200.0 | 3200.0 | 8.0 |
| 5122.0 | 5122.0 | 14.0 |
| 5389.0 | 5389.0 | 14.0 |
| 3888.0 | 3888.0 | 10.0 |
| 4545.0 | 4545.0 | 12.0 |
| 4401.0 | 4401.0 | 12.0 |
| 830.0 | 830 | 2 |
| 611.0 | 611 | 1 |
| 2355.0 | 2355 | 6 |
| 3027.0 | 3027 | 8 |
| 141.0 | 141 | 0 |
| 3.0 | 3 | 0 |
| 909.0 | 909 | 2 |
| 618.0 | 618 | 1 |
| 832.0 | 832 | 2 |
| 130.0 | 130 | 0 |
| 1035.0 | 1035 | 2 |
| 950.0 | 950 | 2 |
| 611.0 | 611 | 1 |
| 808.0 | 808 | 2 |
| 231.0 | 231 | 0 |
| 175.0 | 175 | 0 |
| 3.0 | 3 | 0 |
| 383.0 | 383 | 1 |
| 307.0 | 307 | 0 |
| 346.0 | 346 | 0 |
| 231.0 | 231 | 0 |
| 679.0 | 679 | 1 |
| 1555.0 | 1555 | 4 |
| 527.0 | 527 | 1 |
| 650.0 | 650 | 1 |
| 265.0 | 265 | 0 |
| 265.0 | 265 | 0 |
| 1035.0 | 1035 | 2 |
| 863.0 | 863 | 2 |
| 749.0 | 749 | 2 |

These fields are very similar, but contain slightly different
information, so I decided to keep them

``` r
target_clinical_accessions |>
  dplyr::select(disease_type, name.project, primary_diagnosis.diagnoses) |>
  head()
```

| disease_type | name.project | primary_diagnosis.diagnoses |
|:---|:---|:---|
| Lymphoid Leukemias | Acute Lymphoblastic Leukemia - Phase II | Acute lymphocytic leukemia |
| Lymphoid Leukemias | Acute Lymphoblastic Leukemia - Phase II | Acute lymphocytic leukemia |
| Lymphoid Leukemias | Acute Lymphoblastic Leukemia - Phase II | Acute lymphocytic leukemia |
| Lymphoid Leukemias | Acute Lymphoblastic Leukemia - Phase II | Acute lymphocytic leukemia |
| Lymphoid Leukemias | Acute Lymphoblastic Leukemia - Phase II | Acute lymphocytic leukemia |
| Lymphoid Leukemias | Acute Lymphoblastic Leukemia - Phase II | Acute lymphocytic leukemia |

### Columns excluded from the metadata because they contain uninformative and/or non-clinical information

vital_status (values: Dead, Alive), name.program.project (value =
TARGET),

## Combine target metadata with SRA sheet of samples that were downloaded

``` r
# merge target accessions with clinical metadata 
target_clinical_accessions <- target_clinical_accessions |>
  # select for columns that are relevant to clinical metadata
  dplyr::select(
    # these columns have the same values in other columns
    # I select for the most clearly named column of these, or rename the column with the most values for clarity
    age_at_diagnosis_days = age_at_diagnosis.diagnoses,
    case_id,
    age_at_earliest_diagnosis_in_years.diagnoses.xena_derived, # selected over age_at_index.demographic due to having unrounded numbers
    # the remaining columns are not repeated in the data frame and are selected for clinical relevance
    disease_type,
    primary_site,
    cause_of_death.demographic,
    race.demographic,
    gender.demographic,
    ethnicity.demographic,
    project_id.project,
    name.project,
    tissue_or_organ_of_origin.diagnoses,
    morphology.diagnoses,
    classification_of_tumor.diagnoses,
    icd_10_code.diagnoses,
    primary_diagnosis.diagnoses,
    sample_type_id.samples,
    tumor_descriptor.samples,
    sample_type.samples,
    tumor_code_id.samples,
    sample_type.samples,
    tissue_type.samples,
    inss_stage.diagnoses,
    cog_neuroblastoma_risk_group.diagnoses,
    mitosis_karyorrhexis_index.diagnoses,
    inpc_grade.diagnoses,
    pathology_detail_id.pathology_details.diagnoses,
    necrosis_percent.pathology_details.diagnoses,
    percent_tumor_nuclei.pathology_details.diagnoses,
    pediatric_kidney_staging.diagnoses,
    wilms_tumor_histologic_subtype.diagnoses,
    OS.time,
    OS,
    `_PATIENT`,
    # these columns are from the dbGaP SRA metadata
    Run,
  analyte_type,
  `Assay Type`,
  BioProject,
  BioSample,
  biospecimen_repository,
  biospecimen_repository_sample_id,
  `Center Name`,
  Experiment,
  `dbGaP accession`,
  Instrument,
  `Library Name`,
  LibraryLayout,
  LibrarySelection,
  LibrarySource,
  Organism,
  Platform,
  `Sample Name`,
  `SRA Study`,
  study_design,
  study_name,
  submitted_subject_id,
  Bases,
  AvgSpotLen,
  histological_type,
  body_site
  )

# print head to eyeball contents
target_clinical_accessions |>
  dplyr::slice_sample( n = 5) |>
  head(n = 5)
```

| age_at_diagnosis_days | case_id | age_at_earliest_diagnosis_in_years.diagnoses.xena_derived | disease_type | primary_site | cause_of_death.demographic | race.demographic | gender.demographic | ethnicity.demographic | project_id.project | name.project | tissue_or_organ_of_origin.diagnoses | morphology.diagnoses | classification_of_tumor.diagnoses | icd_10_code.diagnoses | primary_diagnosis.diagnoses | sample_type_id.samples | tumor_descriptor.samples | sample_type.samples | tumor_code_id.samples | tissue_type.samples | inss_stage.diagnoses | cog_neuroblastoma_risk_group.diagnoses | mitosis_karyorrhexis_index.diagnoses | inpc_grade.diagnoses | pathology_detail_id.pathology_details.diagnoses | necrosis_percent.pathology_details.diagnoses | percent_tumor_nuclei.pathology_details.diagnoses | pediatric_kidney_staging.diagnoses | wilms_tumor_histologic_subtype.diagnoses | OS.time | OS | \_PATIENT | Run | analyte_type | Assay Type | BioProject | BioSample | biospecimen_repository | biospecimen_repository_sample_id | Center Name | Experiment | dbGaP accession | Instrument | Library Name | LibraryLayout | LibrarySelection | LibrarySource | Organism | Platform | Sample Name | SRA Study | study_design | study_name | submitted_subject_id | Bases | AvgSpotLen | histological_type | body_site |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 141 | ab0f41dc-eb02-531a-bf58-9eae9bb18e53 | 0.3863013698630137 | Complex Mixed and Stromal Neoplasms | Kidney | NA | white | male | Unknown | TARGET-RT | Rhabdoid Tumor | Kidney, NOS | 8963/3 | primary | C64 | Malignant rhabdoid tumor | 01 | Primary | Primary Tumor | 52 | Tumor | NA | NA | NA | NA | NA | NA | NA | Nephrectomy specimen with tumor that penetrates the renal capsule or involves the renal sinus with negative margins and negative lymph nodes; no distant metastasis | NA | NA | NA | NA | SRR2042838 | RNA | RNA-Seq | PRJNA89539 | SAMN03153988 | NCI_TARGET | TARGET-52-PAJMBW-01A-01R | BCCAGSC | SRX1041134 | phs000470 | Illumina HiSeq 2500 | A37348 | PAIRED | cDNA | TRANSCRIPTOMIC | Homo sapiens | ILLUMINA | TARGET-52-PAJMBW-01A-01R | SRP012005 | Tumor vs. Matched-Normal | TARGET: Kidney, Rhabdoid Tumor (RT) | TARGET-52-PAJMBW | 15931793700 | 150 | Kidney Tumors | Primary Solid Tumor |
| 4782.0 | 0842ae10-4f6b-5c68-9a8a-a24a46acbff3 | 13.101369863013698 | Lymphoid Leukemias | Hematopoietic and reticuloendothelial systems | Not Reported | Unknown | male | hispanic or latino | TARGET-ALL-P2 | Acute Lymphoblastic Leukemia - Phase II | Bone marrow | 9835/3 | primary | C91.0 | Acute lymphocytic leukemia | 04 | Recurrence | Recurrent Blood Derived Cancer - Bone Marrow | 10 | Tumor | NA | NA | NA | NA | NA | NA | NA | NA | NA | 1533.0 | 1 | TARGET-10-PAPSPG | SRR1797122 | RNA | RNA-Seq | PRJNA89529 | SAMN02386053 | NCI_TARGET | TARGET-10-PAPSPG-04A-01R | BCCAGSC | SRX543461 | phs000464 | Illumina HiSeq 2000 | A32678 | PAIRED | cDNA | TRANSCRIPTOMIC | Homo sapiens | ILLUMINA | TARGET-10-PAPSPG-04A-01R | SRP011999 | Tumor vs. Matched-Normal | TARGET: Acute Lymphoblastic Leukemia (ALL) Expansion Phase 2 | TARGET-10-PAPSPG | 24394568550 | 150 | ALL | Recurrent Blood Derived Cancer - Bone Marrow |
| 851.0 | 13b25ed5-cf0d-4a0f-bb1b-a3e200cdaebb | 2.3315068493150686 | Lymphoid Leukemias | Hematopoietic and reticuloendothelial systems | Not Reported | white | female | hispanic or latino | TARGET-ALL-P2 | Acute Lymphoblastic Leukemia - Phase II | Bone marrow | 9835/3 | primary | C91.0 | Acute lymphocytic leukemia | 03 | Primary | Primary Blood Derived Cancer - Peripheral Blood | NA | Tumor | NA | NA | NA | NA | NA | NA | NA | NA | NA | 676.0 | 1 | TARGET-10-PATELK | SRR4416287 | RNA | RNA-Seq | PRJNA89529 | SAMN05784709 | NCI_TARGET | TARGET-10-PATELK-03B-01R | BCCAGSC | SRX2239264 | phs000464 | Illumina HiSeq 2500 | A59330 | PAIRED | cDNA | TRANSCRIPTOMIC | Homo sapiens | ILLUMINA | TARGET-10-PATELK-03B-01R | SRP011999 | Tumor vs. Matched-Normal | TARGET: Acute Lymphoblastic Leukemia (ALL) Expansion Phase 2 | TARGET-10-PATELK | 15350428650 | 150 | ALL | Primary Blood Derived Cancer - Peripheral Blood |
| 1449.0 | bf5addca-ecad-554c-ab71-8f771a88313a | 3.96986301369863 | Lymphoid Leukemias | Hematopoietic and reticuloendothelial systems | NA | black or african american | male | Unknown | TARGET-ALL-P2 | Acute Lymphoblastic Leukemia - Phase II | Bone marrow | 9835/3 | primary | C91.0 | Acute lymphocytic leukemia | 09 | Primary | Primary Blood Derived Cancer - Bone Marrow | 10 | Tumor | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | SRR1797127 | RNA | RNA-Seq | PRJNA89529 | SAMN02386014 | NCI_TARGET | TARGET-10-PARLAF-09A-01R | BCCAGSC | SRX543458 | phs000464 | Illumina HiSeq 2000 | A32627 | PAIRED | cDNA | TRANSCRIPTOMIC | Homo sapiens | ILLUMINA | TARGET-10-PARLAF-09A-01R | SRP011999 | Tumor vs. Matched-Normal | TARGET: Acute Lymphoblastic Leukemia (ALL) Expansion Phase 2 | TARGET-10-PARLAF | 14644995450 | 150 | ALL | Primary Blood Derived Cancer - Bone Marrow |
| 4401.0 | 9d24c3a5-8445-54d8-b87b-c036aad8647f | 12.057534246575342 | Myeloid Leukemias | Hematopoietic and reticuloendothelial systems | NA | white | female | not hispanic or latino | TARGET-AML | Acute Myeloid Leukemia | Bone marrow | 9861/3 | primary | C92.0 | Acute myeloid leukemia, NOS | 03 | Primary | Primary Blood Derived Cancer - Peripheral Blood | 20 | Tumor | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | SRR2239720 | RNA | RNA-Seq | PRJNA89525 | SAMN01778127 | NCI_TARGET | TARGET-20-PATELT-03A-01R | HAIB | SRX1182087 | phs000465 | Illumina HiSeq 2000 | 1170-SL-0058 | PAIRED | cDNA | TRANSCRIPTOMIC | Homo sapiens | ILLUMINA | TARGET-20-PATELT-03A-01R | SRP012000 | Tumor vs. Matched-Normal | TARGET: Acute Myeloid Leukemia (AML) | TARGET-20-PATELT | 2715930336 | 102 | AML | Primary Blood Derived Cancer - Peripheral Blood |

## Write output

``` r
# holding off on this until the code has gone through review
```
