# Add clinical metadata to TARGET compendium v1 sample sheet
Cindy Liang (celiang@ucsc.edu)
2026-06-07

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

- Something with a single number or value per column so it can be
  computed on easily

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
  })

target_survival_metadata <- target_survival_metadata_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default="c"))
  })

# combine clinical and survival metadata 
combined_target_clinical <- purrr::list_rbind(target_clinical_metadata)

combined_target_survival <- purrr::list_rbind(target_survival_metadata)
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
target_metadata <- combined_target_survival |>
  dplyr::full_join(combined_target_clinical, by="sample")

# make a list of target accessions that were downloaded and processed
target_downloaded <- compendium_samples |>
  dplyr::filter(group == "target") |>
  dplyr::pull(sample)

# filter target accessions for accessions that are in the sample sheet
target_accessions <- target_sra |>
  dplyr::filter(Run %in% target_downloaded) |>
  # select for relevant columns
  dplyr::select(
    Run, biospecimen_repository, `DATASTORE filetype`, body_site, Bytes,
    study_name, `DATASTORE provider`, analyte_type, `Assay Type`, `SRA Study`,
    BioProject, BioSample, LibraryLayout, LibrarySelection, `Center Name`,
    biospecimen_repository_sample_id
  ) |>
  # derive patient ID from sample ID 
  dplyr::mutate(
    # remove the hyphen followed by characters that are not hyphens that are at the end of the string
    patient_id = sub("-[^-]+$", "", biospecimen_repository_sample_id)
  )
  
# Check target metadata contents
colnames(target_metadata)
```

     [1] "sample"                                                   
     [2] "OS.time"                                                  
     [3] "OS"                                                       
     [4] "_PATIENT"                                                 
     [5] "id"                                                       
     [6] "disease_type"                                             
     [7] "case_id"                                                  
     [8] "submitter_id"                                             
     [9] "primary_site"                                             
    [10] "cause_of_death.demographic"                               
    [11] "race.demographic"                                         
    [12] "gender.demographic"                                       
    [13] "ethnicity.demographic"                                    
    [14] "vital_status.demographic"                                 
    [15] "age_at_index.demographic"                                 
    [16] "days_to_birth.demographic"                                
    [17] "age_is_obfuscated.demographic"                            
    [18] "days_to_death.demographic"                                
    [19] "primary_site.project"                                     
    [20] "project_id.project"                                       
    [21] "disease_type.project"                                     
    [22] "name.project"                                             
    [23] "name.program.project"                                     
    [24] "entity_submitter_id.annotations"                          
    [25] "notes.annotations"                                        
    [26] "submitter_id.annotations"                                 
    [27] "classification.annotations"                               
    [28] "entity_id.annotations"                                    
    [29] "created_datetime.annotations"                             
    [30] "annotation_id.annotations"                                
    [31] "entity_type.annotations"                                  
    [32] "updated_datetime.annotations"                             
    [33] "case_id.annotations"                                      
    [34] "state.annotations"                                        
    [35] "category.annotations"                                     
    [36] "status.annotations"                                       
    [37] "case_submitter_id.annotations"                            
    [38] "tissue_or_organ_of_origin.diagnoses"                      
    [39] "age_at_diagnosis.diagnoses"                               
    [40] "morphology.diagnoses"                                     
    [41] "classification_of_tumor.diagnoses"                        
    [42] "icd_10_code.diagnoses"                                    
    [43] "days_to_diagnosis.diagnoses"                              
    [44] "primary_diagnosis.diagnoses"                              
    [45] "year_of_diagnosis.diagnoses"                              
    [46] "diagnosis_is_primary_disease.diagnoses"                   
    [47] "site_of_resection_or_biopsy.diagnoses"                    
    [48] "age_at_earliest_diagnosis.diagnoses.xena_derived"         
    [49] "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"
    [50] "protocol_identifier.treatments.diagnoses"                 
    [51] "updated_datetime.treatments.diagnoses"                    
    [52] "treatment_id.treatments.diagnoses"                        
    [53] "submitter_id.treatments.diagnoses"                        
    [54] "state.treatments.diagnoses"                               
    [55] "treatment_or_therapy.treatments.diagnoses"                
    [56] "created_datetime.treatments.diagnoses"                    
    [57] "sample_type_id.samples"                                   
    [58] "tumor_descriptor.samples"                                 
    [59] "sample_id.samples"                                        
    [60] "sample_type.samples"                                      
    [61] "tumor_code.samples"                                       
    [62] "preservation_method.samples"                              
    [63] "freezing_method.samples"                                  
    [64] "tumor_code_id.samples"                                    
    [65] "oct_embedded.samples"                                     
    [66] "specimen_type.samples"                                    
    [67] "tissue_type.samples"                                      
    [68] "last_known_disease_status.diagnoses"                      
    [69] "days_to_last_follow_up.diagnoses"                         
    [70] "tumor_grade.diagnoses"                                    
    [71] "progression_or_recurrence.diagnoses"                      
    [72] "sites_of_involvement.diagnoses"                           
    [73] "timepoint_category.treatments.diagnoses"                  
    [74] "treatment_type.treatments.diagnoses"                      
    [75] "course_number.treatments.diagnoses"                       
    [76] "reason_treatment_ended.treatments.diagnoses"              
    [77] "treatment_outcome.treatments.diagnoses"                   
    [78] "therapeutic_agents.treatments.diagnoses"                  
    [79] "annotations.samples"                                      
    [80] "inss_stage.diagnoses"                                     
    [81] "cog_neuroblastoma_risk_group.diagnoses"                   
    [82] "mitosis_karyorrhexis_index.diagnoses"                     
    [83] "inpc_grade.diagnoses"                                     
    [84] "pathology_detail_id.pathology_details.diagnoses"          
    [85] "updated_datetime.pathology_details.diagnoses"             
    [86] "submitter_id.pathology_details.diagnoses"                 
    [87] "state.pathology_details.diagnoses"                        
    [88] "created_datetime.pathology_details.diagnoses"             
    [89] "necrosis_percent.pathology_details.diagnoses"             
    [90] "percent_tumor_nuclei.pathology_details.diagnoses"         
    [91] "pediatric_kidney_staging.diagnoses"                       
    [92] "wilms_tumor_histologic_subtype.diagnoses"                 

There are more fields in the metadata than we planned to include, so
some discussion on what may or may not be useful as a group will be
helpful. For now, I won’t exclude any metadata and will just focus on
merging the metadata with the accessions sheet

Combine target metadata with SRA sheet of samples that were downloaded

``` r
# merge target accessions with clinical metadata 
target_clinical_accessions <- target_accessions |>
  dplyr::left_join(target_metadata, by = dplyr::join_by(patient_id == sample))

# print head to eyeball contents
target_clinical_accessions |>
  dplyr::slice_sample( n = 5) |>
  head(n = 5)
```

| Run | biospecimen_repository | DATASTORE filetype | body_site | Bytes | study_name | DATASTORE provider | analyte_type | Assay Type | SRA Study | BioProject | BioSample | LibraryLayout | LibrarySelection | Center Name | biospecimen_repository_sample_id | patient_id | OS.time | OS | \_PATIENT | id | disease_type | case_id | submitter_id | primary_site | cause_of_death.demographic | race.demographic | gender.demographic | ethnicity.demographic | vital_status.demographic | age_at_index.demographic | days_to_birth.demographic | age_is_obfuscated.demographic | days_to_death.demographic | primary_site.project | project_id.project | disease_type.project | name.project | name.program.project | entity_submitter_id.annotations | notes.annotations | submitter_id.annotations | classification.annotations | entity_id.annotations | created_datetime.annotations | annotation_id.annotations | entity_type.annotations | updated_datetime.annotations | case_id.annotations | state.annotations | category.annotations | status.annotations | case_submitter_id.annotations | tissue_or_organ_of_origin.diagnoses | age_at_diagnosis.diagnoses | morphology.diagnoses | classification_of_tumor.diagnoses | icd_10_code.diagnoses | days_to_diagnosis.diagnoses | primary_diagnosis.diagnoses | year_of_diagnosis.diagnoses | diagnosis_is_primary_disease.diagnoses | site_of_resection_or_biopsy.diagnoses | age_at_earliest_diagnosis.diagnoses.xena_derived | age_at_earliest_diagnosis_in_years.diagnoses.xena_derived | protocol_identifier.treatments.diagnoses | updated_datetime.treatments.diagnoses | treatment_id.treatments.diagnoses | submitter_id.treatments.diagnoses | state.treatments.diagnoses | treatment_or_therapy.treatments.diagnoses | created_datetime.treatments.diagnoses | sample_type_id.samples | tumor_descriptor.samples | sample_id.samples | sample_type.samples | tumor_code.samples | preservation_method.samples | freezing_method.samples | tumor_code_id.samples | oct_embedded.samples | specimen_type.samples | tissue_type.samples | last_known_disease_status.diagnoses | days_to_last_follow_up.diagnoses | tumor_grade.diagnoses | progression_or_recurrence.diagnoses | sites_of_involvement.diagnoses | timepoint_category.treatments.diagnoses | treatment_type.treatments.diagnoses | course_number.treatments.diagnoses | reason_treatment_ended.treatments.diagnoses | treatment_outcome.treatments.diagnoses | therapeutic_agents.treatments.diagnoses | annotations.samples | inss_stage.diagnoses | cog_neuroblastoma_risk_group.diagnoses | mitosis_karyorrhexis_index.diagnoses | inpc_grade.diagnoses | pathology_detail_id.pathology_details.diagnoses | updated_datetime.pathology_details.diagnoses | submitter_id.pathology_details.diagnoses | state.pathology_details.diagnoses | created_datetime.pathology_details.diagnoses | necrosis_percent.pathology_details.diagnoses | percent_tumor_nuclei.pathology_details.diagnoses | pediatric_kidney_staging.diagnoses | wilms_tumor_histologic_subtype.diagnoses |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| SRR5115239 | NCI_TARGET | fastq,run.zq,sra | Primary Blood Derived Cancer - Bone Marrow | 7979059451 | TARGET: Acute Lymphoblastic Leukemia (ALL) Expansion Phase 2 | gs,ncbi,s3 | RNA | RNA-Seq | SRP011999 | PRJNA89529 | SAMN05784924 | PAIRED | cDNA | BCCAGSC | TARGET-10-PASDYM-09B-01R | TARGET-10-PASDYM-09B | NA | NA | NA | 49f90981-996a-4742-acbc-365382d4851c | Lymphoid Leukemias | 49f90981-996a-4742-acbc-365382d4851c | TARGET-10-PASDYM | Hematopoietic and reticuloendothelial systems | NA | white | male | not hispanic or latino | Alive | 19.0 | -7209.0 | False | NA | Hematopoietic and reticuloendothelial systems | TARGET-ALL-P2 | Lymphoid Leukemias | Acute Lymphoblastic Leukemia - Phase II | TARGET | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | Bone marrow | 7209.0 | 9835/3 | primary | C91.0 | 0.0 | Acute lymphocytic leukemia | 2008 | True | Not Reported | 7209.0 | 19.75068493150685 | AALL0232 | 2024-03-14T10:56:04.835747-05:00 | 7d216b00-b0f7-4bb7-988c-2914c6ee9af0 | TARGET-10-PASDYM_treatment | released | yes | 2023-07-27T15:33:27.373144-05:00 | 09 | Primary | 9b0a02d5-d559-440a-88b9-d65f749e2dfc | Primary Blood Derived Cancer - Bone Marrow | NA | Unknown | NA | NA | NA | Bone Marrow NOS | Tumor | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| SRR1799049 | NCI_TARGET | run.zq,fastq,sra | Primary Solid Tumor | 11925342718 | TARGET: Kidney, Rhabdoid Tumor (RT) | gs,ncbi,s3 | RNA | RNA-Seq | SRP012005 | PRJNA89539 | SAMN01778470 | PAIRED | cDNA | BCCAGSC | TARGET-52-PARPFY-01A-01R | TARGET-52-PARPFY-01A | 102 | 1 | TARGET-52-PARPFY | 5cdd05ea-5285-50b7-971a-8bc005d01669 | Complex Mixed and Stromal Neoplasms | 5cdd05ea-5285-50b7-971a-8bc005d01669 | TARGET-52-PARPFY | Kidney | Toxicity | white | male | not hispanic or latino | Dead | 0 | -361 | False | 102.0 | \[‘Lip’, ‘Liver and intrahepatic bile ducts’, ‘Kidney’\] | TARGET-RT | Complex Mixed and Stromal Neoplasms | Rhabdoid Tumor | TARGET | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | Kidney, NOS | 361 | 8963/3 | primary | C64 | 0 | Malignant rhabdoid tumor | 2007 | True | Not Reported | 361.0 | 0.989041095890411 | AREN03B2 | 2024-03-14T11:20:31.137478-05:00 | d73693c7-8880-4001-8b8e-e9b367d3246a | TARGET-52-PARPFY_treatment | released | yes | 2023-03-28T19:15:24.199187-05:00 | 01 | Primary | ae6b9660-eb5e-55cf-8ef0-282a49904ead | Primary Tumor | Rhabdoid tumor (kidney) (RT) | Unknown | None | 52 | None | Solid Tissue | Tumor | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | Nephrectomy specimen with tumor that is present at the surgical margin of resection or within regional lymph nodes; no distant metastasis | NA |
| SRR1821181 | NCI_TARGET | sra,run.zq,fastq | Primary Blood Derived Cancer - Bone Marrow | 12669745334 | TARGET: Acute Myeloid Leukemia (AML) | s3,ncbi,gs | RNA | RNA-Seq | SRP012000 | PRJNA89525 | SAMN02443372 | PAIRED | cDNA | BCCAGSC | TARGET-20-PANKFG-09A-02R | TARGET-20-PANKFG-09A | NA | NA | NA | 274a8a17-65ce-5e4b-85a1-531cf5bf6191 | Myeloid Leukemias | 274a8a17-65ce-5e4b-85a1-531cf5bf6191 | TARGET-20-PANKFG | Hematopoietic and reticuloendothelial systems | NA | white | male | not hispanic or latino | Alive | 1.0 | -408.0 | False | NA | \[‘Unknown’, ‘Hematopoietic and reticuloendothelial systems’\] | TARGET-AML | \[‘Not Applicable’, ‘Myeloid Leukemias’\] | Acute Myeloid Leukemia | TARGET | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | Bone marrow | 408.0 | 9861/3 | primary | C92.0 | 0.0 | Acute myeloid leukemia, NOS | 2004 | True | Not Reported | 408.0 | 1.1178082191780823 | \[‘AAML03P1’, ’‘,’‘,’‘,’’\] | \[‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’\] | \[‘6ac5706a-fe14-44dd-b29a-6add70768739’, ‘a3150459-acad-4be7-bdd3-a054b9906b5b’, ‘ae790713-cda3-44ba-8318-045ecf572146’, ‘b1bd0986-cec8-4539-ab1d-c579f1a70c5f’, ‘d5eebe41-5bf2-4fea-b166-83913ec6ed81’\] | \[‘TARGET-20-PANKFG_treatment’, ‘TARGET-20-PANKFG_treatment3’, ‘TARGET-20-PANKFG_treatment2’, ‘TARGET-20-PANKFG_treatment4’, ‘TARGET-20-PANKFG_treatment5’\] | \[‘released’, ‘released’, ‘released’, ‘released’, ‘released’\] | \[‘yes’, ‘yes’, ‘no’, ‘yes’, ‘yes’\] | \[‘2023-07-20T17:22:53.780014-05:00’, ‘2023-07-20T17:22:53.780014-05:00’, ‘2023-07-20T17:22:53.780014-05:00’, ‘2023-07-20T17:22:53.780014-05:00’, ‘2023-07-20T17:22:53.780014-05:00’\] | 09 | Primary | 05245329-cd95-50de-972b-c72cbf4dedd1 | Primary Blood Derived Cancer - Bone Marrow | Acute myeloid leukemia (AML) | Unknown | None | 20 | None | Bone Marrow NOS | Tumor | NA | NA | NA | NA | NA | \[’‘, ’End of Treatment Course’, ‘First Complete Response’, ‘End of Treatment Course’, ’’\] | \[’‘,’‘, ’Stem Cell Transplantation, NOS’, ’‘, ’Pharmaceutical Therapy, NOS’\] | \[’‘, ’1.0’, ’‘, ’2.0’, ’’\] | \[’‘, ’Course of Therapy Completed’, ’‘, ’Course of Therapy Completed’, ’’\] | \[’‘, ’Complete Response’, ’‘, ’Complete Response’, ’’\] | \[’‘,’‘,’‘,’‘, ’Gemtuzumab Ozogamicin’\] | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| SRR1797068 | NCI_TARGET | fastq,run.zq,sra | Recurrent Blood Derived Cancer - Bone Marrow | 12679474067 | TARGET: Acute Lymphoblastic Leukemia (ALL) Expansion Phase 2 | s3,gs,ncbi | RNA | RNA-Seq | SRP011999 | PRJNA89529 | SAMN02386114 | PAIRED | cDNA | BCCAGSC | TARGET-10-PAPRCS-04A-01R | TARGET-10-PAPRCS-04A | 460.0 | 1 | TARGET-10-PAPRCS | 75018b77-b741-5c40-9d73-ecf31599e048 | Lymphoid Leukemias | 75018b77-b741-5c40-9d73-ecf31599e048 | TARGET-10-PAPRCS | Hematopoietic and reticuloendothelial systems | Not Reported | white | male | not hispanic or latino | Dead | 14.0 | -5147.0 | False | 460.0 | Hematopoietic and reticuloendothelial systems | TARGET-ALL-P2 | Lymphoid Leukemias | Acute Lymphoblastic Leukemia - Phase II | TARGET | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | Bone marrow | 5147.0 | 9835/3 | primary | C91.0 | 0.0 | Acute lymphocytic leukemia | 2006 | True | Not Reported | 5147.0 | 14.101369863013698 | AALL0232 | 2024-03-14T10:56:04.835747-05:00 | 4fccd31b-8f54-48cc-b24b-1178734c2936 | TARGET-10-PAPRCS_treatment | released | yes | 2023-07-27T13:10:22.144146-05:00 | 04 | Recurrence | 26dae768-2507-5ed0-a5ee-eb21181c3546 | Recurrent Blood Derived Cancer - Bone Marrow | Acute lymphoblastic leukemia (ALL) | Unknown | None | 10 | None | Bone Marrow NOS | Tumor | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| SRR2083153 | NCI_TARGET | fastq,run.zq,sra | Primary Blood Derived Cancer - Peripheral Blood | 7244711912 | TARGET: Acute Myeloid Leukemia (AML) | s3,gs,ncbi | RNA | RNA-Seq | SRP012000 | PRJNA89525 | SAMN01778160 | PAIRED | cDNA | BCCAGSC | TARGET-20-PARBFZ-03A-02R | TARGET-20-PARBFZ-03A | NA | NA | NA | 79a1df12-4b16-58cc-8a81-544993c3877e | Myeloid Leukemias | 79a1df12-4b16-58cc-8a81-544993c3877e | TARGET-20-PARBFZ | Hematopoietic and reticuloendothelial systems | NA | Unknown | female | hispanic or latino | Alive | 14.0 | -5298.0 | False | NA | \[‘Unknown’, ‘Hematopoietic and reticuloendothelial systems’\] | TARGET-AML | \[‘Not Applicable’, ‘Myeloid Leukemias’\] | Acute Myeloid Leukemia | TARGET | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | Bone marrow | 5298.0 | 9861/3 | primary | C92.0 | 0.0 | Acute myeloid leukemia, NOS | 2007 | True | Not Reported | 5298.0 | 14.515068493150684 | \[’‘,’‘,’‘,’‘, ’AAML0531’\] | \[‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’, ‘2024-03-14T11:00:25.965662-05:00’\] | \[‘00e45a3a-eb49-4618-b299-d5fa71e0fdec’, ‘066aaaf5-dce3-469a-b347-a6d6ce6fd0f5’, ‘71ee567d-e514-43b5-a502-4600524a51c3’, ‘e4c13f13-4592-4b82-921a-9d70e2d41d43’, ‘fe16b0e8-65c8-4050-8d4a-ec54a5fbe2a7’\] | \[‘TARGET-20-PARBFZ_treatment3’, ‘TARGET-20-PARBFZ_treatment5’, ‘TARGET-20-PARBFZ_treatment2’, ‘TARGET-20-PARBFZ_treatment4’, ‘TARGET-20-PARBFZ_treatment’\] | \[‘released’, ‘released’, ‘released’, ‘released’, ‘released’\] | \[‘yes’, ‘yes’, ‘no’, ‘yes’, ‘yes’\] | \[‘2023-07-20T18:01:48.268823-05:00’, ‘2023-07-20T18:01:48.268823-05:00’, ‘2023-07-20T18:01:48.268823-05:00’, ‘2023-07-20T18:01:48.268823-05:00’, ‘2023-07-20T18:01:48.268823-05:00’\] | 03 | Primary | 6be0bf6e-1801-5b76-9395-28ab21a9c0b8 | Primary Blood Derived Cancer - Peripheral Blood | Acute myeloid leukemia (AML) | Unknown | None | 20 | None | Peripheral Blood NOS | Tumor | NA | NA | NA | NA | Central nervous system | \[‘End of Treatment Course’, ’‘, ’First Complete Response’, ‘End of Treatment Course’, ’’\] | \[’‘, ’Pharmaceutical Therapy, NOS’, ‘Stem Cell Transplantation, NOS’, ’‘,’’\] | \[‘1.0’, ’‘,’‘, ’2.0’, ’’\] | \[‘Course of Therapy Completed’, ’‘,’‘, ’Course of Therapy Completed’, ’’\] | \[‘Unknown’, ’‘,’‘, ’Complete Response’, ’’\] | \[’‘, ’Gemtuzumab Ozogamicin’, ’‘,’‘,’’\] | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

Not all samples have all clinical/survival metadata fields.

## Write output

``` r
# holding off on this until the code has gone through review
```
