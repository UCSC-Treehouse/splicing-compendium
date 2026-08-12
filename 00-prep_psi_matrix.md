# Prep psi matrix for analysis
Cindy Liang (celiang@ucsc.edu)
2026-08-12

## Introduction

This notebook reads in the PSI tables generated from the Treehouse
Splice Compendium workflow, merges them into one matrix for ease of
downstream analysis, and adds in a gene_name column to associate ENSG
IDs with gene names.

## Setup

Load libraries

``` r
library(rtracklayer)
```

    Loading required package: GenomicRanges

    Loading required package: stats4

    Loading required package: BiocGenerics


    Attaching package: 'BiocGenerics'

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
        table, tapply, union, unique, unsplit, which.max, which.min

    Loading required package: S4Vectors


    Attaching package: 'S4Vectors'

    The following object is masked from 'package:utils':

        findMatches

    The following objects are masked from 'package:base':

        expand.grid, I, unname

    Loading required package: IRanges

    Loading required package: GenomeInfoDb

Define directories and file paths

``` r
# find the root-level repo directory so results files can be accessed
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

## directories ##
# parent directories
exploration_dir <- file.path(repo_root, "exploration")
results_dir <- file.path(repo_root, "results")
references_dir <- file.path(repo_root, "references")

# psi table dir
target_pilot_dir <- file.path(exploration_dir, "merging-psi-tables", "target_pilot")
combined_results_dir <- file.path(target_pilot_dir, "shiba_combined", "results", "splicing")

# output dir of psi table with gene name and splice event annotation status
v1_compendium_resuts_dir <- file.path(results_dir, "v1_compendium")

## files ##
# gtf file for converting ensg id to gene names
gtf_file <- file.path(references_dir, "gencode.v47.primary_assembly.annotation.gtf")

# psi matrix to make nice
psi_matrix_file <- file.path(combined_results_dir, "PSI_matrix_sample.txt")
# output nice psi matrix file
out_file <- file.path(v1_compendium_resuts_dir, "v1_compendium_psi_matrix.txt")

# define list of PSI results files produced by Shiba
psi_files <- c(
  se = "PSI_SE.txt",
  afe = "PSI_AFE.txt",
  ale = "PSI_ALE.txt",
  five = "PSI_FIVE.txt",
  three = "PSI_THREE.txt",
  mse = "PSI_MSE.txt",
  mxe = "PSI_MXE.txt",
  ri = "PSI_RI.txt"
)

# construct psi table paths of combined run shiba results
shiba_psi_paths <-file.path(combined_results_dir, psi_files)
names(shiba_psi_paths) <- names(psi_files)

# Check output dir exists; if not, create it
if (!dir.exists(v1_compendium_resuts_dir)) {
  dir.create(v1_compendium_resuts_dir)
}
```

Read in files

``` r
## splice result files ##
# read merged splice table from v1 workflow on the TARGET pilot samples
psi_results <- shiba_psi_paths |>
  purrr::map(\(path){
    readr::read_tsv(path, col_types=readr::cols(.default = "c")) |>
      # select for columns that will be used in downstream analysis
      dplyr::select(pos_id, gene_id, label, contains("_PSI")) |>
      # convert PSI values to numeric
      dplyr::mutate(across(contains("_PSI"), \(x) as.numeric(x)))
  })

# combine PSI values of samples run together into one dataframe to compare against other PSI tables
psi_table <- purrr::list_rbind(psi_results, names_to = "event_type") |>
  # REMOVE THIS CODE WHEN YOU RUN ON THE V1 FILES! subset to save memory on laptop
  dplyr::slice_sample(n = 10)

## gene annotation file ##
# import gtf
gtf <- import(gtf_file)

# convert gtf to data frame for manipulating fields
gtf_df <- as.data.frame(gtf) |>
  # we only need gene id and gene name columns
  dplyr::select(gene_id,
                gene_name)
```

## Add gene names to psi matrix made from combining all psi tables

Merge PSI matrix with gene info and annotation columns from merged psi
table Each row of this matrix corresponds to a unique splice event.
columns are event_type, pos_id, gene_id, gene_name, label (annotation
status), and the sample-level PSI values (columns start with “SRR”)

``` r
nice_psi_df <- dplyr::left_join(
  psi_table,
  gtf_df,
  by = "gene_id",
  # many splice events are present for each gene
  relationship = "many-to-many"
)

# print column names of non-PSI value columns
colnames(nice_psi_df[,sapply(nice_psi_df,is.character)])
```

    [1] "event_type" "pos_id"     "gene_id"    "label"      "gene_name" 

## Write output

``` r
readr::write_tsv(nice_psi_df, out_file)
```

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
    [1] stats4    stats     graphics  grDevices datasets  utils     methods  
    [8] base     

    other attached packages:
    [1] rtracklayer_1.66.0   GenomicRanges_1.58.0 GenomeInfoDb_1.42.3 
    [4] IRanges_2.40.1       S4Vectors_0.44.0     BiocGenerics_0.52.0 

    loaded via a namespace (and not attached):
     [1] SummarizedExperiment_1.36.0 rjson_0.2.23               
     [3] xfun_0.60                   Biobase_2.66.0             
     [5] lattice_0.22-7              tzdb_0.5.0                 
     [7] vctrs_0.7.3                 tools_4.4.3                
     [9] bitops_1.0-9                generics_0.1.4             
    [11] curl_7.1.0                  parallel_4.4.3             
    [13] tibble_3.3.1                pkgconfig_2.0.3            
    [15] Matrix_1.7-4                lifecycle_1.0.5            
    [17] GenomeInfoDbData_1.2.13     compiler_4.4.3             
    [19] Rsamtools_2.22.0            Biostrings_2.74.1          
    [21] codetools_0.2-20            htmltools_0.5.9            
    [23] RCurl_1.98-1.19             yaml_2.3.12                
    [25] pillar_1.11.1               crayon_1.5.3               
    [27] BiocParallel_1.40.2         DelayedArray_0.32.0        
    [29] abind_1.4-8                 tidyselect_1.2.1           
    [31] digest_0.6.39               dplyr_1.2.1                
    [33] purrr_1.2.2                 restfulr_0.0.17            
    [35] rprojroot_2.1.1             fastmap_1.2.0              
    [37] grid_4.4.3                  cli_3.6.6                  
    [39] SparseArray_1.6.2           magrittr_2.0.5             
    [41] S4Arrays_1.6.0              XML_3.99-0.23              
    [43] withr_3.0.3                 readr_2.2.0                
    [45] UCSC.utils_1.2.0            bit64_4.8.2                
    [47] rmarkdown_2.31              XVector_0.46.0             
    [49] httr_1.4.8                  matrixStats_1.5.0          
    [51] bit_4.6.0                   otel_0.2.0                 
    [53] hms_1.1.4                   evaluate_1.0.5             
    [55] knitr_1.51                  BiocIO_1.16.0              
    [57] rlang_1.3.0                 glue_1.8.1                 
    [59] renv_1.2.3                  vroom_1.7.1                
    [61] jsonlite_2.0.0              R6_2.6.1                   
    [63] MatrixGenerics_1.18.1       GenomicAlignments_1.42.0   
    [65] zlibbioc_1.52.0            
