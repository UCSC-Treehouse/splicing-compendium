# Biology vignette for splice compendium
celiang@ucsc.edu
2026-11-05

This notebook takes Wilcoxon test statistics and dPSI values for query
and reference groups chosen in `scripts/06-diffsplice_calculation.R` as
input. After filtering for splice events with user-defined padj and dPSI
thresholds, this notebook performs over-representation analysis (ORA) on
the parent genes of the filtered splice events.

Caveats:

- The input to this notebook is currently the AML (stranded) vs. Whole
  Blood (unstranded) diffsplice results so doesn’t reflect the
  biological comparison we want to make for the manuscript.

## Setup

Load libraries

``` r
library(ggplot2)
library(msigdbr)
library(clusterProfiler)
```

    clusterProfiler v4.14.6 Learn more at https://yulab-smu.top/contribution-knowledge-mining/

    Please cite:

    Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
    clusterProfiler: an R package for comparing biological themes among
    gene clusters. OMICS: A Journal of Integrative Biology. 2012,
    16(5):284-287


    Attaching package: 'clusterProfiler'

    The following object is masked from 'package:stats':

        filter

``` r
library(enrichplot)
```

    enrichplot v1.26.6 Learn more at https://yulab-smu.top/contribution-knowledge-mining/

    Please cite:

    Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
    clusterProfiler: an R package for comparing biological themes among
    gene clusters. OMICS: A Journal of Integrative Biology. 2012,
    16(5):284-287

Define directory and file paths

``` r
## directories
# repo root
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

## files
# whole blood vs AML splice results with gene name, median PSI, diffsplice dpsi, and p-adj columns
diffsplice_results_file <- file.path(repo_root, params$in_file)

# output ORA results
out_ora_file <- file.path(repo_root, params$out_ora_file)
out_sig_diffsplice_file <- file.path(repo_root, params$out_sig_diffsplice_file)
```

Read in files

``` r
diffsplice_results <- readr::read_tsv(diffsplice_results_file) 
```

    Rows: 440178 Columns: 7
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): pos_id, gene_names
    dbl (5): p, p_adj, median_PSI_query, median_PSI_ref, dPSI

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

## ORA

Obtain msigdb gene sets relevant to humans

``` r
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
```

### Retrieve gene lists from splice results for ORA

Background genes (`detected_genes`) are all genes with non-NA (numeric)
PSI values quantified for at least min_samples in the reference and
query groups. In the test data used here, min_samples is 5.

``` r
detected_genes <- diffsplice_results |>
  # remove NA p values, as they are present in the diffsplice data frame
  dplyr::filter(!is.na(p_adj)) |>
  # extract gene names column for ORA
  dplyr::pull(gene_names) |>
  # take unique set of parent genes
  unique()
```

Filter for differential splicing results within padj and significance
thresholds and obtain the diffsplice gene list

``` r
# filter for significant diffsplice events
filtered_diffsplice_results <- diffsplice_results |>
  dplyr::filter(
    p_adj < params$padj,
    abs(dPSI) > params$dpsi
  ) |>
  # rank significant diff splice events by abs(dPSI)
  dplyr::arrange(desc(abs(dPSI)))

# obtain vector of significant diffsplice parent genes
diffsplice_genes <- unique(filtered_diffsplice_results$gene_names)
```

### ORA on significantly differentially spliced genes

For a first pass, I am not filtering for specific gene sets and applying
a more stringent p filter. In an earlier iteration, I tried to perform
ORA on the splicing gene sets defined by `splice_ids`, but no set was
enriched even when `pvalueCutoff` = 1.

``` r
ora_result <- enricher(
  gene = diffsplice_genes, 
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  universe = detected_genes,
  TERM2GENE = dplyr::select(
    hs_msigdb_df,
    gs_name,
    gene_symbol
  )
)
```

## Plot ORA results

Plot overlapping genes in each gene set with an upset plot

``` r
upsetplot(ora_result) +
  ggtitle("Overlap of genes in ORA gene sets for AML vs Whole Blood")
```

    Warning: `aes_()` was deprecated in ggplot2 3.0.0.
    ℹ Please use tidy evaluation idioms with `aes()`
    ℹ The deprecated feature was likely used in the enrichplot package.
      Please report the issue at
      <https://github.com/GuangchuangYu/enrichplot/issues>.

    Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ℹ Please use `linewidth` instead.
    ℹ The deprecated feature was likely used in the ggupset package.
      Please report the issue at <https://github.com/const-ae/ggupset/issues>.

![](ora_of_diffsplice_genes_files/figure-commonmark/unnamed-chunk-8-1.png)

Plot gene set gene count and significance from ORA: dot plot of gene set
gene abundance (gene ratio) and ORA p values

``` r
dotplot(ora_result)
```

![](ora_of_diffsplice_genes_files/figure-commonmark/unnamed-chunk-9-1.png)

## Write outputs to file

``` r
readr::write_tsv(as.data.frame(ora_result), out_ora_file)
readr::write_tsv(diffsplice_results, out_sig_diffsplice_file)
```
