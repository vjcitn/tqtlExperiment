# PCA exploration of cis-QTL t-statistics with qtlPCABrowser

## Overview

[`qtlPCABrowser()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlPCABrowser.md)
launches a Shiny application that runs PCA on the per-covariate
t-statistics produced by
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/articles/regressionStats.md)
for a set of SNP-gene pairs, and lets you explore the result
interactively:

- a **PCA scatter** tab (`plotly`) with selectable PC axes; hovering
  over a point shows a beeswarm of that pair’s gene expression by
  genotype, in the panel to the right, optionally colored by a
  sample-level covariate,
- a **PC pairs** tab, a base-R
  [`pairs()`](https://rdrr.io/r/graphics/pairs.html) plot of the first
  five PCs.

Because each row of the t-statistic matrix is one SNP-gene pair and each
column is one covariate term, PC1/PC2 of this matrix summarize how pairs
differ in the *pattern* of covariate effects driving their association –
useful for spotting pairs whose signal is dominated by a batch or
population effect rather than genotype.

``` r

library(tQTLExperiment)
library(SummarizedExperiment)
```

## Two ways to launch it

[`qtlPCABrowser()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlPCABrowser.md)
has two modes, selected by which argument you supply.

### Gene-symbol mode (`tqe` only)

Supply just a `tQTLExperiment` with gene symbols already added
([`addGeneSymbols()`](https://github.com/vjcitn/tQTLExperiment/reference/addGeneSymbols.md)).
The app shows a gene-symbol selector; on “Compute regression stats/PCA”
it calls `qtlRegressionStats(tqe, symbol = ...)` for the chosen gene,
then PCA’s the result. Repeat calls for the same gene are memoised
within the Shiny session.

``` r

exdir <- system.file("extdata", package = "tQTLExperiment")
tqe <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
    covFile     = file.path(exdir, "cov-n100-tqtl.tsv"),
    genome      = "hg38"
)
tqe <- addGeneSymbols(tqe)
qtlPCABrowser(tqe)
```

This mode is the right default when you don’t yet know which genes are
interesting: compute on demand, gene by gene.

### Pre-computed mode (`res` supplied)

If you already have a wide
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
result – for instance for every gene in a region, computed ahead of time
and possibly in parallel – pass it as `res`; the gene-symbol selector is
skipped and the app opens directly on the combined PCA.

``` r

symbols <- na.omit(unique(mcols(rowRanges(tqe))[["gene_name"]]))
res_list <- lapply(symbols, function(s)
    tryCatch(qtlRegressionStats(tqe, symbol = s), error = function(e) NULL))
res_all <- do.call(rbind, res_list)
qtlPCABrowser(res_all, tqe = tqe)
```

`res` must contain `phenotype_id` and `variant_id` columns (the wide,
`t_only = TRUE` shape
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
returns by default).

## What happens to the t-statistic columns

Before PCA,
[`qtlPCABrowser()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlPCABrowser.md)
runs
[`collapseFactorTstats()`](https://github.com/vjcitn/tQTLExperiment/reference/collapseFactorTstats.md)
on `res` (or on each gene’s freshly computed result), using the
`collapse_patterns` argument – by default tuned for the MAGE dataset’s
`factor_batch_*` and `population*` dummy columns. This keeps a factor
with many levels from dominating the PCA the way
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
alone would, since it would otherwise contribute one raw t-statistic
column per level. Pass your own `collapse_patterns` list if your
covariates use different names – see
[`?collapseFactorTstats`](https://github.com/vjcitn/tQTLExperiment/reference/collapseFactorTstats.md).

``` r

qtlPCABrowser(
    tqe,
    collapse_patterns = list(
        site = "^t\\.site",
        sex  = "^t\\.sex"
    )
)
```

PCA itself is [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)
on the resulting t-statistic matrix, centered but not scaled
(`center = TRUE, scale. = FALSE`), so columns with naturally larger
t-statistics – typically genotype – dominate the leading PCs, which is
usually the point.

## Interpreting the scatter-plus-beeswarm view

Hovering a point in the PCA scatter plot identifies a SNP-gene pair,
whose genotype effect is then drawn in the beeswarm panel using the same
data extraction as
[`plotGenotypeEffect()`](https://github.com/vjcitn/tQTLExperiment/reference/plotGenotypeEffect.md):
genotype (0/1/2 alternate alleles) on the x-axis, phenotype value on the
y-axis. The “Color beeswarm by” selector reconstructs the original
factor levels from the covariate dummy columns (e.g. `batch`,
`population`) plus any remaining numeric covariates, so you can
immediately check whether an outlying PCA point corresponds to a pair
whose apparent effect is really driven by one covariate level rather
than genotype.

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Golden Gate 27.0.1
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] tQTLExperiment_0.1.50       SummarizedExperiment_1.43.0
#>  [3] Biobase_2.73.2              GenomicRanges_1.65.4       
#>  [5] Seqinfo_1.3.2               IRanges_2.47.5             
#>  [7] S4Vectors_0.51.10           BiocGenerics_0.59.12       
#>  [9] generics_0.1.4              MatrixGenerics_1.25.0      
#> [11] matrixStats_1.5.0           BiocStyle_2.41.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] DBI_1.3.0            httr2_1.3.0          rlang_1.3.0         
#>  [4] magrittr_2.0.5       otel_0.2.0           compiler_4.6.1      
#>  [7] RSQLite_3.53.3       png_0.1-9            systemfonts_1.3.2   
#> [10] vctrs_0.7.3          pkgconfig_2.0.3      crayon_1.5.3        
#> [13] fastmap_1.2.0        dbplyr_2.6.0         XVector_0.53.0      
#> [16] promises_1.5.0       rmarkdown_2.32       UCSC.utils_1.9.0    
#> [19] ggbeeswarm_0.7.3     ragg_1.5.2           purrr_1.2.2         
#> [22] bit_4.6.0            xfun_0.62            cachem_1.1.0        
#> [25] GenomeInfoDb_1.49.1  jsonlite_2.0.0       blob_1.3.0          
#> [28] later_1.4.8          DelayedArray_0.39.8  BiocParallel_1.47.0 
#> [31] parallel_4.6.1       R6_2.6.1             bslib_0.12.0        
#> [34] RColorBrewer_1.1-3   jquerylib_0.1.4      Rcpp_1.1.2          
#> [37] bookdown_0.48        knitr_1.52           BiocBaseUtils_1.15.1
#> [40] httpuv_1.6.17        Matrix_1.7-6         tidyselect_1.2.1    
#> [43] dichromat_2.0-1      abind_1.4-8          yaml_2.3.12         
#> [46] codetools_0.2-20     curl_8.1.0           lattice_0.23-1      
#> [49] tibble_3.3.1         shiny_1.14.0         KEGGREST_1.53.6     
#> [52] S7_0.2.2             evaluate_1.0.5       desc_1.4.3          
#> [55] BiocFileCache_3.3.1  Biostrings_2.81.9    pillar_1.11.1       
#> [58] BiocManager_1.30.27  filelock_1.0.3       DT_0.34.0           
#> [61] plotly_4.12.1        ggplot2_4.0.3        scales_1.4.0        
#> [64] xtable_1.8-8         glue_1.8.1           tools_4.6.1         
#> [67] data.table_1.18.6.1  fs_2.1.0             grid_4.6.1          
#> [70] tidyr_1.3.2          AnnotationDbi_1.75.2 beeswarm_0.4.0      
#> [73] vipor_0.4.7          cli_3.6.6            textshaping_1.0.5   
#> [76] S4Arrays_1.13.2      viridisLite_0.4.3    dplyr_1.2.1         
#> [79] gtable_0.3.6         sass_0.4.10          digest_0.6.39       
#> [82] SparseArray_1.13.4   htmlwidgets_1.6.4    farver_2.1.2        
#> [85] memoise_2.0.1        htmltools_0.5.9      pkgdown_2.2.1       
#> [88] lifecycle_1.0.5      httr_1.4.9           mime_0.13           
#> [91] bit64_4.8.6
```
