# Interactive cis-QTL browsing with run_cis_qtl_browser

## Overview

[`run_cis_qtl_browser()`](https://github.com/vjcitn/tQTLExperiment/reference/run_cis_qtl_browser.md)
launches a Shiny application for browsing `tensorQTL` cis-association
results one gene at a time: pick a gene, pick one of its cis-SNPs, and
see the genotype effect plot, the underlying linear model summary, and a
filterable table of every association tested for that gene. It is the
interactive complement to
[`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/articles/tQTLExperiment.html#readtqtl)
and
[`plotGenotypeEffect()`](https://github.com/vjcitn/tQTLExperiment/reference/plotGenotypeEffect.md):
rather than writing `pairs[pairs$gene_name == ...]` by hand for each
gene of interest, the app puts a searchable gene list and SNP dropdown
in front of the same data.

``` r

library(tQTLExperiment)
```

## Inputs

[`run_cis_qtl_browser()`](https://github.com/vjcitn/tQTLExperiment/reference/run_cis_qtl_browser.md)
needs two things:

- `tqe`: a `tQTLExperiment` with genotype data attached (so it can
  extract phenotype and genotype vectors for the plot and model), and
- `cis_gr`: a `GRanges` of cis-association results, typically the
  `pairs` element returned by `readTQTL(..., mode = "cis_nominal")`,
  with metadata columns identifying the gene, SNP, phenotype, and
  nominal p-value for each row (`gene_col`, `snp_col`, `phenotype_col`,
  `pval_col` – their defaults match
  [`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)’s
  and
  [`annotateTensorQTLResults()`](https://github.com/vjcitn/tQTLExperiment/reference/annotateTensorQTLResults.md)’s
  output).

This vignette reuses the pre-computed MAGE chr17 demo results from the
[main
vignette](https://github.com/vjcitn/tQTLExperiment/articles/tQTLExperiment.html#demo-data),
which need no `tensorQTL` installation:

``` r

cc <- cache_demo_data()
#> [ cache hit   ] covariates.tsv
#> [ cache hit   ] mage_chr17_cis.rda
#> [ cache hit   ] pheno.bed
#> [ cache hit   ] tqtl_out.cis_qtl_pairs.17.parquet
#> [ cache hit   ] tqtl_out.tensorQTL.cis_nominal.log
#> [ complete    ] All demo files cached in: /Users/vincentcarey/Library/Caches/org.R-project.R/R/BiocFileCache/tqtlExperiment_demodir
demodir <- dirname(cc$log)

res <- readTQTL(demodir, mode = "cis_nominal")
res$pairs
#> GRanges object with 5991505 ranges and 10 metadata columns:
#>             seqnames    ranges strand |       phenotype_id          variant_id
#>                <Rle> <IRanges>  <Rle> |        <character>         <character>
#>         [1]       17    114101      * |  ENSG00000280279.1       17:114101:G:A
#>         [2]       17    114226      * |  ENSG00000280279.1       17:114226:A:G
#>         [3]       17    116159      * |  ENSG00000280279.1       17:116159:G:C
#>         [4]       17    116270      * |  ENSG00000280279.1       17:116270:G:C
#>         [5]       17    116354      * |  ENSG00000280279.1       17:116354:C:T
#>         ...      ...       ...    ... .                ...                 ...
#>   [5991501]       17  83232567      * | ENSG00000176845.13     17:83232567:A:T
#>   [5991502]       17  83233377      * | ENSG00000176845.13   17:83233377:GAC:G
#>   [5991503]       17  83234128      * | ENSG00000176845.13     17:83234128:C:T
#>   [5991504]       17  83234421      * | ENSG00000176845.13 17:83234421:TAAAA:T
#>   [5991505]       17  83234453      * | ENSG00000176845.13     17:83234453:C:T
#>             start_distance end_distance        af ma_samples  ma_count
#>                  <integer>    <integer> <numeric>  <integer> <integer>
#>         [1]          37236        37235  0.450068        510       658
#>         [2]          37361        37360  0.197674        247       289
#>         [3]          39294        39293  0.450752        515       659
#>         [4]          39405        39404  0.264706        326       387
#>         [5]          39489        39488  0.486320        515       711
#>         ...            ...          ...       ...        ...       ...
#>   [5991501]         152959       152958 0.1080711        158       158
#>   [5991502]         153769       153768 0.0560875         82        82
#>   [5991503]         154520       154519 0.1374829        200       201
#>   [5991504]         154812       154811 0.0937072        135       137
#>   [5991505]         154845       154844 0.1039672        151       152
#>             pval_nominal       slope  slope_se
#>                <numeric>   <numeric> <numeric>
#>         [1]     0.857064 -0.00633339 0.0351502
#>         [2]     0.441229  0.03399045 0.0441110
#>         [3]     0.883984 -0.00517267 0.0354354
#>         [4]     0.329784  0.03938431 0.0403842
#>         [5]     0.805759 -0.00807438 0.0328232
#>         ...          ...         ...       ...
#>   [5991501]    0.2689163   0.0809385 0.0731511
#>   [5991502]    0.2167746  -0.1158684 0.0937229
#>   [5991503]    0.0934929   0.1111321 0.0661664
#>   [5991504]    0.4333309   0.0579827 0.0739608
#>   [5991505]    0.7809999  -0.0197885 0.0711495
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

`gene_symbol` is added by
[`annotateTensorQTLResults()`](https://github.com/vjcitn/tQTLExperiment/reference/annotateTensorQTLResults.md),
which
[`run_cis_qtl_browser()`](https://github.com/vjcitn/tQTLExperiment/reference/run_cis_qtl_browser.md)’s
default `gene_col` expects:

``` r

newpairs <- annotateTensorQTLResults(res$pairs)
```

For the full `tqe` (genotype-attached) construction used below, see the
main vignette’s [“Building a
tQTLExperiment”](https://github.com/vjcitn/tQTLExperiment/articles/tQTLExperiment.html#building-a-tqtlexperiment)
section, which builds `tqe` from `CSHLvc2026::mageSEfilt` and the cached
MAGE chr17 PLINK files via
[`tQTLExperimentFromRSE()`](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperimentFromRSE.md).

## Launching the browser

``` r

if (!requireNamespace("CSHLvc2026")) BiocManager::install("vjcitn/CSHLvc2026")
library(CSHLvc2026)

plink_paths <- cache_mage_chr17_plink()
plpre <- tools::file_path_sans_ext(plink_paths[["bed"]])

data(mageSEfilt, package = "CSHLvc2026")
cd <- as.data.frame(colData(mageSEfilt))
mm <- model.matrix(~ batch + population + sex, data = cd)
mm <- mm[, -1, drop = FALSE]

tqe <- tQTLExperimentFromRSE(
    se              = mageSEfilt,
    plinkPrefix     = plpre,
    covariateMatrix = mm,
    genome          = "hg38"
)
tqe <- addGeneSymbols(tqe)

newpairs <- annotateTensorQTLResults(res$pairs)

# only show associations at pval_nominal <= 0.05
run_cis_qtl_browser(tqe, newpairs, pval_thresh = 0.05)
```

`pval_thresh` is applied once, up front, to `cis_gr` before the app UI
is built – it is a fixed filter on the whole result set, not an
interactive control in the app itself. Leave it `NULL` (the default) to
browse every tested pair.

## Using the app

The sidebar has:

- **Gene symbol**: a searchable dropdown (`selectizeInput`) over every
  distinct value of `gene_col` in `cis_gr` (after the p-value filter),
- **cis-SNP**: repopulated whenever the gene changes, listing that
  gene’s distinct SNPs (from `snp_col`), sorted by nominal p-value,
- **Plot**: an explicit action button – the genotype-effect plot and
  model summary only (re)render when clicked, not on every dropdown
  change, so that selecting a gene and then a SNP doesn’t trigger two
  redraws.

Three tabs share the current gene/SNP selection:

- **Genotype effect**: the beeswarm plot from
  [`plotGenotypeEffect()`](https://github.com/vjcitn/tQTLExperiment/reference/plotGenotypeEffect.md),
  phenotype value by genotype (0/1/2 alternate alleles), for the
  selected phenotype-SNP pair.
- **Regression model**:
  `summary(lm(phenotype ~ genotype + <covariates>))` for that same pair,
  printed verbatim – the same model
  [`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
  fits per pair, but shown in full rather than reduced to a tidy table.
  Useful for checking covariate significance or residual behavior behind
  a specific point.
- **cis associations**: a
  [`DT::datatable()`](https://rdrr.io/pkg/DT/man/datatable.html) of
  every row of `cis_gr` for the selected gene (ignoring the SNP
  dropdown), sorted by nominal p-value, with per-column filters
  (`filter = "top"`) for further narrowing by hand.

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
#> [34] RColorBrewer_1.1-3   jquerylib_0.1.4      assertthat_0.2.1    
#> [37] Rcpp_1.1.2           bookdown_0.48        knitr_1.52          
#> [40] BiocBaseUtils_1.15.1 httpuv_1.6.17        Matrix_1.7-6        
#> [43] tidyselect_1.2.1     dichromat_2.0-1      abind_1.4-8         
#> [46] yaml_2.3.12          codetools_0.2-20     curl_8.1.0          
#> [49] lattice_0.23-1       tibble_3.3.1         shiny_1.14.0        
#> [52] KEGGREST_1.53.6      S7_0.2.2             evaluate_1.0.5      
#> [55] desc_1.4.3           BiocFileCache_3.3.1  Biostrings_2.81.9   
#> [58] pillar_1.11.1        BiocManager_1.30.27  filelock_1.0.3      
#> [61] DT_0.34.0            plotly_4.12.1        ggplot2_4.0.3       
#> [64] scales_1.4.0         xtable_1.8-8         glue_1.8.1          
#> [67] tools_4.6.1          data.table_1.18.6.1  fs_2.1.0            
#> [70] grid_4.6.1           tidyr_1.3.2          AnnotationDbi_1.75.2
#> [73] beeswarm_0.4.0       vipor_0.4.7          cli_3.6.6           
#> [76] textshaping_1.0.5    S4Arrays_1.13.2      viridisLite_0.4.3   
#> [79] arrow_26.0.0         dplyr_1.2.1          gtable_0.3.6        
#> [82] sass_0.4.10          digest_0.6.39        SparseArray_1.13.4  
#> [85] htmlwidgets_1.6.4    farver_2.1.2         memoise_2.0.1       
#> [88] htmltools_0.5.9      pkgdown_2.2.1        lifecycle_1.0.5     
#> [91] httr_1.4.9           mime_0.13            bit64_4.8.6
```
