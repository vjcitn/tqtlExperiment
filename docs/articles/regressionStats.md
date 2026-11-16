# Per-pair regression statistics with qtlRegressionStats

## Overview

[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
fits, for each requested SNP-gene pair, a linear model of phenotype on
genotype plus all covariates stored in the `tQTLExperiment`’s `colData`:

    phenotype ~ genotype + <covariates>

It exists alongside `tensorQTL`, not as a replacement for it:
`tensorQTL` is the tool for mapping cis-QTLs genome-wide, while
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
gives full [`lm()`](https://rdrr.io/r/stats/lm.html) output
(coefficients, standard errors, t-statistics for every covariate, not
just genotype) for a targeted set of pairs, in-process in R. This is
useful for follow-up on specific genes, for building the input to
[`qtlPCABrowser()`](https://github.com/vjcitn/tQTLExperiment/articles/pcaBrowser.md),
and any time you want the whole regression table rather than tensorQTL’s
genotype-only summary.

Pairs can be processed in parallel across phenotypes (or, in gene mode,
across variants) via
[`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html).

``` r

library(tQTLExperiment)
library(GenomicRanges)
library(S4Vectors)
```

## Example data

The small chr22 example dataset (20 genes, 100 samples) ships with the
package and needs no external downloads.

``` r

exdir <- system.file("extdata", package = "tQTLExperiment")
tqe <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
    covFile     = file.path(exdir, "cov-n100-tqtl.tsv"),
    genome      = "hg38"
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
tqe
#> class: tQTLExperiment
#> features: 20  samples: 100 
#> assays( 1 ): pheno 
#> rowRanges: GRanges with 20 features
#> colData( 11 ) covariates: int, sex, age, expr_pc1 ... 
#> geno: 100 samples x 69638 variants [BEDMatrix - lazy]
#> plinkPrefix: /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpcbUgAS/temp_libpath2147cb89ee3/tQTLExperiment/extdata/chr22-n100 
#> use prepareTQTL() to write inputs and get CLI command
colnames(colData(tqe))
#>  [1] "int"       "sex"       "age"       "expr_pc1"  "expr_pc2"  "geno_pc1" 
#>  [7] "geno_pc2"  "geno_pc3"  "geno_pc4"  "geno_pc5"  "geno_pc6"  "fam_index"
```

The covariates here are a mix of a fixed effect (`sex`), a continuous
covariate (`age`), and expression/genotype PCs – no factor-expanded
columns, so
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)’s
output has one column per covariate rather than one per factor level.

## Two ways to specify pairs

### Explicit pairs

Supply a data frame with `phenotype_id` and `variant_id` columns naming
exactly the pairs to fit.

``` r

pairs <- data.frame(
    phenotype_id = rownames(tqe)[1],
    variant_id   = mcols(tqtlVariantRanges(tqe))[["snp_id"]][1]
)
qtlRegressionStats(tqe, pairs)
#>      phenotype_id     variant_id t._Intercept_ t.genotype     t.sex  t.age
#> 1 ENSG00000100181 22:16849573A-G          1.81      -1.07 -0.000359 -0.785
#>   t.expr_pc1 t.expr_pc2 t.geno_pc1 t.geno_pc2 t.geno_pc3 t.geno_pc4 t.geno_pc5
#> 1      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#>   t.geno_pc6
#> 1      0.972
```

### Gene symbol and cis window

Alternatively, supply `symbol` and a `window` (bp) to fit every variant
within that window of the gene’s start. This requires a `gene_name`
column in `rowRanges(tqe)`, added by
[`addGeneSymbols()`](https://github.com/vjcitn/tQTLExperiment/reference/addGeneSymbols.md).

``` r

if (requireNamespace("EnsDb.Hsapiens.v79", quietly = TRUE)) {
    tqe <- addGeneSymbols(tqe)
    head(mcols(rowRanges(tqe))[, c("phenotype_id", "gene_name")])
}
#> DataFrame with 6 rows and 2 columns
#>                    phenotype_id   gene_name
#>                     <character> <character>
#> ENSG00000100181 ENSG00000100181      TPTEP1
#> ENSG00000237438 ENSG00000237438       CECR7
#> ENSG00000177663 ENSG00000177663      IL17RA
#> ENSG00000069998 ENSG00000069998       CECR5
#> ENSG00000093072 ENSG00000093072       CECR1
#> ENSG00000131100 ENSG00000131100    ATP6V1E1
```

``` r

res_gene <- qtlRegressionStats(tqe, symbol = "TPTEP1", window = 1000000L)
dim(res_gene)
#> [1] 1864   14
head(res_gene)
#>      phenotype_id     variant_id t._Intercept_ t.genotype     t.sex  t.age
#> 1 ENSG00000100181 22:16849573A-G          1.81      -1.07 -0.000359 -0.785
#> 2 ENSG00000100181 22:16849971A-T          1.81      -1.07 -0.000359 -0.785
#> 3 ENSG00000100181 22:16850437G-A          1.81      -1.07 -0.000359 -0.785
#> 4 ENSG00000100181 22:16851225C-T          1.81      -1.07 -0.000359 -0.785
#> 5 ENSG00000100181 22:16851356C-T          1.81      -1.07 -0.000359 -0.785
#> 6 ENSG00000100181 22:16851640C-T          1.81      -1.07 -0.000359 -0.785
#>   t.expr_pc1 t.expr_pc2 t.geno_pc1 t.geno_pc2 t.geno_pc3 t.geno_pc4 t.geno_pc5
#> 1      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 2      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 3      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 4      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 5      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 6      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#>   t.geno_pc6
#> 1      0.972
#> 2      0.972
#> 3      0.972
#> 4      0.972
#> 5      0.972
#> 6      0.972
```

`pairs` and `symbol` are mutually exclusive – supply exactly one.

## Output shape: t_only and wideFormat

Two arguments control the shape of the result:

- `t_only`: `TRUE` (default) returns just the t-statistic
  (`estimate / std.error`) per term; `FALSE` also returns the estimate,
  standard error, raw t-statistic, and p-value.
- `wideFormat`: `TRUE` (default) pivots to one row per pair, with one
  (or two, if `t_only = FALSE`) columns per model term; `FALSE` returns
  the long tidy form, one row per pair x term.

``` r

# long form, full statistics -- one row per pair x covariate
long <- qtlRegressionStats(tqe, symbol = "TPTEP1", t_only = FALSE, wideFormat = FALSE)
head(long)
#>      phenotype_id     variant_id        term  estimate std.error statistic
#> 1 ENSG00000100181 22:16849573A-G (Intercept)  1.79e-02   0.00985  1.812558
#> 2 ENSG00000100181 22:16849573A-G    genotype -2.07e-03   0.00194 -1.067511
#> 3 ENSG00000100181 22:16849573A-G         sex -1.03e-06   0.00285 -0.000359
#> 4 ENSG00000100181 22:16849573A-G         age -1.02e-04   0.00013 -0.785064
#> 5 ENSG00000100181 22:16849573A-G    expr_pc1 -3.05e-03   0.00218 -1.394788
#> 6 ENSG00000100181 22:16849573A-G    expr_pc2 -1.93e-03   0.00260 -0.742653
#>   p.value
#> 1  0.0733
#> 2  0.2887
#> 3  0.9997
#> 4  0.4345
#> 5  0.1666
#> 6  0.4597

# wide form, t-statistics only (the input qtlPCABrowser() expects)
head(qtlRegressionStats(tqe, symbol = "TPTEP1"))
#>      phenotype_id     variant_id t._Intercept_ t.genotype     t.sex  t.age
#> 1 ENSG00000100181 22:16849573A-G          1.81      -1.07 -0.000359 -0.785
#> 2 ENSG00000100181 22:16849971A-T          1.81      -1.07 -0.000359 -0.785
#> 3 ENSG00000100181 22:16850437G-A          1.81      -1.07 -0.000359 -0.785
#> 4 ENSG00000100181 22:16851225C-T          1.81      -1.07 -0.000359 -0.785
#> 5 ENSG00000100181 22:16851356C-T          1.81      -1.07 -0.000359 -0.785
#> 6 ENSG00000100181 22:16851640C-T          1.81      -1.07 -0.000359 -0.785
#>   t.expr_pc1 t.expr_pc2 t.geno_pc1 t.geno_pc2 t.geno_pc3 t.geno_pc4 t.geno_pc5
#> 1      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 2      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 3      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 4      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 5      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 6      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#>   t.geno_pc6
#> 1      0.972
#> 2      0.972
#> 3      0.972
#> 4      0.972
#> 5      0.972
#> 6      0.972
```

The wide, `t_only = TRUE` form is exactly what
[`qtlPCABrowser()`](https://github.com/vjcitn/tQTLExperiment/articles/pcaBrowser.md)
and
[`collapseFactorTstats()`](https://github.com/vjcitn/tQTLExperiment/reference/collapseFactorTstats.md)
expect as input: one row per SNP-gene pair, one `<term>.t` column per
covariate.

## Collapsing factor-expanded columns

When a covariate matrix includes a multi-level factor (e.g. `batch` or
`population`),
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) expands it
into several dummy columns, so the wide
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
output carries one `t.<var><level>` column per level.
[`collapseFactorTstats()`](https://github.com/vjcitn/tQTLExperiment/reference/collapseFactorTstats.md)
replaces each such group with a single root-mean-square summary column,
which keeps a PCA on t-statistics from being dominated by however many
levels a factor happens to have.

``` r

wide <- qtlRegressionStats(tqe, symbol = "TPTEP1")
names(wide)
#>  [1] "phenotype_id"  "variant_id"    "t._Intercept_" "t.genotype"   
#>  [5] "t.sex"         "t.age"         "t.expr_pc1"    "t.expr_pc2"   
#>  [9] "t.geno_pc1"    "t.geno_pc2"    "t.geno_pc3"    "t.geno_pc4"   
#> [13] "t.geno_pc5"    "t.geno_pc6"
# no factor-expanded columns in this toy covariate set, so this is a no-op
head(collapseFactorTstats(wide))
#>      phenotype_id     variant_id t._Intercept_ t.genotype     t.sex  t.age
#> 1 ENSG00000100181 22:16849573A-G          1.81      -1.07 -0.000359 -0.785
#> 2 ENSG00000100181 22:16849971A-T          1.81      -1.07 -0.000359 -0.785
#> 3 ENSG00000100181 22:16850437G-A          1.81      -1.07 -0.000359 -0.785
#> 4 ENSG00000100181 22:16851225C-T          1.81      -1.07 -0.000359 -0.785
#> 5 ENSG00000100181 22:16851356C-T          1.81      -1.07 -0.000359 -0.785
#> 6 ENSG00000100181 22:16851640C-T          1.81      -1.07 -0.000359 -0.785
#>   t.expr_pc1 t.expr_pc2 t.geno_pc1 t.geno_pc2 t.geno_pc3 t.geno_pc4 t.geno_pc5
#> 1      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 2      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 3      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 4      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 5      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#> 6      -1.39     -0.743       1.17       1.02       1.11     -0.612      0.215
#>   t.geno_pc6
#> 1      0.972
#> 2      0.972
#> 3      0.972
#> 4      0.972
#> 5      0.972
#> 6      0.972
```

On the MAGE chr17 dataset used in the main vignette, the covariate
matrix includes `factor_batch_*` and `population*` dummy columns, and
[`collapseFactorTstats()`](https://github.com/vjcitn/tQTLExperiment/reference/collapseFactorTstats.md)
reduces each group to `t.batch` and `t.population`. See
[`?collapseFactorTstats`](https://github.com/vjcitn/tQTLExperiment/reference/collapseFactorTstats.md)
for the pattern arguments.

## Parallelization

`BPPARAM` is passed straight to
[`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html).
Pairs are grouped by phenotype: with multiple phenotypes, groups are
distributed across workers; in gene-symbol mode (a single phenotype,
many variants), the inner loop over variants is parallelized instead.

``` r

qtlRegressionStats(
    tqe,
    symbol  = "TPTEP1",
    BPPARAM = BiocParallel::MulticoreParam(4)
)
```

## Chromosome-naming caveat

`rowRanges(tqe)` built from a `SummarizedExperiment` typically carries
UCSC-style seqnames (`"chr17"`), while `tqtlVariantRanges(tqe)` reflects
the PLINK `.bim` convention of bare chromosome numbers (`"17"`).
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
normalizes both by stripping any leading `"chr"` before comparing when
`symbol` is used to select the cis window. If you write code that
intersects `rowRanges(tqe)` and `tqtlVariantRanges(tqe)` directly, apply
the same normalization:

``` r

sub("^chr", "", as.character(seqnames(gr)))
```

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
#>   [1] DBI_1.3.0                 bitops_1.1-0             
#>   [3] httr2_1.3.0               rlang_1.3.0              
#>   [5] magrittr_2.0.5            otel_0.2.0               
#>   [7] compiler_4.6.1            RSQLite_3.53.3           
#>   [9] GenomicFeatures_1.65.0    png_0.1-9                
#>  [11] systemfonts_1.3.2         vctrs_0.7.3              
#>  [13] ProtGenerics_1.45.0       pkgconfig_2.0.3          
#>  [15] crayon_1.5.3              fastmap_1.2.0            
#>  [17] dbplyr_2.6.0              XVector_0.53.0           
#>  [19] Rsamtools_2.29.0          promises_1.5.0           
#>  [21] rmarkdown_2.32            UCSC.utils_1.9.0         
#>  [23] ggbeeswarm_0.7.3          ragg_1.5.2               
#>  [25] purrr_1.2.2               bit_4.6.0                
#>  [27] xfun_0.62                 cachem_1.1.0             
#>  [29] cigarillo_1.3.1           GenomeInfoDb_1.49.1      
#>  [31] jsonlite_2.0.0            blob_1.3.0               
#>  [33] later_1.4.8               DelayedArray_0.39.8      
#>  [35] BiocParallel_1.47.0       parallel_4.6.1           
#>  [37] R6_2.6.1                  bslib_0.12.0             
#>  [39] RColorBrewer_1.1-3        rtracklayer_1.73.0       
#>  [41] jquerylib_0.1.4           Rcpp_1.1.2               
#>  [43] bookdown_0.48             knitr_1.52               
#>  [45] BiocBaseUtils_1.15.1      httpuv_1.6.17            
#>  [47] Matrix_1.7-6              tidyselect_1.2.1         
#>  [49] dichromat_2.0-1           abind_1.4-8              
#>  [51] yaml_2.3.12               codetools_0.2-20         
#>  [53] curl_8.1.0                lattice_0.23-1           
#>  [55] tibble_3.3.1              shiny_1.14.0             
#>  [57] KEGGREST_1.53.6           S7_0.2.2                 
#>  [59] evaluate_1.0.5            BEDMatrix_2.0.4          
#>  [61] desc_1.4.3                BiocFileCache_3.3.1      
#>  [63] crochet_2.3.0             Biostrings_2.81.9        
#>  [65] filelock_1.0.3            pillar_1.11.1            
#>  [67] BiocManager_1.30.27       DT_0.34.0                
#>  [69] plotly_4.12.1             RCurl_1.98-1.20          
#>  [71] ensembldb_2.37.3          ggplot2_4.0.3            
#>  [73] scales_1.4.0              xtable_1.8-8             
#>  [75] glue_1.8.1                lazyeval_0.2.3           
#>  [77] tools_4.6.1               BiocIO_1.23.3            
#>  [79] data.table_1.18.6.1       GenomicAlignments_1.49.2 
#>  [81] fs_2.1.0                  XML_3.99-0.25            
#>  [83] grid_4.6.1                tidyr_1.3.2              
#>  [85] AnnotationDbi_1.75.2      EnsDb.Hsapiens.v79_2.99.0
#>  [87] beeswarm_0.4.0            restfulr_0.0.17          
#>  [89] vipor_0.4.7               cli_3.6.6                
#>  [91] textshaping_1.0.5         viridisLite_0.4.3        
#>  [93] S4Arrays_1.13.2           dplyr_1.2.1              
#>  [95] AnnotationFilter_1.37.0   gtable_0.3.6             
#>  [97] sass_0.4.10               digest_0.6.39            
#>  [99] SparseArray_1.13.4        rjson_0.2.23             
#> [101] htmlwidgets_1.6.4         farver_2.1.2             
#> [103] memoise_2.0.1             htmltools_0.5.9          
#> [105] pkgdown_2.2.1             lifecycle_1.0.5          
#> [107] httr_1.4.9                mime_0.13                
#> [109] bit64_4.8.6
```
