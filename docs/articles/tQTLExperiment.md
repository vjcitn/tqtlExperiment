<div id="main" class="col-md-9" role="main">

# eQTL mapping with tQTLExperiment and tensorQTL

<div class="section level2">

## Overview

`tQTLExperiment` wraps the
[tensorQTL](https://github.com/broadinstitute/tensorqtl) GPU-accelerated
eQTL mapping tool in a
[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment)-style
container. The key design principle is a **two-step workflow** that
keeps R and Python in separate environments:

1.  `prepareTQTL()` — writes phenotype and covariate files to a
    directory and returns the complete `python -m tensorqtl` command
    string.
2.  The user runs that command in a terminal with tensorQTL installed.
3.  `readTQTL()` — reads the output files back into R as
    [GRanges](https://bioconductor.org/packages/GenomicRanges).

This avoids all R/Python environment conflicts (conda PATH issues,
`rpy2` segfaults, duplicate OpenMP runtimes).

Genotype data are represented **lazily** via
[BEDMatrix](https://CRAN.R-project.org/package=BEDMatrix) — the PLINK
`.bed` file is never loaded into R memory. This makes it practical to
work with whole-genome genotype files containing millions of variants.

<div class="section level3">

### Installation

<div id="cb1" class="sourceCode">

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("vjcitn/tQTLExperiment")
```

</div>

tensorQTL must be installed in a Python environment:

<div id="cb2" class="sourceCode">

``` bash
pip install tensorqtl
```

</div>

</div>

</div>

<div class="section level2">

## MAGE example data

This vignette uses data from the
[MAGE](https://www.science.org/doi/10.1126/science.adi6529) study
(Multi-Ancestry Genomics, Epigenomics, and Transcriptomics), a
multi-ancestry eQTL resource. The `CSHLvc2026` companion package
provides cached access to chromosome 17 genotypes and a pre-filtered
pseudobulk B-cell SummarizedExperiment.

<div id="cb3" class="sourceCode">

``` r
library(tQTLExperiment)
library(CSHLvc2026)

data(mageSEfilt, package = "CSHLvc2026")
mageSEfilt
```

</div>

    #> Warning: multiple methods tables found for 'transform'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'IRanges'
    #> Warning: multiple methods tables found for 'transform'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'Seqinfo'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'GenomicRanges'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'SummarizedExperiment'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'S4Arrays'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'DelayedArray'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'SparseArray'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'XVector'
    #> Warning: multiple methods tables found for 'scale'
    #> Warning: replacing previous import 'BiocGenerics::scale' by
    #> 'DelayedArray::scale' when loading 'SummarizedExperiment'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'AnnotationDbi'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'Biostrings'
    #> Warning: replacing previous import 'utils::data' by 'BiocGenerics::data' when
    #> loading 'GenomeInfoDb'
    #> Warning: replacing previous import 'BiocGenerics::transform' by
    #> 'S4Vectors::transform' when loading 'GenomeInfoDb'

The genotype data live in PLINK format on OSN cloud storage and are
downloaded once via `BiocFileCache`:

<div id="cb5" class="sourceCode">

``` r
plink_paths <- cache_mage_chr17_plink()
plpre <- tools::file_path_sans_ext(plink_paths[["bed"]])
```

</div>

</div>

<div class="section level2">

## Building a tQTLExperiment

<div class="section level3">

### Covariate matrix

tensorQTL requires numeric covariates. We build a model matrix from
`colData` using `model.matrix()`, then drop the intercept column —
tensorQTL does not expect one. Dropping the first column after
construction (rather than using `~ 0 + ...`) preserves reference-level
dummy coding for factor variables.

<div id="cb6" class="sourceCode">

``` r
cd <- as.data.frame(colData(mageSEfilt))
mm <- model.matrix(~ batch + population + sex, data = cd)
mm <- mm[, -1, drop = FALSE]   # remove (Intercept)
```

</div>

</div>

<div class="section level3">

### Constructor

`tQTLExperimentFromRSE()` takes an existing `RangedSummarizedExperiment`
and PLINK prefix, attaches the lazy genotype matrix, and optionally
records the genome build in `seqinfo`.

<div id="cb7" class="sourceCode">

``` r
tqe <- tQTLExperimentFromRSE(
    se              = mageSEfilt,
    plinkPrefix     = plpre,
    covariateMatrix = mm,
    genome          = "hg38"
)
tqe
```

</div>

</div>

<div class="section level3">

### Adding gene symbols

`addGeneSymbols()` looks up HGNC names via `EnsDb.Hsapiens.v79` and adds
a `gene_name` column to `mcols(rowRanges(tqe))`. Ensembl version
suffixes (e.g. `.3`) are stripped automatically.

<div id="cb8" class="sourceCode">

``` r
tqe <- addGeneSymbols(tqe)
head(mcols(rowRanges(tqe))[, c("phenotype_id", "gene_name")])
```

</div>

</div>

</div>

<div class="section level2">

## Running tensorQTL

<div class="section level3">

### Subset to chr17 features

The PLINK file covers chromosome 17 only, so we filter the SE to
matching features before writing the phenotype file.

<div id="cb9" class="sourceCode">

``` r
is_17 <- as.character(seqnames(rowRanges(tqe))) == "chr17"
tqe17 <- tqe[which(is_17), ]
dim(tqe17)
```

</div>

</div>

<div class="section level3">

### prepareTQTL

`prepareTQTL()` writes `pheno.bed` and `covariates.tsv` to `outDir` and
returns the shell command. The built-in `--maf_threshold` flag is one of
tensorQTL’s key advantages — no pre-filtering step is required.

<div id="cb10" class="sourceCode">

``` r
outdir <- "~/mage_tqtl_run"
dir.create(outdir, showWarnings = FALSE)

cmd <- prepareTQTL(
    tqe17,
    outDir       = outdir,
    mode         = "cis_nominal",
    mafThreshold = 0.05
)
cat(cmd, "\n")
```

</div>

The printed command looks like:

    python3 -m tensorqtl /path/to/CCDG_mage_chr17 ~/mage_tqtl_run/pheno.bed \
        ~/mage_tqtl_run/tqtl_out --mode cis_nominal --maf_threshold 0.05 \
        --window 1000000 -o ~/mage_tqtl_run \
        --covariates ~/mage_tqtl_run/covariates.tsv

Paste this into a terminal where your tensorQTL conda environment is
active:

<div id="cb12" class="sourceCode">

``` bash
conda activate plink-env
KMP_DUPLICATE_LIB_OK=TRUE <paste command>
```

</div>

</div>

<div class="section level3">

### readTQTL

Once tensorQTL finishes, load results back into R. Passing `x = tqe17`
propagates the genome build and gene symbols automatically.

<div id="cb13" class="sourceCode">

``` r
res <- readTQTL(outdir, mode = "cis_nominal", x = tqe17)
res$pairs
```

</div>

</div>

</div>

<div class="section level2">

## Exploring results

<div id="cb14" class="sourceCode">

``` r
pairs <- res$pairs

# top associations by nominal p-value
pairs[order(pairs$pval_nominal)][1:10,
      c("gene_name", "variant_id", "af", "slope", "pval_nominal")]

# associations for a specific gene
tp53bp1 <- pairs[pairs$gene_name == "TP53BP1" & !is.na(pairs$gene_name)]
tp53bp1[order(tp53bp1$pval_nominal)][1:5]

# distribution of nominal p-values
hist(pairs$pval_nominal, breaks = 50,
     main = "tensorQTL nominal p-values (chr17, cis)",
     xlab = "p-value")
```

</div>

</div>

<div class="section level2">

## Demo data

The package includes two datasets for testing and exploration.

<div class="section level3">

### Small example (chr22, 20 genes)

For testing object construction without external dependencies:

<div id="cb15" class="sourceCode">

``` r
exdir <- system.file("extdata", package = "tQTLExperiment")

tqe_ex <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
    genome      = "hg38"
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
tqe_ex
#> class: tQTLExperiment
#> features: 20  samples: 100 
#> assays( 1 ): pheno 
#> rowRanges: GRanges with 20 features
#> colData( 0 ) covariates:   
#> geno: 100 samples x 69638 variants [BEDMatrix - lazy]
#> plinkPrefix: /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpuS4nx1/temp_libpath5a4b140770e5/tQTLExperiment/extdata/chr22-n100 
#> use prepareTQTL() to write inputs and get CLI command
```

</div>

</div>

<div class="section level3">

### Realistic demo (chr17, 50 genes, pre-computed results)

Pre-computed tensorQTL `cis_nominal` results that demonstrate the full
workflow:

<div id="cb16" class="sourceCode">

``` r
demodir <- system.file("demodir", package = "tQTLExperiment")

# Load results: 71k+ feature-variant pairs on chr17
res <- readTQTL(demodir, mode = "cis_nominal")
res$pairs
#> GRanges object with 393196 ranges and 10 metadata columns:
#>            seqnames    ranges strand |      phenotype_id     variant_id
#>               <Rle> <IRanges>  <Rle> |       <character>    <character>
#>        [1]       17    114101      * | ENSG00000280279.1  17:114101:G:A
#>        [2]       17    114226      * | ENSG00000280279.1  17:114226:A:G
#>        [3]       17    116159      * | ENSG00000280279.1  17:116159:G:C
#>        [4]       17    116270      * | ENSG00000280279.1  17:116270:G:C
#>        [5]       17    116354      * | ENSG00000280279.1  17:116354:C:T
#>        ...      ...       ...    ... .               ...            ...
#>   [393192]       17   5221648      * | ENSG00000263219.1 17:5221648:C:T
#>   [393193]       17   5221790      * | ENSG00000263219.1 17:5221790:C:A
#>   [393194]       17   5221954      * | ENSG00000263219.1 17:5221954:A:G
#>   [393195]       17   5221997      * | ENSG00000263219.1 17:5221997:C:T
#>   [393196]       17   5222038      * | ENSG00000263219.1 17:5222038:A:T
#>            start_distance end_distance        af ma_samples  ma_count
#>                 <integer>    <integer> <numeric>  <integer> <integer>
#>        [1]          37236        37235  0.450068        510       658
#>        [2]          37361        37360  0.197674        247       289
#>        [3]          39294        39293  0.450752        515       659
#>        [4]          39405        39404  0.264706        326       387
#>        [5]          39489        39488  0.486320        515       711
#>        ...            ...          ...       ...        ...       ...
#>   [393192]         999558       999557 0.0595075         81        87
#>   [393193]         999700       999699 0.3522572        396       515
#>   [393194]         999864       999863 0.4418605        484       646
#>   [393195]         999907       999906 0.0827633        112       121
#>   [393196]         999948       999947 0.4418605        484       646
#>            pval_nominal      slope  slope_se
#>               <numeric>  <numeric> <numeric>
#>        [1]     0.213275  0.0618208 0.0496256
#>        [2]     0.491346  0.0431121 0.0626146
#>        [3]     0.263432  0.0559091 0.0499538
#>        [4]     0.164120  0.0794705 0.0570578
#>        [5]     0.291421  0.0490658 0.0464724
#>        ...          ...        ...       ...
#>   [393192]     0.883109 -0.0172224 0.1170930
#>   [393193]     0.248010  0.0752312 0.0650699
#>   [393194]     0.686174  0.0240168 0.0594152
#>   [393195]     0.762729  0.0299093 0.0990322
#>   [393196]     0.686174  0.0240168 0.0594152
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

</div>

Top associations by nominal p-value:

<div id="cb17" class="sourceCode">

``` r
top_hits <- res$pairs[order(res$pairs$pval_nominal)][1:10]
top_hits[, c("phenotype_id", "variant_id", "pval_nominal", "slope")]
#> GRanges object with 10 ranges and 4 metadata columns:
#>        seqnames    ranges strand |      phenotype_id    variant_id pval_nominal
#>           <Rle> <IRanges>  <Rle> |       <character>   <character>    <numeric>
#>    [1]       17    410391      * | ENSG00000187624.9 17:410391:T:C 3.11053e-163
#>    [2]       17    410351      * | ENSG00000187624.9 17:410351:G:T 6.25116e-163
#>    [3]       17    409513      * | ENSG00000187624.9 17:409513:A:G 3.30974e-162
#>    [4]       17    409857      * | ENSG00000187624.9 17:409857:G:A 1.00181e-159
#>    [5]       17    410971      * | ENSG00000187624.9 17:410971:G:A 8.13960e-113
#>    [6]       17    410956      * | ENSG00000187624.9 17:410956:G:A 1.56592e-110
#>    [7]       17    399065      * | ENSG00000187624.9 17:399065:C:T 4.89702e-107
#>    [8]       17    412764      * | ENSG00000187624.9 17:412764:A:G  5.00737e-73
#>    [9]       17    416189      * | ENSG00000187624.9 17:416189:G:C  5.05483e-72
#>   [10]       17    409016      * | ENSG00000187624.9 17:409016:C:A  5.70008e-72
#>            slope
#>        <numeric>
#>    [1]   1.92088
#>    [2]   1.92108
#>    [3]   1.92089
#>    [4]   1.91371
#>    [5]   1.78995
#>    [6]   1.76924
#>    [7]  -1.71253
#>    [8]   1.48159
#>    [9]   1.48194
#>   [10]   1.47057
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

</div>

<div id="cb18" class="sourceCode">

``` r
data(mageSEfilt, package="CSHLvc2026")
plink_paths <- cache_mage_chr17_plink()
#> [ cache hit   ] CCDG_mage_chr17.fam
#> [ cache hit   ] CCDG_mage_chr17.bim
#> [ cache hit   ] CCDG_mage_chr17.bed
#> [ validate    ] Checking .bed magic bytes ...
#> [ validate    ] .bed magic bytes OK.
plpre <- tools::file_path_sans_ext(plink_paths[["bed"]])
cd <- as.data.frame(colData(mageSEfilt))
mm <- model.matrix(~ batch + population + sex, data = cd)
mm <- mm[, -1, drop = FALSE]   # remove (Intercept)
tqe <- tQTLExperimentFromRSE(
    se              = mageSEfilt,
    plinkPrefix     = plpre,
    covariateMatrix = mm,
    genome          = "hg38"
)
#> Extracting number of samples and rownames from CCDG_mage_chr17.fam...
#> Extracting number of variants and colnames from CCDG_mage_chr17.bim...
tqe
#> class: tQTLExperiment
#> features: 15116  samples: 731 
#> assays( 1 ): pseudocounts 
#> rowRanges: GRanges with 15116 features
#> colData( 27 ) covariates: batch, populationASW, populationBEB, populationCDX ... 
#> geno: 731 samples x 2075523 variants [BEDMatrix - lazy]
#> plinkPrefix: /Users/vincentcarey/Library/Caches/org.R-project.R/R/BiocFileCache/CCDG_mage_chr17_plink/CCDG_mage_chr17 
#> use prepareTQTL() to write inputs and get CLI command
plotGenotypeEffect(tqe, "17:410391:T:C", "ENSG00000187624.9")
```

</div>

![](tQTLExperiment_files/figure-html/shviz-1.png)

</div>

</div>

<div class="section level2">

## Session info

<div id="cb19" class="sourceCode">

``` r
sessionInfo()
#> R version 4.6.0 beta (2026-04-12 r89882)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Sequoia 15.7.7
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
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
#>  [1] tQTLExperiment_0.1.16       SummarizedExperiment_1.42.0
#>  [3] Biobase_2.72.0              GenomicRanges_1.64.0       
#>  [5] Seqinfo_1.2.0               IRanges_2.46.0             
#>  [7] S4Vectors_0.50.1            BiocGenerics_0.59.7        
#>  [9] generics_0.1.4              MatrixGenerics_1.24.0      
#> [11] matrixStats_1.5.0           BiocStyle_2.40.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1     vipor_0.4.7          dplyr_1.2.1         
#>  [4] farver_2.1.2         blob_1.3.0           filelock_1.0.3      
#>  [7] arrow_24.0.0         Biostrings_2.80.1    S7_0.2.2            
#> [10] fastmap_1.2.0        BiocFileCache_3.2.0  digest_0.6.39       
#> [13] BEDMatrix_2.0.4      lifecycle_1.0.5      KEGGREST_1.52.0     
#> [16] RSQLite_3.53.1       magrittr_2.0.5       compiler_4.6.0      
#> [19] rlang_1.2.0          sass_0.4.10          tools_4.6.0         
#> [22] yaml_2.3.12          data.table_1.18.4    knitr_1.51          
#> [25] labeling_0.4.3       S4Arrays_1.12.0      htmlwidgets_1.6.4   
#> [28] bit_4.6.0            curl_7.1.0           DelayedArray_0.38.2 
#> [31] RColorBrewer_1.1-3   abind_1.4-8          withr_3.0.2         
#> [34] purrr_1.2.2          desc_1.4.3           grid_4.6.0          
#> [37] ggplot2_4.0.3        scales_1.4.0         dichromat_2.0-0.1   
#> [40] cli_3.6.6            rmarkdown_2.31       crayon_1.5.3        
#> [43] ragg_1.5.2           otel_0.2.0           httr_1.4.8          
#> [46] ggbeeswarm_0.7.3     DBI_1.3.0            cachem_1.1.0        
#> [49] assertthat_0.2.1     AnnotationDbi_1.74.0 BiocManager_1.30.27 
#> [52] XVector_0.52.0       vctrs_0.7.3          Matrix_1.7-5        
#> [55] jsonlite_2.0.0       bookdown_0.46        bit64_4.8.2         
#> [58] beeswarm_0.4.0       systemfonts_1.3.2    jquerylib_0.1.4     
#> [61] crochet_2.3.0        glue_1.8.1           pkgdown_2.2.0       
#> [64] gtable_0.3.6         GenomeInfoDb_1.48.0  UCSC.utils_1.8.0    
#> [67] tibble_3.3.1         pillar_1.11.1        rappdirs_0.3.4      
#> [70] htmltools_0.5.9      R6_2.6.1             dbplyr_2.5.2        
#> [73] httr2_1.2.2          textshaping_1.0.5    evaluate_1.0.5      
#> [76] lattice_0.22-9       png_0.1-9            memoise_2.0.1       
#> [79] bslib_0.11.0         SparseArray_1.12.2   xfun_0.58           
#> [82] fs_2.1.0             pkgconfig_2.0.3
```

</div>

</div>

</div>
