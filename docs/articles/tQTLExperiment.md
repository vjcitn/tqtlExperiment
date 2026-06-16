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

## Example data (chr22)

The package ships with a small example dataset (20 genes, 100 samples,
chromosome 22) that can be used without external data or a tensorQTL
installation to verify the object construction steps.

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
#> plinkPrefix: /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpotkzK9/temp_libpath10c16232f7b5e/tQTLExperiment/extdata/chr22-n100 
#> use prepareTQTL() to write inputs and get CLI command
rowRanges(tqe_ex)
#> GRanges object with 20 ranges and 1 metadata column:
#>                   seqnames    ranges strand |    phenotype_id
#>                      <Rle> <IRanges>  <Rle> |     <character>
#>   ENSG00000100181       22  17082777      * | ENSG00000100181
#>   ENSG00000237438       22  17517460      * | ENSG00000237438
#>   ENSG00000177663       22  17565844      * | ENSG00000177663
#>   ENSG00000069998       22  17646177      * | ENSG00000069998
#>   ENSG00000093072       22  17702879      * | ENSG00000093072
#>               ...      ...       ...    ... .             ...
#>   ENSG00000260924       22  19158908      * | ENSG00000260924
#>   ENSG00000100075       22  19166343      * | ENSG00000100075
#>   ENSG00000185608       22  19419425      * | ENSG00000185608
#>   ENSG00000100084       22  19435224      * | ENSG00000100084
#>   ENSG00000185065       22  19435416      * | ENSG00000185065
#>   -------
#>   seqinfo: 1 sequence from hg38 genome; no seqlengths
tqtlGeno(tqe_ex)
#> BEDMatrix: 100 x 69638 [/private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpotkzK9/temp_libpath10c16232f7b5e/tQTLExperiment/extdata/chr22-n100.bed]
```

</div>

The tensorQTL command for this example would be obtained with:

<div id="cb16" class="sourceCode">

``` r
od  <- tempdir()
cmd <- prepareTQTL(tqe_ex, outDir = od, mode = "cis_nominal")
#> Run the following command in a terminal with tensorqtl available:
#> 
#> python3 -m tensorqtl /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpotkzK9/temp_libpath10c16232f7b5e/tQTLExperiment/extdata/chr22-n100 /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpwkAvcF/pheno.bed /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpwkAvcF/tqtl_out --mode cis_nominal --maf_threshold 0.05 --window 1000000 -o /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpwkAvcF
#> 
#> Then call readTQTL('/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpwkAvcF', mode = 'cis_nominal') to load results into R.
cat(cmd, "\n")
#> python3 -m tensorqtl /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpotkzK9/temp_libpath10c16232f7b5e/tQTLExperiment/extdata/chr22-n100 /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpwkAvcF/pheno.bed /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpwkAvcF/tqtl_out --mode cis_nominal --maf_threshold 0.05 --window 1000000 -o /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpwkAvcF
```

</div>

</div>

<div class="section level2">

## Session info

<div id="cb17" class="sourceCode">

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
#>  [1] tQTLExperiment_0.1.14       SummarizedExperiment_1.42.0
#>  [3] Biobase_2.72.0              GenomicRanges_1.64.0       
#>  [5] Seqinfo_1.2.0               IRanges_2.46.0             
#>  [7] S4Vectors_0.50.1            BiocGenerics_0.59.7        
#>  [9] generics_0.1.4              MatrixGenerics_1.24.0      
#> [11] matrixStats_1.5.0           BiocStyle_2.40.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] KEGGREST_1.52.0      httr2_1.2.2          xfun_0.58           
#>  [4] bslib_0.11.0         htmlwidgets_1.6.4    lattice_0.22-9      
#>  [7] vctrs_0.7.3          tools_4.6.0          curl_7.1.0          
#> [10] tibble_3.3.1         AnnotationDbi_1.74.0 RSQLite_3.53.1      
#> [13] blob_1.3.0           pkgconfig_2.0.3      Matrix_1.7-5        
#> [16] data.table_1.18.4    dbplyr_2.5.2         desc_1.4.3          
#> [19] BEDMatrix_2.0.4      lifecycle_1.0.5      compiler_4.6.0      
#> [22] textshaping_1.0.5    Biostrings_2.80.1    GenomeInfoDb_1.48.0 
#> [25] htmltools_0.5.9      sass_0.4.10          yaml_2.3.12         
#> [28] pkgdown_2.2.0        pillar_1.11.1        crayon_1.5.3        
#> [31] jquerylib_0.1.4      DelayedArray_0.38.2  cachem_1.1.0        
#> [34] abind_1.4-8          tidyselect_1.2.1     digest_0.6.39       
#> [37] dplyr_1.2.1          bookdown_0.46        fastmap_1.2.0       
#> [40] grid_4.6.0           cli_3.6.6            SparseArray_1.12.2  
#> [43] magrittr_2.0.5       S4Arrays_1.12.0      UCSC.utils_1.8.0    
#> [46] filelock_1.0.3       rappdirs_0.3.4       bit64_4.8.2         
#> [49] rmarkdown_2.31       XVector_0.52.0       httr_1.4.8          
#> [52] bit_4.6.0            otel_0.2.0           ragg_1.5.2          
#> [55] png_0.1-9            memoise_2.0.1        evaluate_1.0.5      
#> [58] knitr_1.51           crochet_2.3.0        BiocFileCache_3.2.0 
#> [61] rlang_1.2.0          glue_1.8.1           DBI_1.3.0           
#> [64] BiocManager_1.30.27  jsonlite_2.0.0       R6_2.6.1            
#> [67] systemfonts_1.3.2    fs_2.1.0
```

</div>

</div>

</div>
