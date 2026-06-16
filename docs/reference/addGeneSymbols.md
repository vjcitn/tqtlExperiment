<div id="main" class="col-md-9" role="main">

# Add gene symbols to rowRanges

<div class="ref-description section level2">

Looks up HGNC gene symbols for Ensembl gene IDs in `rownames(x)` using
an [ensembldb::EnsDb](https://rdrr.io/pkg/ensembldb/man/EnsDb.html)
annotation object and adds them as a `gene_name` column in
`mcols(rowRanges(x))`. Version suffixes (e.g. `.3` in
`ENSG00000100181.3`) are stripped before lookup.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
addGeneSymbols(x, ensdb = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A
    [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
    (or any
    [SummarizedExperiment::RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
    with Ensembl IDs as `rownames`).

-   ensdb:

    An [ensembldb::EnsDb](https://rdrr.io/pkg/ensembldb/man/EnsDb.html)
    object. Defaults to `EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79` if that
    package is installed.

</div>

<div class="section level2">

## Value

`x` with an additional `gene_name` metadata column in `rowRanges(x)`.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
exdir <- system.file("extdata", package = "tQTLExperiment")
tqe <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
    genome      = "hg38"
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
if (requireNamespace("EnsDb.Hsapiens.v79", quietly = TRUE))
    tqe <- addGeneSymbols(tqe)
```

</div>

</div>

</div>
