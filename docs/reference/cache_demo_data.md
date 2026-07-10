<div id="main" class="col-md-9" role="main">

# Cache demo dataset files

<div class="ref-description section level2">

Downloads the complete tQTLExperiment demo dataset (chr17 tensorQTL
results and associated data) from the GitHub release into a
`BiocFileCache` subdirectory. Files are cached locally and only
re-downloaded when absent.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
cache_demo_data(
  release_url = paste0("https://github.com/vjcitn/tQTLExperiment/",
    "releases/download/v0.1.0-demodir"),
  bfc = BiocFileCache::BiocFileCache(),
  verbose = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   release\_url:

    Base URL of the GitHub release assets. Defaults to the
    `v0.1.0-demodir` release of `vjcitn/tQTLExperiment`.

-   bfc:

    A
    [BiocFileCache::BiocFileCache](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
    object. Defaults to the user-level cache.

-   verbose:

    Logical. Emit progress messages? Default `TRUE`.

</div>

<div class="section level2">

## Value

A named list with elements:

-   `covariates_tsv`: Path to covariates.tsv

-   `mage_rda`: Path to mage\_chr17\_cis.rda (pre-computed results)

-   `pheno_bed`: Path to pheno.bed (phenotype data)

-   `parquet`: Path to tqtl\_out.cis\_qtl\_pairs.17.parquet (tensorQTL
    output)

-   `log`: Path to tqtl\_out.tensorQTL.cis\_nominal.log

All files share the same directory, accessible via
`dirname(result[["covariates_tsv"]])`.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
demo_data <- cache_demo_data()
demodir <- dirname(demo_data[["covariates_tsv"]])
} # }
```

</div>

</div>

</div>
