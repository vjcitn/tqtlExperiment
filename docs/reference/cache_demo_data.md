# Cache demo dataset files

Downloads the complete tQTLExperiment demo dataset (chr17 tensorQTL
results and associated data) from the GitHub release into a
`BiocFileCache` subdirectory. Files are cached locally and only
re-downloaded when absent.

## Usage

``` r
cache_demo_data(
  release_url = paste0("https://github.com/vjcitn/tQTLExperiment/",
    "releases/download/v0.1.0-demodir"),
  bfc = BiocFileCache::BiocFileCache(),
  verbose = TRUE
)
```

## Arguments

- release_url:

  Base URL of the GitHub release assets. Defaults to the
  `v0.1.0-demodir` release of `vjcitn/tQTLExperiment`.

- bfc:

  A
  [BiocFileCache::BiocFileCache](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
  object. Defaults to the user-level cache.

- verbose:

  Logical. Emit progress messages? Default `TRUE`.

## Value

A named list with elements:

- `covariates_tsv`: Path to covariates.tsv

- `mage_rda`: Path to mage_chr17_cis.rda (pre-computed results)

- `pheno_bed`: Path to pheno.bed (phenotype data)

- `parquet`: Path to tqtl_out.cis_qtl_pairs.17.parquet (tensorQTL
  output)

- `log`: Path to tqtl_out.tensorQTL.cis_nominal.log

All files share the same directory, accessible via
`dirname(result[["covariates_tsv"]])`.

## Examples

``` r
if (FALSE) { # \dontrun{
demo_data <- cache_demo_data()
demodir <- dirname(demo_data[["covariates_tsv"]])
} # }
```
