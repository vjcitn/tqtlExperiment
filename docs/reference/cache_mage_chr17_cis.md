<div id="main" class="col-md-9" role="main">

# Cache pre-computed MAGE chr17 tensorQTL results

<div class="ref-description section level2">

Downloads the pre-computed `cis_nominal` tensorQTL results for MAGE
chromosome 17 from the tQTLExperiment GitHub release into a
`BiocFileCache` subdirectory. Results are cached as a serialized R
object and loaded into memory as a
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
cache_mage_chr17_cis(
  release_url = paste0("https://github.com/vjcitn/tQTLExperiment/",
    "releases/download/v0.1.0-demo-data"),
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
    `v0.1.0-demo-data` release of `vjcitn/tQTLExperiment`.

-   bfc:

    A
    [BiocFileCache::BiocFileCache](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
    object. Defaults to the user-level cache.

-   verbose:

    Logical. Emit progress messages? Default `TRUE`.

</div>

<div class="section level2">

## Value

A
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object containing 71,000+ feature-variant pairs from tensorQTL
`cis_nominal` results. Metadata columns include `phenotype_id`,
`variant_id`, `tss_distance`, `ma_samples`, `ma_count`, `maf`,
`pval_nominal`, and `slope`.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
gr <- cache_mage_chr17_cis()
head(gr)
} # }
```

</div>

</div>

</div>
