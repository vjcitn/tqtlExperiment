<div id="main" class="col-md-9" role="main">

# Cache MAGE chr17 PLINK files from GitHub releases

<div class="ref-description section level2">

Downloads the three PLINK files (`CCDG_mage_chr17.fam`, `.bim`, `.bed`)
from the `vjcitn/CSHLvc2026` GitHub release into a `BiocFileCache`
subdirectory. Files are only re-downloaded when absent from the cache.
No authentication is required.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
cache_mage_chr17_plink(
  release_url = paste0("https://github.com/vjcitn/CSHLvc2026/",
    "releases/download/v0.1.0-data"),
  bfc = BiocFileCache::BiocFileCache(),
  verbose = TRUE,
  validate_bed = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   release\_url:

    Base URL of the GitHub release assets. Defaults to the `v0.1.0-data`
    release of `vjcitn/CSHLvc2026`.

-   bfc:

    A
    [BiocFileCache::BiocFileCache](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
    object. Defaults to the user-level cache.

-   verbose:

    Logical. Emit progress messages? Default `TRUE`.

-   validate\_bed:

    Logical. Check `.bed` magic bytes after download? Default `TRUE`.

</div>

<div class="section level2">

## Value

A named list with elements `fam`, `bim`, and `bed`, each a character
scalar giving the local path to the cached file. All three share the
same directory and stem, so `tools::file_path_sans_ext(result[["bed"]])`
gives the PLINK prefix for use with `tQTLExperimentFromRSE()`.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
plink <- cache_mage_chr17_plink()
plpre <- tools::file_path_sans_ext(plink[["bed"]])
} # }
```

</div>

</div>

</div>
