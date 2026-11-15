# Cache and load the MAGE chr17 SummarizedExperiment

Downloads `mageSEfilt.rda` from the Open Storage Network bucket into a
`BiocFileCache` subdirectory and returns the contained
[SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object. The file is only downloaded once; subsequent calls return the
cached copy.

## Usage

``` r
get_cached_mage17SE(
  url = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocMAGE17SE/mageSEfilt.rda",
  bfc = BiocFileCache::BiocFileCache(),
  verbose = TRUE
)
```

## Arguments

- url:

  Full URL to `mageSEfilt.rda`. Defaults to the OSN BiocMAGE17SE bucket.

- bfc:

  A
  [BiocFileCache::BiocFileCache](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
  object. Defaults to the user-level cache.

- verbose:

  Logical. Emit progress messages? Default `TRUE`.

## Value

A
[SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object (`mageSEfilt`).

## Examples

``` r
if (FALSE) { # \dontrun{
se <- get_cached_mage17SE()
se
} # }
```
