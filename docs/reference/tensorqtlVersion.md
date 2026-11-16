# Report the tensorQTL version from a Python installation

Report the tensorQTL version from a Python installation

## Usage

``` r
tensorqtlVersion(python = findTQTL())
```

## Arguments

- python:

  Path to the Python executable. Defaults to
  [`findTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/findTQTL.md).

## Value

A character string with the version, or `NA`.

## Examples

``` r
tryCatch(tensorqtlVersion(), error = function(e) NA_character_)
#> [1] "Warning: 'rfunc' cannot be imported. R with the 'qvalue' library and the 'rpy2' Python package are needed to compute q-values."
```
