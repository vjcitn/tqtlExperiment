# Access the bound tensorQTL cis_nominal results directory

Stores the path to a directory of tensorQTL `cis_nominal` parquet
output, without loading any of it into memory – the actual data stays on
disk and is queried lazily via
[`tqtlResults()`](https://github.com/vjcitn/tQTLExperiment/reference/tqtlResults.md).
Set automatically by
[`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)
when `mode = "cis_nominal"` and `x` is supplied, or set directly to bind
a directory from a run done outside R.

## Usage

``` r
tqtlResultsDir(x)

tqtlResultsDir(x) <- value

# S4 method for class 'tQTLExperiment'
tqtlResultsDir(x)

# S4 method for class 'tQTLExperiment'
tqtlResultsDir(x) <- value
```

## Arguments

- x:

  A
  [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

- value:

  A length-one character string giving the directory path.

## Value

The `tqtlResultsDir` slot (a path, or `NA_character_` if unset).

## See also

[`tqtlResults()`](https://github.com/vjcitn/tQTLExperiment/reference/tqtlResults.md)
to query the bound directory lazily.

## Examples

``` r
exdir <- system.file("extdata", package = "tQTLExperiment")
tqe <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
tqtlResultsDir(tqe)  # NA_character_ -- nothing bound yet
#> [1] NA
```
