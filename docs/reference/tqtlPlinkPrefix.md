# Access the PLINK file prefix

Access the PLINK file prefix

## Usage

``` r
tqtlPlinkPrefix(x)

tqtlPlinkPrefix(x) <- value

# S4 method for class 'tQTLExperiment'
tqtlPlinkPrefix(x)

# S4 method for class 'tQTLExperiment'
tqtlPlinkPrefix(x) <- value
```

## Arguments

- x:

  A
  [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

- value:

  A length-one character string.

## Value

The PLINK prefix path.

## Examples

``` r
exdir <- system.file("extdata", package = "tQTLExperiment")
tqe <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
tqtlPlinkPrefix(tqe)
#> [1] "/private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpPoihVR/temp_libpath173973d2e8831/tQTLExperiment/extdata/chr22-n100"
```
