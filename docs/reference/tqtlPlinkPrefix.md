<div id="main" class="col-md-9" role="main">

# Access the PLINK file prefix

<div class="ref-description section level2">

Access the PLINK file prefix

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
tqtlPlinkPrefix(x)

tqtlPlinkPrefix(x) <- value

# S4 method for class 'tQTLExperiment'
tqtlPlinkPrefix(x)

# S4 method for class 'tQTLExperiment'
tqtlPlinkPrefix(x) <- value
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A
    [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

-   value:

    A length-one character string.

</div>

<div class="section level2">

## Value

The PLINK prefix path.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
exdir <- system.file("extdata", package = "tQTLExperiment")
tqe <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
tqtlPlinkPrefix(tqe)
#> [1] "/private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpotkzK9/temp_libpath10c16232f7b5e/tQTLExperiment/extdata/chr22-n100"
```

</div>

</div>

</div>
