<div id="main" class="col-md-9" role="main">

# Access the lazy genotype matrix

<div class="ref-description section level2">

Access the lazy genotype matrix

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
tqtlGeno(x)

tqtlGeno(x) <- value

# S4 method for class 'tQTLExperiment'
tqtlGeno(x)

# S4 method for class 'tQTLExperiment'
tqtlGeno(x) <- value
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A
    [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

-   value:

    A
    [BEDMatrix::BEDMatrix](https://rdrr.io/pkg/BEDMatrix/man/BEDMatrix.html)
    (samples × variants).

</div>

<div class="section level2">

## Value

The `geno` slot (samples × variants BEDMatrix).

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
tqtlGeno(tqe)
#> BEDMatrix: 100 x 69638 [/private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpuS4nx1/temp_libpath5a4b140770e5/tQTLExperiment/extdata/chr22-n100.bed]
```

</div>

</div>

</div>
