<div id="main" class="col-md-9" role="main">

# Subset a tQTLExperiment

<div class="ref-description section level2">

Sample subsetting propagates to the lazy `geno` matrix (BEDMatrix rows).

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S4 method for class 'tQTLExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A
    [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

-   i:

    Row (feature) index.

-   j:

    Column (sample) index.

-   ...:

    Ignored.

-   drop:

    Ignored.

</div>

<div class="section level2">

## Value

A
[tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

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
tqe[1:5, 1:10]
#> class: tQTLExperiment
#> features: 5  samples: 10 
#> assays( 1 ): pheno 
#> rowRanges: GRanges with 5 features
#> colData( 0 ) covariates:   
#> geno: 10 samples x 69638 variants [BEDMatrix - lazy]
#> plinkPrefix: /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpZyaoNi/temp_libpatha318583333e0/tQTLExperiment/extdata/chr22-n100 
#> use prepareTQTL() to write inputs and get CLI command
```

</div>

</div>

</div>
