# Subset a tQTLExperiment

Sample subsetting propagates to the lazy `geno` matrix (BEDMatrix rows).

## Usage

``` r
# S4 method for class 'tQTLExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A
  [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

- i:

  Row (feature) index.

- j:

  Column (sample) index.

- ...:

  Ignored.

- drop:

  Ignored.

## Value

A
[tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

## Examples

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
#> plinkPrefix: /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpPoihVR/temp_libpath173973d2e8831/tQTLExperiment/extdata/chr22-n100 
#> use prepareTQTL() to write inputs and get CLI command
```
