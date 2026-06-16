<div id="main" class="col-md-9" role="main">

# Access variant genomic ranges

<div class="ref-description section level2">

Access variant genomic ranges

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
tqtlVariantRanges(x)

tqtlVariantRanges(x) <- value

# S4 method for class 'tQTLExperiment'
tqtlVariantRanges(x)

# S4 method for class 'tQTLExperiment'
tqtlVariantRanges(x) <- value
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
    [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html).

</div>

<div class="section level2">

## Value

A `GRanges` with metadata columns `snp_id`, `ref`, `alt`.

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
tqtlVariantRanges(tqe)
#> GRanges object with 69638 ranges and 3 metadata columns:
#>                  seqnames    ranges strand |         snp_id         ref
#>                     <Rle> <IRanges>  <Rle> |    <character> <character>
#>   22:16849573A-G       22  16849573      * | 22:16849573A-G           A
#>   22:16849971A-T       22  16849971      * | 22:16849971A-T           A
#>   22:16850437G-A       22  16850437      * | 22:16850437G-A           G
#>   22:16851225C-T       22  16851225      * | 22:16851225C-T           C
#>   22:16851356C-T       22  16851356      * | 22:16851356C-T           C
#>              ...      ...       ...    ... .            ...         ...
#>   22:51202748A-G       22  51202748      * | 22:51202748A-G           A
#>   22:51208568G-T       22  51208568      * | 22:51208568G-T           G
#>   22:51211031A-G       22  51211031      * | 22:51211031A-G           A
#>   22:51213613C-T       22  51213613      * | 22:51213613C-T           C
#>   22:51216564T-C       22  51216564      * | 22:51216564T-C           T
#>                          alt
#>                  <character>
#>   22:16849573A-G           G
#>   22:16849971A-T           T
#>   22:16850437G-A           A
#>   22:16851225C-T           T
#>   22:16851356C-T           T
#>              ...         ...
#>   22:51202748A-G           G
#>   22:51208568G-T           T
#>   22:51211031A-G           G
#>   22:51213613C-T           T
#>   22:51216564T-C           C
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

</div>

</div>

</div>
