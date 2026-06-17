<div id="main" class="col-md-9" role="main">

# Read tensorQTL results into R

<div class="ref-description section level2">

Reads the output files written by a tensorQTL run into
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
objects, matched to the features in the originating
[tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
readTQTL(
  outDir,
  mode = c("cis_nominal", "cis", "cis_independent", "trans"),
  x = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   outDir:

    The same directory passed to `prepareTQTL()`.

-   mode:

    The mode used for the tensorQTL run.

-   x:

    The
    [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
    used to generate the input files, used to attach feature coordinates
    to region-level results.

</div>

<div class="section level2">

## Value

A named list of
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html):

-   `pairs`:

    (`cis_nominal`) Feature-variant pairs with association statistics.

-   `hits`:

    (`cis`) One range per feature with permutation p-value.

-   `output_files`:

    (other modes) Paths to the raw output files.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`prepareTQTL()` to write input files and obtain the CLI command.

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# Load pre-computed tensorQTL results from demo data
demodir <- system.file("demodir", package = "tQTLExperiment")
res <- readTQTL(demodir, mode = "cis_nominal")
head(res$pairs)
#> GRanges object with 6 ranges and 10 metadata columns:
#>       seqnames    ranges strand |      phenotype_id    variant_id
#>          <Rle> <IRanges>  <Rle> |       <character>   <character>
#>   [1]       17    114101      * | ENSG00000280279.1 17:114101:G:A
#>   [2]       17    114226      * | ENSG00000280279.1 17:114226:A:G
#>   [3]       17    116159      * | ENSG00000280279.1 17:116159:G:C
#>   [4]       17    116270      * | ENSG00000280279.1 17:116270:G:C
#>   [5]       17    116354      * | ENSG00000280279.1 17:116354:C:T
#>   [6]       17    118358      * | ENSG00000280279.1 17:118358:C:T
#>       start_distance end_distance        af ma_samples  ma_count pval_nominal
#>            <integer>    <integer> <numeric>  <integer> <integer>    <numeric>
#>   [1]          37236        37235 0.4500684        510       658     0.213275
#>   [2]          37361        37360 0.1976744        247       289     0.491346
#>   [3]          39294        39293 0.4507524        515       659     0.263432
#>   [4]          39405        39404 0.2647059        326       387     0.164120
#>   [5]          39489        39488 0.4863201        515       711     0.291421
#>   [6]          41493        41492 0.0978112        136       143     0.656597
#>           slope  slope_se
#>       <numeric> <numeric>
#>   [1] 0.0618208 0.0496256
#>   [2] 0.0431121 0.0626146
#>   [3] 0.0559091 0.0499538
#>   [4] 0.0794705 0.0570578
#>   [5] 0.0490658 0.0464724
#>   [6] 0.0402160 0.0904124
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

</div>

</div>

</div>
