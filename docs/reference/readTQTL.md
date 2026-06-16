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
exdir <- system.file("extdata", package = "tQTLExperiment")
tqe <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
od <- tempdir()
cmd <- prepareTQTL(tqe, outDir = od, mode = "cis_nominal")
#> Run the following command in a terminal with tensorqtl available:
#> 
#> python3 -m tensorqtl /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpotkzK9/temp_libpath10c16232f7b5e/tQTLExperiment/extdata/chr22-n100 /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpkmYph4/pheno.bed /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpkmYph4/tqtl_out --mode cis_nominal --maf_threshold 0.05 --window 1000000 -o /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpkmYph4
#> 
#> Then call readTQTL('/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpkmYph4', mode = 'cis_nominal') to load results into R.
# (user runs cmd in terminal)
# res <- readTQTL(od, mode = "cis_nominal", x = tqe)
```

</div>

</div>

</div>
