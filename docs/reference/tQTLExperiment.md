<div id="main" class="col-md-9" role="main">

# Construct a tQTLExperiment from input files

<div class="ref-description section level2">

Reads phenotype, covariate, and variant metadata files into a
tQTLExperiment object. Genotype data is represented lazily via
[BEDMatrix::BEDMatrix](https://rdrr.io/pkg/BEDMatrix/man/BEDMatrix.html).

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
tQTLExperiment(
  plinkPrefix,
  phenoFile,
  covFile = NULL,
  assayName = "pheno",
  genome = NA_character_
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   plinkPrefix:

    Path prefix for the PLINK file set (`.bed`, `.bim`, `.fam`).

-   phenoFile:

    Path to the phenotype BED file (tab-separated, columns `#chr`,
    `start`, `end`, `phenotype_id`, then one column per sample; 0-based
    half-open coordinates).

-   covFile:

    Path to the covariate file. tensorQTL expects covariates as rows
    with the first column being the covariate name and remaining columns
    being sample values (covariates × samples). If `NULL`, no covariate
    file is used.

-   assayName:

    Name to assign the phenotype matrix in `assays()`. Defaults to
    `"pheno"`.

-   genome:

    Genome build string (e.g. `"hg38"`) assigned to the `seqinfo` of
    `rowRanges`. Defaults to `NA` (unspecified).

</div>

<div class="section level2">

## Value

A tQTLExperiment object.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
exdir <- system.file("extdata", package = "tQTLExperiment")
tqe <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
    covFile     = file.path(exdir, "cov-n100-tqtl.tsv"),
    genome      = "hg38"
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
tqe
#> class: tQTLExperiment
#> features: 20  samples: 100 
#> assays( 1 ): pheno 
#> rowRanges: GRanges with 20 features
#> colData( 11 ) covariates: int, sex, age, expr_pc1 ... 
#> geno: 100 samples x 69638 variants [BEDMatrix - lazy]
#> plinkPrefix: /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpotkzK9/temp_libpath10c16232f7b5e/tQTLExperiment/extdata/chr22-n100 
#> use prepareTQTL() to write inputs and get CLI command
```

</div>

</div>

</div>
