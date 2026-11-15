# Construct a tQTLExperiment from a RangedSummarizedExperiment

Converts an existing
[SummarizedExperiment::RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
into a
[tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
by attaching PLINK genotype data.

## Usage

``` r
tQTLExperimentFromRSE(
  se,
  plinkPrefix,
  covariateMatrix = NULL,
  assayName = NULL,
  featureIdColumn = NULL,
  genome = NA_character_
)
```

## Arguments

- se:

  A
  [SummarizedExperiment::RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html).

- plinkPrefix:

  Path prefix for the PLINK file set.

- covariateMatrix:

  A numeric matrix (rows = samples, cols = covariates) built e.g. with
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  on `as.data.frame(colData(se))`. See also the
  [ExploreModelMatrix](https://bioconductor.org/packages/ExploreModelMatrix)
  package. Rename `(Intercept)` to `int` if including an intercept.

- assayName:

  Name of the assay in `se` to use. Defaults to first.

- featureIdColumn:

  Optional mcols column name to use as `phenotype_id`; defaults to
  `rownames(se)`.

- genome:

  Genome build string (e.g. `"hg38"`) assigned to the `seqinfo` of
  `rowRanges`. Defaults to `NA` (unspecified).

## Value

A
[tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
object.

## Examples

``` r
exdir <- system.file("extdata", package = "tQTLExperiment")
tqe_ref <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
se <- as(tqe_ref, "RangedSummarizedExperiment")
cd <- as.data.frame(colData(tqe_ref))
mm <- model.matrix(~ 1, data = cd)
colnames(mm)[1] <- "int"
tqe2 <- tQTLExperimentFromRSE(se,
    plinkPrefix     = file.path(exdir, "chr22-n100"),
    covariateMatrix = mm)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
tqe2
#> class: tQTLExperiment
#> features: 20  samples: 100 
#> assays( 1 ): pheno 
#> rowRanges: GRanges with 20 features
#> colData( 1 ) covariates: int  
#> geno: 100 samples x 69638 variants [BEDMatrix - lazy]
#> plinkPrefix: /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpOhRh54/temp_libpath175734e41dcb3/tQTLExperiment/extdata/chr22-n100 
#> use prepareTQTL() to write inputs and get CLI command
```
