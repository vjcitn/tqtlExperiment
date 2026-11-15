# Collapse factor-expanded t-statistic columns into group summaries

When a factor variable (e.g. `batch`) or a multi-level character
variable (e.g. `population`) is included in the regression model, the
wide output of
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
contains one `t.<var><level>` column per level. This function identifies
groups of such columns by regular-expression patterns, replaces each
group with the root-mean-square (RMS) of its t-statistics, and returns a
leaner data frame suitable for PCA.

## Usage

``` r
collapseFactorTstats(
  wide_df,
  patterns = list(batch = "^t\\.factor_batch_", population = "^t\\.population")
)
```

## Arguments

- wide_df:

  A wide data frame from
  [`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
  with `t_only = TRUE`.

- patterns:

  A named list of regular expressions. Each element matches a set of
  columns to collapse; the name becomes the new column name (prefixed
  with `"t."`). Defaults to patterns that detect the MAGE chr17
  dataset's batch and population indicator columns.

## Value

`wide_df` with matched column groups replaced by single `t.<name>` RMS
columns.

## Examples

``` r
exdir <- system.file("extdata", package = "tQTLExperiment")
tqe   <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
    covFile     = file.path(exdir, "cov-n100-tqtl.tsv"),
    genome      = "hg38"
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
tqe <- addGeneSymbols(tqe)
res <- qtlRegressionStats(tqe, symbol = "TPTEP1")
collapsed <- collapseFactorTstats(res)
```
