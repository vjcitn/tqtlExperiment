# Interactive PCA browser for cis-QTL t-statistics

Collapses factor-expanded t-statistic columns via
[`collapseFactorTstats()`](https://github.com/vjcitn/tQTLExperiment/reference/collapseFactorTstats.md),
computes PCA on the result, and launches a Shiny application with a
selectable PC-vs-PC scatter plot (plotly) and a pairs plot of the first
five PCs. Hovering over a SNP-gene point in the scatter plot displays a
beeswarm of gene expression by genotype for that pair, with an optional
colour-by selector for sample-level variables (e.g. sex, batch). If
regression statistics are not supplied, the app provides a gene-symbol
selector and computes/memoises
[`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
results within the Shiny session.

## Usage

``` r
qtlPCABrowser(
  res = NULL,
  tqe = NULL,
  assayName = NULL,
  collapse_patterns = list(batch = "^t\\.factor_batch_", population =
    "^t\\.population"),
  window = 1000000L,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- res:

  Optional wide data frame from
  [`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
  with `t_only = TRUE`. Must contain `phenotype_id` and `variant_id`
  columns. If omitted, the app lists available gene symbols from
  `rowRanges(tqe)$gene_name` and computes
  `qtlRegressionStats(tqe, symbol = ...)` on request.

- tqe:

  A
  [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
  with genotype data.

- assayName:

  Name of the assay to use. Defaults to the first assay.

- collapse_patterns:

  Named list of regex patterns passed to
  [`collapseFactorTstats()`](https://github.com/vjcitn/tQTLExperiment/reference/collapseFactorTstats.md).
  Defaults handle MAGE-style batch and population columns.

- window:

  Integer. Cis-window half-width passed to
  [`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
  when `res` is omitted.

- BPPARAM:

  A
  [BiocParallel::BiocParallelParam](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object passed to
  [`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
  when `res` is omitted.

## Value

A Shiny application object (invisibly).

## Examples

``` r
if (interactive()) {
  exdir <- system.file("extdata", package = "tQTLExperiment")
  tqe   <- tQTLExperiment(
      plinkPrefix = file.path(exdir, "chr22-n100"),
      phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
      covFile     = file.path(exdir, "cov-n100-tqtl.tsv"),
      genome      = "hg38"
  )
  tqe <- addGeneSymbols(tqe)
  qtlPCABrowser(tqe)
}
```
