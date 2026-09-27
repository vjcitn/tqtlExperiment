# Recent developments on the PCA browser branch

This branch collects recent work on `tQTLExperiment` focused on testing support
and extending `qtlPCABrowser()` from a precomputed PCA display into a more
complete interactive workflow for gene-specific cis-QTL regression summaries.

## Package version

The branch currently reports:

```text
Version: 0.1.49
```

## Test infrastructure

- Added a standard `tests/testthat.R` entry point.
- Added unit tests for `collapseFactorTstats()`.
- The tests check that factor-expanded t-statistic columns are collapsed by RMS
  into group-level summaries and that unmatched/singleton groups are left
  unchanged.

## `qtlPCABrowser()` UI and PCA enhancements

The PCA browser was expanded beyond its original fixed PC1-vs-PC2 scatter plot:

- The app title was changed to **"cis-QTL t-stat PCA"**.
- The plotly PCA scatter plot now has selectable X and Y principal components.
- A new **PC pairs** panel displays pairwise plots for the first five PCs.
- PCA hover events are scoped to the scatter plot so linked genotype-effect
  plots remain tied to the selected PCA point.

## SNP-gene pair labeling

The browser now treats each PCA point as a specific SNP-gene regression result,
rather than as a SNP alone.

- Hover text includes both the SNP and gene information.
- The selected point is keyed by row identity to avoid ambiguity when the same
  SNP appears with multiple genes.
- Genotype-effect plots now label the response as gene expression for the
  selected gene.
- User-facing labels use gene symbols when available, while retaining ENSG IDs
  in hover text for disambiguation.

## Gene-driven Shiny workflow

`qtlPCABrowser()` can now be launched directly from a `tQTLExperiment`:

```r
qtlPCABrowser(tqe)
```

In this mode, the app:

- reads available gene symbols from `rowRanges(tqe)$gene_name`;
- presents them in a `selectInput`;
- computes `qtlRegressionStats(tqe, symbol = selected_gene, ...)` on request;
- memoises the regression-statistics result within the Shiny session so the
  same gene/window/assay combination is not recomputed; and
- reuses the same PCA scatter, PC pairs, and genotype-effect plot UI.

The previous precomputed-result workflow is preserved:

```r
res <- qtlRegressionStats(tqe, symbol = "TPTEP1")
qtlPCABrowser(res, tqe)
```

## Documentation updates

The manual page `man/qtlPCABrowser.Rd` was updated to describe:

- optional precomputed regression statistics;
- the gene-symbol selector workflow;
- `window` and `BPPARAM` arguments passed through to `qtlRegressionStats()`;
- selectable PC axes;
- first-five-PC pairs plots; and
- memoised regression-statistics computation within a Shiny session.

## Validation performed

The current test suite was run after the changes:

```text
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 6 ]
```

Both supported browser launch modes were also instantiated on bundled chr22
example data:

```r
qtlPCABrowser(res, tqe)
qtlPCABrowser(tqe)
```
