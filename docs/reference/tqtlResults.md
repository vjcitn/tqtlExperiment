# Query bound tensorQTL cis_nominal results lazily

Queries a `cis_nominal` output directory bound to a
[tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
via
[`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)
or `tqtlResultsDir<-`, without loading the whole genome-wide table into
memory. Filtering is pushed down to the parquet files by
[`arrow::open_dataset()`](https://arrow.apache.org/docs/r/reference/open_dataset.html),
so only the matching rows (e.g. one gene's pairs) are ever materialized
in R.

## Usage

``` r
tqtlResults(
  x,
  phenotype_id = NULL,
  variant_id = NULL,
  as = c("GRanges", "data.frame")
)
```

## Arguments

- x:

  A
  [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
  with `tqtlResultsDir(x)` set.

- phenotype_id:

  Optional character vector of phenotype/gene IDs to filter to. If
  `NULL` (default), no phenotype filter is applied.

- variant_id:

  Optional character vector of variant IDs to filter to. If `NULL`
  (default), no variant filter is applied.

- as:

  One of `"GRanges"` (default) or `"data.frame"`.

## Value

A
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
(or data frame) of the matching `cis_nominal` rows.

## See also

[`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)
to bind a results directory after a run, `tqtlResultsDir<-` to bind one
directly.

## Examples

``` r
if (FALSE) { # \dontrun{
tqe <- readTQTL(outdir, mode = "cis_nominal", x = tqe)$tqe
tqtlResults(tqe, phenotype_id = "ENSG00000187624.9")
} # }
```
