# Linear regression statistics for a collection of SNP-gene pairs

For each (phenotype_id, variant_id) pair, fits a linear model of
phenotype on genotype plus all covariates stored in the
[tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
and returns a tidy data frame of coefficients, standard errors,
t-statistics, and p-values for every model term. Pairs are grouped by
phenotype and processed in parallel via
[`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html).

## Usage

``` r
qtlRegressionStats(
  tqe,
  pairs = NULL,
  symbol = NULL,
  window = 1000000L,
  assayName = NULL,
  t_only = TRUE,
  wideFormat = TRUE,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- tqe:

  A
  [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
  with genotype data.

- pairs:

  A data frame with columns `phenotype_id` and `variant_id` identifying
  the SNP-gene pairs to analyse. Mutually exclusive with `symbol`.

- symbol:

  Character scalar. Gene symbol to look up in the `gene_name` column of
  `rowRanges(tqe)`. All variants within `window` bp of the gene's start
  are used as the SNP set. The `tQTLExperiment` must have a `gene_name`
  column added via
  [`addGeneSymbols()`](https://github.com/vjcitn/tQTLExperiment/reference/addGeneSymbols.md);
  the function stops with an informative error if it is absent. Mutually
  exclusive with `pairs`.

- window:

  Integer. Cis-window half-width in bp used when `symbol` is specified.
  Defaults to `1000000L`.

- assayName:

  Name of the assay to use. Defaults to the first assay.

- t_only:

  Logical. If `TRUE` (default), return only the t-statistic (estimate /
  std.error) for each model term. In wide format columns are named
  `<term>.t`. If `FALSE`, return estimates and standard errors (and
  additionally statistic and p.value in long format).

- wideFormat:

  Logical. If `TRUE` (default), pivot to one row per pair. With
  `t_only = TRUE`: one column per term named `<term>.t`. With
  `t_only = FALSE`: two columns per term, `estimate.<term>` and
  `std.error.<term>`. If `FALSE`, return the long tidy form with one row
  per pair x term.

- BPPARAM:

  A
  [BiocParallel::BiocParallelParam](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object controlling parallelisation. Defaults to
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  (single core). Use e.g.
  [`BiocParallel::MulticoreParam()`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)
  for multicore.

## Value

A data frame. See `t_only` and `wideFormat` for the shape of the output.

## Details

**Chromosome naming discrepancy.** The `rowRanges` of a `tQTLExperiment`
built from a `SummarizedExperiment` typically carries UCSC-style
seqnames (e.g. `"chr17"`), while
[`tqtlVariantRanges()`](https://github.com/vjcitn/tQTLExperiment/reference/tqtlVariantRanges.md)
reflects the PLINK BIM file convention of bare chromosome numbers (e.g.
`"17"`). When `symbol` is used to select cis-variants, both seqname
vectors are normalised by stripping any leading `"chr"` before
comparison, so the window filter works regardless of which convention
each source uses.

If you write code that compares or intersects `rowRanges` and
`tqtlVariantRanges` directly, apply the same normalisation:


      sub("^chr", "", as.character(seqnames(gr)))

## Examples

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
pairs <- data.frame(
    phenotype_id = rownames(tqe)[1],
    variant_id   = S4Vectors::mcols(tqtlVariantRanges(tqe))[["snp_id"]][1]
)
qtlRegressionStats(tqe, pairs)
#>      phenotype_id     variant_id t._Intercept_ t.genotype         t.sex
#> 1 ENSG00000100181 22:16849573A-G      1.812558  -1.067511 -0.0003593606
#>        t.age t.expr_pc1 t.expr_pc2 t.geno_pc1 t.geno_pc2 t.geno_pc3 t.geno_pc4
#> 1 -0.7850642  -1.394788 -0.7426527   1.167902   1.021817   1.112029 -0.6122195
#>   t.geno_pc5 t.geno_pc6
#> 1  0.2154714   0.972025
qtlRegressionStats(tqe, pairs, t_only = FALSE, wideFormat = FALSE)
#>       phenotype_id     variant_id        term      estimate    std.error
#> 1  ENSG00000100181 22:16849573A-G (Intercept)  1.785488e-02 0.0098506513
#> 2  ENSG00000100181 22:16849573A-G    genotype -2.065937e-03 0.0019352838
#> 3  ENSG00000100181 22:16849573A-G         sex -1.025895e-06 0.0028547777
#> 4  ENSG00000100181 22:16849573A-G         age -1.024120e-04 0.0001304505
#> 5  ENSG00000100181 22:16849573A-G    expr_pc1 -3.047032e-03 0.0021845848
#> 6  ENSG00000100181 22:16849573A-G    expr_pc2 -1.933679e-03 0.0026037463
#> 7  ENSG00000100181 22:16849573A-G    geno_pc1  3.377120e-02 0.0289161266
#> 8  ENSG00000100181 22:16849573A-G    geno_pc2  5.047878e-02 0.0494010174
#> 9  ENSG00000100181 22:16849573A-G    geno_pc3  4.345973e-02 0.0390814609
#> 10 ENSG00000100181 22:16849573A-G    geno_pc4 -3.190808e-02 0.0521186941
#> 11 ENSG00000100181 22:16849573A-G    geno_pc5  1.060592e-02 0.0492219665
#> 12 ENSG00000100181 22:16849573A-G    geno_pc6  4.014786e-02 0.0413033217
#>        statistic    p.value
#> 1   1.8125578678 0.07331014
#> 2  -1.0675113501 0.28866080
#> 3  -0.0003593606 0.99971409
#> 4  -0.7850642387 0.43452439
#> 5  -1.3947876087 0.16659070
#> 6  -0.7426527250 0.45966982
#> 7   1.1679018389 0.24600089
#> 8   1.0218166139 0.30966858
#> 9   1.1120292605 0.26915435
#> 10 -0.6122195229 0.54197122
#> 11  0.2154713740 0.82989845
#> 12  0.9720250353 0.33370234
```
