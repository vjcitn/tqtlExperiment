<div id="main" class="col-md-9" role="main">

# Plot phenotype values grouped by SNP genotype

<div class="ref-description section level2">

Creates a beeswarm plot comparing a phenotype across genotype groups (0,
1, 2 copies of the alternate allele) for a given SNP. Useful for
visualizing the effect of a significant eQTL variant.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
plotGenotypeEffect(
  x,
  snp_id,
  phenotype_id,
  assayName = NULL,
  size = 2,
  color = "steelblue",
  title = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A
    [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

-   snp\_id:

    Character. SNP identifier (row name in `variantRanges(x)`).

-   phenotype\_id:

    Character. Phenotype/gene identifier (row name in `rowRanges(x)`).

-   assayName:

    Name of the assay to use. Defaults to the first assay.

-   size:

    Point size for beeswarm. Default 2.

-   color:

    Color for points. Default "steelblue".

-   title:

    Main plot title. If `NULL`, auto-generated from SNP and gene.

</div>

<div class="section level2">

## Value

A ggplot2 plot object.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (!requireNamespace("CSHLvc2026")) BiocManager::install("vjcitn/CSHLvc2026")
#> Loading required namespace: CSHLvc2026
data(mageSEfilt, package="CSHLvc2026")
plink_paths <- cache_mage_chr17_plink()
#> [ cache hit   ] CCDG_mage_chr17.fam
#> [ cache hit   ] CCDG_mage_chr17.bim
#> [ cache hit   ] CCDG_mage_chr17.bed
#> [ validate    ] Checking .bed magic bytes ...
#> [ validate    ] .bed magic bytes OK.
plpre <- tools::file_path_sans_ext(plink_paths[["bed"]])
cd <- as.data.frame(colData(mageSEfilt))
mm <- model.matrix(~ batch + population + sex, data = cd)
mm <- mm[, -1, drop = FALSE]   # remove (Intercept)
tqe <- tQTLExperimentFromRSE(
    se              = mageSEfilt,
    plinkPrefix     = plpre,
    covariateMatrix = mm,
    genome          = "hg38"
)
#> Extracting number of samples and rownames from CCDG_mage_chr17.fam...
#> Extracting number of variants and colnames from CCDG_mage_chr17.bim...
plotGenotypeEffect(tqe, "17:410391:T:C", "ENSG00000187624.9")
```

</div>

</div>

</div>
