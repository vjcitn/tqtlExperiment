<div id="main" class="col-md-9" role="main">

# Annotate tensorQTL results with gene symbols and biotypes

<div class="ref-description section level2">

Adds gene symbols (HGNC names) and gene biotypes to a
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object containing tensorQTL results. Useful for annotating the output
from `readTQTL()`.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
annotateTensorQTLResults(gr, phenotype_col = "phenotype_id")
```

</div>

</div>

<div class="section level2">

## Arguments

-   gr:

    A
    [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
    object with tensorQTL results, typically from `readTQTL()`.

-   phenotype\_col:

    Name of the metadata column in `gr` containing Ensembl gene IDs
    (ENSG...). Default: `"phenotype_id"`.

</div>

<div class="section level2">

## Value

The input `gr` with two additional metadata columns:

-   `gene_symbol`: HGNC gene symbol

-   `gene_biotype`: Ensembl gene biotype (e.g., "protein\_coding",
    "lincRNA")

</div>

<div class="section level2">

## Details

Ensembl version suffixes (e.g., `.3` in `ENSG00000100181.3`) are
automatically stripped before lookup in EnsDb.Hsapiens.v79.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
demo_data <- cache_demo_data()
#> [ cache hit   ] covariates.tsv
#> [ cache hit   ] mage_chr17_cis.rda
#> [ cache hit   ] pheno.bed
#> [ cache hit   ] tqtl_out.cis_qtl_pairs.17.parquet
#> [ cache hit   ] tqtl_out.tensorQTL.cis_nominal.log
#> [ complete    ] All demo files cached in: /Users/vincentcarey/Library/Caches/org.R-project.R/R/BiocFileCache/tqtlExperiment_demodir
demodir = dirname(demo_data[["parquet"]])
res <- readTQTL(demodir, mode = "cis_nominal")
res_anno <- annotateTensorQTLResults(res$pairs)
head(mcols(res_anno)[, c("phenotype_id", "gene_symbol", "gene_biotype")])
#> DataFrame with 6 rows and 3 columns
#>        phenotype_id    gene_symbol gene_biotype
#>         <character>    <character>  <character>
#> 1 ENSG00000280279.1 RP11-1228E12.1      lincRNA
#> 2 ENSG00000280279.1 RP11-1228E12.1      lincRNA
#> 3 ENSG00000280279.1 RP11-1228E12.1      lincRNA
#> 4 ENSG00000280279.1 RP11-1228E12.1      lincRNA
#> 5 ENSG00000280279.1 RP11-1228E12.1      lincRNA
#> 6 ENSG00000280279.1 RP11-1228E12.1      lincRNA
```

</div>

</div>

</div>
