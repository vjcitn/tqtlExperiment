# Annotate tensorQTL results with gene symbols and biotypes

Adds gene symbols (HGNC names) and gene biotypes to a
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object containing tensorQTL results. Useful for annotating the output
from
[`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md).

## Usage

``` r
annotateTensorQTLResults(gr, phenotype_col = "phenotype_id")
```

## Arguments

- gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object with tensorQTL results, typically from
  [`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md).

- phenotype_col:

  Name of the metadata column in `gr` containing Ensembl gene IDs
  (ENSG...). Default: `"phenotype_id"`.

## Value

The input `gr` with two additional metadata columns:

- `gene_symbol`: HGNC gene symbol

- `gene_biotype`: Ensembl gene biotype (e.g., "protein_coding",
  "lincRNA")

## Details

Ensembl version suffixes (e.g., `.3` in `ENSG00000100181.3`) are
automatically stripped before lookup in EnsDb.Hsapiens.v79.

## Examples

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
