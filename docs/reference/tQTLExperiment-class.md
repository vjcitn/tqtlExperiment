<div id="main" class="col-md-9" role="main">

# tQTLExperiment: a SummarizedExperiment for tensorQTL eQTL mapping

<div class="ref-description section level2">

`tQTLExperiment` extends
[SummarizedExperiment::RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
to hold all data needed for the tensorQTL eQTL mapping tool:

</div>

<div class="section level2">

## Details

-   **assays**: phenotype matrix (features × samples)

-   **rowRanges**: genomic coordinates for each phenotype feature (gene)

-   **colData**: per-sample metadata including covariates

-   **geno**: lazy
    [BEDMatrix::BEDMatrix](https://rdrr.io/pkg/BEDMatrix/man/BEDMatrix.html)
    over the PLINK `.bed` file (samples × variants); data are not loaded
    until subscripted

-   **variantRanges**:
    [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
    for each variant, with metadata columns `snp_id`, `ref`, `alt`

-   **plinkPrefix**: path prefix for the PLINK file set

tensorQTL is invoked via its CLI (`python -m tensorqtl`) and supports
GPU acceleration through CUDA or Apple MPS (Metal Performance Shaders)
on Apple Silicon.

</div>

<div class="section level2">

## Slots

-   `geno`:

    A
    [BEDMatrix::BEDMatrix](https://rdrr.io/pkg/BEDMatrix/man/BEDMatrix.html)
    (samples × variants).

-   `variantRanges`:

    A
    [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
    of length equal to `ncol(geno)`.

-   `plinkPrefix`:

    A length-one character string giving the PLINK file prefix.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`tQTLExperiment()` for the constructor, `prepareTQTL()` to write input
files and get the CLI command, `readTQTL()` to load results,
`findTQTL()` to locate Python.

</div>

</div>

</div>
