# tQTLExperiment: a SummarizedExperiment for tensorQTL eQTL mapping. `tQTLExperiment` extends [SummarizedExperiment::RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html) to hold all data needed for the tensorQTL eQTL mapping tool:

- **assays**: phenotype matrix (features × samples)

- **rowRanges**: genomic coordinates for each phenotype feature (gene)

- **colData**: per-sample metadata including covariates

- **geno**: lazy
  [BEDMatrix::BEDMatrix](https://rdrr.io/pkg/BEDMatrix/man/BEDMatrix.html)
  over the PLINK `.bed` file (samples × variants); data are not loaded
  until subscripted

- **variantRanges**:
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  for each variant, with metadata columns `snp_id`, `ref`, `alt`

- **plinkPrefix**: path prefix for the PLINK file set

- **tqtlResultsDir**: path to a directory of tensorQTL `cis_nominal`
  parquet output, bound by
  [`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)
  after a successful run, and queried lazily (without loading the whole
  table) via
  [`tqtlResults()`](https://github.com/vjcitn/tQTLExperiment/reference/tqtlResults.md)

## Details

tensorQTL is invoked via its CLI (`python -m tensorqtl`) and supports
GPU acceleration through CUDA or Apple MPS (Metal Performance Shaders)
on Apple Silicon.

## Slots

- `geno`:

  A
  [BEDMatrix::BEDMatrix](https://rdrr.io/pkg/BEDMatrix/man/BEDMatrix.html)
  (samples × variants).

- `variantRanges`:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of length equal to `ncol(geno)`.

- `plinkPrefix`:

  A length-one character string giving the PLINK file prefix.

- `tqtlResultsDir`:

  A length-one character string giving the directory of bound
  `cis_nominal` tensorQTL output, or `NA_character_` if none has been
  bound yet.

## See also

[`tQTLExperiment()`](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
for the constructor,
[`prepareTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/prepareTQTL.md)
to write input files and get the CLI command,
[`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)
to load results,
[`tqtlResults()`](https://github.com/vjcitn/tQTLExperiment/reference/tqtlResults.md)
to query bound results lazily,
[`findTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/findTQTL.md)
to locate Python.
