# Package index

## All functions

- [`addGeneSymbols()`](https://github.com/vjcitn/tQTLExperiment/reference/addGeneSymbols.md)
  : Add gene symbols to rowRanges

- [`annotateTensorQTLResults()`](https://github.com/vjcitn/tQTLExperiment/reference/annotateTensorQTLResults.md)
  : Annotate tensorQTL results with gene symbols and biotypes

- [`cache_demo_data()`](https://github.com/vjcitn/tQTLExperiment/reference/cache_demo_data.md)
  : Cache demo dataset files

- [`cache_mage_chr17_cis()`](https://github.com/vjcitn/tQTLExperiment/reference/cache_mage_chr17_cis.md)
  : Cache pre-computed MAGE chr17 tensorQTL results

- [`cache_mage_chr17_plink()`](https://github.com/vjcitn/tQTLExperiment/reference/cache_mage_chr17_plink.md)
  : Cache MAGE chr17 PLINK files from GitHub releases

- [`collapseFactorTstats()`](https://github.com/vjcitn/tQTLExperiment/reference/collapseFactorTstats.md)
  : Collapse factor-expanded t-statistic columns into group summaries

- [`findTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/findTQTL.md)
  : Find a Python installation with tensorQTL

- [`get_cached_mage17SE()`](https://github.com/vjcitn/tQTLExperiment/reference/get_cached_mage17SE.md)
  : Cache and load the MAGE chr17 SummarizedExperiment

- [`plotGenotypeEffect()`](https://github.com/vjcitn/tQTLExperiment/reference/plotGenotypeEffect.md)
  : Plot phenotype values grouped by SNP genotype

- [`prepareTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/prepareTQTL.md)
  :

  Prepare input files and return the tensorQTL CLI command. Writes the
  phenotype BED file and (optionally) the covariate file to a
  user-specified directory, then returns the complete
  `python -m tensorqtl` command string. The user runs this command in a
  terminal where their Python/conda environment is properly configured,
  then calls
  [`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)
  to load the results back into R.

- [`qtlPCABrowser()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlPCABrowser.md)
  : Interactive PCA browser for cis-QTL t-statistics

- [`qtlRegressionStats()`](https://github.com/vjcitn/tQTLExperiment/reference/qtlRegressionStats.md)
  : Linear regression statistics for a collection of SNP-gene pairs

- [`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)
  : Read tensorQTL results into R

- [`run_cis_qtl_browser()`](https://github.com/vjcitn/tQTLExperiment/reference/run_cis_qtl_browser.md)
  :

  Interactive cis-QTL browser. Launches a Shiny application for
  exploring cis-QTL associations. Users can select genes, view
  associated SNPs, and visualize the genotype effect on phenotypes using
  [`plotGenotypeEffect()`](https://github.com/vjcitn/tQTLExperiment/reference/plotGenotypeEffect.md).

- [`show(`*`<tQTLExperiment>`*`)`](https://github.com/vjcitn/tQTLExperiment/reference/show-tQTLExperiment-method.md)
  : Show a tQTLExperiment

- [`` `[`( ``*`<tQTLExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://github.com/vjcitn/tQTLExperiment/reference/sub-tQTLExperiment-ANY-ANY-ANY-method.md)
  : Subset a tQTLExperiment

- [`tQTLExperiment-class`](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment-class.md)
  :

  tQTLExperiment: a SummarizedExperiment for tensorQTL eQTL mapping.
  `tQTLExperiment` extends
  [SummarizedExperiment::RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  to hold all data needed for the tensorQTL eQTL mapping tool:

- [`tQTLExperiment()`](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md)
  : Construct a tQTLExperiment from input files

- [`tQTLExperimentFromRSE()`](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperimentFromRSE.md)
  : Construct a tQTLExperiment from a RangedSummarizedExperiment

- [`tensorqtlVersion()`](https://github.com/vjcitn/tQTLExperiment/reference/tensorqtlVersion.md)
  : Report the tensorQTL version from a Python installation

- [`tqtlGeno()`](https://github.com/vjcitn/tQTLExperiment/reference/tqtlGeno.md)
  [`` `tqtlGeno<-`() ``](https://github.com/vjcitn/tQTLExperiment/reference/tqtlGeno.md)
  : Access the lazy genotype matrix

- [`tqtlPlinkPrefix()`](https://github.com/vjcitn/tQTLExperiment/reference/tqtlPlinkPrefix.md)
  [`` `tqtlPlinkPrefix<-`() ``](https://github.com/vjcitn/tQTLExperiment/reference/tqtlPlinkPrefix.md)
  : Access the PLINK file prefix

- [`tqtlVariantRanges()`](https://github.com/vjcitn/tQTLExperiment/reference/tqtlVariantRanges.md)
  [`` `tqtlVariantRanges<-`() ``](https://github.com/vjcitn/tQTLExperiment/reference/tqtlVariantRanges.md)
  : Access variant genomic ranges
