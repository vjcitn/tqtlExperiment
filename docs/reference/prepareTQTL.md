# Prepare input files and return the tensorQTL CLI command. Writes the phenotype BED file and (optionally) the covariate file to a user-specified directory, then returns the complete `python -m tensorqtl` command string. The user runs this command in a terminal where their Python/conda environment is properly configured, then calls [`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md) to load the results back into R.

Prepare input files and return the tensorQTL CLI command. Writes the
phenotype BED file and (optionally) the covariate file to a
user-specified directory, then returns the complete
`python -m tensorqtl` command string. The user runs this command in a
terminal where their Python/conda environment is properly configured,
then calls
[`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)
to load the results back into R.

## Usage

``` r
prepareTQTL(
  x,
  outDir,
  mode = c("cis_nominal", "cis", "cis_independent", "trans"),
  assayName = NULL,
  mafThreshold = 0.05,
  window = 1000000L,
  permutations = 1000L,
  python = .detect_python(),
  ...
)
```

## Arguments

- x:

  A
  [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

- outDir:

  Path to an existing directory that will hold both the intermediate
  input files and the tensorQTL output. **Must be specified by the
  caller** - no default is provided.

- mode:

  One of `"cis_nominal"`, `"cis"`, `"cis_independent"`, `"trans"`.
  Passed to `--mode`.

- assayName:

  Name of the assay to use as the phenotype matrix. Defaults to the
  first assay.

- mafThreshold:

  MAF threshold. Defaults to `0.05`.

- window:

  Cis-window in base pairs. Defaults to `1000000`.

- permutations:

  Permutations for `"cis"` mode. Defaults to `1000`.

- python:

  Path to the Python executable to embed in the returned command string.
  Defaults to auto-detection: uses the Python inside the active conda
  environment (`CONDA_PREFIX`) or virtualenv (`VIRTUAL_ENV`) if set,
  otherwise the first `python3` on `PATH`.

- ...:

  Additional `--key value` flags passed to tensorQTL verbatim.

## Value

A character string containing the complete shell command to run. The
string is also printed via
[`message()`](https://rdrr.io/r/base/message.html) for easy copy-paste.

## See also

[`readTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/readTQTL.md)
to load results after the command has been run.

## Examples

``` r
exdir <- system.file("extdata", package = "tQTLExperiment")
tqe <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
)
#> Extracting number of samples and rownames from chr22-n100.fam...
#> Extracting number of variants and colnames from chr22-n100.bim...
cmd <- prepareTQTL(tqe, outDir = tempdir(), mode = "cis_nominal")
#> Run the following command in a terminal with tensorqtl available:
#> 
#> /Users/vincentcarey/miniforge3/bin/python3 -m tensorqtl /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpcbUgAS/temp_libpath2147cb89ee3/tQTLExperiment/extdata/chr22-n100 /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmppIQXj2/pheno.bed /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmppIQXj2/tqtl_out --mode cis_nominal --maf_threshold 0.05 --window 1000000 -o /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmppIQXj2
#> 
#> Then call readTQTL('/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmppIQXj2', mode = 'cis_nominal') to load results into R.
cat(cmd, "\n")
#> /Users/vincentcarey/miniforge3/bin/python3 -m tensorqtl /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpcbUgAS/temp_libpath2147cb89ee3/tQTLExperiment/extdata/chr22-n100 /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmppIQXj2/pheno.bed /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmppIQXj2/tqtl_out --mode cis_nominal --maf_threshold 0.05 --window 1000000 -o /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmppIQXj2 
```
