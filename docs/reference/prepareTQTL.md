<div id="main" class="col-md-9" role="main">

# Prepare input files and return the tensorQTL CLI command

<div class="ref-description section level2">

Writes the phenotype BED file and (optionally) the covariate file to a
user-specified directory, then returns the complete
`python -m tensorqtl` command string. The user runs this command in a
terminal where their Python/conda environment is properly configured,
then calls `readTQTL()` to load the results back into R.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
prepareTQTL(
  x,
  outDir,
  mode = c("cis_nominal", "cis", "cis_independent", "trans"),
  assayName = NULL,
  mafThreshold = 0.05,
  window = 1000000L,
  permutations = 1000L,
  python = "python3",
  ...
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A
    [tQTLExperiment](https://github.com/vjcitn/tQTLExperiment/reference/tQTLExperiment.md).

-   outDir:

    Path to an existing directory that will hold both the intermediate
    input files and the tensorQTL output. **Must be specified by the
    caller** - no default is provided.

-   mode:

    One of `"cis_nominal"`, `"cis"`, `"cis_independent"`, `"trans"`.
    Passed to `--mode`.

-   assayName:

    Name of the assay to use as the phenotype matrix. Defaults to the
    first assay.

-   mafThreshold:

    MAF threshold. Defaults to `0.05`.

-   window:

    Cis-window in base pairs. Defaults to `1000000`.

-   permutations:

    Permutations for `"cis"` mode. Defaults to `1000`.

-   python:

    Path to the Python executable to embed in the returned command
    string. Defaults to `"python3"`.

-   ...:

    Additional `--key value` flags passed to tensorQTL verbatim.

</div>

<div class="section level2">

## Value

A character string containing the complete shell command to run. The
string is also printed via `message()` for easy copy-paste.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`readTQTL()` to load results after the command has been run.

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

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
#> python3 -m tensorqtl /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpZyaoNi/temp_libpatha318583333e0/tQTLExperiment/extdata/chr22-n100 /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmphg7Wmk/pheno.bed /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmphg7Wmk/tqtl_out --mode cis_nominal --maf_threshold 0.05 --window 1000000 -o /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmphg7Wmk
#> 
#> Then call readTQTL('/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmphg7Wmk', mode = 'cis_nominal') to load results into R.
cat(cmd, "\n")
#> python3 -m tensorqtl /private/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T/RtmpZyaoNi/temp_libpatha318583333e0/tQTLExperiment/extdata/chr22-n100 /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmphg7Wmk/pheno.bed /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmphg7Wmk/tqtl_out --mode cis_nominal --maf_threshold 0.05 --window 1000000 -o /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmphg7Wmk 
```

</div>

</div>

</div>
