<div id="main" class="col-md-9" role="main">

# Report the tensorQTL version from a Python installation

<div class="ref-description section level2">

Report the tensorQTL version from a Python installation

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
tensorqtlVersion(python = findTQTL())
```

</div>

</div>

<div class="section level2">

## Arguments

-   python:

    Path to the Python executable. Defaults to `findTQTL()`.

</div>

<div class="section level2">

## Value

A character string with the version, or `NA`.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
tryCatch(tensorqtlVersion(), error = function(e) NA_character_)
#> Warning: running command ''/Users/vincentcarey/miniforge3/bin/python3' 2>/dev/null < '/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//RtmpkmYph4/file10c4b11fbde70'' had status 1
#> [1] "Warning: 'rfunc' cannot be imported. R with the 'qvalue' library and the 'rpy2' Python package are needed to compute q-values."
```

</div>

</div>

</div>
