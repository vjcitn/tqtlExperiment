<div id="main" class="col-md-9" role="main">

# Find a Python installation with tensorQTL

<div class="ref-description section level2">

Searches for a Python executable that has tensorQTL installed, for use
in constructing the command string returned by `prepareTQTL()`. The
search checks common conda/mamba locations before falling back to PATH.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
findTQTL()
```

</div>

</div>

<div class="section level2">

## Value

Path to a Python executable as a character string.

</div>

<div class="section level2">

## Details

Because `prepareTQTL()` returns a shell command for the user to run
separately, this function is only needed to populate the `python`
argument of `prepareTQTL()`. Set `options(tQTLExperiment.python)` to
skip the search entirely.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
tryCatch(findTQTL(), error = function(e) invisible(NULL))
#> [1] "/Users/vincentcarey/miniforge3/bin/python3"
```

</div>

</div>

</div>
