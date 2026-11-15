# Find a Python installation with tensorQTL

Searches for a Python executable that has tensorQTL installed, for use
in constructing the command string returned by
[`prepareTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/prepareTQTL.md).
The search checks common conda/mamba locations before falling back to
PATH.

## Usage

``` r
findTQTL()
```

## Value

Path to a Python executable as a character string.

## Details

Because
[`prepareTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/prepareTQTL.md)
returns a shell command for the user to run separately, this function is
only needed to populate the `python` argument of
[`prepareTQTL()`](https://github.com/vjcitn/tQTLExperiment/reference/prepareTQTL.md).
Set `options(tQTLExperiment.python)` to skip the search entirely.

## Examples

``` r
tryCatch(findTQTL(), error = function(e) invisible(NULL))
#> [1] "/Users/vincentcarey/miniforge3/bin/python3"
```
