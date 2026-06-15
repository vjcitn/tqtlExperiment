#' Find a Python installation with tensorQTL
#'
#' Searches for a Python executable that has tensorQTL installed, for use
#' in constructing the command string returned by [prepareTQTL()].  The
#' search checks common conda/mamba locations before falling back to PATH.
#'
#' Because [prepareTQTL()] returns a shell command for the user to run
#' separately, this function is only needed to populate the `python`
#' argument of [prepareTQTL()].  Set `options(tQTLExperiment.python)` to
#' skip the search entirely.
#'
#' @return Path to a Python executable as a character string.
#'
#' @examples
#' tryCatch(findTQTL(), error = function(e) invisible(NULL))
#'
#' @export
findTQTL <- function() {
    opt <- getOption("tQTLExperiment.python")
    if (!is.null(opt)) return(opt)

    for (py in .python_candidates()) {
        if (file.exists(py)) return(py)
    }

    # fall back to PATH names
    for (nm in c("python3", "python")) {
        py <- Sys.which(nm)
        if (nzchar(py)) return(py)
    }

    stop(
        "No Python executable found. Set the path with:\n",
        "  options(tQTLExperiment.python = '/path/to/python3')"
    )
}

.python_candidates <- function() {
    home <- path.expand("~")
    roots <- file.path(home, c("miniforge3", "mambaforge",
                                "miniconda3",  "anaconda3",
                                "opt/miniconda3", "opt/anaconda3"))
    base_pythons <- file.path(roots, "bin", "python3")
    env_pythons  <- unlist(lapply(roots, function(r) {
        envs <- list.dirs(file.path(r, "envs"),
                          full.names = TRUE, recursive = FALSE)
        file.path(envs, "bin", "python3")
    }))
    unique(c(base_pythons, env_pythons,
             "/usr/local/bin/python3", "/opt/homebrew/bin/python3"))
}

#' Report the tensorQTL version from a Python installation
#'
#' @param python Path to the Python executable. Defaults to [findTQTL()].
#' @return A character string with the version, or `NA`.
#' @examples
#' tryCatch(tensorqtlVersion(), error = function(e) NA_character_)
#' @export
tensorqtlVersion <- function(python = findTQTL()) {
    out <- tryCatch(
        system2(python, stdout = TRUE, stderr = FALSE,
                input = "import tensorqtl; print(tensorqtl.__version__)"),
        error = function(e) character(0)
    )
    if (length(out) == 0L || !nzchar(out[1L])) NA_character_ else out[1L]
}
