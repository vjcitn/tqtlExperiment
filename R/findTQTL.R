#' Find a Python installation with tensorQTL
#'
#' Searches for a Python executable that has the tensorQTL package installed.
#' The search order is:
#' 1. `options(tQTLExperiment.python)` — explicit user preference
#' 2. `python3` on PATH
#' 3. `python` on PATH
#'
#' @return The path to the Python executable as a character string.
#'
#' @examples
#' tryCatch(findTQTL(), error = function(e) invisible(NULL))
#'
#' @export
findTQTL <- function() {
    opt <- getOption("tQTLExperiment.python")
    if (!is.null(opt)) {
        .check_tensorqtl(opt)
        return(opt)
    }
    for (py in c("python3", "python")) {
        path <- Sys.which(py)
        if (nzchar(path)) {
            if (isTRUE(.check_tensorqtl(path, error = FALSE)))
                return(path)
        }
    }
    stop(
        "No Python installation with tensorqtl found. ",
        "Install with: pip install tensorqtl\n",
        "Then set options(tQTLExperiment.python = '/path/to/python3') ",
        "if it is not on PATH."
    )
}

.check_tensorqtl <- function(python, error = TRUE) {
    res <- system2(python, args = c("-c", "import tensorqtl"),
                   stdout = FALSE, stderr = FALSE)
    if (res != 0L) {
        if (error)
            stop("tensorqtl not importable from '", python, "'. ",
                 "Install with: pip install tensorqtl")
        return(FALSE)
    }
    TRUE
}

#' Report the installed tensorQTL version
#'
#' @param python Path to the Python executable. Defaults to the result of
#'   [findTQTL()].
#' @return A character string with the tensorQTL version, or `NA` if not
#'   found.
#' @examples
#' tryCatch(tensorqtlVersion(), error = function(e) NA_character_)
#' @export
tensorqtlVersion <- function(python = findTQTL()) {
    out <- system2(python,
                   args = c("-c",
                       "import tensorqtl; print(tensorqtl.__version__)"),
                   stdout = TRUE, stderr = FALSE)
    if (length(out) == 0L) NA_character_ else out[[1L]]
}
