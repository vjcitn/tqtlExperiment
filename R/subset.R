#' Subset a tQTLExperiment
#'
#' Sample subsetting propagates to the lazy `geno` matrix (BEDMatrix rows).
#'
#' @param x A [tQTLExperiment].
#' @param i Row (feature) index.
#' @param j Column (sample) index.
#' @param ... Ignored.
#' @param drop Ignored.
#' @return A [tQTLExperiment].
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
#' )
#' tqe[1:5, 1:10]
#' @export
setMethod("[", "tQTLExperiment", function(x, i, j, ..., drop = FALSE) {
    se_sub <- if (missing(i) && missing(j))
        callNextMethod()
    else if (missing(i))
        callNextMethod(x, , j, ..., drop = drop)
    else if (missing(j))
        callNextMethod(x, i, , ..., drop = drop)
    else
        callNextMethod(x, i, j, ..., drop = drop)

    j_idx <- if (missing(j)) seq_len(ncol(x)) else j
    se_sub@geno <- x@geno[j_idx, , drop = FALSE]
    se_sub
})
