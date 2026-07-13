#' Collapse factor-expanded t-statistic columns into group summaries
#'
#' When a factor variable (e.g. \code{batch}) or a multi-level character
#' variable (e.g. \code{population}) is included in the regression model, the
#' wide output of [qtlRegressionStats()] contains one \code{t.<var><level>}
#' column per level.  This function identifies groups of such columns by
#' regular-expression patterns, replaces each group with the root-mean-square
#' (RMS) of its t-statistics, and returns a leaner data frame suitable for
#' PCA.
#'
#' @param wide_df A wide data frame from [qtlRegressionStats()] with
#'   \code{t_only = TRUE}.
#' @param patterns A named list of regular expressions.  Each element matches
#'   a set of columns to collapse; the name becomes the new column name
#'   (prefixed with \code{"t."}).  Defaults to patterns that detect the MAGE
#'   chr17 dataset's batch and population indicator columns.
#'
#' @return \code{wide_df} with matched column groups replaced by single
#'   \code{t.<name>} RMS columns.
#'
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe   <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100-tqtl.tsv"),
#'     genome      = "hg38"
#' )
#' tqe <- addGeneSymbols(tqe)
#' res <- qtlRegressionStats(tqe, symbol = "TPTEP1")
#' collapsed <- collapseFactorTstats(res)
#' @export
collapseFactorTstats <- function(
        wide_df,
        patterns = list(
            batch      = "^t\\.factor_batch_",
            population = "^t\\.population"
        )
) {
    result <- wide_df
    all_cols <- names(result)

    for (grp_name in names(patterns)) {
        pat  <- patterns[[grp_name]]
        cols <- grep(pat, all_cols, value = TRUE)
        if (length(cols) < 2L) next

        mat <- as.matrix(result[, cols, drop = FALSE])
        rms <- sqrt(rowMeans(mat^2, na.rm = TRUE))

        new_col <- paste0("t.", grp_name)
        result[[new_col]] <- rms
        result[, cols] <- NULL
    }

    result
}
