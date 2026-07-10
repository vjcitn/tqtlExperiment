#' Linear regression statistics for a collection of SNP-gene pairs
#'
#' For each (phenotype_id, variant_id) pair, fits a linear model of phenotype
#' on genotype plus all covariates stored in the [tQTLExperiment] and returns
#' a tidy data frame of coefficients, standard errors, t-statistics, and
#' p-values for every model term.
#'
#' @param tqe A [tQTLExperiment] with genotype data.
#' @param pairs A data frame with columns \code{phenotype_id} and
#'   \code{variant_id} identifying the SNP-gene pairs to analyse.
#' @param assayName Name of the assay to use. Defaults to the first assay.
#' @param wideFormat Logical. If \code{TRUE} (default), pivot to one row per
#'   pair with two columns per model term: \code{estimate.<term>} and
#'   \code{std.error.<term>} (e.g. \code{estimate.genotype},
#'   \code{std.error.genotype}), suitable for multivariate analysis. If
#'   \code{FALSE}, return the long tidy form with one row per pair x term
#'   including \code{statistic} and \code{p.value} as well.
#'
#' @return A data frame. In wide format (default): one row per pair, two
#'   columns per model term (\code{estimate.<term>} and
#'   \code{std.error.<term>}). In long format: one row per pair x model term
#'   with columns \code{phenotype_id}, \code{variant_id}, \code{term},
#'   \code{estimate}, \code{std.error}, \code{statistic}, \code{p.value}.
#'
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100-tqtl.tsv"),
#'     genome      = "hg38"
#' )
#' pairs <- data.frame(
#'     phenotype_id = rownames(tqe)[1],
#'     variant_id   = S4Vectors::mcols(tqtlVariantRanges(tqe))[["snp_id"]][1]
#' )
#' qtlRegressionStats(tqe, pairs)
#' qtlRegressionStats(tqe, pairs, wideFormat = FALSE)
#' @importFrom stats lm reshape
#' @export
qtlRegressionStats <- function(tqe, pairs, assayName = NULL, wideFormat = TRUE) {
    stopifnot(is.data.frame(pairs))
    if (!all(c("phenotype_id", "variant_id") %in% colnames(pairs)))
        stop("'pairs' must have columns 'phenotype_id' and 'variant_id'")

    if (is.null(assayName))
        assayName <- SummarizedExperiment::assayNames(tqe)[1L]

    # Build lookup tables once
    row_gr     <- SummarizedExperiment::rowRanges(tqe, use.names = TRUE)
    pheno_names <- names(row_gr)

    var_gr    <- tqtlVariantRanges(tqe)
    var_names <- S4Vectors::mcols(var_gr)[["snp_id"]]

    cov_cols <- setdiff(colnames(SummarizedExperiment::colData(tqe)), "fam_index")
    cov_data <- as.data.frame(SummarizedExperiment::colData(tqe))[, cov_cols, drop = FALSE]

    bed       <- tqtlGeno(tqe)
    assay_mat <- SummarizedExperiment::assay(tqe, assayName)

    results <- vector("list", nrow(pairs))

    for (i in seq_len(nrow(pairs))) {
        pid <- pairs$phenotype_id[i]
        vid <- pairs$variant_id[i]

        pheno_idx <- match(pid, pheno_names)
        var_idx   <- match(vid, var_names)

        if (is.na(pheno_idx) || is.na(var_idx)) {
            warning("pair not found - skipped: ", pid, " / ", vid)
            next
        }

        pheno <- as.numeric(assay_mat[pheno_idx, ])
        geno  <- as.integer(bed[, var_idx])

        reg_data <- data.frame(
            phenotype = pheno,
            genotype  = geno,
            cov_data,
            stringsAsFactors = FALSE
        )

        cf <- summary(lm(phenotype ~ ., data = reg_data))$coefficients

        results[[i]] <- data.frame(
            phenotype_id = pid,
            variant_id   = vid,
            term         = rownames(cf),
            estimate     = cf[, "Estimate"],
            std.error    = cf[, "Std. Error"],
            statistic    = cf[, "t value"],
            p.value      = cf[, "Pr(>|t|)"],
            row.names    = NULL,
            stringsAsFactors = FALSE
        )
    }

    long <- do.call(rbind, results)

    if (!wideFormat)
        return(long)

    # Sanitise term names so they make valid column name components
    long$term <- gsub("[^A-Za-z0-9_]", "_", long$term)

    wide <- reshape(long,
                    idvar     = c("phenotype_id", "variant_id"),
                    timevar   = "term",
                    v.names   = c("estimate", "std.error"),
                    direction = "wide",
                    drop      = c("statistic", "p.value"))
    rownames(wide) <- NULL
    wide
}
