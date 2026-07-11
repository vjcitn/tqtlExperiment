#' Linear regression statistics for a collection of SNP-gene pairs
#'
#' For each (phenotype_id, variant_id) pair, fits a linear model of phenotype
#' on genotype plus all covariates stored in the [tQTLExperiment] and returns
#' a tidy data frame of coefficients, standard errors, t-statistics, and
#' p-values for every model term.  Pairs are grouped by phenotype and processed
#' in parallel via [BiocParallel::bplapply()].
#'
#' @param tqe A [tQTLExperiment] with genotype data.
#' @param pairs A data frame with columns \code{phenotype_id} and
#'   \code{variant_id} identifying the SNP-gene pairs to analyse.
#' @param assayName Name of the assay to use. Defaults to the first assay.
#' @param t_only Logical. If \code{TRUE} (default), return only the
#'   t-statistic (estimate / std.error) for each model term. In wide format
#'   columns are named \code{<term>.t}. If \code{FALSE}, return estimates and
#'   standard errors (and additionally statistic and p.value in long format).
#' @param wideFormat Logical. If \code{TRUE} (default), pivot to one row per
#'   pair. With \code{t_only = TRUE}: one column per term named
#'   \code{<term>.t}. With \code{t_only = FALSE}: two columns per term,
#'   \code{estimate.<term>} and \code{std.error.<term>}. If \code{FALSE},
#'   return the long tidy form with one row per pair x term.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling
#'   parallelisation. Defaults to [BiocParallel::SerialParam()] (single
#'   core). Use e.g. [BiocParallel::MulticoreParam()] for multicore.
#'
#' @return A data frame. See \code{t_only} and \code{wideFormat} for the
#'   shape of the output.
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
#' qtlRegressionStats(tqe, pairs, t_only = FALSE, wideFormat = FALSE)
#' @importFrom stats lm reshape
#' @importFrom BiocParallel bplapply SerialParam
#' @export
qtlRegressionStats <- function(tqe, pairs, assayName = NULL,
                               t_only = TRUE, wideFormat = TRUE,
                               BPPARAM = BiocParallel::SerialParam()) {
    stopifnot(is.data.frame(pairs))
    if (!all(c("phenotype_id", "variant_id") %in% colnames(pairs)))
        stop("'pairs' must have columns 'phenotype_id' and 'variant_id'")

    if (is.null(assayName))
        assayName <- SummarizedExperiment::assayNames(tqe)[1L]

    # Build lookup tables once
    row_gr      <- SummarizedExperiment::rowRanges(tqe, use.names = TRUE)
    pheno_names <- names(row_gr)

    var_gr    <- tqtlVariantRanges(tqe)
    var_names <- S4Vectors::mcols(var_gr)[["snp_id"]]

    cov_cols <- setdiff(colnames(SummarizedExperiment::colData(tqe)), "fam_index")
    cov_data <- as.data.frame(SummarizedExperiment::colData(tqe))[, cov_cols, drop = FALSE]

    bed       <- tqtlGeno(tqe)
    assay_mat <- SummarizedExperiment::assay(tqe, assayName)

    # Group by phenotype so each phenotype vector is extracted only once
    groups <- split(seq_len(nrow(pairs)), pairs$phenotype_id)

    process_group <- function(idx) {
        pid       <- pairs$phenotype_id[idx[1L]]
        pheno_idx <- match(pid, pheno_names)

        if (is.na(pheno_idx)) {
            warning("phenotype not found - skipped: ", pid)
            return(NULL)
        }

        pheno <- as.numeric(assay_mat[pheno_idx, ])

        lapply(idx, function(i) {
            vid     <- pairs$variant_id[i]
            var_idx <- match(vid, var_names)

            if (is.na(var_idx)) {
                warning("variant not found - skipped: ", vid)
                return(NULL)
            }

            geno <- as.integer(bed[, var_idx])

            reg_data <- data.frame(
                phenotype = pheno,
                genotype  = geno,
                cov_data,
                stringsAsFactors = FALSE
            )

            cf <- summary(lm(phenotype ~ ., data = reg_data))$coefficients

            if (t_only) {
                tval <- cf[, "Estimate"] / cf[, "Std. Error"]
                data.frame(
                    phenotype_id = pid,
                    variant_id   = vid,
                    term         = rownames(cf),
                    t            = tval,
                    row.names    = NULL,
                    stringsAsFactors = FALSE
                )
            } else {
                data.frame(
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
        })
    }

    group_results <- BiocParallel::bplapply(groups, process_group,
                                            BPPARAM = BPPARAM)

    long <- do.call(rbind, unlist(group_results, recursive = FALSE))

    # Restore original row order
    long <- long[order(match(paste(long$phenotype_id, long$variant_id),
                             paste(pairs$phenotype_id, pairs$variant_id))), ]
    rownames(long) <- NULL

    if (!wideFormat)
        return(long)

    long$term <- gsub("[^A-Za-z0-9_]", "_", long$term)

    if (t_only) {
        wide <- reshape(long,
                        idvar     = c("phenotype_id", "variant_id"),
                        timevar   = "term",
                        v.names   = "t",
                        direction = "wide")
    } else {
        wide <- reshape(long,
                        idvar     = c("phenotype_id", "variant_id"),
                        timevar   = "term",
                        v.names   = c("estimate", "std.error"),
                        direction = "wide",
                        drop      = c("statistic", "p.value"))
    }
    rownames(wide) <- NULL
    wide
}
