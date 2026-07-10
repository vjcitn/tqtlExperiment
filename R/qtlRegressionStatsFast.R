#' Fast linear regression statistics for a collection of SNP-gene pairs
#'
#' Optimised version of [qtlRegressionStats()] that uses QR rank-1 updates
#' from the \pkg{fastQR} package.  For each phenotype group the covariate
#' design matrix is factorised once via \code{qr_lm_init()} (O(np\eqn{^2}));
#' each SNP within the group then adds its genotype column to the existing QR
#' via \code{qr_lm_update()} (O(np)) rather than re-factorising from scratch.
#' Phenotype groups are processed in parallel via [BiocParallel::bplapply()].
#'
#' \pkg{fastQR} is not on CRAN or Bioconductor; install it from the local
#' source or GitHub before calling this function.
#'
#' @param tqe A [tQTLExperiment] with genotype data.
#' @param pairs A data frame with columns \code{phenotype_id} and
#'   \code{variant_id} identifying the SNP-gene pairs to analyse.
#' @param assayName Name of the assay to use. Defaults to the first assay.
#' @param wideFormat Logical. If \code{TRUE} (default), pivot to one row per
#'   pair with two columns per model term (\code{estimate.<term>} and
#'   \code{std.error.<term>}) suitable for multivariate analysis. If
#'   \code{FALSE}, return the long tidy form with one row per pair x term
#'   including \code{statistic} and \code{p.value}.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling
#'   parallelisation. Defaults to [BiocParallel::SerialParam()] (single
#'   core). Use e.g. [BiocParallel::MulticoreParam()] for multicore.
#'
#' @return Same structure as [qtlRegressionStats()]: wide format (default)
#'   with two columns per model term, or long format with one row per pair x
#'   term.
#'
#' @importFrom stats pt reshape
#' @importFrom BiocParallel bplapply SerialParam
#' @export
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
#'     phenotype_id = rownames(tqe)[1:3],
#'     variant_id   = S4Vectors::mcols(tqtlVariantRanges(tqe))[["snp_id"]][1:3]
#' )
#' if (requireNamespace("fastQR", quietly = TRUE))
#'     qtlRegressionStatsFast(tqe, pairs)
qtlRegressionStatsFast <- function(tqe, pairs, assayName = NULL,
                                   wideFormat = TRUE,
                                   BPPARAM = BiocParallel::SerialParam()) {
    if (!requireNamespace("fastQR", quietly = TRUE))
        stop("Package 'fastQR' is required. Install it from the local source or GitHub.")

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

    # Remove constant columns (e.g. an explicit intercept column from tensorQTL
    # covariates). lm() drops these silently via QR pivoting; qr_lm_init does not.
    const_cols <- vapply(cov_data,
                         function(x) length(unique(x[!is.na(x)])) < 2L,
                         logical(1L))
    if (any(const_cols))
        cov_data <- cov_data[, !const_cols, drop = FALSE]

    bed       <- tqtlGeno(tqe)
    assay_mat <- SummarizedExperiment::assay(tqe, assayName)

    # Split pairs into per-phenotype groups
    groups <- split(seq_len(nrow(pairs)), pairs$phenotype_id)

    process_group <- function(idx) {
        pid       <- pairs$phenotype_id[idx[1L]]
        pheno_idx <- match(pid, pheno_names)

        if (is.na(pheno_idx)) {
            warning("phenotype not found - skipped: ", pid)
            return(NULL)
        }

        pheno     <- as.numeric(assay_mat[pheno_idx, ])
        base_data <- data.frame(phenotype = pheno, cov_data, stringsAsFactors = FALSE)

        # Full QR factorisation of covariate design matrix — done once per phenotype
        base_fit <- fastQR::qr_lm_init(phenotype ~ ., data = base_data)

        lapply(idx, function(i) {
            vid     <- pairs$variant_id[i]
            var_idx <- match(vid, var_names)

            if (is.na(var_idx)) {
                warning("variant not found - skipped: ", vid)
                return(NULL)
            }

            geno <- as.integer(bed[, var_idx])

            # O(np) column-add via Givens rotations — no re-factorisation
            fit <- fastQR::qr_lm_update(base_fit,
                                         add      = matrix(geno, ncol = 1L),
                                         add_name = "genotype")

            cf   <- fit$coef
            se   <- fit$coef.se
            tval <- cf / se
            pval <- 2 * stats::pt(abs(tval), df = fit$df, lower.tail = FALSE)

            data.frame(
                phenotype_id = pid,
                variant_id   = vid,
                term         = names(cf),
                estimate     = cf,
                std.error    = se,
                statistic    = tval,
                p.value      = pval,
                row.names    = NULL,
                stringsAsFactors = FALSE
            )
        })
    }

    group_results <- BiocParallel::bplapply(groups, process_group, BPPARAM = BPPARAM)

    long <- do.call(rbind, unlist(group_results, recursive = FALSE))

    # Restore original row order
    long <- long[order(match(paste(long$phenotype_id, long$variant_id),
                             paste(pairs$phenotype_id, pairs$variant_id))), ]
    rownames(long) <- NULL

    if (!wideFormat)
        return(long)

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
