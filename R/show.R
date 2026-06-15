#' Show a tQTLExperiment
#' @param object A [tQTLExperiment].
#' @importFrom methods show
setMethod("show", "tQTLExperiment", function(object) {
    cat("class: tQTLExperiment\n")
    cat("features:", nrow(object), " samples:", ncol(object), "\n")
    cat("assays(", length(SummarizedExperiment::assayNames(object)), "):",
        paste(SummarizedExperiment::assayNames(object), collapse = ", "), "\n")
    cat("rowRanges: GRanges with", nrow(object), "features\n")
    cov_cols <- setdiff(colnames(SummarizedExperiment::colData(object)),
                        "fam_index")
    cat("colData(", length(cov_cols), ") covariates:",
        paste(head(cov_cols, 4), collapse = ", "),
        if (length(cov_cols) > 4) "..." else "", "\n")
    cat("geno:", nrow(object@geno), "samples x", ncol(object@geno),  # direct slot access
        "variants [BEDMatrix - lazy]\n")
    cat("plinkPrefix:", object@plinkPrefix, "\n")
    py <- tryCatch(findTQTL(), error = function(e) NA_character_)
    ver <- if (!is.na(py)) tryCatch(tensorqtlVersion(py),
                                    error = function(e) "?") else "not found"
    cat("tensorqtl:", ver, "\n")
    invisible(object)
})
