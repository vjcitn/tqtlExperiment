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
    cat("use prepareTQTL() to write inputs and get CLI command\n")
    invisible(object)
})
