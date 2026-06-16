#' @import methods
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @importFrom S4Vectors DataFrame SimpleList mcols
#' @importFrom Matrix Matrix forceSymmetric
#' @importFrom GenomeInfoDb genome
#' @importFrom AnnotationDbi select
#' @importFrom utils head read.table write.table
#' @importFrom stats setNames model.matrix
NULL

#' tQTLExperiment: a SummarizedExperiment for tensorQTL eQTL mapping
#'
#' `tQTLExperiment` extends [SummarizedExperiment::RangedSummarizedExperiment]
#' to hold all data needed for the tensorQTL eQTL mapping tool:
#'
#' - **assays**: phenotype matrix (features × samples)
#' - **rowRanges**: genomic coordinates for each phenotype feature (gene)
#' - **colData**: per-sample metadata including covariates
#' - **geno**: lazy [BEDMatrix::BEDMatrix] over the PLINK `.bed` file
#'   (samples × variants); data are not loaded until subscripted
#' - **variantRanges**: [GenomicRanges::GRanges] for each variant, with
#'   metadata columns `snp_id`, `ref`, `alt`
#' - **plinkPrefix**: path prefix for the PLINK file set
#'
#' tensorQTL is invoked via its CLI (`python -m tensorqtl`) and supports
#' GPU acceleration through CUDA or Apple MPS (Metal Performance Shaders)
#' on Apple Silicon.
#'
#' @slot geno A [BEDMatrix::BEDMatrix] (samples × variants).
#' @slot variantRanges A [GenomicRanges::GRanges] of length equal to
#'   `ncol(geno)`.
#' @slot plinkPrefix A length-one character string giving the PLINK file
#'   prefix.
#'
#' @seealso [tQTLExperiment()] for the constructor,
#'   [prepareTQTL()] to write input files and get the CLI command,
#'   [readTQTL()] to load results, [findTQTL()] to locate Python.
#'
#' @exportClass tQTLExperiment
setClass("tQTLExperiment",
    contains = "RangedSummarizedExperiment",
    representation(
        geno          = "ANY",
        variantRanges = "GRanges",
        plinkPrefix   = "character"
    )
)

setValidity("tQTLExperiment", function(object) {
    msg <- character(0)

    # BEDMatrix is samples x variants
    ns_geno <- nrow(object@geno)
    ns_se   <- ncol(object)
    if (ns_geno != ns_se)
        msg <- c(msg, sprintf(
            "nrow(geno) (%d) must equal ncol(experiment) (%d)",
            ns_geno, ns_se))

    nv_geno <- ncol(object@geno)
    nv      <- length(object@variantRanges)
    if (nv_geno != nv)
        msg <- c(msg, sprintf(
            "ncol(geno) (%d) must equal length(variantRanges) (%d)",
            nv_geno, nv))

    if (length(object@plinkPrefix) != 1L || is.na(object@plinkPrefix))
        msg <- c(msg, "'plinkPrefix' must be a single non-NA character string")

    if (length(msg)) msg else TRUE
})
