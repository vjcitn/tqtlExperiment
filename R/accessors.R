#' Access the lazy genotype matrix
#'
#' @param x A [tQTLExperiment].
#' @param value A [BEDMatrix::BEDMatrix] (samples × variants).
#' @return The `geno` slot (samples × variants BEDMatrix).
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
#' )
#' tqtlGeno(tqe)
#' @export
setGeneric("tqtlGeno", function(x) standardGeneric("tqtlGeno"))

#' @rdname tqtlGeno
#' @export
setGeneric("tqtlGeno<-", function(x, value) standardGeneric("tqtlGeno<-"))

#' @rdname tqtlGeno
setMethod("tqtlGeno", "tQTLExperiment", function(x) x@geno)

#' @rdname tqtlGeno
setReplaceMethod("tqtlGeno", "tQTLExperiment", function(x, value) {
    x@geno <- value; validObject(x); x
})

#' Access variant genomic ranges
#'
#' @param x A [tQTLExperiment].
#' @param value A [GenomicRanges::GRanges].
#' @return A `GRanges` with metadata columns `snp_id`, `ref`, `alt`.
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
#' )
#' tqtlVariantRanges(tqe)
#' @export
setGeneric("tqtlVariantRanges", function(x) standardGeneric("tqtlVariantRanges"))

#' @rdname tqtlVariantRanges
#' @export
setGeneric("tqtlVariantRanges<-", function(x, value) standardGeneric("tqtlVariantRanges<-"))

#' @rdname tqtlVariantRanges
setMethod("tqtlVariantRanges", "tQTLExperiment", function(x) x@variantRanges)

#' @rdname tqtlVariantRanges
setReplaceMethod("tqtlVariantRanges", "tQTLExperiment", function(x, value) {
    x@variantRanges <- value; validObject(x); x
})

#' Access the PLINK file prefix
#'
#' @param x A [tQTLExperiment].
#' @param value A length-one character string.
#' @return The PLINK prefix path.
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
#' )
#' tqtlPlinkPrefix(tqe)
#' @export
setGeneric("tqtlPlinkPrefix", function(x) standardGeneric("tqtlPlinkPrefix"))

#' @rdname tqtlPlinkPrefix
#' @export
setGeneric("tqtlPlinkPrefix<-", function(x, value) standardGeneric("tqtlPlinkPrefix<-"))

#' @rdname tqtlPlinkPrefix
setMethod("tqtlPlinkPrefix", "tQTLExperiment", function(x) x@plinkPrefix)

#' @rdname tqtlPlinkPrefix
setReplaceMethod("tqtlPlinkPrefix", "tQTLExperiment", function(x, value) {
    x@plinkPrefix <- value; validObject(x); x
})
