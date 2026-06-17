#' Plot phenotype values grouped by SNP genotype
#'
#' Creates a beeswarm plot comparing a phenotype across genotype groups
#' (0, 1, 2 copies of the alternate allele) for a given SNP.  Useful for
#' visualizing the effect of a significant eQTL variant.
#'
#' @param x A [tQTLExperiment].
#' @param snp_id Character. SNP identifier (row name in `variantRanges(x)`).
#' @param phenotype_id Character. Phenotype/gene identifier (row name in
#'   `rowRanges(x)`).
#' @param assayName Name of the assay to use. Defaults to the first assay.
#' @param size Point size for beeswarm. Default 2.
#' @param color Color for points. Default "steelblue".
#' @param title Main plot title. If `NULL`, auto-generated from SNP and gene.
#'
#' @return A ggplot2 plot object.
#'
#' @import rlang
#' @import ggplot2 ggbeeswarm
#'
#' @examples
#' if (!requireNamespace("CSHLvc2026")) BiocManager::install("vjcitn/CSHLvc2026")
#' data(mageSEfilt, package="CSHLvc2026")
#' plink_paths <- cache_mage_chr17_plink()
#' plpre <- tools::file_path_sans_ext(plink_paths[["bed"]])
#' cd <- as.data.frame(colData(mageSEfilt))
#' mm <- model.matrix(~ batch + population + sex, data = cd)
#' mm <- mm[, -1, drop = FALSE]   # remove (Intercept)
#' tqe <- tQTLExperimentFromRSE(
#'     se              = mageSEfilt,
#'     plinkPrefix     = plpre,
#'     covariateMatrix = mm,
#'     genome          = "hg38"
#' )
#' plotGenotypeEffect(tqe, "17:410391:T:C", "ENSG00000187624.9")
#' @export
plotGenotypeEffect <- function(x, snp_id, phenotype_id,
                               assayName = NULL,
                               size = 2,
                               color = "steelblue",
                               title = NULL) {

    if (is.null(assayName))
        assayName <- SummarizedExperiment::assayNames(x)[1L]

    # Get variant row (search by snp_id metadata)
    var_gr <- tqtlVariantRanges(x)
    if (!("snp_id" %in% names(S4Vectors::mcols(var_gr))))
        stop("variantRanges missing snp_id metadata")
    var_idx <- which(S4Vectors::mcols(var_gr)[["snp_id"]] == snp_id)
    if (length(var_idx) == 0L)
        stop("SNP not found: ", snp_id)

    # Get phenotype row
    row_gr <- SummarizedExperiment::rowRanges(x, use.names = TRUE)
    if (!phenotype_id %in% names(row_gr))
        stop("Phenotype not found: ", phenotype_id)
    pheno_idx <- which(names(row_gr) == phenotype_id)

    # Extract genotypes for this SNP (from BEDMatrix)
    bed <- tqtlGeno(x)
    geno <- as.integer(bed[, var_idx])   # samples x variants; extract column

    # Extract phenotype values
    pheno <- as.numeric(SummarizedExperiment::assay(x, assayName)[pheno_idx, ])

    # Build data frame
    df <- data.frame(
        genotype = factor(geno, levels = 0:2,
                         labels = c("0", "1", "2")),
        phenotype = pheno,
        stringsAsFactors = FALSE
    )

    # Remove NAs
    df <- df[!is.na(df$genotype) & !is.na(df$phenotype), ]

    if (is.null(title)) {
        snp_info <- var_gr[var_idx]
        gene_name <- if ("gene_name" %in% names(S4Vectors::mcols(row_gr)))
            S4Vectors::mcols(row_gr)[["gene_name"]][pheno_idx] else phenotype_id
        title <- paste0(gene_name, " (", snp_id, ")")
    }

    # Create plot
    ggplot2::ggplot(df, ggplot2::aes(x = .data$genotype, y = .data$phenotype)) +
        ggbeeswarm::geom_beeswarm(size = size, color = color, alpha = 0.7) +
        ggplot2::geom_boxplot(alpha = 0.3, width = 0.3) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Genotype (# alternate alleles)") +
        ggplot2::ylab("Phenotype value") +
        ggplot2::ggtitle(title)
}
