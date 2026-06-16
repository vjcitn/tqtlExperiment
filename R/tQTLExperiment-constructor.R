#' Construct a tQTLExperiment from input files
#'
#' Reads phenotype, covariate, and variant metadata files into a
#' [tQTLExperiment] object. Genotype data is represented lazily via
#' [BEDMatrix::BEDMatrix].
#'
#' @param plinkPrefix Path prefix for the PLINK file set (`.bed`, `.bim`,
#'   `.fam`).
#' @param phenoFile Path to the phenotype BED file (tab-separated, columns
#'   `#chr`, `start`, `end`, `phenotype_id`, then one column per sample;
#'   0-based half-open coordinates).
#' @param covFile Path to the covariate file. tensorQTL expects covariates
#'   as rows with the first column being the covariate name and remaining
#'   columns being sample values (covariates × samples). If `NULL`, no
#'   covariate file is used.
#' @param assayName Name to assign the phenotype matrix in `assays()`.
#'   Defaults to `"pheno"`.
#'
#' @return A [tQTLExperiment] object.
#'
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100-tqtl.tsv")
#' )
#' tqe
#'
#' @export
tQTLExperiment <- function(plinkPrefix, phenoFile, covFile = NULL,
                           assayName = "pheno") {
    # ---- phenotype BED --------------------------------------------------
    pheno_raw <- read.table(phenoFile, header = TRUE, sep = "\t",
                            check.names = FALSE, comment.char = "")
    colnames(pheno_raw)[1] <- "chr"

    sample_ids  <- colnames(pheno_raw)[-(1:4)]
    feature_ids <- pheno_raw[["phenotype_id"]]

    pheno_mat <- as.matrix(pheno_raw[, -(1:4), drop = FALSE])
    rownames(pheno_mat) <- feature_ids
    mode(pheno_mat) <- "double"

    # rowRanges: BED is 0-based half-open → 1-based closed
    row_gr <- GRanges(
        seqnames = pheno_raw[["chr"]],
        ranges   = IRanges::IRanges(
            start = pheno_raw[["start"]] + 1L,
            end   = pheno_raw[["end"]]
        ),
        phenotype_id = feature_ids
    )
    names(row_gr) <- feature_ids

    # ---- covariates (optional) ------------------------------------------
    # tensorQTL format: covariates x samples (rows = covariates)
    cov_df <- S4Vectors::DataFrame(row.names = sample_ids)
    if (!is.null(covFile)) {
        cov_raw <- read.table(covFile, header = TRUE, sep = "\t",
                              check.names = FALSE, row.names = 1)
        # cov_raw is covariates x samples; transpose to samples x covariates
        # for colData storage
        cov_t <- as.data.frame(t(cov_raw[, sample_ids, drop = FALSE]))
        cov_df <- S4Vectors::DataFrame(cov_t, row.names = sample_ids)
    }

    # ---- PLINK / BEDMatrix ----------------------------------------------
    bim_file <- paste0(plinkPrefix, ".bim")
    fam_file <- paste0(plinkPrefix, ".fam")
    for (f in c(paste0(plinkPrefix, ".bed"), bim_file, fam_file))
        if (!file.exists(f))
            stop("PLINK file not found: ", f)

    geno_mat <- BEDMatrix::BEDMatrix(plinkPrefix)

    fam     <- read.table(fam_file, header = FALSE,
                          col.names = c("fid", "iid", "pat", "mat", "sex", "phen"))
    fam_idx <- .match_fam_ids(sample_ids, fam[["iid"]])
    missing_geno <- sample_ids[is.na(fam_idx)]
    if (length(missing_geno))
        warning("Samples in phenotype file not found in .fam file: ",
                paste(missing_geno, collapse = ", "))
    cov_df[["fam_index"]] <- fam_idx

    # ---- variant ranges (from .bim) -------------------------------------
    bim <- read.table(bim_file, header = FALSE,
                      col.names = c("chrom", "snp_id", "cm", "pos", "alt", "ref"),
                      colClasses = c("character", "character", "numeric",
                                     "integer",   "character", "character"))
    var_gr <- GRanges(
        seqnames = bim[["chrom"]],
        ranges   = IRanges::IRanges(start = bim[["pos"]], width = 1L),
        snp_id   = bim[["snp_id"]],
        ref      = bim[["ref"]],
        alt      = bim[["alt"]]
    )
    names(var_gr) <- bim[["snp_id"]]

    # ---- assemble -------------------------------------------------------
    se <- SummarizedExperiment::SummarizedExperiment(
        assays    = S4Vectors::SimpleList(setNames(list(pheno_mat), assayName)),
        rowRanges = row_gr,
        colData   = cov_df
    )

    new("tQTLExperiment",
        se,
        geno          = geno_mat,
        variantRanges = var_gr,
        plinkPrefix   = as.character(plinkPrefix)
    )
}

#' Construct a tQTLExperiment from a RangedSummarizedExperiment
#'
#' Converts an existing [SummarizedExperiment::RangedSummarizedExperiment]
#' into a [tQTLExperiment] by attaching PLINK genotype data.
#'
#' @param se A [SummarizedExperiment::RangedSummarizedExperiment].
#' @param plinkPrefix Path prefix for the PLINK file set.
#' @param covariateMatrix A numeric matrix (rows = samples, cols =
#'   covariates) built e.g. with [stats::model.matrix()] on
#'   `as.data.frame(colData(se))`. See also the
#'   [ExploreModelMatrix](https://bioconductor.org/packages/ExploreModelMatrix)
#'   package. Rename `(Intercept)` to `int` if including an intercept.
#' @param assayName Name of the assay in `se` to use. Defaults to first.
#' @param featureIdColumn Optional mcols column name to use as
#'   `phenotype_id`; defaults to `rownames(se)`.
#'
#' @return A [tQTLExperiment] object.
#'
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe_ref <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
#' )
#' se <- as(tqe_ref, "RangedSummarizedExperiment")
#' cd <- as.data.frame(colData(tqe_ref))
#' mm <- model.matrix(~ 1, data = cd)
#' colnames(mm)[1] <- "int"
#' tqe2 <- tQTLExperimentFromRSE(se,
#'     plinkPrefix     = file.path(exdir, "chr22-n100"),
#'     covariateMatrix = mm)
#' tqe2
#'
#' @export
tQTLExperimentFromRSE <- function(se, plinkPrefix,
                                   covariateMatrix  = NULL,
                                   assayName        = NULL,
                                   featureIdColumn  = NULL) {
    if (!is(se, "RangedSummarizedExperiment"))
        stop("'se' must be a RangedSummarizedExperiment")

    if (is.null(assayName))
        assayName <- SummarizedExperiment::assayNames(se)[1]

    pheno_mat <- as.matrix(SummarizedExperiment::assay(se, assayName))
    mode(pheno_mat) <- "double"

    sample_ids  <- colnames(se)
    feature_ids <- if (!is.null(featureIdColumn))
        mcols(SummarizedExperiment::rowRanges(se))[[featureIdColumn]]
    else
        rownames(se)
    if (is.null(feature_ids) || anyNA(feature_ids))
        stop("Could not determine feature IDs: set rownames(se) or supply 'featureIdColumn'")
    rownames(pheno_mat) <- feature_ids

    rr <- SummarizedExperiment::rowRanges(se)
    names(rr) <- feature_ids
    if (!"phenotype_id" %in% names(S4Vectors::mcols(rr)))
        rr$phenotype_id <- feature_ids

    # covariates
    if (!is.null(covariateMatrix)) {
        if (!is.matrix(covariateMatrix) || !is.numeric(covariateMatrix))
            stop("'covariateMatrix' must be a numeric matrix")
        if (nrow(covariateMatrix) != ncol(se))
            stop("nrow(covariateMatrix) must equal ncol(se)")
        cov_df <- S4Vectors::DataFrame(covariateMatrix, row.names = sample_ids)
    } else {
        cov_df <- S4Vectors::DataFrame(row.names = sample_ids)
    }

    bim_file <- paste0(plinkPrefix, ".bim")
    fam_file <- paste0(plinkPrefix, ".fam")
    for (f in c(paste0(plinkPrefix, ".bed"), bim_file, fam_file))
        if (!file.exists(f)) stop("PLINK file not found: ", f)

    geno_mat <- BEDMatrix::BEDMatrix(plinkPrefix)

    fam     <- read.table(fam_file, header = FALSE,
                          col.names = c("fid", "iid", "pat", "mat", "sex", "phen"))
    fam_idx <- .match_fam_ids(sample_ids, fam[["iid"]])
    if (anyNA(fam_idx))
        warning("Samples in se not found in .fam file: ",
                paste(sample_ids[is.na(fam_idx)], collapse = ", "))
    cov_df[["fam_index"]] <- fam_idx

    bim <- read.table(bim_file, header = FALSE,
                      col.names = c("chrom", "snp_id", "cm", "pos", "alt", "ref"),
                      colClasses = c("character", "character", "numeric",
                                     "integer",   "character", "character"))
    var_gr <- GRanges(
        seqnames = bim[["chrom"]],
        ranges   = IRanges::IRanges(start = bim[["pos"]], width = 1L),
        snp_id   = bim[["snp_id"]],
        ref      = bim[["ref"]],
        alt      = bim[["alt"]]
    )
    names(var_gr) <- bim[["snp_id"]]

    se_new <- SummarizedExperiment::SummarizedExperiment(
        assays    = S4Vectors::SimpleList(setNames(list(pheno_mat), assayName)),
        rowRanges = rr,
        colData   = cov_df
    )

    new("tQTLExperiment",
        se_new,
        geno          = geno_mat,
        variantRanges = var_gr,
        plinkPrefix   = as.character(plinkPrefix)
    )
}
