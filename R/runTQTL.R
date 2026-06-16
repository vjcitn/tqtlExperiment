#' Prepare input files and return the tensorQTL CLI command
#'
#' Writes the phenotype BED file and (optionally) the covariate file to a
#' user-specified directory, then returns the complete `python -m tensorqtl`
#' command string.  The user runs this command in a terminal where their
#' Python/conda environment is properly configured, then calls [readTQTL()]
#' to load the results back into R.
#'
#' @param x A [tQTLExperiment].
#' @param outDir Path to an existing directory that will hold both the
#'   intermediate input files and the tensorQTL output.  **Must be specified
#'   by the caller** - no default is provided.
#' @param mode One of `"cis_nominal"`, `"cis"`, `"cis_independent"`,
#'   `"trans"`. Passed to `--mode`.
#' @param assayName Name of the assay to use as the phenotype matrix.
#'   Defaults to the first assay.
#' @param mafThreshold MAF threshold. Defaults to `0.05`.
#' @param window Cis-window in base pairs. Defaults to `1000000`.
#' @param permutations Permutations for `"cis"` mode. Defaults to `1000`.
#' @param python Path to the Python executable to embed in the returned
#'   command string. Defaults to `"python3"`.
#' @param ... Additional `--key value` flags passed to tensorQTL verbatim.
#'
#' @return A character string containing the complete shell command to run.
#'   The string is also printed via [message()] for easy copy-paste.
#'
#' @seealso [readTQTL()] to load results after the command has been run.
#'
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
#' )
#' cmd <- prepareTQTL(tqe, outDir = tempdir(), mode = "cis_nominal")
#' cat(cmd, "\n")
#'
#' @export
prepareTQTL <- function(x, outDir,
                        mode      = c("cis_nominal", "cis",
                                      "cis_independent", "trans"),
                        assayName    = NULL,
                        mafThreshold = 0.05,
                        window       = 1000000L,
                        permutations = 1000L,
                        python       = "python3",
                        ...) {
    mode <- match.arg(mode)

    if (missing(outDir) || !nzchar(outDir))
        stop("'outDir' must be specified - it will hold input and output files")
    outDir <- path.expand(outDir)
    if (!dir.exists(outDir))
        stop("'outDir' does not exist: ", outDir)

    if (is.null(assayName))
        assayName <- SummarizedExperiment::assayNames(x)[1]

    pheno_file <- file.path(outDir, "pheno.bed")
    prefix     <- file.path(outDir, "tqtl_out")

    .tqtl_write_pheno_bed(x, assayName, pheno_file)

    cov_cols <- setdiff(colnames(SummarizedExperiment::colData(x)), "fam_index")
    cov_arg  <- character(0)
    if (length(cov_cols) > 0L) {
        cov_file <- file.path(outDir, "covariates.tsv")
        .tqtl_write_covariates(x, cov_cols, cov_file)
        cov_arg <- c("--covariates", shQuote(cov_file))
    }

    extra <- list(...)
    extra_args <- character(0)
    for (nm in names(extra)) {
        val <- extra[[nm]]
        if (isTRUE(val))
            extra_args <- c(extra_args, paste0("--", nm))
        else
            extra_args <- c(extra_args, paste0("--", nm), as.character(val))
    }

    perm_arg <- if (mode == "cis")
        c("--permutations", permutations) else character(0)

    parts <- c(
        python, "-m", "tensorqtl",
        x@plinkPrefix,
        pheno_file,
        prefix,
        "--mode",          mode,
        "--maf_threshold", mafThreshold,
        "--window",        window,
        "-o",              outDir,
        cov_arg,
        perm_arg,
        extra_args
    )

    cmd <- paste(parts, collapse = " ")
    message("Run the following command in a terminal with tensorqtl available:\n\n",
            cmd, "\n\nThen call readTQTL('", outDir, "', mode = '", mode,
            "') to load results into R.")
    invisible(cmd)
}

#' Read tensorQTL results into R
#'
#' Reads the output files written by a tensorQTL run into
#' [GenomicRanges::GRanges] objects, matched to the features in the
#' originating [tQTLExperiment].
#'
#' @param outDir The same directory passed to [prepareTQTL()].
#' @param mode The mode used for the tensorQTL run.
#' @param x The [tQTLExperiment] used to generate the input files, used to
#'   attach feature coordinates to region-level results.
#'
#' @return A named list of [GenomicRanges::GRanges]:
#'   \describe{
#'     \item{`pairs`}{(`cis_nominal`) Feature-variant pairs with association
#'       statistics.}
#'     \item{`hits`}{(`cis`) One range per feature with permutation p-value.}
#'     \item{`output_files`}{(other modes) Paths to the raw output files.}
#'   }
#'
#' @seealso [prepareTQTL()] to write input files and obtain the CLI command.
#'
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
#' )
#' od <- tempdir()
#' cmd <- prepareTQTL(tqe, outDir = od, mode = "cis_nominal")
#' # (user runs cmd in terminal)
#' # res <- readTQTL(od, mode = "cis_nominal", x = tqe)
#'
#' @export
readTQTL <- function(outDir,
                     mode = c("cis_nominal", "cis",
                              "cis_independent", "trans"),
                     x    = NULL) {
    mode <- match.arg(mode)
    outDir <- path.expand(outDir)
    if (!dir.exists(outDir))
        stop("'outDir' not found: ", outDir)

    outfiles <- list.files(outDir, pattern = "tqtl_out",
                           full.names = TRUE)
    if (length(outfiles) == 0L)
        stop("No tensorQTL output files found in '", outDir,
             "'. Has the CLI command been run yet?")

    if (mode == "cis_nominal") {
        parquet_files <- outfiles[grepl("parquet", outfiles)]
        if (length(parquet_files) == 0L)
            stop("No .parquet files found. Has the CLI command been run yet?")
        if (!requireNamespace("arrow", quietly = TRUE))
            stop("Package 'arrow' needed to read parquet output: ",
                 "install.packages('arrow')")
        df <- do.call(rbind, lapply(parquet_files, arrow::read_parquet))

        # Parse GRanges coordinates from variant_id.
        # Handles common formats:
        #   "17:114101:G:A"   (chrom:pos:ref:alt, 4 colon-separated fields)
        #   "22:16849573A-G"  (chrom:posRef-Alt,  2 fields, pos is numeric prefix)
        #   "17_114101_G_A"   (underscore-separated)
        vid   <- df[["variant_id"]]
        chrom <- pos <- NULL

        # try colon-separated
        parts <- strsplit(vid, ":", fixed = TRUE)
        if (all(lengths(parts) >= 2L)) {
            chrom <- vapply(parts, `[[`, character(1), 1L)
            pos2  <- vapply(parts, `[[`, character(1), 2L)
            pos   <- suppressWarnings(as.integer(pos2))
            if (anyNA(pos))   # pos2 may be "16849573A-G" — take numeric prefix
                pos <- suppressWarnings(
                    as.integer(sub("^([0-9]+).*", "\\1", pos2)))
        }

        # try underscore-separated if colon parse failed
        if (is.null(pos) || anyNA(pos)) {
            parts2 <- strsplit(vid, "_", fixed = TRUE)
            if (all(lengths(parts2) >= 2L)) {
                chrom <- vapply(parts2, `[[`, character(1), 1L)
                pos   <- suppressWarnings(
                    as.integer(vapply(parts2, `[[`, character(1), 2L)))
            }
        }

        if (!is.null(pos) && !anyNA(pos)) {
            gr <- GRanges(
                seqnames = chrom,
                ranges   = IRanges::IRanges(start = pos, width = 1L)
            )
            S4Vectors::mcols(gr) <- S4Vectors::DataFrame(df)
            return(list(pairs = gr))
        }

        message("Could not parse chrom/pos from variant_id (e.g. '",
                vid[1], "'); returning DataFrame")
        return(list(pairs = S4Vectors::DataFrame(df)))
    }

    if (mode == "cis") {
        hit_files <- outfiles[grepl("cis_qtl\\.txt", outfiles)]
        if (length(hit_files) == 0L)
            stop("No cis_qtl output found. Has the CLI command been run yet?")
        df <- read.table(hit_files[1L], header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
        if (!is.null(x)) {
            feat_gr <- SummarizedExperiment::rowRanges(x)
            idx     <- match(df[["phenotype_id"]], names(feat_gr))
            gr      <- feat_gr[idx[!is.na(idx)]]
            S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
                df[!is.na(idx), , drop = FALSE])
        } else {
            gr <- GRanges()
            S4Vectors::mcols(gr) <- S4Vectors::DataFrame(df)
        }
        return(list(hits = gr))
    }

    list(output_files = outfiles)
}

# ---- internal file writers ----------------------------------------------

.tqtl_write_pheno_bed <- function(x, assayName, path) {
    rr  <- SummarizedExperiment::rowRanges(x)
    mat <- SummarizedExperiment::assay(x, assayName)

    chr_raw <- sub("^chr", "", as.character(GenomicRanges::seqnames(rr)))
    chr_raw[chr_raw == "X"]            <- "23"
    chr_raw[chr_raw == "Y"]            <- "24"
    chr_raw[chr_raw == "XY"]           <- "25"
    chr_raw[chr_raw %in% c("M","MT")] <- "26"
    chr <- suppressWarnings(as.integer(chr_raw))
    if (anyNA(chr)) {
        bad  <- unique(as.character(GenomicRanges::seqnames(rr))[is.na(chr)])
        warning("Dropping ", sum(is.na(chr)), " features with unmappable ",
                "seqnames: ", paste(bad, collapse = ", "))
        keep <- !is.na(chr)
        rr   <- rr[keep]; mat <- mat[keep, , drop = FALSE]; chr <- chr[keep]
    }

    df <- data.frame(
        `#chr`       = chr,
        start        = GenomicRanges::start(rr) - 1L,
        end          = GenomicRanges::end(rr),
        phenotype_id = names(rr),
        check.names  = FALSE, stringsAsFactors = FALSE
    )
    out <- cbind(df, mat[, colnames(x), drop = FALSE])
    write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

.tqtl_write_covariates <- function(x, cov_cols, path) {
    cd <- as.data.frame(SummarizedExperiment::colData(x))[, cov_cols,
                                                           drop = FALSE]
    non_num <- names(cd)[!vapply(cd, is.numeric, logical(1))]
    if (length(non_num))
        stop("Non-numeric covariate columns: ",
             paste(non_num, collapse = ", "),
             ". Use tQTLExperimentFromRSE() with a covariateMatrix.")
    # tensorQTL format: covariates x samples
    out <- cbind(data.frame(ID = colnames(cd), stringsAsFactors = FALSE),
                 as.data.frame(t(cd)))
    write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
