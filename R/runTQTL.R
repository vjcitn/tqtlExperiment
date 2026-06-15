#' Run the tensorQTL eQTL mapping CLI
#'
#' Writes temporary input files from a [tQTLExperiment], invokes
#' `python -m tensorqtl`, reads the results, and returns them as a named
#' list of [GenomicRanges::GRanges] objects.
#'
#' @param x A [tQTLExperiment].
#' @param python Path to the Python executable with tensorqtl installed.
#'   Defaults to [findTQTL()].
#' @param mode One of `"cis"`, `"cis_nominal"`, `"cis_independent"`,
#'   `"trans"`. Passed to `--mode`. Defaults to `"cis_nominal"`.
#' @param assayName Name of the assay to use as phenotype. Defaults to
#'   first assay.
#' @param mafThreshold MAF threshold for variant filtering. Passed to
#'   `--maf_threshold`. Defaults to `0.05`.
#' @param window Cis-window size in base pairs. Defaults to `1000000`.
#' @param permutations Number of permutations for `"cis"` mode. Ignored for
#'   `"cis_nominal"`.
#' @param debug Logical; if `TRUE`, the temporary directory is kept after
#'   the run and its path is printed along with the first lines of each
#'   input file.
#' @param ... Additional arguments passed to the tensorQTL CLI as
#'   `--key value` pairs.
#'
#' @return A named list of [GenomicRanges::GRanges]:
#'   \describe{
#'     \item{`pairs`}{For `cis_nominal`: all tested feature–variant pairs.}
#'     \item{`hits`}{For `cis`: one range per feature with permutation
#'       p-value.}
#'   }
#'   The exact contents depend on mode; raw output files are preserved when
#'   `debug = TRUE`.
#'
#' @examples
#' exdir <- system.file("extdata", package = "tQTLExperiment")
#' tqe <- tQTLExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed")
#' )
#' tryCatch({
#'     res <- runTQTL(tqe, mode = "cis_nominal")
#'     res$pairs
#' }, error = function(e) invisible(NULL))
#'
#' @export
setGeneric("runTQTL",
    function(x, python = NULL,
             mode = c("cis_nominal", "cis", "cis_independent", "trans"),
             assayName = NULL, mafThreshold = 0.05, window = 1000000L,
             permutations = 1000L, debug = FALSE, ...)
        standardGeneric("runTQTL")
)

#' @rdname runTQTL
setMethod("runTQTL", "tQTLExperiment",
    function(x, python = NULL,
             mode = c("cis_nominal", "cis", "cis_independent", "trans"),
             assayName = NULL, mafThreshold = 0.05, window = 1000000L,
             permutations = 1000L, debug = FALSE, ...) {

        mode <- match.arg(mode)
        if (is.null(python)) python <- findTQTL()
        if (is.null(assayName))
            assayName <- SummarizedExperiment::assayNames(x)[1]

        tmpdir <- tempfile("tqtl_")
        dir.create(tmpdir)
        if (debug) {
            message("debug=TRUE: intermediate files kept in ", tmpdir)
        } else {
            on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
        }

        pheno_file <- file.path(tmpdir, "pheno.bed")
        prefix     <- file.path(tmpdir, "tqtl_out")

        .tqtl_write_pheno_bed(x, assayName, pheno_file)

        if (debug) {
            message("--- pheno.bed (first 2 lines) ---")
            message(paste(readLines(pheno_file, n = 2L), collapse = "\n"))
        }

        args <- c("-m", "tensorqtl",
                  x@plinkPrefix,
                  pheno_file,
                  prefix,
                  "--mode", mode,
                  "--maf_threshold", mafThreshold,
                  "--window", window,
                  "-o", tmpdir)

        # add covariates if any (excluding fam_index)
        cov_cols <- setdiff(colnames(SummarizedExperiment::colData(x)),
                            "fam_index")
        if (length(cov_cols) > 0L) {
            cov_file <- file.path(tmpdir, "covariates.tsv")
            .tqtl_write_covariates(x, cov_cols, cov_file)
            args <- c(args, "--covariates", cov_file)
            if (debug) {
                message("--- covariates.tsv (first 2 lines) ---")
                message(paste(readLines(cov_file, n = 2L), collapse = "\n"))
            }
        }

        if (mode == "cis")
            args <- c(args, "--permutations", permutations)

        extra <- list(...)
        for (nm in names(extra)) {
            val <- extra[[nm]]
            if (isTRUE(val))
                args <- c(args, paste0("--", nm))
            else
                args <- c(args, paste0("--", nm), as.character(val))
        }

        # clear KMP/OMP thread-limit vars so tensorQTL's torch sees all cores
        omp_limit_vars <- c("OMP_THREAD_LIMIT", "KMP_DEVICE_THREAD_LIMIT",
                            "KMP_TEAMS_THREAD_LIMIT", "KMP_ALL_THREADS")
        saved_env <- Sys.getenv(omp_limit_vars, names = TRUE,
                                unset = NA_character_)
        vars_to_unset <- names(saved_env)[!is.na(saved_env)]
        if (length(vars_to_unset)) {
            Sys.unsetenv(vars_to_unset)
            on.exit(
                do.call(Sys.setenv, as.list(saved_env[vars_to_unset])),
                add = TRUE
            )
        }

        status <- system2(python, args = args)
        if (status != 0L)
            stop("tensorqtl exited with status ", status)

        .tqtl_read_results(tmpdir, basename(prefix), mode, x)
    }
)

# ---- internal helpers ---------------------------------------------------

.tqtl_write_pheno_bed <- function(x, assayName, path) {
    rr  <- SummarizedExperiment::rowRanges(x)
    mat <- SummarizedExperiment::assay(x, assayName)

    # tensorQTL expects seqnames without 'chr' prefix as integers
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
    # tensorQTL format: covariates x samples (transpose of colData)
    out <- cbind(data.frame(ID = rownames(t(cd)), stringsAsFactors = FALSE),
                 as.data.frame(t(cd)))
    write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

.tqtl_read_results <- function(tmpdir, prefix, mode, x) {
    outfiles <- list.files(tmpdir, pattern = prefix, full.names = TRUE)

    if (mode == "cis_nominal") {
        parquet_files <- outfiles[grepl("\\.parquet$", outfiles)]
        if (length(parquet_files) == 0L)
            stop("No parquet output files found in ", tmpdir,
                 ". Is the 'arrow' package installed?")
        if (!requireNamespace("arrow", quietly = TRUE))
            stop("Package 'arrow' is required to read tensorQTL parquet output. ",
                 "Install with: install.packages('arrow')")
        dfs <- lapply(parquet_files, arrow::read_parquet)
        df  <- do.call(rbind, dfs)
        gr  <- GRanges(
            seqnames = df[["variant_id"]],   # placeholder; parse if needed
            ranges   = IRanges::IRanges(start = df[["pos"]], width = 1L)
        )
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(df)
        return(list(pairs = gr))
    }

    if (mode == "cis") {
        hit_files <- outfiles[grepl("cis_qtl\\.txt", outfiles)]
        if (length(hit_files) == 0L)
            stop("No cis_qtl output file found in ", tmpdir)
        df  <- read.table(hit_files[1L], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
        feat_gr <- SummarizedExperiment::rowRanges(x)
        idx     <- match(df[["phenotype_id"]], names(feat_gr))
        gr      <- feat_gr[idx[!is.na(idx)]]
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(df[!is.na(idx), ,
                                                        drop = FALSE])
        return(list(hits = gr))
    }

    # other modes: return raw file paths for user inspection
    list(output_files = outfiles)
}
