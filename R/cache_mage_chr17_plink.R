#' Validate PLINK .bed magic bytes
#'
#' @param path Character. Local path to the .bed file.
#' @return Invisibly returns `TRUE`; stops on failure.
#' @keywords internal
.validate_bed_magic <- function(path) {
    con <- file(path, open = "rb")
    on.exit(close(con))
    magic    <- readBin(con, what = "raw", n = 3L)
    expected <- as.raw(c(0x6c, 0x1b, 0x01))
    if (!identical(magic, expected))
        stop(
            "PLINK .bed magic byte check failed for: ", path, "\n",
            "  Expected: ", paste(expected, collapse = " "), "\n",
            "  Got:      ", paste(magic,    collapse = " "), "\n",
            "The file may be corrupt or a PLINK 1 (pre-1.9) binary."
        )
    invisible(TRUE)
}

#' Cache MAGE chr17 PLINK files from GitHub releases
#'
#' Downloads the three PLINK files (`CCDG_mage_chr17.fam`, `.bim`, `.bed`)
#' from the `vjcitn/CSHLvc2026` GitHub release into a `BiocFileCache`
#' subdirectory.  Files are only re-downloaded when absent from the cache.
#' No authentication is required.
#'
#' @param release_url Base URL of the GitHub release assets.  Defaults to
#'   the `v0.1.0-data` release of `vjcitn/CSHLvc2026`.
#' @param bfc A [BiocFileCache::BiocFileCache] object.  Defaults to the
#'   user-level cache.
#' @param verbose Logical.  Emit progress messages?  Default `TRUE`.
#' @param validate_bed Logical.  Check `.bed` magic bytes after download?
#'   Default `TRUE`.
#'
#' @return A named list with elements `fam`, `bim`, and `bed`, each a
#'   character scalar giving the local path to the cached file.  All three
#'   share the same directory and stem, so
#'   `tools::file_path_sans_ext(result[["bed"]])` gives the PLINK prefix
#'   for use with [tQTLExperimentFromRSE()].
#'
#' @importFrom BiocFileCache BiocFileCache bfccache
#' @importFrom curl curl_download new_handle
#' @importFrom utils download.file
#' @export
#'
#' @examples
#' \dontrun{
#' plink <- cache_mage_chr17_plink()
#' plpre <- tools::file_path_sans_ext(plink[["bed"]])
#' }
cache_mage_chr17_plink <- function(
    release_url  = paste0("https://github.com/vjcitn/CSHLvc2026/",
                          "releases/download/v0.1.0-data"),
    bfc          = BiocFileCache::BiocFileCache(),
    verbose      = TRUE,
    validate_bed = TRUE
) {
    stem  <- "CCDG_mage_chr17"
    exts  <- c("fam", "bim", "bed")

    plink_dir <- file.path(BiocFileCache::bfccache(bfc),
                           "CCDG_mage_chr17_plink")
    dir.create(plink_dir, showWarnings = FALSE, recursive = TRUE)

    paths <- setNames(
        file.path(plink_dir, paste0(stem, ".", exts)),
        exts
    )

    for (ext in exts) {
        dest <- paths[[ext]]
        if (file.exists(dest)) {
            if (verbose) message("[ cache hit   ] ", basename(dest))
        } else {
            url <- paste0(release_url, "/", basename(dest))
            if (verbose) message("[ downloading ] ", basename(dest))
            tryCatch({
                curl::curl_download(url, destfile = dest,
                                    handle = curl::new_handle(followlocation = TRUE),
                                    quiet = !verbose)
            }, error = function(e) {
                if (verbose) message("  (curl unavailable; using standard download)")
                utils::download.file(url, destfile = dest, mode = "wb",
                                     quiet = !verbose)
            })
        }
    }

    if (validate_bed) {
        if (verbose) message("[ validate    ] Checking .bed magic bytes ...")
        .validate_bed_magic(paths[["bed"]])
        if (verbose) message("[ validate    ] .bed magic bytes OK.")
    }

    paths
}
