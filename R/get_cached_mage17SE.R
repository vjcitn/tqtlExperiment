#' Cache and load the MAGE chr17 SummarizedExperiment
#'
#' Downloads `mageSEfilt.rda` from the Open Storage Network bucket into a
#' `BiocFileCache` subdirectory and returns the contained
#' [SummarizedExperiment::SummarizedExperiment] object.  The file is only
#' downloaded once; subsequent calls return the cached copy.
#'
#' @param url Full URL to `mageSEfilt.rda`.  Defaults to the OSN BiocMAGE17SE
#'   bucket.
#' @param bfc A [BiocFileCache::BiocFileCache] object.  Defaults to the
#'   user-level cache.
#' @param verbose Logical.  Emit progress messages?  Default `TRUE`.
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] object (`mageSEfilt`).
#'
#' @importFrom BiocFileCache BiocFileCache bfccache
#' @importFrom curl curl_download new_handle
#' @importFrom utils download.file
#' @export
#'
#' @examples
#' \dontrun{
#' se <- get_cached_mage17SE()
#' se
#' }
get_cached_mage17SE <- function(
    url     = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocMAGE17SE/mageSEfilt.rda",
    bfc     = BiocFileCache::BiocFileCache(),
    verbose = TRUE
) {
    fname     <- "mageSEfilt.rda"
    cache_dir <- file.path(BiocFileCache::bfccache(bfc), "mage17SE")
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

    dest <- file.path(cache_dir, fname)

    if (file.exists(dest)) {
        if (verbose) message("[ cache hit   ] ", fname)
    } else {
        if (verbose) message("[ downloading ] ", fname)
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

    if (verbose) message("[ loading     ] Restoring ", fname)
    env <- new.env()
    load(dest, envir = env)
    obj <- get(ls(env)[1L], envir = env)

    if (verbose) message("[ loaded      ] ", class(obj)[1L],
                         " with ", nrow(obj), " features x ", ncol(obj), " samples")
    obj
}
