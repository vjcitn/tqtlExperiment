#' Cache pre-computed MAGE chr17 tensorQTL results
#'
#' Downloads the pre-computed `cis_nominal` tensorQTL results for MAGE
#' chromosome 17 from the tQTLExperiment GitHub release into a
#' `BiocFileCache` subdirectory. Results are cached as a serialized R object
#' and loaded into memory as a [GenomicRanges::GRanges] object.
#'
#' @param release_url Base URL of the GitHub release assets. Defaults to
#'   the `v0.1.0-demo-data` release of `vjcitn/tQTLExperiment`.
#' @param bfc A [BiocFileCache::BiocFileCache] object. Defaults to the
#'   user-level cache.
#' @param verbose Logical. Emit progress messages? Default `TRUE`.
#'
#' @return A [GenomicRanges::GRanges] object containing 71,000+ feature-variant
#'   pairs from tensorQTL `cis_nominal` results. Metadata columns include
#'   `phenotype_id`, `variant_id`, `tss_distance`, `ma_samples`, `ma_count`,
#'   `maf`, `pval_nominal`, and `slope`.
#'
#' @importFrom BiocFileCache BiocFileCache bfccache
#' @importFrom curl curl_download new_handle
#' @importFrom utils download.file
#' @export
#'
#' @examples
#' \dontrun{
#' gr <- cache_mage_chr17_cis()
#' head(gr)
#' }
cache_mage_chr17_cis <- function(
    release_url = paste0("https://github.com/vjcitn/tQTLExperiment/",
                         "releases/download/v0.1.0-demo-data"),
    bfc         = BiocFileCache::BiocFileCache(),
    verbose     = TRUE
) {
    fname <- "magec17cis.rda"
    cache_dir <- file.path(BiocFileCache::bfccache(bfc),
                           "mage_chr17_cis")
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

    dest <- file.path(cache_dir, fname)

    if (file.exists(dest)) {
        if (verbose) message("[ cache hit   ] ", fname)
    } else {
        url <- paste0(release_url, "/", fname)
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
                        " with ", length(obj), " records")
    obj
}
