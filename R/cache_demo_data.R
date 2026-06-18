#' Cache demo dataset files
#'
#' Downloads the complete tQTLExperiment demo dataset (chr17 tensorQTL results
#' and associated data) from the GitHub release into a
#' `BiocFileCache` subdirectory. Files are cached locally and only
#' re-downloaded when absent.
#'
#' @param release_url Base URL of the GitHub release assets. Defaults to
#'   the `v0.1.0-demodir` release of `vjcitn/tQTLExperiment`.
#' @param bfc A [BiocFileCache::BiocFileCache] object. Defaults to the
#'   user-level cache.
#' @param verbose Logical. Emit progress messages? Default `TRUE`.
#'
#' @return A named list with elements:
#'   - `covariates_tsv`: Path to covariates.tsv
#'   - `mage_rda`: Path to mage_chr17_cis.rda (pre-computed results)
#'   - `pheno_bed`: Path to pheno.bed (phenotype data)
#'   - `parquet`: Path to tqtl_out.cis_qtl_pairs.17.parquet (tensorQTL output)
#'   - `log`: Path to tqtl_out.tensorQTL.cis_nominal.log
#'
#'   All files share the same directory, accessible via
#'   `dirname(result[["covariates_tsv"]])`.
#'
#' @importFrom BiocFileCache BiocFileCache bfccache
#' @importFrom curl curl_download new_handle
#' @importFrom utils download.file
#' @export
#'
#' @examples
#' \dontrun{
#' demo_data <- cache_demo_data()
#' demodir <- dirname(demo_data[["covariates_tsv"]])
#' }
cache_demo_data <- function(
    release_url = paste0("https://github.com/vjcitn/tQTLExperiment/",
                         "releases/download/v0.1.0-demodir"),
    bfc         = BiocFileCache::BiocFileCache(),
    verbose     = TRUE
) {
    fnames <- c(
        "covariates.tsv",
        "mage_chr17_cis.rda",
        "pheno.bed",
        "tqtl_out.cis_qtl_pairs.17.parquet",
        "tqtl_out.tensorQTL.cis_nominal.log"
    )

    cache_dir <- file.path(BiocFileCache::bfccache(bfc), "tqtlExperiment_demodir")
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

    paths <- setNames(file.path(cache_dir, fnames), fnames)

    for (fname in fnames) {
        dest <- paths[[fname]]
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
    }

    if (verbose) message("[ complete    ] All demo files cached in: ", cache_dir)

    # Return with simplified names
    list(
        covariates_tsv = paths[["covariates.tsv"]],
        mage_rda       = paths[["mage_chr17_cis.rda"]],
        pheno_bed      = paths[["pheno.bed"]],
        parquet        = paths[["tqtl_out.cis_qtl_pairs.17.parquet"]],
        log            = paths[["tqtl_out.tensorQTL.cis_nominal.log"]]
    )
}
