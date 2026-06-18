#' Annotate tensorQTL results with gene symbols and biotypes
#'
#' Adds gene symbols (HGNC names) and gene biotypes to a [GenomicRanges::GRanges]
#' object containing tensorQTL results. Useful for annotating the output from
#' [readTQTL()].
#'
#' @param gr A [GenomicRanges::GRanges] object with tensorQTL results, typically
#'   from [readTQTL()].
#' @param phenotype_col Name of the metadata column in `gr` containing Ensembl
#'   gene IDs (ENSG...). Default: `"phenotype_id"`.
#'
#' @return The input `gr` with two additional metadata columns:
#'   - `gene_symbol`: HGNC gene symbol
#'   - `gene_biotype`: Ensembl gene biotype (e.g., "protein_coding", "lincRNA")
#'
#' @details
#' Ensembl version suffixes (e.g., `.3` in `ENSG00000100181.3`) are
#' automatically stripped before lookup in EnsDb.Hsapiens.v79.
#'
#' @export
#'
#' @examples
#' demo_data <- cache_demo_data()
#' demodir = dirname(demo_data[["parquet"]])
#' res <- readTQTL(demodir, mode = "cis_nominal")
#' res_anno <- annotateTensorQTLResults(res$pairs)
#' head(mcols(res_anno)[, c("phenotype_id", "gene_symbol", "gene_biotype")])
annotateTensorQTLResults <- function(gr, phenotype_col = "phenotype_id") {
  # Extract ENSG IDs, stripping version suffixes
  raw_ids <- S4Vectors::mcols(gr)[[phenotype_col]]
  ensg_ids <- sub("\\..*$", "", raw_ids)

  # Query EnsDb for gene_name (symbol) and gene_biotype
  anno_gr <- ensembldb::genes(EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79)
  anno <- as.data.frame(S4Vectors::mcols(anno_gr))

  # Match back to the (possibly non-unique) vector of stripped IDs
  idx <- match(ensg_ids, anno$gene_id)

  S4Vectors::mcols(gr)$gene_symbol  <- anno$gene_name[idx]
  S4Vectors::mcols(gr)$gene_biotype <- anno$gene_biotype[idx]

  gr
}
