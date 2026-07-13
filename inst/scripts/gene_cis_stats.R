# gene_cis_stats.R
#
# Demonstrates qtlRegressionStats() with the symbol= interface:
# retrieve full regression t-statistics for all cis-variants of one gene.

library(tQTLExperiment)

exdir <- system.file("extdata", package = "tQTLExperiment")
tqe   <- tQTLExperiment(
    plinkPrefix = file.path(exdir, "chr22-n100"),
    phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
    covFile     = file.path(exdir, "cov-n100-tqtl.tsv"),
    genome      = "hg38"
)

# Gene symbols must be present - add them if not already done.
tqe <- addGeneSymbols(tqe)

# Show available symbols in this small demo dataset.
syms <- S4Vectors::mcols(SummarizedExperiment::rowRanges(tqe))[["gene_name"]]
cat("Available gene symbols:", paste(head(syms[!is.na(syms)], 10), collapse = ", "), "\n")

# Collect t-statistics for all cis-variants of one gene (default window = 1 Mb).
res <- qtlRegressionStats(tqe, symbol = "TPTEP1")

cat("Dimensions:", dim(res), "\n")
cat("Columns:   ", names(res), "\n")
print(head(res[, 1:6]))
