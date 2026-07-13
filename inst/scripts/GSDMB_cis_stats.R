# GSDMB_cis_stats.R
#
# Retrieve full cis-regression t-statistics for the gene GSDMB using the
# MAGE chromosome 17 dataset. Requires internet access on first run to
# cache the SE and PLINK files via BiocFileCache.

library(tQTLExperiment)
library(SummarizedExperiment)

# ---- Build tQTLExperiment --------------------------------------------------

se          <- get_cached_mage17SE(verbose = FALSE)
plink_paths <- cache_mage_chr17_plink(verbose = FALSE)
plpre       <- tools::file_path_sans_ext(plink_paths[["bed"]])

cd <- as.data.frame(colData(se))
mm <- stats::model.matrix(~ factor(batch) + population + sex, data = cd)
mm <- mm[, -1, drop = FALSE]

tqe <- tQTLExperimentFromRSE(
    se              = se,
    plinkPrefix     = plpre,
    covariateMatrix = mm,
    genome          = "hg38"
)
tqe <- addGeneSymbols(tqe)
message("tQTLExperiment ready: ", nrow(tqe), " features x ", ncol(tqe), " samples")

# ---- Regression t-statistics for GSDMB ------------------------------------

t0  <- proc.time()
res <- qtlRegressionStats(tqe, symbol = "GSDMB")
message("elapsed: ", round((proc.time() - t0)[["elapsed"]], 1), " s")

message("dim: ", paste(dim(res), collapse = " x "))
message("columns: ", paste(names(res), collapse = ", "))
print(head(res[, 1:6]))
