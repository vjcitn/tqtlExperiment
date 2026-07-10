library(tQTLExperiment)

library(BEDMatrix)

pl = cache_mage_chr17_plink()
mageSEfilt = get_cached_mage17SE()

bb = BEDMatrix(pl[["bed"]])
plpre = file.path(dirname(pl[["bed"]]), gsub(".bed", "", basename(pl[["bed"]])))

#plpre = "/Users/vincentcarey/Library/Caches/org.R-project.R/R/BiocFileCache/CCDG_mage_chr17_plink/CCDG_mage_chr17"

  # Build covariate matrix
  cd <- as.data.frame(colData(mageSEfilt))
  mm <- model.matrix(~ batch + population + sex, data = cd)
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "int"

  # Construct (reuse plpre from before)
  tqe <- tQTLExperimentFromRSE(mageSEfilt, plpre, covariateMatrix = mm)

  # Subset to chr17 features
  is_17 <- as.character(seqnames(rowRanges(tqe))) == "chr17"
  tqe17 <- tqe[which(is_17),]  # [1:50], ]

  # Prepare files and get command
  outdir <- "~/mage_tqtl_runJulybig2"
  dir.create(outdir, showWarnings = FALSE)
  cmd <- prepareTQTL(tqe17, outDir = outdir, mode = "cis_nominal",
                     mafThreshold = 0.05)

