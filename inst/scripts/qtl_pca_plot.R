# qtl_pca_plot.R
#
# PCA on the t-statistic matrix from qtlRegressionStats(), with an
# interactive PC1 vs PC2 plot (plotly) where hovering shows the SNP name.
#
# Expects a wide data frame from qtlRegressionStats(t_only=TRUE) as input.
# Run GSDMB_cis_stats.R first to produce GSDMB_cis_tstat.rds, then:
#
#   res <- readRDS("GSDMB_cis_tstat.rds")
#   source("qtl_pca_plot.R")

library(ggplot2)
library(plotly)

# ---- load results ----------------------------------------------------------

if (!exists("res")) {
    rds <- "GSDMB_cis_tstat.rds"
    if (!file.exists(rds))
        stop("Run GSDMB_cis_stats.R first to produce ", rds,
             ", or assign a qtlRegressionStats() result to 'res'.")
    res <- readRDS(rds)
}

# ---- clean and prepare matrix ----------------------------------------------

res_clean <- na.omit(res)
message("Rows after na.omit: ", nrow(res_clean),
        " (removed ", nrow(res) - nrow(res_clean), " with NA genotype t-stat)")

# Numeric t-statistic columns only (drop phenotype_id and variant_id)
t_cols  <- setdiff(names(res_clean), c("phenotype_id", "variant_id"))
t_mat   <- as.matrix(res_clean[, t_cols])
rownames(t_mat) <- res_clean[["variant_id"]]

# ---- PCA -------------------------------------------------------------------

pca     <- prcomp(t_mat, center = TRUE, scale. = TRUE)
pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
scores  <- as.data.frame(pca$x[, 1:2])
scores$variant_id <- res_clean[["variant_id"]]

message("PC1: ", pct_var[1], "%  PC2: ", pct_var[2], "%")

# ---- interactive plot ------------------------------------------------------

p <- ggplot(scores, aes(x = PC1, y = PC2, text = variant_id)) +
    geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
    labs(
        x = paste0("PC1 (", pct_var[1], "%)"),
        y = paste0("PC2 (", pct_var[2], "%)"),
        title = paste0("cis-QTL coefficient PCA - ",
                       unique(res_clean$phenotype_id))
    ) +
    theme_minimal()

ggplotly(p, tooltip = "text")
