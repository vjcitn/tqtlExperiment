#' Interactive PCA browser for cis-QTL t-statistics
#'
#' Computes a principal component analysis on the t-statistic matrix returned
#' by [qtlRegressionStats()] (with \code{t_only = TRUE}) and launches a Shiny
#' application showing a PC1 vs PC2 scatter plot (plotly).  Hovering over a
#' point displays a beeswarm plot of phenotype expression by genotype for that
#' SNP, rendered via [plotGenotypeEffect()].
#'
#' @param res A wide data frame from [qtlRegressionStats()] with
#'   \code{t_only = TRUE}.  Must contain \code{phenotype_id} and
#'   \code{variant_id} columns.
#' @param tqe A [tQTLExperiment] with genotype data, used to render the
#'   beeswarm plots on hover.
#' @param assayName Name of the assay to use in [plotGenotypeEffect()].
#'   Defaults to the first assay.
#'
#' @return A Shiny application object (invisibly).
#'
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @importFrom plotly plot_ly layout event_data plotlyOutput renderPlotly
#' @importFrom bslib bs_theme
#' @export
#'
#' @examples
#' if (interactive()) {
#'   exdir <- system.file("extdata", package = "tQTLExperiment")
#'   tqe   <- tQTLExperiment(
#'       plinkPrefix = file.path(exdir, "chr22-n100"),
#'       phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'       covFile     = file.path(exdir, "cov-n100-tqtl.tsv"),
#'       genome      = "hg38"
#'   )
#'   tqe <- addGeneSymbols(tqe)
#'   res <- qtlRegressionStats(tqe, symbol = "TPTEP1")
#'   qtlPCABrowser(res, tqe)
#' }
qtlPCABrowser <- function(res, tqe, assayName = NULL) {
    if (!requireNamespace("plotly", quietly = TRUE))
        stop("Package 'plotly' is required.")

    # ---- PCA ---------------------------------------------------------------
    res_clean <- stats::na.omit(res)
    t_cols    <- setdiff(names(res_clean), c("phenotype_id", "variant_id"))
    t_mat     <- as.matrix(res_clean[, t_cols])

    pca     <- stats::prcomp(t_mat, center = TRUE, scale. = FALSE)
    pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

    scores <- data.frame(
        PC1        = pca$x[, 1],
        PC2        = pca$x[, 2],
        variant_id = res_clean[["variant_id"]],
        stringsAsFactors = FALSE
    )
    phenotype_id <- res_clean[["phenotype_id"]][1L]

    # ---- Shiny UI ----------------------------------------------------------
    ui <- fluidPage(
        theme = bslib::bs_theme(bootswatch = "flatly"),
        titlePanel("cis-QTL coefficient PCA"),
        fluidRow(
            column(7, plotly::plotlyOutput("pca", height = "550px")),
            column(5, plotOutput("beeswarm", height = "550px"))
        ),
        fluidRow(
            column(12,
                   p(style = "color:grey; font-size:0.85em;",
                     "Hover over a point to see the genotype effect plot."))
        )
    )

    # ---- Shiny server ------------------------------------------------------
    server <- function(input, output, session) {

        output$pca <- plotly::renderPlotly({
            plotly::plot_ly(
                scores,
                x       = ~PC1,
                y       = ~PC2,
                key     = ~variant_id,
                type    = "scatter",
                mode    = "markers",
                text    = ~variant_id,
                hoverinfo = "text",
                marker  = list(size = 5, color = "steelblue", opacity = 0.6)
            ) |>
                plotly::layout(
                    xaxis = list(title = paste0("PC1 (", pct_var[1], "%)")),
                    yaxis = list(title = paste0("PC2 (", pct_var[2], "%)")),
                    hoverlabel = list(bgcolor = "white")
                )
        })

        output$beeswarm <- renderPlot({
            hover <- plotly::event_data("plotly_hover")
            if (is.null(hover)) {
                plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
                text(0, 0, "Hover over a point\nto see genotype effect",
                     cex = 1.2, col = "grey50")
                return(invisible(NULL))
            }
            vid <- hover$key
            tryCatch(
                plotGenotypeEffect(tqe, snp_id = vid,
                                   phenotype_id = phenotype_id,
                                   assayName    = assayName,
                                   title        = vid),
                error = function(e) {
                    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
                    text(0, 0, paste("No plot:\n", e$message),
                         cex = 0.9, col = "grey50")
                }
            )
        })
    }

    shiny::shinyApp(ui, server)
}
