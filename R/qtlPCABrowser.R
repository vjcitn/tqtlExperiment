#' Interactive PCA browser for cis-QTL t-statistics
#'
#' Collapses factor-expanded t-statistic columns via [collapseFactorTstats()],
#' computes PCA on the result, and launches a Shiny application with a
#' selectable PC-vs-PC scatter plot (plotly) and a pairs plot of the first five
#' PCs.  Hovering over a SNP-phenotype point in the scatter plot displays a
#' beeswarm of phenotype expression by genotype for that pair, with an optional
#' colour-by selector for sample-level variables (e.g. sex, batch).
#'
#' @param res A wide data frame from [qtlRegressionStats()] with
#'   \code{t_only = TRUE}.  Must contain \code{phenotype_id} and
#'   \code{variant_id} columns.
#' @param tqe A [tQTLExperiment] with genotype data.
#' @param assayName Name of the assay to use. Defaults to the first assay.
#' @param collapse_patterns Named list of regex patterns passed to
#'   [collapseFactorTstats()].  Defaults handle MAGE-style batch and
#'   population columns.
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
qtlPCABrowser <- function(res, tqe, assayName = NULL,
                          collapse_patterns = list(
                              batch      = "^t\\.factor_batch_",
                              population = "^t\\.population"
                          )) {
    if (!requireNamespace("plotly", quietly = TRUE))
        stop("Package 'plotly' is required.")

    if (is.null(assayName))
        assayName <- SummarizedExperiment::assayNames(tqe)[1L]

    # ---- Pre-extract data for beeswarm rendering ---------------------------
    row_gr      <- SummarizedExperiment::rowRanges(tqe, use.names = TRUE)
    pheno_names <- names(row_gr)
    var_gr      <- tqtlVariantRanges(tqe)
    var_names   <- S4Vectors::mcols(var_gr)[["snp_id"]]
    bed         <- tqtlGeno(tqe)
    assay_mat   <- SummarizedExperiment::assay(tqe, assayName)

    # colData - reconstruct original categoricals from dummy columns
    cov_cols <- setdiff(colnames(SummarizedExperiment::colData(tqe)), "fam_index")
    cd       <- as.data.frame(SummarizedExperiment::colData(tqe))[, cov_cols, drop = FALSE]

    # Reverse dummy-coding: for each group of indicator columns sharing a prefix,
    # reconstruct the original factor (ref level = all indicators 0).
    dummy_groups <- list(
        batch      = grep("^factor_batch_",  names(cd), value = TRUE),
        population = grep("^population",     names(cd), value = TRUE)
    )
    for (grp in names(dummy_groups)) {
        cols <- dummy_groups[[grp]]
        if (length(cols) < 2L) next
        mat    <- as.matrix(cd[, cols, drop = FALSE])
        # strip prefix to get level labels; ref level gets label "ref"
        pfx    <- if (grp == "batch") "factor_batch_" else "population"
        labels <- sub(pfx, "", cols)
        level_vec <- apply(mat, 1L, function(x) {
            hit <- which(x > 0.5)
            if (length(hit) == 0L) "ref" else labels[hit[1L]]
        })
        cd[[grp]] <- factor(level_vec)
    }
    # colour choices: reconstructed categoricals + remaining numeric covariates
    color_choices <- c("none",
                       intersect(names(dummy_groups), names(cd)),
                       setdiff(cov_cols, unlist(dummy_groups)))

    # ---- Collapse factor columns then PCA ----------------------------------
    res_clean <- stats::na.omit(res)
    res_coll  <- collapseFactorTstats(res_clean, patterns = collapse_patterns)
    t_cols    <- setdiff(names(res_coll), c("phenotype_id", "variant_id"))
    t_mat     <- as.matrix(res_coll[, t_cols])

    pca     <- stats::prcomp(t_mat, center = TRUE, scale. = FALSE)
    pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

    pc_choices <- colnames(pca$x)
    if (is.null(pc_choices))
        pc_choices <- paste0("PC", seq_len(ncol(pca$x)))

    scores <- data.frame(pca$x, check.names = FALSE)
    scores[["variant_id"]] <- res_clean[["variant_id"]]
    scores[["phenotype_id"]] <- res_clean[["phenotype_id"]]
    scores[["row_id"]] <- seq_len(nrow(scores))
    scores[["hover_text"]] <- paste0(
        "SNP: ", scores[["variant_id"]],
        "<br>Phenotype: ", scores[["phenotype_id"]]
    )
    pc_axis_title <- function(pc) {
        idx <- match(pc, pc_choices)
        paste0(pc, " (", pct_var[idx], "%)")
    }
    pair_cols   <- head(pc_choices, 5L)
    pair_labels <- pc_axis_title(pair_cols)

    # ---- Shiny UI ----------------------------------------------------------
    ui <- fluidPage(
        theme = bslib::bs_theme(bootswatch = "flatly"),
        titlePanel("cis-QTL t-stat PCA"),
        fluidRow(
            column(7,
                   tabsetPanel(
                       tabPanel(
                           "PCA scatter",
                           plotly::plotlyOutput("pca", height = "550px")
                       ),
                       tabPanel(
                           "PC pairs",
                           plotOutput("pc_pairs", height = "550px")
                       )
                   )),
            column(5,
                   fluidRow(
                       column(6,
                              selectInput("x_pc", "X-axis PC:",
                                          choices = pc_choices,
                                          selected = pc_choices[1L])),
                       column(6,
                              selectInput("y_pc", "Y-axis PC:",
                                          choices = pc_choices,
                                          selected = pc_choices[min(2L, length(pc_choices))]))
                   ),
                   selectInput("color_var", "Color beeswarm by:",
                               choices  = color_choices,
                               selected = if ("sexXY" %in% color_choices) "sexXY" else "none"),
                   plotOutput("beeswarm", height = "480px"))
        )
    )

    # ---- Shiny server ------------------------------------------------------
    server <- function(input, output, session) {

        output$pca <- plotly::renderPlotly({
            req(input$x_pc, input$y_pc)
            x_pc <- input$x_pc
            y_pc <- input$y_pc
            plotly::plot_ly(
                scores,
                x         = scores[[x_pc]],
                y         = scores[[y_pc]],
                key       = scores[["row_id"]],
                type      = "scatter",
                mode      = "markers",
                text      = scores[["hover_text"]],
                hoverinfo = "text",
                source    = "pca_scores",
                marker    = list(size = 5, color = "steelblue", opacity = 0.6)
            ) |>
                plotly::layout(
                    xaxis = list(title = pc_axis_title(x_pc)),
                    yaxis = list(title = pc_axis_title(y_pc)),
                    hoverlabel = list(bgcolor = "white")
                )
        })

        output$pc_pairs <- renderPlot({
            if (length(pair_cols) < 2L) {
                plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
                text(0, 0, "At least two PCs are required for a pairs plot",
                     cex = 1.1, col = "grey50")
                return(invisible(NULL))
            }
            graphics::pairs(
                scores[, pair_cols, drop = FALSE],
                labels = pair_labels,
                pch    = 16,
                cex    = 0.5,
                col    = grDevices::adjustcolor("steelblue", alpha.f = 0.5),
                main   = "First five PC score pairs"
            )
        })

        output$beeswarm <- renderPlot({
            hover <- plotly::event_data("plotly_hover", source = "pca_scores")
            if (is.null(hover)) {
                plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
                text(0, 0, "Hover over a point\nto see genotype effect",
                     cex = 1.2, col = "grey50")
                return(invisible(NULL))
            }

            row_idx <- suppressWarnings(as.integer(hover$key))
            if (is.na(row_idx) && !is.null(hover$pointNumber))
                row_idx <- hover$pointNumber + 1L
            if (is.na(row_idx) || row_idx < 1L || row_idx > nrow(scores)) {
                plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
                text(0, 0, "Selected PCA point was not found",
                     cex = 0.9, col = "grey50")
                return(invisible(NULL))
            }

            vid <- scores[["variant_id"]][row_idx]
            pid <- scores[["phenotype_id"]][row_idx]

            var_idx   <- match(vid, var_names)
            pheno_idx <- match(pid, pheno_names)

            if (is.na(var_idx) || is.na(pheno_idx)) {
                plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
                text(0, 0, paste("not found:", vid, pid),
                     cex = 0.9, col = "grey50")
                return(invisible(NULL))
            }

            geno  <- as.integer(bed[, var_idx])
            pheno <- as.numeric(assay_mat[pheno_idx, ])

            df <- data.frame(
                genotype  = factor(geno, levels = 0:2, labels = c("0", "1", "2")),
                phenotype = pheno,
                stringsAsFactors = FALSE
            )

            color_var <- input$color_var
            if (color_var != "none" && color_var %in% names(cd)) {
                df[["color_by"]] <- as.factor(cd[[color_var]])
                p <- ggplot2::ggplot(df, ggplot2::aes(
                        x = .data$genotype, y = .data$phenotype,
                        color = .data$color_by)) +
                    ggbeeswarm::geom_beeswarm(size = 2, alpha = 0.7) +
                    ggplot2::geom_boxplot(alpha = 0.2, width = 0.3,
                                         color = "grey40", outlier.shape = NA) +
                    ggplot2::labs(color = color_var)
            } else {
                p <- ggplot2::ggplot(df, ggplot2::aes(
                        x = .data$genotype, y = .data$phenotype)) +
                    ggbeeswarm::geom_beeswarm(size = 2, color = "steelblue",
                                              alpha = 0.7) +
                    ggplot2::geom_boxplot(alpha = 0.2, width = 0.3,
                                         color = "grey40", outlier.shape = NA)
            }

            p + ggplot2::theme_minimal() +
                ggplot2::xlab("Genotype (# alt alleles)") +
                ggplot2::ylab(paste0("Phenotype value (", pid, ")")) +
                ggplot2::ggtitle(paste0("SNP: ", vid, "\nPhenotype: ", pid))
        })
    }

    shiny::shinyApp(ui, server)
}
