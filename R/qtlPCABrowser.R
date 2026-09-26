#' Interactive PCA browser for cis-QTL t-statistics
#'
#' Collapses factor-expanded t-statistic columns via [collapseFactorTstats()],
#' computes PCA on the result, and launches a Shiny application with a
#' selectable PC-vs-PC scatter plot (plotly) and a pairs plot of the first five
#' PCs.  Hovering over a SNP-gene point in the scatter plot displays a
#' beeswarm of gene expression by genotype for that pair, with an optional
#' colour-by selector for sample-level variables (e.g. sex, batch).
#' If regression statistics are not supplied, the app provides a gene-symbol
#' selector and computes/memoises [qtlRegressionStats()] results within the
#' Shiny session.
#'
#' @param res Optional wide data frame from [qtlRegressionStats()] with
#'   \code{t_only = TRUE}.  Must contain \code{phenotype_id} and
#'   \code{variant_id} columns.  If omitted, the app lists available gene
#'   symbols from \code{rowRanges(tqe)$gene_name} and computes
#'   \code{qtlRegressionStats(tqe, symbol = ...)} on request.
#' @param tqe A [tQTLExperiment] with genotype data.
#' @param assayName Name of the assay to use. Defaults to the first assay.
#' @param collapse_patterns Named list of regex patterns passed to
#'   [collapseFactorTstats()].  Defaults handle MAGE-style batch and
#'   population columns.
#' @param window Integer. Cis-window half-width passed to
#'   [qtlRegressionStats()] when \code{res} is omitted.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object passed to
#'   [qtlRegressionStats()] when \code{res} is omitted.
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
#'   qtlPCABrowser(tqe)
#' }
qtlPCABrowser <- function(res = NULL, tqe = NULL, assayName = NULL,
                          collapse_patterns = list(
                              batch      = "^t\\.factor_batch_",
                              population = "^t\\.population"
                          ),
                          window = 1000000L,
                          BPPARAM = BiocParallel::SerialParam()) {
    if (!requireNamespace("plotly", quietly = TRUE))
        stop("Package 'plotly' is required.")

    if (is.null(tqe) && methods::is(res, "tQTLExperiment")) {
        tqe <- res
        res <- NULL
    }
    if (is.null(tqe))
        stop("'tqe' must be supplied, or passed as the first argument")

    if (is.null(assayName))
        assayName <- SummarizedExperiment::assayNames(tqe)[1L]

    # ---- Pre-extract data for beeswarm rendering ---------------------------
    row_gr      <- SummarizedExperiment::rowRanges(tqe, use.names = TRUE)
    pheno_names <- names(row_gr)
    gene_names  <- S4Vectors::mcols(row_gr)[["gene_name"]]
    if (is.null(gene_names))
        gene_names <- rep(NA_character_, length(row_gr))
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

    gene_label <- function(phenotype_id) {
        idx <- match(phenotype_id, pheno_names)
        sym <- gene_names[idx]
        ifelse(!is.na(sym) & nzchar(sym), sym, phenotype_id)
    }

    gene_mode <- is.null(res)
    gene_symbols <- sort(unique(gene_names[!is.na(gene_names) &
                                           nzchar(gene_names)]))
    if (gene_mode && length(gene_symbols) == 0L)
        stop("'tqe' has no gene symbols in rowRanges gene_name - ",
             "run addGeneSymbols() first")

    make_pca_data <- function(stat_res, selected_symbol = NULL) {
        if (!is.data.frame(stat_res))
            stop("'res' must be a data.frame from qtlRegressionStats()")
        required <- c("phenotype_id", "variant_id")
        if (!all(required %in% names(stat_res)))
            stop("'res' must contain columns 'phenotype_id' and 'variant_id'")

        res_clean <- stats::na.omit(stat_res)
        res_coll  <- collapseFactorTstats(res_clean, patterns = collapse_patterns)
        t_cols    <- setdiff(names(res_coll), c("phenotype_id", "variant_id"))
        if (length(t_cols) < 2L)
            stop("At least two regression statistic columns are needed for PCA")
        t_mat <- as.matrix(res_coll[, t_cols, drop = FALSE])

        pca     <- stats::prcomp(t_mat, center = TRUE, scale. = FALSE)
        pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

        pc_choices <- colnames(pca$x)
        if (is.null(pc_choices))
            pc_choices <- paste0("PC", seq_len(ncol(pca$x)))

        scores <- data.frame(pca$x, check.names = FALSE)
        scores[["variant_id"]]   <- res_clean[["variant_id"]]
        scores[["phenotype_id"]] <- res_clean[["phenotype_id"]]
        scores[["gene_symbol"]]  <- if (!is.null(selected_symbol))
            selected_symbol else gene_label(scores[["phenotype_id"]])
        scores[["row_id"]] <- seq_len(nrow(scores))
        ensg_line <- ifelse(scores[["gene_symbol"]] != scores[["phenotype_id"]],
                            paste0("<br>ENSG: ", scores[["phenotype_id"]]),
                            "")
        scores[["hover_text"]] <- paste0(
            "SNP: ", scores[["variant_id"]],
            "<br>Gene: ", scores[["gene_symbol"]],
            ensg_line
        )

        pc_axis_title <- function(pc) {
            idx <- match(pc, pc_choices)
            paste0(pc, " (", pct_var[idx], "%)")
        }
        pair_cols <- head(pc_choices, 5L)
        list(
            scores        = scores,
            pc_choices    = pc_choices,
            pair_cols     = pair_cols,
            pair_labels   = pc_axis_title(pair_cols),
            pc_axis_title = pc_axis_title,
            selected_gene = selected_symbol
        )
    }

    initial_pca <- if (!gene_mode) make_pca_data(res) else NULL
    initial_pc_choices <- if (!is.null(initial_pca))
        initial_pca[["pc_choices"]] else character(0)
    initial_x <- if (length(initial_pc_choices))
        initial_pc_choices[1L] else character(0)
    initial_y <- if (length(initial_pc_choices))
        initial_pc_choices[min(2L, length(initial_pc_choices))] else character(0)

    gene_controls <- if (gene_mode) {
        tagList(
            selectInput("gene_symbol", "Gene symbol:",
                        choices = gene_symbols,
                        selected = gene_symbols[1L]),
            actionButton("compute_gene", "Compute regression stats/PCA",
                         class = "btn-primary mb-3 w-100")
        )
    } else NULL

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
                   gene_controls,
                   fluidRow(
                       column(6,
                              selectInput("x_pc", "X-axis PC:",
                                          choices = initial_pc_choices,
                                          selected = initial_x)),
                       column(6,
                              selectInput("y_pc", "Y-axis PC:",
                                          choices = initial_pc_choices,
                                          selected = initial_y))
                   ),
                   selectInput("color_var", "Color beeswarm by:",
                               choices  = color_choices,
                               selected = if ("sexXY" %in% color_choices) "sexXY" else "none"),
                   plotOutput("beeswarm", height = "480px"))
        )
    )

    # ---- Shiny server ------------------------------------------------------
    server <- function(input, output, session) {
        stat_cache <- new.env(parent = emptyenv())

        get_gene_stats <- function(symbol) {
            key <- paste(symbol, window, assayName, sep = "\r")
            if (exists(key, envir = stat_cache, inherits = FALSE))
                return(get(key, envir = stat_cache, inherits = FALSE))

            value <- withProgress(
                message = paste("Computing regression statistics for", symbol),
                value = 0,
                {
                    qtlRegressionStats(
                        tqe,
                        symbol    = symbol,
                        window    = window,
                        assayName = assayName,
                        BPPARAM   = BPPARAM
                    )
                }
            )
            assign(key, value, envir = stat_cache)
            value
        }

        pca_data <- if (gene_mode) {
            eventReactive(input$compute_gene, {
                req(input$gene_symbol)
                make_pca_data(get_gene_stats(input$gene_symbol),
                              selected_symbol = input$gene_symbol)
            }, ignoreNULL = TRUE)
        } else {
            reactive(initial_pca)
        }

        observeEvent(pca_data(), {
            pd <- pca_data()
            updateSelectInput(session, "x_pc",
                              choices = pd[["pc_choices"]],
                              selected = pd[["pc_choices"]][1L])
            updateSelectInput(session, "y_pc",
                              choices = pd[["pc_choices"]],
                              selected = pd[["pc_choices"]][
                                  min(2L, length(pd[["pc_choices"]]))])
        })

        output$pca <- plotly::renderPlotly({
            pd <- pca_data()
            req(pd, input$x_pc, input$y_pc)
            x_pc <- input$x_pc
            y_pc <- input$y_pc
            req(x_pc %in% pd[["pc_choices"]], y_pc %in% pd[["pc_choices"]])
            scores <- pd[["scores"]]
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
                    xaxis = list(title = pd[["pc_axis_title"]](x_pc)),
                    yaxis = list(title = pd[["pc_axis_title"]](y_pc)),
                    hoverlabel = list(bgcolor = "white")
                )
        })

        output$pc_pairs <- renderPlot({
            pd <- pca_data()
            req(pd)
            pair_cols <- pd[["pair_cols"]]
            if (length(pair_cols) < 2L) {
                plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
                text(0, 0, "At least two PCs are required for a pairs plot",
                     cex = 1.1, col = "grey50")
                return(invisible(NULL))
            }
            graphics::pairs(
                pd[["scores"]][, pair_cols, drop = FALSE],
                labels = pd[["pair_labels"]],
                pch    = 16,
                cex    = 0.5,
                col    = grDevices::adjustcolor("steelblue", alpha.f = 0.5),
                main   = "First five PC score pairs"
            )
        })

        output$beeswarm <- renderPlot({
            pd <- pca_data()
            req(pd)
            scores <- pd[["scores"]]
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
            sym <- scores[["gene_symbol"]][row_idx]

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
                ggplot2::ylab(paste0("Gene expression value (", sym, ")")) +
                ggplot2::ggtitle(paste0("SNP: ", vid, "\nGene: ", sym))
        })
    }

    shiny::shinyApp(ui, server)
}
