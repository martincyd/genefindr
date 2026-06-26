#' Get gene characterization data across multiple cancer sites as a table
#'
#' @param gene Character. Gene symbol (e.g. "TP53")
#' @param sites Character vector of cancer sites (e.g. c("breast", "lung", "colon"))
#' @return A data frame with results for each site
#' @export
#' @examples
#' \donttest{
#' results <- findr_sites_table("TP53", sites = c("breast", "lung", "colon"))
#' findr_plot_sites(results)
#' }
findr_sites_table <- function(gene, sites) {
  results <- list()
  for (s in sites) {
    res <- suppressMessages(findr(gene, site = s))
    if (!is.null(res)) results[[s]] <- res
  }
  do.call(rbind, results)
}


#' Plot gene association scores for multiple genes in a single disease context
#'
#' @param results A data frame returned by findr_multi() with output = "table"
#' @param title Optional plot title
#' @return A ggplot2 object
#' @export
#' @examples
#' \donttest{
#' results <- findr_multi(c("TP53", "BRCA1", "MYC"), site = "breast", output = "table")
#' findr_plot_genes(results)
#' }
findr_plot_genes <- function(results, title = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Install with install.packages('ggplot2')")
  }

  # check required columns exist
  if (!all(c("gene", "association_score") %in% names(results))) {
    stop("Results must contain 'gene' and 'association_score' columns. Use findr_multi(output = 'table')")
  }

  # remove rows with no association score
  plot_data <- results[!is.na(results$association_score), ]

  if (nrow(plot_data) == 0) {
    stop("No association scores found in results.")
  }

  # sort by association score
  plot_data$gene <- factor(plot_data$gene,
                           levels = plot_data$gene[order(plot_data$association_score)])

  # build plot
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = association_score,
                                    y = gene,
                                    fill = mutation_frequency_pct)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                        color = "gray40", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = "#FFE4E1", high = "#C71585",
                                 name = "Mutation\nfrequency (%)",
                                 na.value = "gray80") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = if (!is.null(title)) title else paste("Gene Association Scores -", unique(results$context)),
      x = "Open Targets association score",
      y = NULL,
      caption = "Dashed line = 0.5 threshold | Color = TCGA mutation frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "right"
    )

  return(p)
}


#' Plot association scores for a single gene across multiple disease sites
#'
#' @param results A data frame returned by findr_multi() with output = "table"
#' @param title Optional plot title
#' @return A ggplot2 object
#' @export
#' @examples
#' \donttest{
#' results <- findr_multi("TP53", site = c("breast", "lung", "colon"), output = "table")
#' findr_plot_sites(results)
#' }
findr_plot_sites <- function(results, title = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Install with install.packages('ggplot2')")
  }

  if (!all(c("gene", "context", "association_score") %in% names(results))) {
    stop("Results must contain 'gene', 'context' and 'association_score' columns.")
  }

  plot_data <- results[!is.na(results$association_score), ]

  if (nrow(plot_data) == 0) {
    stop("No association scores found in results.")
  }

  plot_data$context <- factor(plot_data$context,
                              levels = plot_data$context[order(plot_data$association_score)])

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = association_score,
                                    y = context,
                                    fill = mutation_frequency_pct)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                        color = "gray40", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = "#FFE4E1", high = "#C71585",
                                 name = "Mutation\nfrequency (%)",
                                 na.value = "gray80") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = if (!is.null(title)) title else paste(unique(results$gene), "- Association Scores by Site"),
      x = "Open Targets association score",
      y = NULL,
      caption = "Dashed line = 0.5 threshold | Color = TCGA mutation frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "right"
    )

  return(p)
}


#' Plot association scores for a single gene across multiple disease sites
#'
#' @param results A data frame returned by findr_multi() with output = "table"
#' @param title Optional plot title
#' @return A ggplot2 object
#' @export
#' @examples
#' \donttest{
#' results <- findr_multi("TP53", site = c("breast", "lung", "colon"), output = "table")
#' findr_plot_sites(results)
#' }
findr_plot_sites <- function(results, title = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Install with install.packages('ggplot2')")
  }

  if (!all(c("gene", "context", "association_score") %in% names(results))) {
    stop("Results must contain 'gene', 'context' and 'association_score' columns.")
  }

  plot_data <- results[!is.na(results$association_score), ]

  if (nrow(plot_data) == 0) {
    stop("No association scores found in results.")
  }

  plot_data$context <- factor(plot_data$context,
                              levels = plot_data$context[order(plot_data$association_score)])

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = association_score,
                                    y = context,
                                    fill = mutation_frequency_pct)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                        color = "gray40", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = "#FFE4E1", high = "#C71585",
                                 name = "Mutation\nfrequency (%)",
                                 na.value = "gray80") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = if (!is.null(title)) title else paste(unique(results$gene), "- Association Scores by Site"),
      x = "Open Targets association score",
      y = NULL,
      caption = "Dashed line = 0.5 threshold | Color = TCGA mutation frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "right"
    )

  return(p)
}
