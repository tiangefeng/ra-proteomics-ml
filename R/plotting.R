# -----------------------------------------------------------------------------
# Plotting utilities
# -----------------------------------------------------------------------------
# Reusable plotting functions for sample distribution, ROC curves, top-protein
# bar plots, and NPX expression boxplots across subgroups.
# -----------------------------------------------------------------------------


#' Plot sample distribution (disease / sex / age)
#'
#' @param data The full dataset.
#' @return A patchwork composite plot.
plot_sample_distribution <- function(data) {
  p1 <- ggplot2::ggplot(data, ggplot2::aes(x = Disease, fill = Subdiagnosis)) +
    ggplot2::geom_bar() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Disease Distribution with Subdiagnosis",
      x = "Disease", y = "count", fill = "Subdiagnosis"
    )

  p2 <- ggplot2::ggplot(data, ggplot2::aes(Sex)) +
    ggplot2::geom_bar(fill = "#4E79A7") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Sex Distribution", x = "", y = "Count")

  p3 <- ggplot2::ggplot(data, ggplot2::aes(Age)) +
    ggplot2::geom_histogram(fill = "#F28E2B", bins = 20, color = "white") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Age Distribution", x = "Age", y = "Count")

  (p1 | p2) / p3
}


#' Plot an ROC curve with AUC annotation
#'
#' @param roc_data Output of `compute_roc()`.
#' @param auc_value Numeric AUC estimate(s).
#' @param title Plot title.
#' @return A ggplot object.
plot_roc <- function(roc_data, auc_value, title) {
  auc_label <- paste0("AUC = ", round(mean(auc_value, na.rm = TRUE), 3))

  ggplot2::autoplot(roc_data) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(title) +
    ggplot2::annotate(
      "text",
      x = 0.6, y = 0.2,
      label = auc_label,
      size = 3
    )
}


#' Bar plot of top proteins by signed coefficient, faceted by class for
#' multiclass models
#'
#' @param top_proteins A tibble from `extract_top_proteins()`.
#' @param title Plot title.
#' @param facet_by_class If TRUE, create one panel per outcome class.
#' @return A ggplot object.
plot_top_proteins <- function(top_proteins,
                              title = "Top Proteins",
                              facet_by_class = FALSE) {
  p <- ggplot2::ggplot(
    top_proteins,
    ggplot2::aes(x = stats::reorder(term, abs(estimate)),
                 y = estimate,
                 fill = class)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      x = "Protein",
      y = "Coefficient Estimate",
      title = title
    )

  if (facet_by_class) {
    p <- p + ggplot2::facet_wrap(~class, scales = "free_y")
  }

  p
}


#' Boxplots of NPX expression for selected proteins across subgroups
#'
#' @param data Data containing the grouping column and protein columns.
#' @param proteins Character vector of protein columns to plot.
#' @param group_col Name of the grouping column as a string. Default
#'   "Subdiagnosis".
#' @param title Plot title.
#' @return A ggplot object.
plot_expression_boxplots <- function(data,
                                     proteins,
                                     group_col = "Subdiagnosis",
                                     title = "Biological Expression (NPX)") {
  plot_df <- data %>%
    dplyr::select(dplyr::all_of(c(group_col, proteins))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(proteins),
      names_to = "Protein",
      values_to = "NPX"
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[group_col]], y = NPX, fill = .data[[group_col]])
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
    ggplot2::facet_wrap(~Protein, scales = "free_y") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::ggtitle(title)
}
