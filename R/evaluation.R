# -----------------------------------------------------------------------------
# Model evaluation utilities
# -----------------------------------------------------------------------------
# Helpers for computing and extracting classification metrics (multiclass
# ROC AUC, accuracy) and confusion matrices from a fitted model.
# -----------------------------------------------------------------------------


#' Compute multiclass ROC curve data from a last_fit() result
#'
#' Works for both binary and multiclass outcomes. For binary problems supply
#' `event_level = "second"` when the positive class is the second factor
#' level.
#'
#' @param final_fit A `last_fit()` result.
#' @param outcome Unquoted name of the outcome column.
#' @param event_level For binary outcomes only. Either "first" or "second".
#' @return A tibble ready for `autoplot()`.
compute_roc <- function(final_fit, outcome, event_level = NULL) {
  outcome_sym <- rlang::ensym(outcome)
  preds <- tune::collect_predictions(final_fit)

  n_classes <- length(levels(preds[[rlang::as_name(outcome_sym)]]))

  if (n_classes == 2 && !is.null(event_level)) {
    # Binary: a single positive-class probability column
    pos_level <- levels(preds[[rlang::as_name(outcome_sym)]])[
      if (event_level == "second") 2 else 1
    ]
    prob_col <- paste0(".pred_", pos_level)
    yardstick::roc_curve(
      preds,
      truth = !!outcome_sym,
      !!rlang::sym(prob_col),
      event_level = event_level
    )
  } else {
    # Multiclass: use all probability columns
    yardstick::roc_curve(
      preds,
      truth = !!outcome_sym,
      dplyr::starts_with(".pred_") & !".pred_class"
    )
  }
}


#' Compute ROC AUC from a last_fit() result
#'
#' @inheritParams compute_roc
#' @return A tibble with .estimate holding the AUC value(s).
compute_auc <- function(final_fit, outcome, event_level = NULL) {
  outcome_sym <- rlang::ensym(outcome)
  preds <- tune::collect_predictions(final_fit)

  n_classes <- length(levels(preds[[rlang::as_name(outcome_sym)]]))

  if (n_classes == 2 && !is.null(event_level)) {
    pos_level <- levels(preds[[rlang::as_name(outcome_sym)]])[
      if (event_level == "second") 2 else 1
    ]
    prob_col <- paste0(".pred_", pos_level)
    yardstick::roc_auc(
      preds,
      truth = !!outcome_sym,
      !!rlang::sym(prob_col),
      event_level = event_level
    )
  } else {
    yardstick::roc_auc(
      preds,
      truth = !!outcome_sym,
      dplyr::starts_with(".pred_") & !".pred_class"
    )
  }
}


#' Compute a confusion matrix from a last_fit() result
#'
#' @param final_fit A `last_fit()` result.
#' @param outcome Unquoted name of the outcome column.
#' @return A yardstick `conf_mat` object.
compute_confmat <- function(final_fit, outcome) {
  outcome_sym <- rlang::ensym(outcome)
  tune::collect_predictions(final_fit) %>%
    yardstick::conf_mat(truth = !!outcome_sym, estimate = .pred_class)
}
