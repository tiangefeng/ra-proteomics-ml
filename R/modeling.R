# -----------------------------------------------------------------------------
# Elastic net modeling
# -----------------------------------------------------------------------------
# Generic functions that build, tune, and fit multinomial/binary elastic net
# classifiers. All four models in the project (Model 0, 1, 2, 3) are special
# cases of `fit_elastic_net()` with different formulas and covariates.
#
# `fit_elastic_net()` accepts the outcome as a bare name (e.g. `Disease`);
# internally it is captured as a string and passed on to `build_workflow()`,
# which works exclusively with strings.
# -----------------------------------------------------------------------------


#' Build the standard elastic net tuning grid
#'
#' Penalty (lambda) varies across 10^-5 to 10^1; mixture (alpha) from 0 (ridge)
#' to 1 (lasso). A 10x10 regular grid is used throughout the project.
#'
#' @return A tibble grid usable by `tune::tune_grid()`.
make_tune_grid <- function() {
  grid_regular(
    penalty(range = c(-5, 1)),
    mixture(range = c(0, 1)),
    levels = 10
  )
}


#' Build a tidymodels workflow for elastic net classification
#'
#' @param train Training data, containing the outcome column plus predictors.
#' @param outcome_name Name of the outcome column as a **string**.
#' @param include_dummies If TRUE, `step_dummy()` is applied to nominal
#'   predictors before normalization. Use this when covariates like Sex are
#'   part of the feature set.
#' @return A tidymodels `workflow` with a tunable elastic net spec.
build_workflow <- function(train, outcome_name, include_dummies = FALSE) {
  formula_obj <- stats::as.formula(paste(outcome_name, "~ ."))

  if (include_dummies) {
    rec <- recipes::recipe(formula_obj, data = train) %>%
      recipes::step_dummy(recipes::all_nominal_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric_predictors())
  } else {
    rec <- recipes::recipe(formula_obj, data = train) %>%
      recipes::step_normalize(recipes::all_predictors())
  }

  spec <- parsnip::multinom_reg(
    penalty = tune(),
    mixture = tune()
  ) %>%
    parsnip::set_engine("glmnet") %>%
    parsnip::set_mode("classification")

  workflows::workflow(rec, spec)
}


#' Fit an elastic net classifier end-to-end
#'
#' Splits data, creates 5-fold CV folds, tunes penalty/mixture, selects the
#' one-standard-error model by ROC AUC, and fits the final model on the full
#' training set.
#'
#' @param df Data frame containing the outcome and all predictors.
#' @param outcome Unquoted name of the outcome column (e.g. `Disease`).
#' @param include_dummies Pass TRUE when nominal covariates are present.
#' @param prop Proportion of samples in the training split. Default 0.75.
#' @param v Number of CV folds. Default 5.
#' @param seed Random seed for reproducibility.
#' @return A list with the split, CV folds, tuning results, best params,
#'   final fit, and the fitted parsnip model object.
fit_elastic_net <- function(df,
                            outcome,
                            include_dummies = FALSE,
                            prop = 0.75,
                            v = 5,
                            seed = 1122) {
  set.seed(seed)
  outcome_name <- rlang::as_name(rlang::ensym(outcome))

  split <- rsample::initial_split(df, prop = prop, strata = !!outcome_name)
  train <- rsample::training(split)
  folds <- rsample::vfold_cv(train, v = v, strata = !!outcome_name)

  wf <- build_workflow(train, outcome_name, include_dummies = include_dummies)
  grid <- make_tune_grid()

  tuned <- tune::tune_grid(
    wf,
    folds,
    grid = grid,
    metrics = yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy)
  )

  best <- tune::select_by_one_std_err(
    tuned,
    metric = "roc_auc",
    desc(penalty)
  )

  final_fit <- wf %>%
    tune::finalize_workflow(best) %>%
    tune::last_fit(split)

  fit_parsnip <- tune::extract_fit_parsnip(final_fit$.workflow[[1]])

  list(
    split = split,
    folds = folds,
    tuned = tuned,
    best = best,
    final_fit = final_fit,
    fit_parsnip = fit_parsnip
  )
}


#' Extract top proteins by absolute coefficient magnitude
#'
#' @param fit_parsnip The `fit_parsnip` element returned by fit_elastic_net().
#' @param n Number of top predictors to return. Default 20.
#' @return A tibble with one row per (term, class) pair, sorted by |estimate|.
extract_top_proteins <- function(fit_parsnip, n = 20) {
  broom::tidy(fit_parsnip) %>%
    dplyr::filter(term != "(Intercept)", estimate != 0) %>%
    dplyr::mutate(abs_est = abs(estimate)) %>%
    dplyr::arrange(dplyr::desc(abs_est)) %>%
    utils::head(n)
}
