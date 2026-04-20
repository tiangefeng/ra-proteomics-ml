# -----------------------------------------------------------------------------
# Data loading and preprocessing
# -----------------------------------------------------------------------------
# Utilities for reading the Olink Explore HT dataset and preparing it for
# downstream modeling.
#
# Expected input: a TSV file with columns
#     DAid, Age, Sex, Disease, Subdiagnosis, <protein_1>, ..., <protein_N>
# where the protein columns contain NPX values.
# -----------------------------------------------------------------------------


#' Load the raw Olink proteomics dataset
#'
#' @param path Path to the TSV file.
#' @return A tibble with all columns as-read.
load_proteomics <- function(path) {
  if (!file.exists(path)) {
    stop("Data file not found: ", path,
         "\nPlace the raw TSV under the `data/` directory, or update the ",
         "`data_path` variable in scripts/run_all.R.")
  }
  readr::read_tsv(path, show_col_types = FALSE)
}


#' Extract the vector of protein column names
#'
#' All columns other than the known metadata columns are treated as proteins.
#'
#' @param data The full dataset.
#' @return Character vector of protein column names.
get_protein_cols <- function(data) {
  meta_cols <- c("DAid", "Age", "Sex", "Disease", "Subdiagnosis")
  setdiff(colnames(data), meta_cols)
}


#' Identify rows with missing protein measurements
#'
#' Two samples (DA19171, DA19162) in the original cohort have missing NPX
#' values. This helper returns their row indices so they can be excluded
#' consistently across all models.
#'
#' @param data The full dataset.
#' @param protein_cols Character vector of protein columns.
#' @return Integer vector of row indices with any NA in protein columns.
find_na_rows <- function(data, protein_cols) {
  which(rowSums(is.na(data[, protein_cols])) > 0)
}


#' Recode Subdiagnosis so that healthy controls get their own factor level
#'
#' In the raw data, Subdiagnosis is NA for healthy controls. For multinomial
#' classification we replace those NAs with "Healthy" and return the column
#' as a factor.
#'
#' @param data The full dataset.
#' @return A data frame with Subdiagnosis as a three-level factor.
recode_subdiagnosis <- function(data) {
  data %>%
    dplyr::mutate(
      Subdiagnosis = ifelse(is.na(Subdiagnosis), "Healthy", Subdiagnosis),
      Subdiagnosis = factor(Subdiagnosis)
    )
}
