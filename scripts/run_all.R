# =============================================================================
# Machine Learning-Based Protein Biomarker Discovery for Rheumatoid Arthritis
# =============================================================================
# Main analysis pipeline. Fits four elastic net classification models on
# Olink Explore HT plasma proteomics data and extracts biomarker candidates.
#
# Models:
#   Model 0 : RA vs Healthy (binary),  proteins only
#   Model 1 : Healthy vs ACPA+ vs ACPA- (multinomial), proteins only
#   Model 2 : Same as Model 1 + Age and Sex as covariates
#   Model 3 : Sex-stratified version of Model 1 (Female / Male separately)
#
# Usage:
#   Set `data_path` below to point at the Olink TSV file, then run:
#     Rscript scripts/run_all.R
#   or source() this file from an interactive R session.
# =============================================================================

# ---- Setup ------------------------------------------------------------------
# Resolve project root so paths work from anywhere
here_script <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
project_root <- normalizePath(file.path(here_script, ".."), mustWork = FALSE)
if (!dir.exists(file.path(project_root, "R"))) {
  project_root <- getwd()
}

# Source shared functions
source(file.path(project_root, "R", "dependencies.R"))
source(file.path(project_root, "R", "data_prep.R"))
source(file.path(project_root, "R", "modeling.R"))
source(file.path(project_root, "R", "evaluation.R"))
source(file.path(project_root, "R", "plotting.R"))

# Path configuration (edit as needed)
data_path   <- file.path(project_root, "data", "RHE1_20251117.tsv")
figures_dir <- file.path(project_root, "results", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1122)

# ---- Load data --------------------------------------------------------------
data <- load_proteomics(data_path)
protein_cols <- get_protein_cols(data)
na_rows <- find_na_rows(data, protein_cols)
cat("Removing", length(na_rows), "sample(s) with missing NPX values:",
    paste(na_rows, collapse = ", "), "\n")

# Sample distribution figure
p_sample <- plot_sample_distribution(data)
ggsave(
  file.path(figures_dir, "sample_distribution.png"),
  plot = p_sample, width = 6, height = 4, units = "in"
)

# =============================================================================
# Model 0 : RA vs Healthy (binary)
# =============================================================================
cat("\n=== Model 0: RA vs Healthy ===\n")

model0_df <- data %>%
  select(Disease, all_of(protein_cols))
model0_df <- model0_df[-na_rows, ]

model0 <- fit_elastic_net(model0_df, Disease, include_dummies = FALSE)

print(collect_metrics(model0$final_fit))

# ROC
roc0 <- compute_roc(model0$final_fit, Disease, event_level = "second")
auc0 <- compute_auc(model0$final_fit, Disease, event_level = "second")
p_roc0 <- plot_roc(roc0, auc0$.estimate, "Model0 ROC Curve: RA vs Healthy")
ggsave(file.path(figures_dir, "model0_roc.png"),
       plot = p_roc0, width = 5, height = 4)

# Top proteins
top0 <- extract_top_proteins(model0$fit_parsnip, n = 20)
p_top0 <- plot_top_proteins(top0, title = "Model 0 Top Proteins")
ggsave(file.path(figures_dir, "model0_top_proteins.png"),
       plot = p_top0, width = 6, height = 5)

# =============================================================================
# Model 1 : Subdiagnosis (Healthy / ACPA+ / ACPA-)
# =============================================================================
cat("\n=== Model 1: Healthy vs ACPA+ vs ACPA- ===\n")

ra_data <- recode_subdiagnosis(data)

model1_df <- ra_data %>%
  select(Subdiagnosis, all_of(protein_cols))
model1_df <- model1_df[-na_rows, ]

model1 <- fit_elastic_net(model1_df, Subdiagnosis, include_dummies = FALSE)

print(collect_metrics(model1$final_fit))

roc1 <- compute_roc(model1$final_fit, Subdiagnosis)
auc1 <- compute_auc(model1$final_fit, Subdiagnosis)
p_roc1 <- plot_roc(roc1, auc1$.estimate,
                   "Model1 ROC Curve: Healthy vs ACPA+ vs ACPA-")
ggsave(file.path(figures_dir, "model1_roc.png"),
       plot = p_roc1, width = 8, height = 4)

# Confusion matrix
cm1 <- compute_confmat(model1$final_fit, Subdiagnosis)
p_cm1 <- autoplot(cm1, type = "heatmap") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()
ggsave(file.path(figures_dir, "model1_confmat.png"),
       plot = p_cm1, width = 5, height = 4)

# Top proteins
top1 <- extract_top_proteins(model1$fit_parsnip, n = 20)
p_top1 <- plot_top_proteins(top1, title = "Model 1 Top Proteins")
ggsave(file.path(figures_dir, "model1_top_proteins.png"),
       plot = p_top1, width = 6, height = 6)

# NPX expression boxplots of the top proteins
top1_names <- unique(top1$term)
p_npx1 <- plot_expression_boxplots(
  ra_data[-na_rows, ],
  proteins = top1_names,
  title = "Model1 Biological Expression (NPX) of Top Proteins"
)
ggsave(file.path(figures_dir, "model1_expression_boxplots.png"),
       plot = p_npx1, width = 10, height = 8)

# =============================================================================
# Model 2 : Subdiagnosis + Age + Sex
# =============================================================================
cat("\n=== Model 2: Subdiagnosis + Age + Sex ===\n")

model2_df <- ra_data %>%
  select(Subdiagnosis, Sex, Age, all_of(protein_cols))
model2_df <- model2_df[-na_rows, ]
model2_df$Age <- as.numeric(model2_df$Age)

model2 <- fit_elastic_net(model2_df, Subdiagnosis, include_dummies = TRUE)

print(collect_metrics(model2$final_fit))

roc2 <- compute_roc(model2$final_fit, Subdiagnosis)
auc2 <- compute_auc(model2$final_fit, Subdiagnosis)
p_roc2 <- plot_roc(roc2, auc2$.estimate,
                   "Model2 ROC Curve: Subdiagnosis + Age + Sex")
ggsave(file.path(figures_dir, "model2_roc.png"),
       plot = p_roc2, width = 8, height = 4)

top2 <- extract_top_proteins(model2$fit_parsnip, n = 20)
p_top2 <- plot_top_proteins(top2, title = "Model 2 Top Proteins")
ggsave(file.path(figures_dir, "model2_top_proteins.png"),
       plot = p_top2, width = 6, height = 6)

# =============================================================================
# Model 3 : Sex-stratified
# =============================================================================
cat("\n=== Model 3: Sex-stratified ===\n")

#' Run Model 1 separately within one sex
run_model3 <- function(sex_value, full_df, protein_cols, na_rows) {
  df_sex <- full_df %>%
    mutate(Sex = factor(Sex)) %>%
    select(Subdiagnosis, Sex, all_of(protein_cols))
  df_sex <- df_sex[-na_rows, ]
  df_sex <- df_sex %>%
    filter(Sex == sex_value) %>%
    select(-Sex)

  cat("  Sex =", sex_value, "| class counts:\n")
  print(table(df_sex$Subdiagnosis))

  fit <- fit_elastic_net(df_sex, Subdiagnosis, include_dummies = FALSE)
  roc <- compute_roc(fit$final_fit, Subdiagnosis)
  auc <- compute_auc(fit$final_fit, Subdiagnosis)
  top <- extract_top_proteins(fit$fit_parsnip, n = 10)

  list(sex = sex_value, fit = fit, roc = roc, auc = auc, top = top)
}

res3_F <- run_model3("F", ra_data, protein_cols, na_rows)
res3_M <- run_model3("M", ra_data, protein_cols, na_rows)

for (res in list(res3_F, res3_M)) {
  p_roc <- plot_roc(
    res$roc, res$auc$.estimate,
    paste0("Model3 ROC Curve (Sex = ", res$sex, ")")
  )
  ggsave(
    file.path(figures_dir, paste0("model3_roc_", res$sex, ".png")),
    plot = p_roc, width = 8, height = 4
  )

  p_top <- plot_top_proteins(
    res$top,
    title = paste0("Model3 Top 10 Proteins (Sex = ", res$sex, ")")
  )
  ggsave(
    file.path(figures_dir, paste0("model3_top_proteins_", res$sex, ".png")),
    plot = p_top, width = 6, height = 5
  )
}

# =============================================================================
# Functional enrichment on Model 1 top proteins
# =============================================================================
cat("\n=== GO Enrichment (Model 1 top proteins) ===\n")

if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("enrichplot", quietly = TRUE)) {

  source(file.path(project_root, "R", "enrichment.R"))

  top50_model1 <- extract_top_proteins(model1$fit_parsnip, n = 50)
  ego1 <- run_go_enrichment(
    target_symbols = unique(top50_model1$term),
    background_symbols = protein_cols
  )
  print(head(ego1))

  p_ego1 <- plot_go_dotplot(ego1, show_n = 15,
                            title = "Model1 GO Pathway Enrichment")
  ggsave(file.path(figures_dir, "model1_go_enrichment.png"),
         plot = p_ego1, width = 7, height = 5)
} else {
  cat("  Skipping enrichment: install 'clusterProfiler', 'org.Hs.eg.db', ",
      "'enrichplot' via BiocManager.\n")
}

cat("\nAll done. Figures written to:", figures_dir, "\n")
