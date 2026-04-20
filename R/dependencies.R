# -----------------------------------------------------------------------------
# Package dependencies
# -----------------------------------------------------------------------------
# All packages used across the project are loaded here. Sourcing this file at
# the top of any script guarantees a consistent environment.
#
# If any of these are missing, install them with:
#   install.packages(c("tidyverse", "tidymodels", "glmnet", "themis",
#                      "yardstick", "patchwork", "ggplot2"))
#   BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
# -----------------------------------------------------------------------------

# Core data manipulation & modeling
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(tidymodels)
  library(glmnet)
  library(themis)
  library(yardstick)

  # Plotting
  library(patchwork)
  library(ggplot2)
})

# Enrichment analysis (Bioconductor). Loaded lazily inside functions that
# need them, so that the rest of the pipeline can run without Bioconductor
# installed. See R/enrichment.R.
