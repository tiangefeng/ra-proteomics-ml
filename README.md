# Machine Learning-Based Protein Biomarker Discovery for Rheumatoid Arthritis

Elastic net classification models applied to Olink Explore HT plasma
proteomics data to identify circulating protein biomarkers that distinguish
rheumatoid arthritis (RA) subtypes (ACPA-positive vs ACPA-negative) from
healthy controls.

> Course project for **CB2050 Project in Molecular Life Science**, KTH /
> SciLifeLab. 

---

## Overview

The pipeline fits four regularized multinomial/binary classifiers on 5,415
plasma proteins measured across 50 healthy controls and 200 RA patients
(100 ACPA+, 100 ACPA-):

| Model   | Task                                       | Covariates          |
|---------|--------------------------------------------|---------------------|
| Model 0 | RA vs Healthy                              | Proteins only       |
| Model 1 | Healthy vs ACPA+ vs ACPA-                  | Proteins only       |
| Model 2 | Healthy vs ACPA+ vs ACPA-                  | Proteins + Age + Sex |
| Model 3 | Sex-stratified Model 1 (Female / Male)     | Proteins only       |

All models use elastic net regularization (`glmnet`) with 5-fold
cross-validated hyperparameter tuning over a 10×10 grid of
`penalty ∈ [10⁻⁵, 10¹]` and `mixture ∈ [0, 1]`. Feature importance is
extracted from the signed model coefficients, and Gene Ontology Biological
Process enrichment is run on the top predictors from Model 1.

Consistently prioritized proteins include **PDCD1, IL6, PADI4, CXCL13**, and
**MSMP** 

---

## Repository layout

```
ra-proteomics-ml/
├── R/                      # Reusable function modules (source this)
│   ├── dependencies.R      # Package loading
│   ├── data_prep.R         # Loading & preprocessing
│   ├── modeling.R          # Elastic net training & tuning
│   ├── evaluation.R        # ROC, AUC, confusion matrices
│   ├── plotting.R          # All plotting helpers
│   └── enrichment.R        # GO enrichment (clusterProfiler)
├── scripts/
│   └── run_all.R           # Main driver; runs all four models end-to-end
├── data/                   # Input TSV goes here (not tracked by git)
│   └── README.md           # How to obtain the data
├── results/
│   └── figures/            # PNG outputs written by run_all.R
├── docs/
└── README.md
```

---

## Getting started

### Requirements

- R ≥ 4.2
- CRAN packages:
  `tidyverse`, `tidymodels`, `glmnet`, `themis`, `yardstick`, `patchwork`,
  `ggplot2`, `rlang`, `broom`
- Bioconductor packages (for the enrichment step only):
  `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`

```r
install.packages(c("tidyverse", "tidymodels", "glmnet", "themis",
                   "yardstick", "patchwork"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
```

### Data

The Olink Explore HT dataset (`RHE1_20251117.tsv`) is **not** included in
this repository. It was obtained from the Human Disease Blood Atlas
(<https://disease.proteinatlas.org/>) under data use terms that do not
permit redistribution.

Place the TSV under `data/` before running any script:

```
data/RHE1_20251117.tsv
```

The file is expected to contain the columns:
`DAid, Age, Sex, Disease, Subdiagnosis, <5,415 protein NPX columns>`.

See [`data/README.md`](data/README.md) for details.

### Running the pipeline

From the project root:

```bash
Rscript scripts/run_all.R
```

or, interactively:

```r
setwd("ra-proteomics-ml")
source("scripts/run_all.R")
```

All figures are written to `results/figures/`.

---

## Reusing the functions

Each file in `R/` is self-contained and can be sourced independently. A
minimal end-to-end example:

```r
source("R/dependencies.R")
source("R/data_prep.R")
source("R/modeling.R")
source("R/evaluation.R")

data         <- load_proteomics("data/RHE1_20251117.tsv")
protein_cols <- get_protein_cols(data)
na_rows      <- find_na_rows(data, protein_cols)

df  <- data %>% select(Disease, all_of(protein_cols))
df  <- df[-na_rows, ]

fit <- fit_elastic_net(df, Disease)
compute_auc(fit$final_fit, Disease, event_level = "second")
```

---

## Methods summary

- **Preprocessing**: two samples (DA19171, DA19162) with missing NPX values
  are removed. Protein expression is z-score normalized within each
  training split via `recipes::step_normalize()`.
- **Model fitting**: `parsnip::multinom_reg(engine = "glmnet")` with
  stratified 75/25 train/test split and 5-fold CV tuning. The final
  hyperparameters are selected by the one-standard-error rule on ROC AUC.
- **Feature selection**: non-zero coefficients from the final glmnet fit
  are ranked by `|coefficient|`.
- **Enrichment**: `clusterProfiler::enrichGO()` on GO Biological Process,
  with all measured proteins as the background universe and a
  Benjamini-Hochberg adjusted p-value threshold of 0.05.

---

## Results highlights

| Model | Task | Test ROC AUC |
|-------|------|--------------|
| 0 | RA vs Healthy | 0.842 |
| 1 | Healthy / ACPA+ / ACPA- | 0.888 (one-vs-rest, all classes) |
| 2 | + Age + Sex | 0.888 |
| 3 | Female-only | 0.767 |
| 3 | Male-only | 0.692 |

Demographic covariates added negligible discriminative power over plasma
proteins alone. Sex-stratified models identified partially distinct
feature sets, though male-specific results are limited by the small
subgroup size.

Enriched GO processes centered on **regulation of cell activation**,
**adaptive immune response**, and innate-immune / bacterial-response
pathways — consistent with the roles of PDCD1, IL6, CXCL13, and PADI4 in
RA pathogenesis.

See [`docs/Project_Report.pdf`](docs/Project_Report.pdf) for full figures
and biological interpretation.

---

## Acknowledgements

Data: Human Disease Blood Atlas.
Proteomics platform: Olink Explore HT (PEA + NGS readout).
