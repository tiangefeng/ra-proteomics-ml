# Data

The proteomics dataset used in this project is **not** included in this
repository. It contains human-subject plasma proteomics measurements
obtained under restricted data-use terms that do not allow redistribution.

## Expected file

Place the following file in this directory before running the pipeline:

```
data/RHE1_20251117.tsv
```

### Format

A tab-separated file with the columns:

| Column          | Type        | Description                               |
|-----------------|-------------|-------------------------------------------|
| `DAid`          | string      | Sample identifier                         |
| `Age`           | numeric     | Age in years                              |
| `Sex`           | categorical | `F` or `M`                                |
| `Disease`       | categorical | `Healthy` or `Rheumatoid arthritis`       |
| `Subdiagnosis`  | categorical | `ACCP_positive`, `ACCP_negative`, or `NA` |
| `<protein_1>` … | numeric     | NPX values for 5,415 proteins             |

`NA` values in `Subdiagnosis` correspond to healthy controls; the
pipeline recodes them to `Healthy` in `R/data_prep.R::recode_subdiagnosis()`.

## Source

Olink Explore HT plasma proteomics data from the
[Human Disease Blood Atlas](https://disease.proteinatlas.org/),
accessed November 2025.

If you want to reproduce this analysis, request access through the Atlas
and save the resulting TSV under this directory.
