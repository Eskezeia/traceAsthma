# traceAsthma

An R package implementing **TRACE-Asthma**, an end-to-end, mechanistically
informed multi-omics asthma risk pipeline:

```
SNP -> MPRS -> CpG methylation -> eQTM score -> TF activity -> TRACE score -> Asthma risk
```


## Installation

```r
#  GitHub:
remotes::install_github("Eskezeia/traceAsthma")

# to build core (Imports-only) functionality locally from source, no
# CRAN/Bioconductor packages required beyond base R:
# R CMD build traceAsthma && R CMD INSTALL traceAsthma_0.2.0.tar.gz
```

Most pipeline stages beyond basic QC and eQTM discovery require
additional packages, listed in `Suggests` rather than `Imports` so the
package installs and the lightweight functions work without them.
Install what you need:

```r
install.packages(c("glmnet", "pROC", "PRROC", "mediation", "lavaan",
                    "xgboost", "randomForest", "rms", "dcurves", "shiny"))

install.packages("BiocManager")
BiocManager::install(c("viper", "minfi", "sva", "dorothea", "decoupleR"))
```

## Quick start

```r
library(traceAsthma)

# Recommended first command after any install: verify everything works
# end to end against the bundled benchmark dataset
run_package_diagnostics()

# No real data yet? Two options with a known ground truth:

# Option A: a fixed, shared community benchmark -- same data for every user,
# useful for comparing results/reporting issues (300 subjects, 500 CpGs, 200 genes)
data(simulated_benchmark_cohort)

# Option B: generate your own, at any size/seed (see vignette("simulating-data"))
sim <- simulate_trace_asthma_cohort(n_subjects = 300, n_cpgs = 500,
                                     n_genes = 200, n_tfs = 8, seed = 42)

result <- run_trace_asthma_auto(
  methylation = sim$methylation, expression = sim$expression,
  clinical = sim$clinical, eqtm_pairs = sim$eqtm_pairs, regulon_db = sim$regulon_db
)

# --- Once you have real data (matching the same object shapes) ---
data(example_cohort)   # tiny, fully synthetic 5-patient demo dataset

# Adaptive: works even if you only have CpG + expression + asthma + age/sex
# (no genotype, no BMI/smoking, no regulon database required)
result <- run_trace_asthma_auto(
  methylation = example_cohort$methylation,
  expression  = example_cohort$expression,
  clinical    = example_cohort$clinical[, c("patient_id","asthma","age","sex")],
  eqtm_pairs  = example_cohort$eqtm_pairs
)

# See what packages your planned analysis actually needs
list_dependencies()

# Publication-ready outputs, and plain-language clinical interpretation
export_table(format_coefficient_table(result$risk_model), "table2_coefficients.csv")
plot_trace_weights(result$trace_fit, path = "figure3_trace_weights.tiff", dpi = 600)
print(interpret_clinical_relevance(result$risk_model, result$model_data))

# --- Moving to a real dataset: same functions, real files ---
data_dictionary("clinical")             # what fields are expected, units, coding
create_data_templates("my_data_prep/")  # blank CSVs to fill in and hand to collaborators

methylation <- import_data("real_methylation.csv")  # or .xlsx/.rds/.parquet/a DBI connection
clinical    <- import_data("real_clinical.csv", as_matrix = FALSE)

validate_matrix(methylation, type = "methylation")
validate_clinical(clinical)
validate_cohort_alignment(methylation = methylation, clinical = clinical)
```

See `vignette("trace-asthma-pipeline")` for the full walkthrough (including
the fully-specified pipeline with genotype/MPRS and TRACE layers),
`vignette("statistical-methods")` for every model formula paired with its
function and a step-by-step worked example, `vignette("simulating-data")`
for testing the pipeline before real data arrives, and
`vignette("clinical-deployment-notes")` for validation, regulatory, and
publishing guidance.

A standalone PDF user manual (purpose, objective, novelty, data
description, every method's formula, step-by-step usage guide, and
conclusion) is also included:

```r
system.file("doc-manual", "traceAsthma_Manual.pdf", package = "traceAsthma")
# open it directly, e.g.:
utils::browseURL(system.file("doc-manual", "traceAsthma_Manual.pdf", package = "traceAsthma"))
```

## Package structure

| File 
| `R/01_qc.R` 
| `R/02_prs.R` 
| `R/03_eqtm.R`  
| `R/04_regulatory_network.R` 
| `R/05_trace_score.R` 
| `R/06_sem.R` 
| `R/07_risk_model.R` 
| `R/08_validation.R` 
| `R/09_predict.R` 
| `R/10_pipeline.R` 
| `R/11_figures.R` 
| `R/12_tables.R` 
| `R/13_auto_pipeline.R` 
| `R/14_clinical_interpretation.R` 
| `R/15_simulate.R` 
| `R/16_data_import.R` 
| `R/17_data_validation.R` 
| `R/18_data_dictionary.R` 
| `inst/shiny/app.R` 



## License

MIT
