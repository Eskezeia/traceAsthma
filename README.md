# traceAsthma

An R package implementing **TRACE-Asthma**, an end-to-end, mechanistically
informed multi-omics asthma risk pipeline:

```
SNP -> MPRS -> CpG methylation -> eQTM score -> TF activity -> TRACE score -> Asthma risk
```

The package covers the full path from raw-ish QC summary statistics
through a clinician-facing single-patient prediction interface, and
includes an optional Shiny point-of-care app.

**Status: research software, not a validated clinical tool.** See
`vignette("clinical-deployment-notes")` before using this for anything
beyond research.

**Author:** Eskezeia Yihunie Dessie, PhD — Postdoctoral Research Fellow,
Division of Pulmonary, Critical Care, Sleep & Occupational Medicine,
Indiana University School of Medicine.

## Project goals

**Scientific goals**
- Provide a reproducible, auditable implementation of the TRACE-Asthma
  mechanistic architecture (SNP -> MPRS -> CpG methylation -> eQTM ->
  TF activity -> TRACE -> asthma risk), so every score in the companion
  manuscript is traceable to exact, versioned code rather than a
  one-off analysis script.
- Formally quantify how much of inherited genetic risk for asthma is
  statistically mediated through methylation-dependent transcriptional
  regulation, via the structural equation modeling layer
  (`fit_mediation_sem()`), rather than only reporting association.
- Extend the sex-stratified eQTM discovery/replication framework from
  prior work into a disease-interaction, TF-mediated risk-prediction
  pipeline with data-driven (not literature-picked) candidate
  transcription factors (`test_regulon_enrichment()`).

**Engineering goals**
- Degrade gracefully to whatever data a site actually has
  (`run_trace_asthma_auto()`) rather than requiring the full genotype +
  methylation + expression + complete clinical panel specification.
- Produce submission-ready output directly from the pipeline: 600 dpi
  Nature-style figures, fully labeled and footnoted tables, and a
  plain-language clinical-relevance translation layer
  (`interpret_clinical_relevance()`), so statistical output and
  manuscript output stay in sync.
- Ship with a realistic data simulator with an injected, known ground
  truth (`simulate_trace_asthma_cohort()`), so pipeline correctness and
  statistical power can be verified before any real cohort data is used.
- Reach a genuinely installable, documented, tested state -- GitHub +
  r-universe first, Bioconductor as the realistic longer-term home given
  the genomics dependency stack (see `vignette("clinical-deployment-notes")`
  for the current gap list toward that).

**Portfolio goal**
- This package is also a working demonstration of an end-to-end
  reproducible bioinformatics/biostatistics pipeline -- QC, statistical
  modeling, mechanistic causal inference, validation, and deployment --
  intended as a public, citable artifact of applied multi-omics and
  causal-ML engineering capability.

## Installation

```r
# once published to r-universe / GitHub:
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

## Package structure

| File | Stage |
|---|---|
| `R/01_qc.R` | Genotype / methylation / RNA-seq QC |
| `R/02_prs.R` | Multi-trait polygenic risk score (MPRS) |
| `R/03_eqtm.R` | Cis/trans eQTM discovery, mediation filtering, eQTM score |
| `R/04_regulatory_network.R` | ARACNe/DoRothEA regulons, VIPER TF activity |
| `R/05_trace_score.R` | TF filtering, stability selection, TRACE score |
| `R/06_sem.R` | MPRS -> eQTM -> TRACE -> Asthma mediation SEM |
| `R/07_risk_model.R` | Final logistic model, nested model comparison, nested CV |
| `R/08_validation.R` | AUROC/AUPRC, calibration, NRI/IDI, decision curve analysis |
| `R/09_predict.R` | Deployable model object + clinician-facing prediction |
| `R/10_pipeline.R` | End-to-end orchestrator (full data specification) |
| `R/11_figures.R` | 600 dpi publication figures (Nature-style labeling) |
| `R/12_tables.R` | Publication tables with footnotes, CSV/xlsx export |
| `R/13_auto_pipeline.R` | Adaptive pipeline for partial/incomplete data + dependency listing |
| `R/14_clinical_interpretation.R` | Plain-language absolute risk, calibration, and mediation narrative |
| `R/15_simulate.R` | Realistic genetic/epigenetic/clinical data simulation with ground truth |
| `R/16_data_import.R` | Multi-format data import/export (CSV, Excel, RDS, Parquet, database) |
| `R/17_data_validation.R` | Structural/range/alignment validation before analysis |
| `R/18_data_dictionary.R` | Data dictionary and blank-template generation for real-data mapping |
| `inst/shiny/app.R` | Point-of-care Shiny app |

## Data requirements

No proprietary, confidential, or patient-identifiable data is included
anywhere in this package. All bundled examples (`example_cohort`), and
every example/test/vignette, use either the tiny fully-synthetic
`example_cohort` or `simulate_trace_asthma_cohort()`'s realistic simulated
data -- see `vignette("simulating-data")`. Every analysis function expects
the same object shapes regardless of data source (simulated or real); see
`data_dictionary()` for the full field-level specification and
`vignette("simulating-data")` section "Moving from simulated to real data"
for the import/validation workflow.

## Development notes

This package was authored and unit-tested (core logic via base-R
equivalents of the `testthat` suite in `tests/testthat/`, since `testthat`
itself requires CRAN access unavailable in the authoring sandbox) in a
sandboxed environment without CRAN/Bioconductor network access. `R CMD
check` confirms 0 errors and full documentation/example coverage; the
remaining warnings (unbuilt vignettes needing `knitr`, a sandbox locale
artifact, and Suggests-package availability for checking) are
environmental to the authoring sandbox, not package defects -- see
`vignette("clinical-deployment-notes")` for the current, itemized status
toward full CRAN/Bioconductor submission readiness.

## License

MIT
