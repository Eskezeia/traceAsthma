# traceAsthma 0.3.0

CRAN-readiness and real-world data deployment:

- **Full documentation coverage**: all 58 exported functions now have
  `.Rd` pages (13 new grouped man pages covering QC, PRS, eQTM, regulatory
  network, TRACE, mediation, risk model, validation, deployment, and
  utility functions that were previously undocumented). `R CMD check`
  confirms "checking for missing documentation entries ... OK" and
  "checking examples ... OK" -- every example runs successfully.
- **Real author/maintainer information**: `DESCRIPTION`, `LICENSE`, and
  `inst/CITATION` updated with real author (Eskezeia Yihunie Dessie, PhD,
  Indiana University School of Medicine), email, and GitHub URL,
  replacing all placeholder values.
- **`Additional_repositories`** added to `DESCRIPTION` for the
  Bioconductor-only `Suggests` packages (`viper`, `minfi`, `sva`,
  `dorothea`, `decoupleR`), resolving a CRAN check gap.
- **Data import/export** (`R/16_data_import.R`): `import_data()` /
  `export_data()` support CSV, TSV, Excel, RDS, Parquet, and (import only)
  live database connections via `DBI`, with format auto-detected from the
  file extension -- so moving from simulated/example data to a real file
  requires changing only the loading call, not any analysis code.
- **Data validation** (`R/17_data_validation.R`): `validate_matrix()`,
  `validate_clinical()`, and `validate_cohort_alignment()` check
  structure, value ranges, missingness, and cross-object subject
  alignment before analysis, with specific, actionable error messages
  (offending row/column/subject IDs named directly) rather than a
  downstream failure or a silently wrong result.
- **Data dictionary and templates** (`R/18_data_dictionary.R`):
  `data_dictionary()` documents every expected field, type, unit, and
  coding scheme across all pipeline inputs; `create_data_templates()`
  writes ready-to-fill CSV templates for real-data mapping.
- **Bug fix**: `simulate_trace_asthma_cohort()`/`simulate_methylation()`
  crashed ("subscript out of bounds") when `n_causal_cpgs` exceeded
  `n_cpgs` (e.g. small test datasets using the function's own default
  `n_causal_cpgs`) -- found via `R CMD check`'s example run, fixed at the
  root cause (index capping in the caller) and defensively (bounds
  checking with an informative warning/error inside `simulate_methylation()`
  itself for any future caller).
- New vignette section: "Moving from simulated to real data" in
  `vignette("simulating-data")`, covering the data dictionary, templates,
  import/export, and validation workflow end to end.

# traceAsthma 0.2.0

New features:

- **Adaptive pipeline**: `detect_available_data()` and
  `run_trace_asthma_auto()` inspect whichever data modalities a user
  actually supplies (methylation + expression + asthma status is the only
  hard requirement) and build the best-supported model, automatically
  skipping the MPRS layer without genotype data and the TRACE layer
  without a regulon database, and using only whichever clinical
  covariates (age/sex/bmi/smoking/ancestry PCs) are present.
- **Publication-quality figures** (`R/11_figures.R`): `plot_nested_auc()`,
  `plot_calibration()`, `plot_roc_curves()`, `plot_trace_weights()`,
  `plot_mediation_path()`, built on `open_nature_device()` /
  `with_nature_device()`. Default 600 dpi, physical sizing in mm matching
  Nature single/double-column widths, sans-serif typography, a
  colorblind-safe categorical palette (`nature_palette`), and direct,
  non-generic labeling throughout (gene symbols, named model predictors --
  never "TF1"/"Model A"). PNG, TIFF (LZW), and PDF (vector) output
  supported, all via base `grDevices` (no new hard dependency).
- **Publication-quality tables** (`R/12_tables.R`):
  `format_nested_model_table()`, `format_coefficient_table()`,
  `format_mediation_table()`, each returning a `traceAsthma_table` object
  with fully spelled-out column headers and a mandatory footnote
  explaining every abbreviation/statistical convention. `export_table()`
  writes CSV + companion footnote file (dependency-free default) or
  `.xlsx` (if `writexl` is installed), with guaranteed UTF-8 fidelity via
  `useBytes = TRUE` regardless of system locale.
- **`list_dependencies()`**: prints every optional package `traceAsthma`
  uses, which pipeline stage needs it, and confirms which stages need no
  dependencies beyond base R (QC, eQTM discovery, regulon
  filtering/enrichment, the final logistic risk model, and all
  table/figure export functions).
- **`interpret_clinical_relevance()`**: translates a fitted risk model
  into plain-language clinical relevance -- absolute predicted risk at
  low/median/high values of a score (not just an odds ratio), a
  calibration trustworthiness statement (from `validate_model()`), and a
  mechanistic mediation statement naming the modifiable pathway (from
  `fit_mediation_sem()`), assembled into one narrative paragraph plus the
  underlying table. Base R only, no dependencies.
- **Realistic data simulation** (`R/15_simulate.R`):
  `simulate_trace_asthma_cohort()` and its component functions
  (`simulate_genotypes()`, `simulate_methylation()`, `simulate_expression()`,
  `simulate_regulon()`, `simulate_clinical()`) generate synthetic cohorts
  matching real distributional properties -- Hardy-Weinberg genotype
  sampling with a realistic MAF spectrum, bimodal Beta-mixture methylation
  matching real EPIC/450K arrays, negative-binomial RNA-seq counts matching
  the DESeq2/edgeR generative model, and realistic clinical covariate
  distributions -- with a known, injected causal architecture
  (SNP -> methylation -> eQTM -> TF activity -> TRACE -> asthma) and a
  returned `ground_truth` element for evaluating discovery precision/recall
  and effect-size recovery. Verified end-to-end: `discover_eqtm()` run on
  simulated data correctly ranked truly-causal CpG-gene pairs far ahead of
  null pairs (median P = 0.011 vs. 0.57 in a representative n=300 run), and
  `fit_risk_model()` recovered injected true effect sizes (0.40/0.50/0.70)
  as approximately 0.37/0.61/0.79 with correct significance patterns.

Bug fixes:

- `export_table()` previously could mangle non-ASCII characters (Greek
  letters, en-dashes, arrows) on non-UTF-8 system locales via
  `write.csv()`'s locale-dependent re-encoding; rewritten to write raw
  UTF-8 bytes directly via `writeLines(..., useBytes = TRUE)`, verified
  correct independent of session locale.

# traceAsthma 0.1.0

Initial release. Implements the full TRACE-Asthma pipeline: QC,
multi-trait PRS, eQTM discovery, regulatory network / TF activity
inference, TRACE score, mechanistic mediation SEM, integrated risk model,
validation, and clinician-facing prediction (including a Shiny
point-of-care app).

Known limitations carried forward to 0.2.0 (see
`vignette("clinical-deployment-notes")`):

- `NAMESPACE` and most `man/*.Rd` pages for lower-level functions were
  authored by hand rather than generated by `roxygen2::roxygenize()`
  (unavailable in the authoring environment); regenerate with
  `devtools::document()` before release. User-facing entry points
  (`run_trace_asthma_pipeline()`, `run_trace_asthma_auto()`,
  `predict_asthma_risk()`, `list_dependencies()`, all `plot_*()`/
  `format_*_table()` functions, `example_cohort`) do have hand-written
  `.Rd` pages.
- `glmnet`/`mediation`/`lavaan`/`viper`/`decoupleR`-dependent code paths
  were verified by code review and by testing their surrounding logic,
  not executed end-to-end (no CRAN/Bioconductor network access in the
  authoring sandbox). Figure/table export and the adaptive-pipeline
  detection logic *were* executed end-to-end against the bundled example
  dataset and verified correct, including at actual 600 dpi.
- No real-world validation cohort has been run through the pipeline.
