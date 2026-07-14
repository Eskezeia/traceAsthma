# traceAsthma 0.3.6

- **Bug fix (found via real-world community testing)**: `infer_tf_activity()`'s
  `"decoupleR"` fallback previously called `decoupleR::run_viper()`, which
  -- despite the name suggesting it was decoupleR's own implementation --
  internally requires the `viper` Bioconductor package to be installed as
  a backend. This meant the documented "fallback if viper is unavailable"
  behavior did not actually work: a user without `viper` installed but
  with `decoupleR` installed would still hit `"there is no package called
  'viper'"`, surfacing as a confusing `run_package_diagnostics()` failure
  several stages downstream. Fixed by switching the decoupleR path to
  `decoupleR::run_ulm()`, decoupleR's own univariate linear model method,
  which is genuinely independent of the `viper` package. `method = "viper"`
  now also automatically falls back to `"decoupleR"` with an informative
  message if `viper` is requested but not installed, rather than
  proceeding to a `viper`-dependent decoupleR call.
- **Diagnostics robustness**: `run_package_diagnostics()`'s SKIPPED-vs-FAIL
  classification now also recognizes generic "there is no package called"
  and "could not find function" errors as SKIPPED (missing optional
  dependency), not just errors matching `requirePkg()`'s exact message
  format -- defense in depth against this class of bug recurring in any
  Suggests-dependent code path, found and fixed as a direct consequence of
  the bug above.

# traceAsthma 0.3.5

- **PDF user manual bundled with the package**: `inst/doc-manual/traceAsthma_Manual.pdf`
  covers purpose, objective, novelty, data description (including the
  simulator's verified fidelity), every method's formula alongside its
  implementing function, a step-by-step usage guide, and a conclusion.
  Access via
  `system.file("doc-manual", "traceAsthma_Manual.pdf", package = "traceAsthma")`.
  Source is `data-raw/manual/traceAsthma_Manual.md`, rendered via
  `pandoc --pdf-engine=xelatex`; regenerate after substantive changes to
  the package's methods or scope.

# traceAsthma 0.3.4

- **New: `run_package_diagnostics()`** -- the recommended first command for
  any new install, a community member verifying a bug report, or checking
  a new release. Runs the full pipeline (QC/validation, eQTM discovery
  with ground-truth recovery checking, mediation, eQTM scoring, TF
  regulon discovery and activity inference with the correct
  full-regulon-for-inference pattern, TRACE scoring, the integrated risk
  model, and validation) against the bundled `simulated_benchmark_cohort`
  (or any user-supplied cohort), and prints a clear PASS/FAIL/SKIPPED
  report per stage. Stages that fail purely because an optional package
  is missing -- including *downstream* stages that only fail as a
  cascading consequence of an earlier skip -- are correctly reported as
  `SKIPPED`, not `FAIL`, so the final summary never wrongly tells a user
  to file a bug report for an uninstalled dependency. Verified: 10 PASS /
  0 FAIL / 5 SKIPPED in this package's own dependency-limited CI
  environment, with all 5 skips correctly and specifically attributed to
  `glmnet`/`decoupleR` being absent.
- This function directly institutionalizes a community-reported diagnostic
  script (see v0.3.2/0.3.3 entries below) so future users get the
  corrected, hardened version by default rather than re-deriving it.

# traceAsthma 0.3.3

- **New vignette, `vignette("statistical-methods")`**: a single reference
  pairing every model formula in the pipeline (trait PRS, MPRS, eQTM
  linear model, mediation indirect effect, eQTM score, hypergeometric TF
  enrichment, VIPER TF activity, TRACE score, mediation SEM path
  equations, integrated logistic risk model, AUROC/calibration/NRI,
  clinical-relevance absolute risk) with the exact function that
  implements it and a runnable step-by-step worked example against the
  bundled `simulated_benchmark_cohort`, plus a one-page quick-reference
  formula table at the end.
- **Bug fix**: `test_regulon_enrichment()` crashed with a confusing
  "incorrect number of dimensions" error when `target_genes` was empty or
  no candidate TF's regulon overlapped `background_genes` (e.g. when an
  upstream trans-eQTM discovery step legitimately found zero significant
  associations, which is realistic behavior at modest sample sizes, not
  necessarily a bug). Found while verifying the new vignette's worked
  example against `simulated_benchmark_cohort` at a strict FDR threshold.
  Now returns a well-formed, zero-row result with the correct columns and
  an informative warning explaining the likely cause, instead of
  erroring.

# traceAsthma 0.3.2

- **Bug fix / robustness**: `infer_tf_activity()` now validates regulon
  size (per TF, after intersecting with `rownames(expression)`) *before*
  delegating to viper/decoupleR, and fails with a specific, actionable
  `traceAsthma`-level error naming exactly which TF(s) are short and by
  how much, rather than surfacing decoupleR's internal
  "Network is empty after intersecting..." message. If some but not all
  supplied TFs are eligible, the ineligible ones are dropped with an
  informative warning rather than the whole call failing. This was found
  via real community usage: a diagnostic pipeline run correctly built a
  target-universe-restricted regulon for TF *discovery*
  (`build_regulons(..., target_universe = <trans-eQTM genes>)`), then
  mistakenly reused that same restricted regulon for activity
  *inference*, where each candidate TF's regulon should be its full
  target set. The new error message names this exact failure mode and
  gives the corrected code pattern directly in the message.

# traceAsthma 0.3.1

- **Bundled community benchmark dataset**: `simulated_benchmark_cohort`
  ships with the package -- a fixed-seed (`20260714`), realistic-scale
  (300 subjects, 500 CpGs, 200 genes, 10 TFs) output of
  `simulate_trace_asthma_cohort()`, identical for every user/install.
  Intended as a shared reference for comparing pipeline behavior across
  systems/versions and for reporting reproducible issues. Verified: loads
  correctly, and `discover_eqtm()` run on it correctly separates causal
  from null CpG-gene pairs (median P ~0 vs. 0.56 in a representative run).
  Regenerate via `data-raw/make_simulated_benchmark.R` if the simulator
  changes; documented in `vignette("simulating-data")`.


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
