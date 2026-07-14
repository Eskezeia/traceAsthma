---
title: "traceAsthma: User Manual"
subtitle: "A Mechanistic Multi-Omics Asthma Risk Model — R Package Reference"
author: "Eskezeia Yihunie Dessie, PhD --- Division of Pulmonary, Critical Care, Sleep & Occupational Medicine, Indiana University School of Medicine"
date: "Version 0.3.4 — July 2026"
toc: true
toc-depth: 3
numbersections: true
geometry: margin=1in
fontsize: 11pt
colorlinks: true
linkcolor: blue
urlcolor: blue
---

\newpage

# Purpose

`traceAsthma` is an R package implementing **TRACE-Asthma**
(Transcriptional Regulatory Activity and CpG-methylation Epigenetic
Asthma score), an end-to-end computational framework that integrates
genetic, epigenetic, and transcriptional regulatory information into a
single, mechanistically interpretable asthma risk model. The package
exists to make a specific scientific hypothesis — that inherited genetic
risk for asthma is transmitted, in substantial part, through
methylation-dependent transcriptional regulation rather than through
fixed DNA sequence alone — reproducible, auditable, and directly usable
by other researchers, rather than existing only as a private analysis
script behind a single manuscript.

Concretely, `traceAsthma` provides:

- A documented, versioned implementation of every analytical step in the
  TRACE-Asthma pipeline, from raw quality control through a
  clinician-facing risk prediction interface.
- An **adaptive** entry point that works with whatever data modalities a
  user actually has (methylation and expression alone, with or without
  genotype, BMI, or smoking status), rather than requiring the full data
  specification.
- Realistic **simulated data**, bundled with the package, so that the
  entire pipeline can be exercised and verified before any real patient
  data is available or shared.
- Publication-standard **figures and tables**, generated directly from
  analysis objects, matching Nature-family journal conventions.
- A single-command **self-diagnostic** (`run_package_diagnostics()`) that
  verifies a fresh installation end to end.

# Objective

The specific objectives of the package are to:

1. **Operationalize the TRACE-Asthma mechanistic architecture** —
   $$\text{SNP} \rightarrow \text{MPRS} \rightarrow \text{CpG methylation} \rightarrow \text{eQTM score} \rightarrow \text{TF activity} \rightarrow \text{TRACE score} \rightarrow \text{Asthma risk}$$
   — as a sequence of independently callable, independently testable R
   functions, so that each stage's assumptions and outputs are inspectable
   rather than hidden inside a single monolithic analysis.
2. **Quantify, not merely assume, mediation.** Provide a structural
   equation modeling layer (`fit_mediation_sem()`) that formally
   decomposes the total effect of inherited genetic risk into a direct
   effect and an indirect effect transmitted through the epigenetic and
   transcriptional layers, with bootstrap confidence intervals.
3. **Discover, rather than presuppose, candidate regulatory mediators.**
   Identify candidate transcription factors from a study's own trans-eQTM
   target genes via regulon enrichment testing
   (`test_regulon_enrichment()`), rather than selecting well-known
   literature TFs a priori — a design choice intended to reduce
   confirmation bias in the TF-discovery step.
4. **Degrade gracefully to real-world data availability.** Provide
   `run_trace_asthma_auto()`, which inspects whichever data modalities a
   user supplies and builds the best-supported model automatically,
   because most real-world users do not have the full genotype +
   methylation + expression + complete clinical panel this framework was
   originally specified against.
5. **Translate statistical output into actionable, plain-language
   information.** Provide `interpret_clinical_relevance()`, which
   converts a fitted model into absolute predicted risk at defined score
   percentiles, a calibration trustworthiness statement, and a
   mechanistic mediation statement — closing the gap between "the
   coefficient is significant" and "here is what this means for a given
   patient."
6. **Support reproducible, regulatory-quality workflows.** Provide data
   validation, a documented data dictionary, multi-format data import, and
   a realistic data simulator, so that the transition from simulated to
   real clinical data requires no change to analysis code — only a change
   in which files are loaded.

# Novelty

As a matter of intellectual honesty, `traceAsthma`'s novelty should be
understood at the level of **integration and translation**, not at the
level of individual statistical methods. Every individual algorithm
wrapped by this package — elastic-net regression, structural equation
modeling, hypergeometric enrichment testing, VIPER-based transcription
factor activity inference — is independently established and available in
other R and Bioconductor packages.

The package's genuine contributions are:

1. **A single, mechanistically ordered pipeline**, from SNP to
   transcriptional regulatory activity to disease risk, implemented as one
   coherent, testable software artifact rather than a collection of
   independently maintained scripts. To the authors' knowledge, no
   existing R package combines multi-trait polygenic scoring,
   mediation-filtered eQTM discovery, data-driven transcription factor
   discovery, structural equation modeling of the full mediation pathway,
   and a deployable clinical prediction interface in one package.
2. **A data-driven, not literature-driven, transcription factor discovery
   step.** `test_regulon_enrichment()` nominates candidate transcription
   factors from a study's own trans-eQTM target genes, which is a more
   defensible design for a hypothesis-generating pipeline than selecting
   well-known disease-associated transcription factors a priori.
3. **An adaptive pipeline architecture** (`run_trace_asthma_auto()`,
   `detect_available_data()`) that is, to the authors' knowledge, unusual
   among domain-specific genomics pipeline packages, most of which require
   the full specified data set and leave graceful degradation to the user.
4. **An automated clinical-translation layer**
   (`interpret_clinical_relevance()`) that converts a fitted statistical
   model directly into a plain-language absolute-risk, calibration, and
   mechanistic narrative, rather than leaving that translation as a manual
   reporting exercise for every user.
5. **Publication output as a first-class function**, not an
   afterthought — Nature-journal-convention figures (600 dpi, physical
   sizing, colorblind-safe palettes, direct non-generic labeling) and
   fully footnoted tables are generated directly by the package.

Readers and reviewers should evaluate the package's scientific claims —
in particular, the magnitude of genetic-effect mediation through the
eQTM/TRACE pathway — on the basis of results obtained from real cohort
data, not on the basis of the software architecture described here.

\newpage

# Data Description

## Bundled data (no real, proprietary, or patient-identifiable data)

`traceAsthma` includes no proprietary, confidential, or
patient-identifiable data of any kind. Three data objects ship with or
are generated by the package, all fully synthetic:

### `example_cohort`

A fixed, minimal dataset (5 subjects; 2 SNPs; 2 CpGs; 4 genes; a 3-TF
regulon) used for instant-run function examples and unit tests throughout
the package documentation. Not sized for meaningful statistical analysis.

### `simulated_benchmark_cohort`

A fixed-seed (`20260714`), realistic-scale simulated cohort bundled with
the package: 300 subjects, 500 CpGs, 200 genes, 10 transcription factors,
and a mild two-cluster ancestry structure. Because the seed is fixed,
every user who installs the package receives an *identical* dataset,
making it suitable as a shared community reference for verifying
installations, comparing pipeline behavior across systems or package
versions, and reporting reproducible issues.

### `simulate_trace_asthma_cohort()`

An on-demand simulator generating a new cohort, at any user-specified
size and random seed, with:

- **Genotypes** sampled under Hardy-Weinberg equilibrium with a realistic
  minor allele frequency spectrum (Uniform(0.05, 0.5)).
- **DNA methylation** drawn from a bimodal Beta-mixture distribution
  matching the marginal distribution of real Illumina EPIC/450K array
  data (most CpGs concentrated near 0 or 1; a minority variably
  methylated).
- **Gene expression** simulated as negative-binomial RNA-sequencing
  counts with sample-level library-size factors, matching the generative
  model underlying DESeq2/edgeR.
- **Clinical covariates** (age, sex, BMI, smoking status, ancestry
  principal components) with realistic marginal distributions and
  cross-correlations.
- A **known, injected causal architecture**
  (SNP $\rightarrow$ methylation $\rightarrow$ eQTM $\rightarrow$ TF
  activity $\rightarrow$ TRACE $\rightarrow$ asthma status), returned as
  a `ground_truth` object, enabling direct evaluation of discovery
  precision/recall and effect-size recovery against a known answer.

The simulator's fidelity was verified empirically during development: at
realistic sample sizes, `discover_eqtm()` run on simulated data correctly
ranked injected causal CpG-gene pairs far ahead of null pairs (median
$P \approx 0$ for causal pairs versus $P \approx 0.56$ for null pairs in a
representative run), and `fit_risk_model()` recovered injected true
effect sizes (0.40, 0.50, 0.70 for the MPRS, eQTM, and TRACE layers
respectively) as approximately 0.37, 0.61, and 0.79, each correctly and
strongly statistically significant.

## Real-data requirements

For deployment on real cohort data, the package provides:

- **`data_dictionary()`** — a structured reference documenting every
  expected data object, field, type, unit, and coding scheme (e.g.
  methylation beta values in $[0,1]$; asthma status coded 0/1; sex coded
  0/1 with automatic detection of F/M text labels).
- **`create_data_templates()`** — generates ready-to-fill CSV templates
  matching the data dictionary, for preparing real data or communicating
  the required format to collaborators.
- **`import_data()` / `export_data()`** — read and write CSV, TSV, Excel,
  RDS, and Parquet files, and read from live database connections via
  `DBI`, so that the transition from simulated to real data requires
  changing only the data-loading call.
- **`validate_matrix()`, `validate_clinical()`,
  `validate_cohort_alignment()`** — structural, value-range,
  missingness, and cross-object subject-alignment checks, producing
  specific, actionable error messages (naming offending rows, columns, or
  subjects directly) before any downstream analysis is attempted.

No component of the analysis engine (the `R/*.R` source files implementing
each pipeline stage) makes any reference to the specific bundled example
data; every function operates purely on the shapes and types documented in
the data dictionary, so real data of the correct shape is a drop-in
replacement for simulated data throughout.

\newpage

# Package Methods

This section documents every analytical stage, its statistical formula,
and the function implementing it. All formulas are also available with
runnable worked examples in `vignette("statistical-methods")`.

## Quality control

| Check | Rule | Function |
|---|---|---|
| Variant call rate | retain if $\geq 0.95$ | `qc_genotype()` |
| Hardy-Weinberg equilibrium | retain if $P_{HWE} \geq 10^{-6}$ | `qc_genotype()` |
| Minor allele frequency | retain if $\geq 0.01$ | `qc_genotype()` |
| Methylation probe detection | retain if fail rate $\leq 0.05$ | `qc_methylation()` |
| Gene expression | retain if count $\geq 6$ in $\geq 3$ samples | `qc_rnaseq()` |

## Multi-trait polygenic risk score (MPRS)

Trait-specific polygenic score, for individual $i$ and trait $t$:
$$PRS_{i,t} = \sum_{j=1}^{m_t} \beta_{j,t} G_{ij}$$
where $G_{ij}$ is genotype dosage and $\beta_{j,t}$ is the posterior SNP
effect size for trait $t$ (obtained upstream, e.g. via PRS-CSx/SBayesRC).
Implemented in `compute_prs()`.

Multi-trait combination, via elastic-net regression against asthma
status:
$$MPRS_i = \sum_t w_t \, PRS_{i,t}$$
Implemented in `compute_mprs()`.

## Cis/trans eQTM discovery

For each candidate CpG-gene pair:
$$Y_g = \alpha_0 + \alpha_1 M_c + \mathbf{C}\gamma + \epsilon$$
where $Y_g$ is gene expression, $M_c$ is CpG methylation, and $\mathbf{C}$
is a covariate matrix (age, sex, ancestry principal components, cell-type
composition, ...). Pairs are classified as cis if within 1 Mb of the
gene's transcription start site, else trans; Benjamini-Hochberg false
discovery rate control is applied separately within each class given
their markedly different multiple-testing burdens. Implemented in
`discover_eqtm()`.

**Asthma-relevant filtering** requires independent CpG- and gene-level
association with asthma status:
$$\text{logit}(P(\text{Asthma})) = \beta_0 + \beta_1 \text{CpG} + \mathbf{C}\gamma$$
Implemented in `filter_asthma_relevant()`.

**Mediation testing** estimates the indirect effect along the pathway
CpG $\rightarrow$ Expression $\rightarrow$ Asthma:
$$IE = a \times b$$
where $a$ is the CpG-expression effect and $b$ is the expression-asthma
effect conditional on CpG, with bootstrap confidence intervals.
Implemented in `test_mediation()`.

**eQTM score construction** applies correlation pruning
($|r| > 0.90$ collapsed to a single representative) followed by
elastic-net feature selection:
$$eQTM_i = \sum_{j=1}^{k} w_j M_{ij}$$
Implemented in `compute_eqtm_score()`.

## Regulatory network and transcription factor activity

**Data-driven candidate transcription factor discovery** tests each
candidate TF's regulon for hypergeometric enrichment among trans-eQTM
target genes:
$$P(X \geq n_{\text{overlap}}) = 1 - F_{\text{hyper}}\left(n_{\text{overlap}} - 1;\ n_{\text{target}},\ n_{\text{bg}} - n_{\text{target}},\ \text{regulon size}\right)$$
Implemented in `build_regulons()` and `test_regulon_enrichment()`.

**Transcription factor activity inference** estimates per-subject TF
activity as the normalized enrichment score of the TF's full regulon
among differentially expressed genes (via VIPER or decoupleR), rather
than the TF's own transcript abundance:
$$A_t = NES_t$$
Implemented in `infer_tf_activity()`. **Important implementation note**:
activity inference must use each candidate TF's complete regulon, not the
`target_universe`-restricted regulon used for the discovery step above;
reusing the restricted regulon for inference is a documented failure mode
(see package `NEWS.md`, versions 0.3.2–0.3.3) that the current
implementation validates against and reports with an informative,
actionable error message.

## TRACE score

Constructed via a four-step procedure: univariate false-discovery-rate
screening, correlation pruning, bootstrap stability selection (retaining
TFs selected in a high proportion of resamples), and final elastic-net
weighting on the stable TF set:
$$TRACE_i = \sum_{t=1}^{m} \theta_t A_{it}$$
Implemented in `compute_trace_score()`.

## Mechanistic mediation (structural equation model)

$$eQTM = a_0 + a_1 MPRS + \epsilon_1$$
$$TRACE = b_0 + b_1 eQTM + \epsilon_2$$
$$Asthma = c_0 + c_1 MPRS + c_2 eQTM + c_3 TRACE + \epsilon_3$$

Direct effect $= c_1$; indirect (mediated) effect
$= a_1 \times b_1 \times c_3$; total effect $=$ direct $+$ indirect, each
with bootstrap confidence intervals. Implemented in
`fit_mediation_sem()`.

## Integrated asthma risk model

$$\text{logit}(P(\text{Asthma})) = \beta_0 + \beta_1 MPRS + \beta_2 eQTM + \beta_3 TRACE + \beta_4 Age + \beta_5 Sex + \beta_6 BMI + \beta_7 Smoking$$

Any subset of predictors may be supplied, accommodating sites without
genotype or complete clinical data. Implemented in `fit_risk_model()`,
with nested model comparison (`compare_nested_models()`) and
cross-validated performance estimation (`nested_cv()`).

## Validation

$$\text{AUROC} = P(\hat{p}_{\text{case}} > \hat{p}_{\text{control}})$$
$$\text{Calibration slope} = \beta_1 \text{ in } \text{logit}(\text{observed}) = \beta_0 + \beta_1 \, \text{logit}(\hat{p})$$
$$\text{NRI} = \left[P(\hat{p}\uparrow \mid \text{case}) - P(\hat{p}\downarrow \mid \text{case})\right] + \left[P(\hat{p}\downarrow \mid \text{control}) - P(\hat{p}\uparrow \mid \text{control})\right]$$

Implemented in `validate_model()` (AUROC, AUPRC, calibration,
Hosmer-Lemeshow test, NRI, IDI) and `decision_curve()` (net benefit
across risk thresholds).

## Clinical interpretation

Absolute predicted risk at defined percentiles of a chosen score, holding
other covariates at their sample median:
$$\Delta \hat{p} = \hat{p}(\text{score}_{90th}) - \hat{p}(\text{score}_{10th})$$
combined with a calibration trustworthiness statement and, where
available, a plain-language mediation statement naming the modifiable
regulatory pathway. Implemented in `interpret_clinical_relevance()`.

\newpage

# Step-by-Step Usage Guide

## 1. Installation

```r
install.packages("remotes")
remotes::install_github("Eskezeia/traceAsthma")
library(traceAsthma)
```

## 2. Verify the installation

```r
run_package_diagnostics()
```

Runs the full pipeline against the bundled `simulated_benchmark_cohort`
and reports a PASS/FAIL/SKIPPED status per stage; stages requiring an
uninstalled optional package are reported as SKIPPED with the exact
install command, not as a failure.

## 3. Load data

```r
# Bundled fixed-seed benchmark
data(simulated_benchmark_cohort)
sim <- simulated_benchmark_cohort

# Or generate a fresh simulated cohort at any size/seed
sim <- simulate_trace_asthma_cohort(n_subjects = 300, n_cpgs = 500,
                                     n_genes = 200, n_tfs = 8, seed = 42)

# Or import real data
methylation <- import_data("real_methylation.csv")
clinical    <- import_data("real_clinical.csv", as_matrix = FALSE)
validate_matrix(methylation, type = "methylation")
validate_clinical(clinical)
```

## 4. Run the pipeline

**Adaptive (recommended for most users)** — automatically detects which
layers your data supports:

```r
result <- run_trace_asthma_auto(
  methylation = sim$methylation, expression = sim$expression,
  clinical = sim$clinical, eqtm_pairs = sim$eqtm_pairs,
  regulon_db = sim$regulon_db
)
```

**Fully specified** — every layer explicit:

```r
result <- run_trace_asthma_pipeline(
  genotype_dosage = sim$dosage, methylation = sim$methylation,
  expression = sim$expression, eqtm_pairs = sim$eqtm_pairs,
  regulon_db = sim$regulon_db, clinical = sim$clinical,
  covariates = sim$covariates, seed = 42
)
```

## 5. Validate and interpret

```r
pred <- predict(result$risk_model, type = "response")
val <- validate_model(result$model_data$asthma, pred)
print(interpret_clinical_relevance(result$risk_model, result$model_data, validation = val))
```

## 6. Generate publication output

```r
export_table(format_coefficient_table(result$risk_model), "table2_coefficients.csv")
plot_trace_weights(result$trace_fit, path = "figure3_trace_weights.tiff", dpi = 600)
```

## 7. Deploy for prediction on new patients

```r
model <- build_deployable_model(
  eqtm_weights = result$eqtm_fit$weights, tf_weights = result$trace_fit$weights,
  regulons = result$deployable_model$regulons, risk_model = result$risk_model
)
predict_asthma_risk(model, methylation = new_patient_methylation,
                     expression = new_patient_expression,
                     clinical = data.frame(age = 34, sex = 0))
```

## Further reading

```r
vignette("trace-asthma-pipeline")     # full walkthrough
vignette("statistical-methods")       # every formula, step by step
vignette("simulating-data")           # simulation and real-data transition
vignette("clinical-deployment-notes") # validation, regulatory, publishing
```

\newpage

# Conclusion

`traceAsthma` packages the TRACE-Asthma mechanistic framework into a
single, documented, versioned, and testable R package, spanning quality
control, multi-layer molecular risk scoring, formal mediation modeling,
prediction, validation, and publication-ready reporting. Its contribution
is architectural and translational: a mechanistically ordered pipeline
with a data-driven transcription factor discovery step, an adaptive
interface that degrades gracefully to real-world data availability, and
an automated bridge from statistical output to plain-language clinical
interpretation — rather than any single novel statistical algorithm.

Every function in the package operates identically on simulated and real
data, and the package ships with a realistic data simulator, a bundled
fixed-seed community benchmark dataset, and a one-command self-diagnostic,
so that the full pipeline can be verified end to end — by the author, by
collaborators, and by any interested member of the community — before it
is ever run on real patient data. Reported associations, mediated effect
sizes, and predictive performance should be evaluated on the basis of
results obtained from real cohort data; the software itself is offered as
a transparent, reproducible instrument for obtaining and reporting those
results, not as a substitute for them.

The package is under active development. Bug reports, feature requests,
and contributions are welcome at
`https://github.com/Eskezeia/traceAsthma`.
