# Generates a bundled, realistic-scale simulated cohort with a known,
# fixed random seed -- so every user who installs traceAsthma gets the
# *exact same* dataset to test against, and results are directly
# comparable and reportable (e.g. in GitHub issues: "on
# simulated_benchmark_cohort, discover_eqtm() gives me X").
#
# This is deliberately larger and more realistic than `example_cohort`
# (which is a tiny 5-patient toy dataset for instant-run examples/tests),
# but still small enough to keep the package's installed size reasonable.

devtools::load_all(".")  # use the in-development simulate_* functions

simulated_benchmark_cohort <- simulate_trace_asthma_cohort(
  n_subjects       = 300,
  n_snps           = 30,
  n_cpgs           = 500,
  n_genes          = 200,
  n_tfs            = 10,
  n_causal_cpgs    = 20,
  n_causal_genes   = 12,
  n_active_tfs     = 2,
  true_effects     = list(mprs = 0.4, eqtm = 0.5, trace = 0.7, age = 0.01, sex = 0.3),
  prevalence       = 0.35,
  ancestry_groups  = 2,     # a mild multi-ancestry structure, for realism
  seed             = 20260714  # fixed: same dataset for every user, every install
)

# sanity checks before shipping this as package data
stopifnot(
  ncol(simulated_benchmark_cohort$methylation) == 300,
  nrow(simulated_benchmark_cohort$methylation) == 500,
  nrow(simulated_benchmark_cohort$expression) == 200,
  length(unique(simulated_benchmark_cohort$regulon_db$tf)) == 10,
  all(simulated_benchmark_cohort$clinical$asthma %in% c(0, 1)),
  length(simulated_benchmark_cohort$ground_truth$causal_cpgs) > 0
)

cat("Realized asthma prevalence:", round(mean(simulated_benchmark_cohort$clinical$asthma), 3), "\n")
cat("Object size:", format(object.size(simulated_benchmark_cohort), units = "MB"), "\n")

usethis::use_data(simulated_benchmark_cohort, overwrite = TRUE, compress = "xz")
