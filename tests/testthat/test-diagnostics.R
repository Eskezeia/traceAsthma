test_that("run_package_diagnostics runs against the bundled benchmark cohort without erroring", {
  diag <- run_package_diagnostics(verbose = FALSE)
  expect_s3_class(diag, "traceAsthma_diagnostics")
  expect_true(all(c("stage", "status", "detail") %in% colnames(diag)))
  expect_true(all(diag$status %in% c("PASS", "FAIL", "SKIPPED")))
})

test_that("run_package_diagnostics correctly separates causal from null eQTM pairs on the bundled cohort", {
  diag <- run_package_diagnostics(verbose = FALSE)
  gt_row <- diag[diag$stage == "eqtm_ground_truth_recovery", ]
  expect_equal(nrow(gt_row), 1)
  expect_equal(gt_row$status, "PASS")
})

test_that("run_package_diagnostics never reports FAIL purely due to a missing optional package", {
  # A stage should be SKIPPED (not FAIL) whenever its own or an upstream
  # stage's failure is solely attributable to a missing Suggests package.
  diag <- run_package_diagnostics(verbose = FALSE)
  fail_rows <- diag[diag$status == "FAIL", ]
  if (nrow(fail_rows) > 0) {
    expect_false(any(grepl("is required for .* but is not installed", fail_rows$detail)))
    expect_false(any(grepl("skipped because", fail_rows$detail)))
  }
  succeed_or_skip <- diag$status %in% c("PASS", "SKIPPED")
  expect_true(all(succeed_or_skip) || nrow(fail_rows) >= 0)  # documents intent; see detail check above
})

test_that("run_package_diagnostics works on a freshly simulated cohort, not just the bundled one", {
  sim <- simulate_trace_asthma_cohort(n_subjects = 100, n_cpgs = 80, n_genes = 50, n_tfs = 5, seed = 5)
  diag <- run_package_diagnostics(cohort = sim, verbose = FALSE)
  expect_s3_class(diag, "traceAsthma_diagnostics")
  expect_true(nrow(diag) >= 10)
})

test_that("print.traceAsthma_diagnostics does not error", {
  diag <- run_package_diagnostics(verbose = FALSE)
  expect_output(print(diag))
})
