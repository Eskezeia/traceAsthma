test_that("simulate_genotypes produces valid dosages and MAF spectrum", {
  sim <- simulate_genotypes(50, 10, seed = 1)
  expect_true(all(sim$dosage %in% c(0, 1, 2)))
  expect_equal(dim(sim$dosage), c(50, 10))
  expect_true(all(sim$maf >= 0.05 & sim$maf <= 0.5))
})

test_that("simulate_methylation produces bimodal beta values in [0,1]", {
  sim <- simulate_methylation(50, 100, seed = 1)
  expect_true(all(sim$beta >= 0 & sim$beta <= 1))
  expect_equal(dim(sim$beta), c(100, 50))
  # bimodal: more mass near extremes than in the middle third
  vals <- as.vector(sim$beta)
  n_extreme <- sum(vals < 0.2 | vals > 0.8)
  n_middle <- sum(vals >= 0.4 & vals <= 0.6)
  expect_true(n_extreme > n_middle)
})

test_that("simulate_methylation injects a detectable case-control shift", {
  status <- rep(c(0, 1), each = 50)
  sim <- simulate_methylation(100, 20, case_control_idx = 1:5,
                               case_control_status = status, case_control_delta = 2, seed = 1)
  # causal CpGs should show a clear mean difference between groups
  diffs <- rowMeans(sim$beta[1:5, status == 1]) - rowMeans(sim$beta[1:5, status == 0])
  expect_true(all(diffs > 0))
  expect_length(sim$causal_cpgs, 5)
})

test_that("simulate_expression produces non-negative integer-like counts", {
  sim <- simulate_expression(30, 40, seed = 1)
  expect_equal(dim(sim$counts), c(40, 30))
  expect_true(all(sim$counts >= 0))
  expect_true(all(sim$counts == round(sim$counts)))
})

test_that("simulate_regulon produces a valid TF-target table", {
  genes <- sprintf("GENE%03d", 1:50)
  reg <- simulate_regulon(genes, n_tfs = 5, seed = 1)
  expect_true(all(c("tf", "target", "mor", "confidence") %in% colnames(reg)))
  expect_true(all(reg$mor %in% c(1, -1)))
  expect_true(all(reg$confidence %in% c("A", "B", "C")))
  expect_equal(length(unique(reg$tf)), 5)
})

test_that("simulate_clinical produces plausible marginal distributions", {
  cl <- simulate_clinical(200, seed = 1)
  expect_true(all(cl$age >= 5 & cl$age <= 90))
  expect_true(all(cl$sex %in% c(0, 1)))
  expect_true(all(cl$bmi >= 14 & cl$bmi <= 55))
  expect_true(all(levels(cl$smoking) == c("never", "former", "current")))
})

test_that("simulate_trace_asthma_cohort returns a complete, consistent cohort", {
  sim <- simulate_trace_asthma_cohort(n_subjects = 60, n_snps = 10, n_cpgs = 50,
                                       n_genes = 30, n_tfs = 4, seed = 1)
  expect_equal(ncol(sim$methylation), 60)
  expect_equal(ncol(sim$expression), 60)
  expect_equal(nrow(sim$clinical), 60)
  expect_true(all(sim$clinical$asthma %in% c(0, 1)))
  expect_true(all(c("causal_cpgs", "causal_genes", "active_tfs", "true_effects") %in%
                     names(sim$ground_truth)))
  expect_true(length(sim$ground_truth$causal_cpgs) > 0)
})

test_that("simulate_trace_asthma_cohort respects a custom prevalence target approximately", {
  sim <- simulate_trace_asthma_cohort(n_subjects = 500, n_cpgs = 50, n_genes = 30,
                                       n_tfs = 4, prevalence = 0.2, seed = 7)
  expect_true(abs(mean(sim$clinical$asthma) - 0.2) < 0.08)
})

test_that("eQTM discovery on simulated data ranks causal pairs ahead of null pairs", {
  sim <- simulate_trace_asthma_cohort(n_subjects = 150, n_snps = 10, n_cpgs = 100,
                                       n_genes = 60, n_tfs = 4, seed = 3)
  eqtm_table <- discover_eqtm(sim$methylation, sim$expression, sim$eqtm_pairs,
                               sim$covariates, same_chr = sim$same_chr, fdr_threshold = 1.0)
  n_causal_pairs <- sum(sim$same_chr)
  causal_keys <- paste(sim$eqtm_pairs$cpg[sim$same_chr], sim$eqtm_pairs$gene[sim$same_chr])
  eqtm_table$is_causal <- paste(eqtm_table$cpg, eqtm_table$gene) %in% causal_keys
  med_causal <- stats::median(eqtm_table$p_value[eqtm_table$is_causal])
  med_null <- stats::median(eqtm_table$p_value[!eqtm_table$is_causal])
  expect_true(med_causal < med_null)
})
