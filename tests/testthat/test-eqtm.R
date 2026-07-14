data(example_cohort, package = "traceAsthma")

test_that("discover_eqtm classifies cis and trans pairs correctly", {
  cov_simple <- data.frame(age = example_cohort$covariates$age)
  rownames(cov_simple) <- rownames(example_cohort$covariates)

  eqtm <- discover_eqtm(
    example_cohort$methylation, example_cohort$expression,
    example_cohort$eqtm_pairs, cov_simple,
    same_chr = example_cohort$same_chr, fdr_threshold = 1.0
  )

  expect_true(all(c("cpg", "gene", "type", "alpha1", "p_value", "fdr") %in% colnames(eqtm)))
  expect_true(all(eqtm$type %in% c("cis", "trans")))
  # the only same-chromosome, <1Mb pair in the example data is cg_IL13-IL13
  expect_equal(eqtm$type[eqtm$cpg == "cg_IL13" & eqtm$gene == "IL13"], "cis")
})

test_that("discover_eqtm returns zero rows when FDR threshold is stringent", {
  cov_simple <- data.frame(age = example_cohort$covariates$age)
  rownames(cov_simple) <- rownames(example_cohort$covariates)
  eqtm <- discover_eqtm(
    example_cohort$methylation, example_cohort$expression,
    example_cohort$eqtm_pairs, cov_simple,
    same_chr = example_cohort$same_chr, fdr_threshold = 1e-10
  )
  expect_equal(nrow(eqtm), 0)
})

test_that("qc_rnaseq filters low-expression genes", {
  keep <- qc_rnaseq(example_cohort$expression, min_count = 1000, min_samples = 1)
  expect_length(keep, 0)
  keep_all <- qc_rnaseq(example_cohort$expression, min_count = 1, min_samples = 1)
  expect_equal(sort(keep_all), sort(rownames(example_cohort$expression)))
})

test_that("build_regulons filters by confidence class", {
  regs <- build_regulons(example_cohort$regulon_db, confidence_classes = c("A"))
  expect_true(all(regs$confidence == "A"))
  expect_true(nrow(regs) < nrow(example_cohort$regulon_db))
})

test_that("test_regulon_enrichment returns valid hypergeometric p-values", {
  enr <- test_regulon_enrichment(
    example_cohort$regulon_db,
    target_genes = c("IL4", "IL5", "IL13"),
    background_genes = rownames(example_cohort$expression),
    fdr_threshold = 1.0
  )
  expect_true(all(enr$p_value >= 0 & enr$p_value <= 1))
  expect_true("GATA3" %in% enr$tf)
})
