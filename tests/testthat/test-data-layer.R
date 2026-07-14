test_that("import_data / export_data round-trip CSV correctly", {
  m <- matrix(runif(12), nrow = 3, dimnames = list(c("a","b","c"), c("P1","P2","P3","P4")))
  tmp <- tempfile(fileext = ".csv")
  export_data(m, tmp)
  reloaded <- import_data(tmp)
  expect_equal(dim(reloaded), dim(m))
  expect_equal(rownames(reloaded), rownames(m))
  expect_equal(unname(reloaded), unname(m), tolerance = 1e-6)
})

test_that("import_data / export_data round-trip RDS correctly", {
  m <- matrix(1:6, nrow = 2, dimnames = list(c("a","b"), c("P1","P2","P3")))
  tmp <- tempfile(fileext = ".rds")
  export_data(m, tmp)
  reloaded <- import_data(tmp)
  expect_true(is.matrix(reloaded))
  expect_equal(dim(reloaded), dim(m))
})

test_that("import_data errors informatively on unrecognized extension", {
  expect_error(import_data("file.xyz"), "Could not infer format")
})

test_that("validate_matrix accepts well-formed methylation data", {
  m <- matrix(runif(20, 0, 1), nrow = 4, dimnames = list(paste0("cg", 1:4), paste0("P", 1:5)))
  expect_true(validate_matrix(m, type = "methylation"))
})

test_that("validate_matrix rejects out-of-range methylation values", {
  m <- matrix(c(0.5, 1.5, 0.3, 0.2), nrow = 2, dimnames = list(c("a","b"), c("P1","P2")))
  expect_error(validate_matrix(m, type = "methylation"), "outside \\[0,1\\]")
})

test_that("validate_matrix rejects negative expression values", {
  m <- matrix(c(5, -2, 3, 8), nrow = 2, dimnames = list(c("a","b"), c("P1","P2")))
  expect_error(validate_matrix(m, type = "expression"), "non-negative")
})

test_that("validate_matrix rejects matrices without dimnames", {
  m <- matrix(runif(6), nrow = 2)
  expect_error(validate_matrix(m, type = "expression"), "row names")
})

test_that("validate_matrix rejects too few samples", {
  m <- matrix(runif(2), nrow = 2, dimnames = list(c("a","b"), "P1"))
  expect_error(validate_matrix(m, type = "expression", min_samples = 3), "at least")
})

test_that("validate_matrix rejects excessive missingness", {
  m <- matrix(c(NA, NA, NA, 0.5, 0.5, 0.5), nrow = 2,
              dimnames = list(c("a","b"), c("P1","P2","P3")))
  expect_error(validate_matrix(m, type = "methylation", max_missing_frac = 0.2), "missing")
})

test_that("validate_clinical accepts a well-formed clinical table", {
  cl <- data.frame(patient_id = paste0("P", 1:5), asthma = c(0,1,0,1,1),
                    age = c(30,40,25,55,33), sex = c(0,1,0,1,0))
  expect_true(validate_clinical(cl))
})

test_that("validate_clinical rejects missing required column", {
  cl <- data.frame(patient_id = paste0("P", 1:3), age = c(30,40,25))
  expect_error(validate_clinical(cl), "missing required")
})

test_that("validate_clinical rejects invalid asthma coding", {
  cl <- data.frame(patient_id = paste0("P", 1:3), asthma = c(0,1,2))
  expect_error(validate_clinical(cl), "0/1")
})

test_that("validate_cohort_alignment passes for aligned inputs", {
  m <- matrix(1, 2, 3, dimnames = list(c("a","b"), c("P1","P2","P3")))
  e <- matrix(1, 2, 3, dimnames = list(c("g1","g2"), c("P1","P2","P3")))
  cl <- data.frame(patient_id = c("P1","P2","P3"), asthma = c(0,1,0))
  expect_true(validate_cohort_alignment(methylation = m, expression = e, clinical = cl))
})

test_that("validate_cohort_alignment detects reordered subjects", {
  m <- matrix(1, 2, 3, dimnames = list(c("a","b"), c("P1","P2","P3")))
  e <- matrix(1, 2, 3, dimnames = list(c("g1","g2"), c("P2","P1","P3")))
  expect_error(validate_cohort_alignment(methylation = m, expression = e), "different order")
})

test_that("validate_cohort_alignment detects missing subjects", {
  m <- matrix(1, 2, 3, dimnames = list(c("a","b"), c("P1","P2","P3")))
  e <- matrix(1, 2, 2, dimnames = list(c("g1","g2"), c("P1","P2")))
  expect_error(validate_cohort_alignment(methylation = m, expression = e), "missing from")
})

test_that("data_dictionary returns the expected structure and filters by object", {
  dd <- data_dictionary()
  expect_true(all(c("object","field","type","units","coding_scheme","required","description") %in% colnames(dd)))
  dd_clin <- data_dictionary("clinical")
  expect_true(all(dd_clin$object == "clinical"))
  expect_true(nrow(dd_clin) < nrow(dd))
})

test_that("create_data_templates writes readable, correctly-shaped files", {
  tmp <- tempfile()
  paths <- create_data_templates(tmp)
  expect_true(all(file.exists(paths)))
  cl <- read.csv(paths["clinical"])
  expect_true(all(c("patient_id","asthma","age","sex") %in% colnames(cl)))
})
