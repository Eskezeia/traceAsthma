test_that("z_score standardizes to mean 0, sd 1", {
  x <- c(1, 2, 3, 4, 5)
  z <- z_score(x)
  expect_equal(mean(z), 0, tolerance = 1e-10)
  expect_equal(sd(z), 1, tolerance = 1e-10)
})

test_that("z_score standardizes within strata independently", {
  x <- c(1, 2, 3, 10, 11, 12)
  strata <- c("A", "A", "A", "B", "B", "B")
  z <- z_score(x, strata = strata)
  expect_equal(mean(z[strata == "A"]), 0, tolerance = 1e-10)
  expect_equal(mean(z[strata == "B"]), 0, tolerance = 1e-10)
})

test_that("weighted_score matches hand-computed example", {
  fm <- matrix(c(0.72, 0.35, 0.61, 0.42), nrow = 2,
               dimnames = list(c("P001", "P002"), c("cg_IL13", "cg_GATA3")))
  w <- c(cg_IL13 = 0.55, cg_GATA3 = 0.45)
  raw <- weighted_score(fm, w, standardize = FALSE)
  expect_equal(unname(raw["P001"]), 0.55 * 0.72 + 0.45 * 0.61, tolerance = 1e-8)
  expect_equal(unname(raw["P002"]), 0.55 * 0.35 + 0.45 * 0.42, tolerance = 1e-8)
})

test_that("weighted_score warns and drops weights with no matching column", {
  fm <- matrix(1:4, nrow = 2, dimnames = list(c("P1", "P2"), c("a", "b")))
  w <- c(a = 1, c = 2)
  expect_warning(weighted_score(fm, w, standardize = FALSE), "not found")
})

test_that("prune_correlated_features drops perfectly correlated duplicates", {
  x <- matrix(c(1, 2, 3, 4, 5,
                2, 4, 6, 8, 10,   # perfectly correlated with column 1
                5, 1, 4, 2, 3),  # independent
              nrow = 5, dimnames = list(NULL, c("a", "b", "c")))
  assoc <- c(a = 1.0, b = 0.5, c = 3.0)
  retained <- prune_correlated_features(x, assoc, threshold = 0.9)
  expect_true("c" %in% retained)
  expect_true(xor("a" %in% retained, "b" %in% retained))  # only one of the pair kept
  expect_false(all(c("a", "b") %in% retained))
})

test_that("requirePkg errors informatively for a missing package", {
  expect_error(requirePkg("thisPackageDoesNotExist12345", "a test"), "not installed")
})
