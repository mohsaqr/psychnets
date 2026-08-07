# Golden regression anchors for the estimators downstream consumers pin against.
#
# Nestimate delegates its cross-sectional estimators to psychnets and holds its
# own frozen baselines at 1e-12, and it has a CRAN reverse dependency, so a
# weight that moves here surfaces as a revdep failure there with no code change
# on either side. These tests exist so that movement is caught in psychnets
# first. Regenerate with tests/testthat/fixtures/make-golden.R only when the
# change of numbers is intentional, and list it in NEWS.md.

# The data generators come from helper-golden.R, which testthat loads
# automatically; make-golden.R is NOT sourced here (it would rewrite the very
# fixtures these tests exist to guard).
golden <- readRDS(test_path("fixtures", "golden-estimators.rds"))

test_that("ebic_glasso default output is unchanged", {
  fit <- ebic_glasso(cor_matrix = golden_cor(), n = 300)
  expect_equal(fit$weights, golden$glasso_base, tolerance = 1e-12)
  expect_equal(fit$lambda, golden$glasso_lambda, tolerance = 1e-12)
})

test_that("mgm_fit base output is unchanged", {
  fit <- mgm_fit(golden_data_mixed(), native = TRUE)
  expect_equal(fit$weights, golden$mgm_base, tolerance = 1e-12)
})

test_that("ising_fit base output is unchanged", {
  fit <- ising_fit(golden_data_binary(), native = TRUE)
  expect_equal(fit$weights, golden$ising_base, tolerance = 1e-12)
})

# The glmnet paths are what Nestimate uses for bit-exact mgm::mgm / IsingFit
# parity, so they are the ones a k-class generalisation of the nodewise
# aggregation could silently regress.
test_that("mgm_fit glmnet output is unchanged", {
  skip_if_not_installed("glmnet")
  fit <- mgm_fit(golden_data_mixed(), native = FALSE)
  expect_equal(fit$weights, golden$mgm_glmnet, tolerance = 1e-12)
})

test_that("ising_fit glmnet output is unchanged", {
  skip_if_not_installed("glmnet")
  fit <- ising_fit(golden_data_binary(), native = FALSE)
  expect_equal(fit$weights, golden$ising_glmnet, tolerance = 1e-12)
})
