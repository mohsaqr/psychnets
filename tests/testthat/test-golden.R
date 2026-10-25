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

# The platform-determinism contract behind the fixture above: inactive lasso
# coordinates are EXACT zeros. Coordinate descent can leave a BLAS-rounding
# crumb (~1e-17) on a coordinate whose gradient sits exactly at the penalty
# boundary -- this data produces one on the g3 node -- and whether the crumb is
# 0.0 or 4e-17 differs between BLAS implementations. If it survives into the
# returned support, the AND rule fabricates a platform-dependent edge.
test_that("base-kernel nodewise support contains no sub-tolerance crumbs", {
  fit <- mgm_fit(golden_data_mixed(), native = TRUE)
  b <- fit$nodewise$beta_std
  expect_true(all(b == 0 | abs(b) > 1e-10))
  fit2 <- ising_fit(golden_data_binary(), native = TRUE)
  b2 <- fit2$nodewise$beta_std
  expect_true(all(b2 == 0 | abs(b2) > 1e-10))
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

# panel_gvar(): a separate fixture so adding the estimator did not touch the
# frozen estimator baselines above. 1e-10 is the policy's own threshold; the ML
# fit is polished to a score of ~1e-13, so the estimate is reproducible far
# below it. Data from helper-panelgvar.R; regenerate with make-panelgvar.R.
golden_pgv <- readRDS(test_path("fixtures", "golden-panelgvar.rds"))

test_that("panel_gvar output is unchanged", {
  cases <- list(default = list(n = 400, p = 3, waves = 4, seed = 1,
                               alpha = 0.01, adjust = "none"),
                large = list(n = 600, p = 5, waves = 5, seed = 8,
                             alpha = 0.01, adjust = "none"),
                bh = list(n = 150, p = 4, waves = 4, seed = 10,
                          alpha = 0.05, adjust = "BH"),
                stepup_small = list(n = 45, p = 4, waves = 4, seed = 21,
                                    alpha = 0.01, adjust = "none"),
                stepup_waves = list(n = 300, p = 4, waves = 5, seed = 12,
                                    alpha = 0.01, adjust = "none"))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    w <- panel_sim_wide(n = cs$n, p = cs$p, waves = cs$waves, seed = cs$seed)
    # The n = 45 design is a known boundary / partly unidentified solution;
    # its boundary and identification warnings are not what this test checks.
    quiet <- function(expr) withCallingHandlers(expr,
      psychnet_boundary = function(w) invokeRestart("muffleWarning"),
      psychnet_singular_information = function(w) invokeRestart("muffleWarning"),
      psychnet_approximate_se = function(w) invokeRestart("muffleWarning"))
    sat <- quiet(panel_gvar(w$data, vars = w$design, select = "none"))
    pru <- quiet(panel_gvar(w$data, vars = w$design, select = "prune",
                            alpha = cs$alpha, adjust = cs$adjust,
                            stepup = FALSE))
    up <- quiet(panel_gvar(w$data, vars = w$design, select = "prune",
                           alpha = cs$alpha, adjust = cs$adjust))
    expect_equal(as.data.frame(sat, what = "parameters")$estimate,
                 golden_pgv[[nm]]$saturated, tolerance = 1e-10, label = nm)
    expect_equal(as.data.frame(pru, what = "parameters")$estimate,
                 golden_pgv[[nm]]$pruned, tolerance = 1e-10, label = nm)
    expect_equal(as.data.frame(up, what = "parameters")$estimate,
                 golden_pgv[[nm]]$stepup, tolerance = 1e-10, label = nm)
  }
})

# select = "ebic": the certified lasso optimum (KKT ~1e-10) on three designs.
# 1e-8 rather than 1e-10: a lasso estimate's precision is its KKT residual
# divided by the local curvature, ~1e-9 here.
test_that("panel_gvar EBIC output is unchanged", {
  skip_on_cran()
  designs <- list(default = list(400, 3, 4, 1), bh = list(150, 4, 4, 10),
                  stepup_waves = list(300, 4, 5, 12))
  for (nm in names(designs)) {
    d <- designs[[nm]]
    w <- panel_sim_wide(n = d[[1]], p = d[[2]], waves = d[[3]], seed = d[[4]])
    fit <- panel_gvar(w$data, vars = w$design, select = "ebic")
    expect_equal(as.data.frame(fit, what = "parameters")$estimate,
                 golden_pgv$ebic[[nm]]$estimate, tolerance = 1e-8, label = nm)
    expect_equal(attr(fit, "settings")$lambda, golden_pgv$ebic[[nm]]$lambda,
                 tolerance = 1e-10, label = nm)
  }
})
