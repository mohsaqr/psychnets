# Regression tests for the pre-release algorithm audit.

test_that("cor_auto refuses pairs with no estimable joint information", {
  x <- cbind(a = c(1, 2, NA, NA), b = c(NA, NA, 1, 2),
             c = c(1, 2, 1, 2))
  expect_error(cor_auto(x), "fewer than two joint observations")
})

test_that("pairwise marginal inference uses edge-specific sample sizes", {
  set.seed(1)
  x <- cbind(a = rnorm(100), b = rnorm(100), c = rnorm(100))
  x[1:80, 2] <- NA
  fit <- cor_network(x, alpha = 0.05)
  r <- fit$cor_matrix[1, 2]
  expected <- 2 * pt(-abs(r * sqrt(18 / (1 - r^2))), df = 18)
  expect_equal(fit$p_values[1, 2], expected)
  expect_error(pcor_network(x, alpha = 0.05), "pairwise-complete")
})

test_that("node betweenness uses undirected-pair normalization", {
  A <- matrix(0, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
  A[1, 2] <- A[2, 1] <- A[2, 3] <- A[3, 2] <- 1
  expect_equal(net_centralities(A, measures = "betweenness")$betweenness,
               c(0, 1, 0))
})

test_that("bootstrap difference p-values treat exact ties as null evidence", {
  z <- matrix(0, 20, 2)
  expect_equal(psychnets:::.psn_diff_pmat(z, c("a", "b"))[1, 2], 1)
  boot <- structure(list(
    edge_boot = z, centrality_boot = list(), edge_labels = c("a", "b"),
    node_labels = character(), measures = character(), ci = 0.95,
    edges = data.frame(observed = c(0, 0)), centrality = data.frame()
  ), class = "psychnet_bootstrap")
  expect_equal(difference_test(boot)$p_value, 1)
})

test_that("rank dichotomization never splits ties by row position", {
  expect_error(dichotomize(matrix(1, 10, 2), "rank"), "constant")
  x <- matrix(c(rep(1, 3), rep(2, 4), rep(3, 3)), ncol = 1)
  b <- dichotomize(x, "rank")
  expect_equal(length(unique(b[x == 2])), 1L)
})

test_that("average aggregation preserves a directed estimator", {
  set.seed(2)
  x <- matrix(rnorm(300 * 4), 300, 4)
  colnames(x) <- paste0("V", 1:4)
  a <- net_aggregate(x, c(1, 1, 2, 2), method = "average",
                     estimator = "relimp")
  expect_true(a$directed)
})

test_that("unavailable and post-threshold certificates are not successes", {
  S <- 0.5^abs(outer(1:5, 1:5, "-"))
  thresholded <- ebic_glasso(cor_matrix = S, n = 300, threshold = 1)
  expect_true(is.na(certificate(thresholded)$certified))
  expect_true(is.finite(thresholded$fit_kkt))
  conditioned <- condition(mgm_fit(
    data.frame(x = rnorm(120), y = rnorm(120), z = rnorm(120), m = rep(0:1, 60)),
    types = c("g", "g", "g", "c"), moderators = 4,
    threshold = "none"), 0)
  expect_true(is.na(certificate(conditioned)$certified))
})

test_that("MGM exposes levels and factor data remain predictable", {
  skip_if_not_installed("glmnet")
  set.seed(3)
  d <- data.frame(g = rnorm(160), b = factor(sample(c("no", "yes"), 160, TRUE)),
                  z = rnorm(160))
  fit <- mgm_fit(d, types = c("g", "c", "g"), native = FALSE,
                 threshold = "none", nlambda = 12, lambda_min_ratio = 0.2)
  expect_equal(unname(fit$levels), c(1L, 2L, 1L))
  expect_equal(fit$nlambda, 12)
  expect_equal(fit$lambda_min_ratio, 0.2)
  expect_s3_class(net_predict(fit, d), "data.frame")
})

test_that("moderated MGM preserves all tuning choices", {
  set.seed(4)
  d <- data.frame(x = rnorm(160), y = rnorm(160), z = rnorm(160),
                  m = rep(0:1, 80))
  fit <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = 4,
                 threshold = "HW", nlambda = 13, lambda_min_ratio = 0.3)
  expect_identical(fit$params$threshold, "HW")
  expect_equal(fit$params$nlambda, 13)
  expect_equal(fit$params$lambda_min_ratio, 0.3)
})

test_that("the unified MGM front door forwards its threshold rule", {
  set.seed(41)
  d <- data.frame(x = rnorm(120), y = rnorm(120), z = rnorm(120),
                  b = rbinom(120, 1, 0.5))
  via_front <- mgm_fit(d, threshold = "none")
  via_router <- psychnet(d, method = "mgm", threshold = "none")
  expect_identical(via_router$threshold, "none")
  expect_equal(via_router$weights, via_front$weights)
})

test_that("resampling aligns rare nodes and translates engine", {
  skip_slow()
  set.seed(5)
  b <- cbind(x = rbinom(50, 1, 0.5), rare = c(rep(0, 49), 1),
             z = rbinom(50, 1, 0.5))
  bs <- net_boot(b, method = "ising", n_boot = 20, cores = 1,
                 engine = "base")
  expect_equal(dim(bs$edge_boot), c(20L, 3L))
  expect_equal(bs$n_success, 20L)
  expect_false(anyNA(bs$edge_boot))
  st <- net_stability(b, method = "ising", iter = 5, drop_prop = 0.3)
  expect_equal(st$n_success, 5L)
})

test_that("a dropped bootstrap node does not erase other predictabilities", {
  set.seed(51)
  x <- cbind(a = rnorm(60), rare = 0, b = rnorm(60))
  fit <- psychnets:::.psn_refit_network(x, "glasso", colnames(x))
  pred <- psychnets:::.psn_predictability(fit)
  expect_true(is.na(pred["rare"]))
  expect_true(all(is.finite(pred[c("a", "b")])))
})

test_that("estimator_args carries estimator thresholds through resampling", {
  set.seed(6)
  x <- matrix(rnorm(120 * 4), 120, 4)
  bs <- net_boot(x, n_boot = 2, cores = 1,
                 estimator_args = list(threshold = 1))
  expect_equal(sum(abs(bs$observed$weights)), 0)
  expect_error(net_split_reliability(x[1:3, ], split = 0.9, iter = 2),
               "at least two")
})

test_that("penalized-diagonal GGM predictability uses fitted variances", {
  S <- 0.5^abs(outer(1:5, 1:5, "-"))
  fit <- ebic_glasso(cor_matrix = S, n = 300, penalize_diagonal = TRUE)
  expected <- 1 - 1 / (diag(fit$precision) * diag(solve(fit$precision)))
  expect_equal(net_predict(fit)$predictability, unname(expected),
               tolerance = 1e-10)
})

test_that("undirected graph diagnostics reject malformed inputs", {
  A <- matrix(c(0, 1, 0, 0), 2, 2)
  expect_error(net_smallworld(A), "symmetric")
  expect_error(net_clustering(A), "symmetric")
  expect_error(net_bridge(A, c(1, 2)), "symmetric")
  expect_error(net_smallworld(diag(3), n_rand = 0), "positive whole")
  expect_error(psychnets:::.new_psychnet(diag(2), c("x", "x"), "test",
                                          FALSE, 10), "unique")
})
