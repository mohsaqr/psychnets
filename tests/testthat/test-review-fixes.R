# Regression tests for the adversarial-review fixes (2026-06-22).

mk <- function(seed, n = 200, p = 5) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(0.4^abs(outer(1:p, 1:p, "-")))
  colnames(X) <- paste0("V", seq_len(p))
  X
}

test_that(".new_psychnet stores a zero-diagonal graph", {
  # a correlation-like matrix with unit diagonal must not leak the diagonal
  g <- matrix(0.5, 4, 4); diag(g) <- 1
  obj <- psychnet:::.new_psychnet(g, paste0("V", 1:4), "test", FALSE, 10)
  expect_true(all(diag(obj$graph) == 0))
  expect_equal(centrality(obj)$strength, rep(1.5, 4))   # 3 * 0.5, no diagonal
})

test_that(".new_psychnet rejects a label/dimension mismatch", {
  expect_error(psychnet:::.new_psychnet(matrix(0, 3, 3), c("a", "b"), "t", FALSE, 5),
               "labels length")
})

test_that("glasso_kkt is infinite for an indefinite precision matrix", {
  theta <- matrix(c(1, 2, 2, 1), 2, 2)                  # eigenvalues 3, -1
  expect_true(is.infinite(glasso_kkt(theta, diag(2), 0.1)))
  # a genuine optimum still certifies near zero
  S <- 0.5^abs(outer(1:5, 1:5, "-"))
  fit <- ebic_glasso(cor_matrix = S, n = 250)
  expect_lt(fit$kkt, 1e-7)
})

test_that("glm_lasso_kkt penalizes a non-stationary intercept", {
  set.seed(1)
  # zero predictors isolate the intercept: the slope gradient is exactly 0, so
  # the certificate reflects only mean(y - mu), the intercept stationarity.
  X <- matrix(0, 100, 2)
  y <- as.numeric(stats::rnorm(100) > 0)
  good <- glm_lasso_kkt(X, y, stats::qlogis(mean(y)), c(0, 0), 0, "binomial")
  bad  <- glm_lasso_kkt(X, y, stats::qlogis(mean(y)) + 1, c(0, 0), 0, "binomial")
  expect_lt(good, 1e-8)
  expect_gt(bad, 0.1)
})

test_that("an isolated binary node gets a logit-scale threshold", {
  set.seed(3)
  # V1 independent of the rest; V2..V5 a chain
  base <- mk(3, n = 600, p = 4)
  b <- cbind(V0 = as.numeric(stats::rnorm(600) > 0.4), (base > 0) * 1)
  colnames(b) <- paste0("V", 0:4)
  fit <- ising_fit(b)
  # threshold for the (near) isolated node should track qlogis(prevalence),
  # i.e. be negative for a low-prevalence node, never a probability in (0,1)
  expect_true(fit$thresholds[["V0"]] < 0)
})

test_that("centrality_stability CS is invariant to drop_prop ordering", {
  set.seed(7)
  x <- mk(8, n = 250, p = 5)
  a <- centrality_stability(x, drop_prop = c(0.1, 0.3, 0.5), iter = 25)
  b <- centrality_stability(x, drop_prop = c(0.5, 0.1, 0.3), iter = 25)
  expect_equal(sort(a$cs), sort(b$cs))
})

test_that("bootstrap keeps all directed edges for relimp", {
  set.seed(9)
  x <- mk(9, n = 150, p = 4)
  bs <- bootstrap_network(x, method = "relimp", n_boot = 20, cores = 1)
  expect_equal(nrow(bs$edges), 4 * 3)                   # p*(p-1) directed edges
})

test_that("paired nct requires equal group sizes", {
  expect_error(nct(mk(1, 10, 3), mk(2, 12, 3), iter = 5, paired = TRUE),
               "equal size")
})

test_that("correlation significance does not fabricate edges when df <= 0", {
  # n - 2 - (p-2) = n - p <= 0 for the partial test: nothing can be significant
  set.seed(2)
  x <- mk(2, n = 6, p = 6)
  fit <- pcor_network(x, alpha = 0.05)
  expect_equal(fit$n_edges, 0L)
})

test_that("centrality rejects a non-square bare matrix", {
  expect_error(centrality(matrix(1, 2, 3)), "square")
})

test_that("relimp rejects an indefinite cor_matrix", {
  S <- matrix(c(1, .9, .9, .9, 1, -.9, .9, -.9, 1), 3, 3)
  expect_error(relimp_network(cor_matrix = S), "positive-semidefinite")
})
