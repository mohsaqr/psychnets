# ggmModSelect refits the UNREGULARIZED MLE on a selected graph; correctness is
# the constrained-MLE stationarity residual (ggm_support_kkt).

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

test_that("ggm_modselect certifies the constrained-MLE optimum", {
  for (S in list(ar1(6, 0.5), ar1(7, 0.6))) {
    fit <- ggm_modselect(cor_matrix = S, n = 300)
    expect_s3_class(fit, "psychnet")
    expect_lt(fit$kkt, 1e-8)
    expect_equal(ggm_support_kkt(fit$precision, S, fit$support), fit$kkt)
  }
})

test_that("the selected support is a symmetric edge set", {
  S <- ar1(6, 0.5)
  fit <- ggm_modselect(cor_matrix = S, n = 300)
  expect_true(is.logical(fit$support))
  expect_equal(fit$support, t(fit$support))
  expect_equal(fit$graph, t(fit$graph))
})

test_that("the full graph recovers the exact Gaussian MLE (W = S)", {
  S <- ar1(5, 0.5)
  full <- matrix(TRUE, 5, 5); diag(full) <- FALSE
  theta <- psychnet:::.ggm_fit_support(S, full)
  expect_lt(max(abs(solve(theta) - S)), 1e-8)        # saturated model = S^{-1}
  expect_lt(max(abs(theta - solve(S))), 1e-8)
})

test_that("ggm_modselect dispatches via psychnet", {
  S <- ar1(5, 0.5)
  X <- matrix(stats::rnorm(300 * 5), 300, 5) %*% chol(S)
  expect_equal(psychnet(X, "ggmModSelect")$method, "ggm")
  expect_equal(psychnet(X, "modselect")$method, "ggm")
})
