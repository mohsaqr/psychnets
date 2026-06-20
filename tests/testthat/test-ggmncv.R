# GGMncv is solved by one-step LLA as a weighted glasso; correctness is the
# weighted KKT residual (glasso_kkt_weighted), no reference solver.

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

test_that("ggmncv_network certifies the weighted optimum for every penalty", {
  S <- ar1(7, 0.5)
  for (pen in c("atan", "scad", "mcp", "lasso")) {
    fit <- ggmncv_network(cor_matrix = S, n = 300, penalty = pen)
    expect_s3_class(fit, "psychnet")
    expect_lt(fit$kkt, 1e-7)
    expect_equal(glasso_kkt_weighted(fit$precision, S, fit$penalty_matrix),
                 fit$kkt)
  }
})

test_that("ggmncv_network graph is a symmetric partial-correlation matrix", {
  S <- ar1(6, 0.5)
  fit <- ggmncv_network(cor_matrix = S, n = 300, penalty = "scad")
  expect_equal(fit$graph, t(fit$graph))
  expect_true(all(diag(fit$graph) == 0))
  expect_true(all(abs(fit$graph) <= 1 + 1e-8))
})

test_that("non-convex penalties shrink strong edges less than the lasso", {
  S <- ar1(6, 0.6)
  lasso <- ggmncv_network(cor_matrix = S, n = 400, penalty = "lasso")
  scad  <- ggmncv_network(cor_matrix = S, n = 400, penalty = "scad")
  # the largest retained partial correlation is at least as large under SCAD
  expect_gte(max(abs(scad$graph)) + 1e-8, max(abs(lasso$graph)))
})

test_that("ggmncv_network dispatches via estimate_network", {
  S <- ar1(5, 0.5)
  X <- matrix(stats::rnorm(300 * 5), 300, 5) %*% chol(S)
  expect_equal(estimate_network(X, "GGMncv", penalty = "mcp")$method, "GGMncv")
  expect_equal(estimate_network(X, "ncv")$method, "GGMncv")
})
