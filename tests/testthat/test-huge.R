# Self-certified against the convex glasso objective on the transformed
# correlation; no reference solver needed.

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

test_that("huge_network certifies the transformed-correlation optimum", {
  set.seed(1)
  X <- matrix(stats::rnorm(400 * 6), 400, 6) %*% chol(ar1(6, 0.5))
  Xn <- exp(X)                                   # break multivariate normality
  colnames(Xn) <- paste0("V", 1:6)
  for (npn in c("shrinkage", "truncation", "skeptic")) {
    fit <- huge_network(Xn, npn = npn)
    expect_s3_class(fit, "psychnet")
    expect_lt(fit$kkt, 1e-6)
    expect_equal(glasso_kkt(fit$precision, fit$cor_matrix, fit$lambda), fit$kkt)
  }
})

test_that("huge_network graph is a symmetric partial-correlation matrix", {
  set.seed(2)
  X <- exp(matrix(stats::rnorm(300 * 5), 300, 5) %*% chol(ar1(5, 0.5)))
  fit <- huge_network(X)
  expect_equal(fit$graph, t(fit$graph))
  expect_true(all(diag(fit$graph) == 0))
  expect_true(all(abs(fit$graph) <= 1 + 1e-8))
})

test_that("huge_network dispatches via estimate_network", {
  set.seed(3)
  X <- exp(matrix(stats::rnorm(300 * 5), 300, 5))
  expect_equal(estimate_network(X, "huge")$method, "huge")
  expect_equal(estimate_network(X, "npn")$method, "huge")
})
