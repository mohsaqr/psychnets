# Nodewise lasso correctness is certified by the GLM stationarity (KKT) residual
# (glm_lasso_kkt), so these tests need no reference solver.

# Two latent factors drive two pairs of binary indicators -> the Ising network
# should connect within-pair nodes and (largely) separate across pairs.
.gen_ising <- function(seed, n = 500) {
  set.seed(seed)
  f1 <- stats::rnorm(n); f2 <- stats::rnorm(n)
  lin <- cbind(f1, f1, f2, f2) + matrix(stats::rnorm(n * 4, sd = 0.5), n, 4)
  b <- (lin > 0) * 1L
  colnames(b) <- paste0("V", 1:4)
  b
}

test_that("ising_fit nodewise regressions reach the penalized-likelihood optimum", {
  b <- .gen_ising(1)
  fit <- ising_fit(b, gamma = 0.25)
  expect_s3_class(fit, "psychnet")
  expect_false(fit$directed)
  expect_lt(fit$kkt, 1e-6)
})

test_that("ising_fit recovers the within-pair structure", {
  b <- .gen_ising(2)
  fit <- ising_fit(b, gamma = 0.25)
  g <- fit$graph
  # the two within-pair edges should be the strongest present edges
  expect_gt(abs(g["V1", "V2"]), 0)
  expect_gt(abs(g["V3", "V4"]), 0)
  expect_gt(abs(g["V1", "V2"]) + abs(g["V3", "V4"]),
            abs(g["V1", "V3"]) + abs(g["V2", "V4"]))
})

test_that("ising_fit is symmetric and respects the AND/OR rule", {
  b <- .gen_ising(3)
  a <- ising_fit(b, rule = "AND")
  o <- ising_fit(b, rule = "OR")
  expect_equal(a$graph, t(a$graph))
  expect_equal(o$graph, t(o$graph))
  # OR keeps at least as many edges as AND
  expect_gte(o$n_edges, a$n_edges)
  expect_length(a$thresholds, 4)
})

test_that("glm_lasso_kkt flags a non-stationary coefficient vector", {
  set.seed(4)
  n <- 300; p <- 5
  X <- matrix(stats::rnorm(n * p), n, p)
  X <- sweep(X, 2, colMeans(X)); X <- sweep(X, 2, sqrt(colMeans(X^2)), "/")
  y <- (X[, 1] + stats::rnorm(n) > 0) * 1L
  fit <- .glm_lasso_fit(X, y, "binomial", lambda = 0.05)
  bad <- fit$beta; bad[2] <- bad[2] + 0.5
  expect_gt(glm_lasso_kkt(X, y, fit$b0, bad, 0.05, "binomial"),
            glm_lasso_kkt(X, y, fit$b0, fit$beta, 0.05, "binomial"))
})

test_that("a slow fit with a finite optimum is not reported as separated", {
  # Nearly collinear predictors; every row has both outcomes, so the
  # unpenalized optimum is finite (glm() converges near +-230) but the IRLS
  # continuation approaches it slowly. 0.8.2 stopped it early as separated.
  g <- expand.grid(z = 0:1, d = c(-1, 1), a = c(-1, 1))
  X <- cbind(g$a, (g$a + 0.02 * g$d) / sqrt(1 + 0.02^2))
  w <- ifelse(g$z == (g$d + 1) / 2, 1, 0.01)
  fit <- .glm_lasso_fit(X, g$z, "binomial", lambda = 0, weights = w)
  expect_false(fit$separation)
  expect_true(fit$converged)
  ref <- suppressWarnings(stats::glm(g$z ~ X, family = stats::quasibinomial,
                                     weights = w, control = stats::glm.control(epsilon = 1e-14, maxit = 1000)))
  expect_equal(unname(fit$beta), unname(stats::coef(ref)[-1L]), tolerance = 1e-4)
})

test_that("a separated fit is still reported as separated", {
  X <- cbind(a = c(-1, 1, -1, 1, -1, 1), b = c(-1, -1, 1, 1, 0, 0))
  fit <- .glm_lasso_fit(X, c(0, 1, 0, 1, 0, 1), "binomial", lambda = 0)
  expect_true(fit$separation)
  expect_false(fit$converged)
})
