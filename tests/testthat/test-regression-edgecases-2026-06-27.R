# Regression tests for the 2026-06-27 exhaustive edge-case audit. Each test
# pins a previously silent-wrong / crash / cryptic-error path to an explicit,
# clear rejection (or, for the integer-coercion fixes, a clean print).

mat_g <- function(seed = 1, n = 200, p = 4) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))
  X
}
mat_b <- function(seed = 1, n = 200, p = 4) {
  (mat_g(seed, n, p) > 0) * 1L
}

# --- HIGH: mgm type/coding validation (was silent garbage) -------------------

test_that("mgm_fit rejects a continuous column declared binary via types", {
  d <- as.data.frame(mat_g(1))                       # all gaussian
  expect_error(mgm_fit(d, types = c("c", "c", "c", "c")),
               "not coded 0/1")
})

test_that("mgm_fit rejects a factor/character column instead of dropping it", {
  d <- as.data.frame(mat_g(2))
  d$fc <- factor(sample(letters[1:3], nrow(d), replace = TRUE))
  expect_error(mgm_fit(d), "non-numeric")
})

test_that("mgm_fit gives an accurate level-count message for a {1,2} column", {
  d <- as.data.frame(mat_g(3, p = 2))
  d$x12 <- sample(1:2, nrow(d), replace = TRUE)
  expect_error(mgm_fit(d), "2 levels not in \\{0, 1\\}")
})

# --- MEDIUM: cor_matrix normalization + PSD ----------------------------------

test_that("a covariance cor_matrix is normalized to unit scale (not raw weights)", {
  S   <- 0.4^abs(outer(1:5, 1:5, "-"))
  Cov <- diag(2:6) %*% S %*% diag(2:6)               # covariance, entries >> 1
  w <- cor_network(cor_matrix = Cov)$weights
  expect_lte(max(abs(w)), 1 + 1e-8)                  # was ~160 before the fix
  # and a proper correlation matrix passes through unchanged
  expect_equal(unname(psychnets:::.check_cor_matrix(S)), unname(S))
})

test_that("an indefinite cor_matrix is rejected, not silently NaN", {
  Ind <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  expect_error(pcor_network(cor_matrix = Ind), "positive semi-definite")
})

test_that("the EBIC GGMs validate n as a single positive scalar", {
  S <- 0.4^abs(outer(1:5, 1:5, "-"))
  expect_error(ebic_glasso(cor_matrix = S, n = 0),          "n > 0")
  expect_error(ebic_glasso(cor_matrix = S, n = c(100, 200)), "length\\(n\\) == 1")
  expect_error(huge_network(cor_matrix = S, n = -5),         "n > 0")
})

# --- MEDIUM: exported certificate validation ---------------------------------

test_that("glm_lasso_kkt validates dimensions and family", {
  X <- scale(matrix(stats::rnorm(200 * 3), 200, 3))
  y <- as.numeric(X %*% c(1, 0, -1) + stats::rnorm(200))
  b <- stats::lm.fit(cbind(1, X), y)$coefficients
  expect_error(glm_lasso_kkt(X, y[1:150], b[1], b[-1], 0), "length\\(y\\) == nrow\\(X\\)")
  expect_error(glm_lasso_kkt(X, y, b[1], b[-1], 0, family = "poisson"), "should be one of")
})

test_that("glasso_kkt / ggm_support_kkt reject mismatched or non-finite inputs", {
  theta <- diag(4); S <- diag(4)
  expect_error(glasso_kkt(theta, diag(3), 0.1), "dim\\(theta\\) == dim")
  expect_error(ggm_support_kkt(theta, S, matrix(TRUE, 3, 3)), "dim\\(support\\)")
  expect_error(ggm_support_kkt(theta, S, matrix(1, 4, 4)),    "is.logical")
})

test_that("structural certificates reject a wrong-type network", {
  cn <- cor_network(mat_g(4))
  expect_error(lmg_certificate(cn),  "relimp network")
  expect_error(tmfg_certificate(cn), "TMFG network")
})

# --- MEDIUM/LOW: nodewise arg validation -------------------------------------

test_that("ising_fit / mgm_fit reject wrong-length labels with a clear message", {
  b <- mat_b(5)
  expect_error(ising_fit(b, labels = c("a", "b")),    "length\\(labels\\) == p")
  expect_error(mgm_fit(mat_g(5), labels = c("a", "b")), "length\\(labels\\) == p")
})

test_that("ising_sampler validates alpha and ising_fit validates the path grid", {
  b <- mat_b(6)
  expect_error(ising_sampler(b, alpha = 2),  "alpha")
  expect_error(ising_fit(b, nlambda = 0),    "nlambda")
  expect_error(ising_fit(b, lambda_min_ratio = 1), "lambda_min_ratio")
})

# --- MEDIUM: integer coercion (was a print crash) ----------------------------

test_that("a fractional n_boot / iter does not corrupt the printed object", {
  skip_slow()
  X <- mat_g(7)
  expect_no_error(utils::capture.output(print(net_boot(X, n_boot = 2.7, cores = 1))))
  expect_no_error(utils::capture.output(print(net_stability(X, drop_prop = 0.5, iter = 2.7))))
})

# --- LOW: centrality input guard ---------------------------------------------

test_that("centrality rejects a non-matrix, non-psychnet argument clearly", {
  expect_error(net_centralities(list(a = 1)), "psychnet object or a square")
  expect_error(net_centralities(NULL),        "psychnet object or a square")
})
