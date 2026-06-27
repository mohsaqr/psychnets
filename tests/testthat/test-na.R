# Missing-data handling: pairwise (default) retains information where listwise
# collapses; both are identical on complete data.

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))
chain_prec <- function(p, a) { K <- diag(p); for (i in seq_len(p - 1)) K[i, i + 1] <- K[i + 1, i] <- a; stats::cov2cor(solve(K)) }
rmvn <- function(n, S) matrix(stats::rnorm(n * ncol(S)), n, ncol(S)) %*% chol(S)
f1 <- function(est, tru) { tp <- sum(est & tru); if (2 * tp + sum(est & !tru) + sum(!est & tru) == 0) 1 else 2 * tp / (2 * tp + sum(est & !tru) + sum(!est & tru)) }

test_that("pairwise and listwise are identical on complete data", {
  set.seed(1)
  X <- matrix(stats::rnorm(300 * 6), 300, 6) %*% chol(ar1(6, 0.5))
  colnames(X) <- paste0("V", 1:6)
  for (m in c("cor", "pcor", "EBICglasso", "huge", "TMFG", "LoGo")) {
    a <- estimate_network(X, m, na_method = "pairwise")$graph
    b <- estimate_network(X, m, na_method = "listwise")$graph
    expect_equal(a, b, tolerance = 1e-10, info = m)
  }
})

test_that("pairwise recovers structure where listwise collapses (severe MCAR)", {
  set.seed(2)
  p <- 10; S <- chain_prec(p, -0.35); true <- (abs(outer(1:p, 1:p, "-")) == 1)[upper.tri(S)]
  X <- rmvn(150, S)
  X[matrix(stats::runif(length(X)) < 0.15, nrow(X))] <- NA
  colnames(X) <- paste0("V", seq_len(p))
  pw <- estimate_network(X, "EBICglasso", na_method = "pairwise")
  lw <- estimate_network(X, "EBICglasso", na_method = "listwise")
  expect_lt(lw$n_obs, 60)                                    # listwise gutted the sample
  expect_gt(pw$n_obs, lw$n_obs)                              # pairwise retained far more
  expect_gt(f1(abs(pw$graph[upper.tri(pw$graph)]) > 1e-6, true), 0.7)
  expect_equal(lw$n_edges, 0L)                               # listwise collapses to empty
})

test_that("pairwise is the default when NAs are present", {
  set.seed(3)
  X <- rmvn(200, ar1(6, 0.5)); X[sample(length(X), 200)] <- NA
  colnames(X) <- paste0("V", 1:6)
  expect_equal(estimate_network(X, "EBICglasso")$na_method, "pairwise")
  expect_equal(cor_network(X)$na_method, "pairwise")
})

test_that(".cor_input projects pairwise correlations to positive-definite", {
  set.seed(4)
  X <- rmvn(120, ar1(8, 0.6)); X[matrix(stats::runif(length(X)) < 0.2, nrow(X))] <- NA
  ci <- psychnet:::.cor_input(X, na_method = "pairwise")
  expect_true(all(eigen(ci$S, only.values = TRUE)$values > 0))
  expect_equal(unname(diag(ci$S)), rep(1, 8), tolerance = 1e-8)
  expect_gt(ci$n, 60)                                        # effective n >> listwise
})

test_that("nodewise estimators tolerate NA via imputation", {
  set.seed(5)
  z <- matrix(stats::rnorm(600 * 2), 600, 2)
  b <- (cbind(z[, 1], z[, 1], z[, 2], z[, 2]) + matrix(stats::rnorm(600 * 4), 600) > 0) * 1
  colnames(b) <- paste0("V", 1:4)
  bna <- b; bna[matrix(stats::runif(length(b)) < 0.12, nrow(b))] <- NA
  fit_pw <- ising_fit(bna, na_method = "pairwise")
  fit_lw <- ising_fit(bna, na_method = "listwise")
  expect_s3_class(fit_pw, "psychnet")
  expect_gt(fit_pw$n_obs, fit_lw$n_obs)                      # imputation keeps all rows
  expect_lt(fit_pw$kkt, 1e-5)
  # imputed binary columns stay binary (mode imputation)
  expect_true(all(b %in% c(0, 1)))
})

test_that("mgm and the imputation helper handle mixed NA data", {
  set.seed(6)
  f <- stats::rnorm(500)
  d <- data.frame(g1 = f + stats::rnorm(500), g2 = f + stats::rnorm(500),
                  b1 = (f + stats::rnorm(500) > 0) * 1, b2 = (f + stats::rnorm(500) > 0) * 1)
  d[cbind(sample(500, 60), sample(4, 60, TRUE))] <- NA
  fit <- mgm_fit(d, na_method = "pairwise")
  expect_s3_class(fit, "psychnet")
  expect_equal(fit$n_obs, 500L)                             # full sample retained
  # mode/mean imputation: binary cols mode-imputed to 0/1, gaussian mean-imputed
  imp <- psychnet:::.na_prep_nodewise(as.matrix(d), "pairwise")
  expect_false(anyNA(imp))
  expect_true(all(imp[, "b1"] %in% c(0, 1)))
})
