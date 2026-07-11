# cor_auto: polychoric / polyserial automatic correlation (qgraph parity).

test_that("the bivariate-normal CDF satisfies its boundary identities", {
  gl <- psychnets:::.gauss_legendre(60L)
  pb <- function(h, k, rho) psychnets:::.pbivnorm(h, k, rho, gl)
  # rho = 0 factorises; +/-Inf margins collapse to the univariate normal
  expect_equal(pb(0.4, -0.7, 0), stats::pnorm(0.4) * stats::pnorm(-0.7))
  expect_equal(pb(Inf, 1.1, 0.6), stats::pnorm(1.1))
  expect_equal(pb(0.3, Inf, -0.5), stats::pnorm(0.3))
  expect_equal(pb(-Inf, 0.2, 0.8), 0)
  # symmetric in its arguments
  expect_equal(pb(0.5, -0.3, 0.4), pb(-0.3, 0.5, 0.4))
})

test_that("polychoric recovers a known latent correlation from Likert data", {
  set.seed(1)
  n <- 3000
  rho <- 0.6
  z <- matrix(stats::rnorm(n * 2), n, 2) %*% chol(matrix(c(1, rho, rho, 1), 2))
  cut5 <- function(v) as.integer(cut(v, c(-Inf, -1, -0.3, 0.3, 1, Inf)))
  X <- cbind(cut5(z[, 1]), cut5(z[, 2]))
  R <- cor_auto(X)
  expect_lt(abs(R[1, 2] - rho), 0.05)          # near the latent truth
  expect_gt(R[1, 2], stats::cor(X[, 1], X[, 2]))  # and above the attenuated Pearson
})

test_that("cor_auto reduces to Pearson on continuous data", {
  set.seed(2)
  X <- matrix(stats::rnorm(500 * 4), 500, 4)    # continuous -> not ordinal
  expect_equal(cor_auto(X), stats::cor(X), ignore_attr = TRUE, tolerance = 1e-8)
})

test_that("cor_method = 'auto' is a valid PD correlation and changes the fit", {
  set.seed(3)
  z <- matrix(stats::rnorm(600 * 5), 600, 5) %*% chol(0.5^abs(outer(1:5, 1:5, "-")))
  X <- apply(z, 2, function(c) as.integer(cut(c, 5)))
  colnames(X) <- paste0("V", 1:5)
  auto <- ebic_glasso(X, cor_method = "auto")
  expect_s3_class(auto, "psychnet")
  expect_gt(min(eigen(auto$cor_matrix, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_false(isTRUE(all.equal(auto$weights, ebic_glasso(X, cor_method = "pearson")$weights)))
})

test_that("cor_auto matches psych::polychoric and qgraph::cor_auto", {
  if (!nzchar(Sys.getenv("PSYCHNET_EQUIV_TESTS"))) {
    skip("PSYCHNET_EQUIV_TESTS not set; gated external-equivalence suite skipped.")
  }
  skip_if_not_installed("psych")
  set.seed(5)
  n <- 1000
  Sig <- 0.5^abs(outer(1:4, 1:4, "-"))
  z <- matrix(stats::rnorm(n * 4), n, 4) %*% chol(Sig)
  X <- apply(z, 2, function(c) as.integer(cut(c, c(-Inf, -1, 0, 1, Inf))))
  ours <- cor_auto(X)
  ref <- suppressWarnings(psych::polychoric(X)$rho)
  expect_lt(max(abs(ours - ref)[upper.tri(ours)]), 1e-3)

  if (requireNamespace("qgraph", quietly = TRUE)) {
    qg <- tryCatch(
      suppressWarnings(suppressMessages(qgraph::cor_auto(X, verbose = FALSE))),
      error = function(e) {
        skip(paste("qgraph::cor_auto unavailable with the installed lavaan:",
                   conditionMessage(e)))
      })
    expect_lt(max(abs(ours - qg)[upper.tri(ours)]), 1e-3)
  }
})
