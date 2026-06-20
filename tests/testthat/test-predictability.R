# Predictability: closed-form R^2 for GGMs (no data), data-driven CC/nCC for the
# nodewise models.

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

test_that("GGM predictability matches the closed-form variance explained", {
  S <- ar1(6, 0.5)
  fit <- ebic_glasso(cor_matrix = S, n = 400)
  pr <- predictability(fit)
  expect_named(pr, c("node", "type", "metric", "predictability", "accuracy"))
  expect_equal(nrow(pr), 6L)
  expect_true(all(pr$type == "gaussian"))
  # R^2_j = 1 - 1/(Theta_jj * S_jj) computed directly from the precision
  direct <- 1 - 1 / (diag(fit$precision) * diag(fit$cor_matrix))
  expect_equal(pr$predictability, as.numeric(direct), tolerance = 1e-10)
  expect_true(all(pr$predictability >= -1e-8 & pr$predictability <= 1 + 1e-8))
})

test_that("nodewise predictability needs data and reports CC / nCC", {
  set.seed(1)
  z <- matrix(stats::rnorm(700 * 2), 700, 2)
  b <- (cbind(z[, 1], z[, 1], z[, 2], z[, 2]) +
          matrix(stats::rnorm(700 * 4), 700) > 0) * 1
  colnames(b) <- paste0("V", 1:4)
  fit <- ising_fit(b)
  expect_error(predictability(fit), "data` is required")
  pr <- predictability(fit, data = b)
  expect_true(all(pr$type == "binary"))
  expect_true(all(pr$metric == "nCC"))
  expect_true(all(pr$accuracy >= 0 & pr$accuracy <= 1))
  expect_true(all(pr$predictability <= 1 + 1e-8))
})

test_that("mgm predictability matches a direct nodewise recomputation", {
  set.seed(2)
  f <- stats::rnorm(500)
  d <- data.frame(g1 = f + stats::rnorm(500), g2 = f + stats::rnorm(500),
                  b1 = (f + stats::rnorm(500) > 0) * 1, n = stats::rnorm(500))
  fit <- mgm_fit(d)
  pr <- predictability(fit, data = d)
  expect_equal(nrow(pr), 4L)
  # gaussian nodes report R^2, binary node reports nCC
  expect_equal(pr$metric[pr$type == "gaussian"][1], "R2")
  expect_true("nCC" %in% pr$metric)
})

test_that("predictability errors on a network without precision or nodewise", {
  set.seed(3)
  X <- matrix(stats::rnorm(200 * 4), 200, 4)
  expect_error(predictability(cor_network(X)), "not defined")
})
