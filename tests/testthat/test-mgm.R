# Mixed graphical model: gaussian + binary nodes on a shared latent factor.
# Correctness is certified per node by the GLM stationarity residual.

.gen_mixed <- function(seed, n = 500) {
  set.seed(seed)
  f <- stats::rnorm(n)
  data.frame(
    g1 = f + stats::rnorm(n, sd = 0.6),
    g2 = f + stats::rnorm(n, sd = 0.6),
    b1 = ((f + stats::rnorm(n, sd = 0.6)) > 0) * 1L,
    noise = stats::rnorm(n)
  )
}

test_that("mgm_fit detects node types and reaches the nodewise optimum", {
  d <- .gen_mixed(1)
  fit <- mgm_fit(d)
  expect_s3_class(fit, "psychnet")
  expect_equal(unname(fit$types), c("g", "g", "c", "g"))
  expect_lt(fit$kkt, 1e-6)
})

test_that("mgm_fit recovers the latent-factor cluster across types", {
  d <- .gen_mixed(2)
  fit <- mgm_fit(d)
  g <- fit$graph
  # g1, g2, b1 load on the same factor -> connected; noise stays isolated-ish
  expect_gt(abs(g["g1", "g2"]), 0)
  expect_equal(unname(fit$graph), unname(t(fit$graph)))
  expect_gt(abs(g["g1", "g2"]) + abs(g["g1", "b1"]),
            abs(g["g1", "noise"]) + abs(g["g2", "noise"]))
})

test_that("mgm_fit errors on a >2-level categorical column", {
  set.seed(3)
  d <- data.frame(g = stats::rnorm(200), cat3 = sample(0:2, 200, replace = TRUE))
  expect_error(mgm_fit(d), "categorical")
})

test_that("psychnet routes to ising and mgm", {
  b <- (matrix(stats::rnorm(400 * 4), 400, 4) > 0) * 1L
  colnames(b) <- paste0("V", 1:4)
  expect_equal(psychnet(b, "ising")$method, "ising")
  expect_equal(psychnet(.gen_mixed(4), "mgm")$method, "mgm")
})
