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

test_that("moderated MGM fits and conditions, recovering a moderated edge", {
  skip_if_not_installed("glmnet")
  set.seed(1)
  n <- 400
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  mod <- rep(0:1, each = n / 2)
  y <- x1 * (mod == 1) + stats::rnorm(n) * 0.5      # x1-y edge only when mod == 1
  d <- data.frame(x1 = x1, x2 = x2, y = y, mod = mod)
  fit <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = 4)
  expect_s3_class(fit, "psychnet_moderated")
  net0 <- condition(fit, 0)
  net1 <- condition(fit, 1)
  expect_s3_class(net1, "psychnet")
  expect_true(inherits(net1, "cograph_network"))
  # the x1-y edge (nodes 1,3) is present at mod = 1 and (near) absent at mod = 0
  expect_gt(abs(net1$weights[1, 3]), abs(net0$weights[1, 3]))
  # the moderator node carries no edges
  expect_equal(unname(net1$weights[4, ]), rep(0, 4))
})
