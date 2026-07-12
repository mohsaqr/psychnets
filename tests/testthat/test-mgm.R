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

.gen_moderated <- function(n = 400, seed = 1) {
  set.seed(seed)
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  mod <- rep(0:1, each = n / 2)
  y <- x1 * (mod == 1) + stats::rnorm(n) * 0.5      # x1-y edge only when mod == 1
  data.frame(x1 = x1, x2 = x2, y = y, mod = mod)
}

test_that("moderated MGM fits and conditions, recovering a moderated edge", {
  d <- .gen_moderated()
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

test_that("the moderated MGM runs on the base kernel with no glmnet, and certifies", {
  d <- .gen_moderated()
  fit <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = 4)  # native = TRUE
  expect_identical(fit$engine, "base")
  # the base engine must never reach glmnet: every nodewise fit is KKT-certified
  # against its own convex objective instead.
  expect_true(is.finite(fit$kkt))
  expect_lt(fit$kkt, 1e-6)
})

test_that("a moderated edge scales with the value of a continuous moderator", {
  set.seed(2)
  n <- 400
  mod <- stats::rnorm(n); x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  y <- 0.7 * x1 * mod + stats::rnorm(n)             # pure interaction, no main effect
  d <- data.frame(x1 = x1, x2 = x2, y = y, mod = mod)
  fit <- mgm_fit(d, types = c("g", "g", "g", "g"), moderators = 4)
  w <- vapply(c(0, 1, 2), function(v) condition(fit, v)$weights["x1", "y"],
              numeric(1))
  expect_equal(w[1], 0, tolerance = 1e-8)           # zero at mod = 0
  expect_gt(w[2], 0)                                # grows with |mod|
  expect_gt(w[3], w[2])
})

test_that("the base moderated kernel rejects a multi-level categorical", {
  d <- .gen_moderated()
  d$x2 <- sample(0:2, nrow(d), replace = TRUE)     # 3-level categorical
  expect_error(
    mgm_fit(d, types = c("g", "c", "g", "c"), moderators = 4),
    "more than 2 levels"
  )
})

test_that("the base moderated kernel matches the glmnet reference in structure", {
  skip_if_not_installed("glmnet")
  d <- .gen_moderated()
  base <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = 4)
  glmn <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = 4,
                  native = FALSE)
  expect_identical(glmn$engine, "glmnet")
  Wb <- condition(base, 1)$weights
  Wg <- condition(glmn, 1)$weights
  expect_identical(Wb != 0, Wg != 0)                # same edge set
  expect_equal(Wb, Wg, tolerance = 0.05)            # independent solver, so magnitudes only close
})
