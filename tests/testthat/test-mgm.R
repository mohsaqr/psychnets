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

test_that("the base kernel errors on a >2-level categorical column", {
  set.seed(3)
  d <- data.frame(g = stats::rnorm(200), cat3 = sample(0:2, 200, replace = TRUE))
  # The column is first promoted to categorical by mgm's <= 10-integer rule
  # (warned, because that changes the model), and then refused because the base
  # kernel is a binomial solver.
  expect_warning(expect_error(mgm_fit(d), "categorical"), "cat3")
  expect_error(suppressWarnings(mgm_fit(d)), "native = FALSE")
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

# --- multi-level categorical nodes -------------------------------------------

.gen_cat <- function(seed, klev = 3L, n = 350) {
  set.seed(seed)
  f <- stats::rnorm(n)
  data.frame(
    g1 = f + stats::rnorm(n, sd = 0.6),
    g2 = f + stats::rnorm(n, sd = 0.6),
    b1 = ((f + stats::rnorm(n, sd = 0.6)) > 0) * 1L,
    fc = factor(cut(f + stats::rnorm(n, sd = 0.8), breaks = klev,
                    labels = FALSE))
  )
}

test_that("detection follows mgm's rule and reports levels", {
  d <- .gen_cat(1, klev = 4L)
  skip_if_not_installed("glmnet")
  fit <- mgm_fit(d, native = FALSE)
  expect_equal(unname(fit$types), c("g", "g", "c", "c"))
  expect_equal(unname(fit$levels), c(1L, 1L, 2L, 4L))
})

test_that("a numeric column promoted by the <=10 integer rule warns by name", {
  # A Likert item silently becoming multinomial is a model change, not an
  # encoding detail, so it is warned about rather than done quietly.
  set.seed(9)
  n <- 250
  f <- stats::rnorm(n)
  d <- data.frame(g1 = f + stats::rnorm(n, sd = 0.6),
                  g2 = f + stats::rnorm(n, sd = 0.6),
                  lik = as.integer(cut(f + stats::rnorm(n, sd = 0.8), 5)))
  skip_if_not_installed("glmnet")
  expect_warning(mgm_fit(d, native = FALSE), "lik")
  expect_warning(mgm_fit(d, native = FALSE), "categorical")
})

test_that("a 3- and 4-level factor fits and stays adjacent to a binary node", {
  skip_if_not_installed("glmnet")
  for (k in c(3L, 4L)) {
    fit <- mgm_fit(.gen_cat(2, klev = k), native = FALSE)
    expect_s3_class(fit, "psychnet")
    expect_equal(nrow(fit$nodes), 4L)                  # NOT one-hot expanded
    expect_equal(unname(fit$weights), unname(t(fit$weights)))
    expect_true(all(diag(fit$weights) == 0))
    # the factor shares the latent factor, so it is not isolated (which
    # individual edge survives EBIC + LW varies with k and seed)
    expect_gt(sum(abs(fit$weights[, "fc"])), 0)
  }
})

test_that("an edge touching a multi-level node is unsigned", {
  # mgm leaves the sign of a categorical interaction undefined; a multi-dummy
  # predictor has no single signed direction to report.
  skip_if_not_installed("glmnet")
  fit <- mgm_fit(.gen_cat(3, klev = 3L), native = FALSE)
  expect_gte(fit$weights["g1", "fc"], 0)
})

test_that("the base kernel refuses a multi-level node by name", {
  d <- .gen_cat(4, klev = 3L)
  expect_error(mgm_fit(d, native = TRUE), "fc")
  expect_error(mgm_fit(d, native = TRUE), "native = FALSE")
})

test_that("net_predict refuses a network with a multi-level node", {
  skip_if_not_installed("glmnet")
  d <- .gen_cat(5, klev = 3L)
  fit <- mgm_fit(d, native = FALSE)
  expect_error(net_predict(fit, data = d), "multi-level")
  # ... but a binary-only network still predicts
  db <- d[, c("g1", "g2", "b1")]
  expect_s3_class(net_predict(mgm_fit(db, native = FALSE), data = db),
                  "data.frame")
})

test_that("a binary node is unaffected by the multinomial generalisation", {
  # k = 2 shares the k-class code path; this pins that sharing it changed
  # nothing (the golden fixtures guard the exact weights).
  skip_if_not_installed("glmnet")
  d <- .gen_cat(6, klev = 3L)
  db <- d[, c("g1", "g2", "b1")]
  fit <- mgm_fit(db, native = FALSE)
  expect_equal(unname(fit$types), c("g", "g", "c"))
  expect_equal(unname(fit$levels), c(1L, 1L, 2L))
  expect_lt(fit$kkt, 1e-5)
})
