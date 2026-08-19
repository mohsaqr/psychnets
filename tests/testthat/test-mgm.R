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

# --- review fixes 2026-08-17: NA, constant columns, certificates, overrides ---

test_that("glmnet path handles missing data in every column type and na_method", {
  skip_if_not_installed("glmnet")
  d <- .gen_cat(7, klev = 3L)
  d$fc[c(3, 10, 40)] <- NA          # categorical NA
  d$g1[c(5, 60)]     <- NA          # gaussian NA
  d$b1[c(8)]         <- NA          # binary NA
  fit_pw <- mgm_fit(d, native = FALSE)                          # pairwise
  expect_equal(fit_pw$n, nrow(d))
  expect_equal(fit_pw$nodes$label, c("g1", "g2", "b1", "fc"))
  fit_lw <- mgm_fit(d, native = FALSE, na_method = "listwise")  # listwise
  expect_equal(fit_lw$n, sum(stats::complete.cases(d)))
  expect_equal(fit_lw$nodes$label, c("g1", "g2", "b1", "fc"))
  # complete data is byte-identical under both, as for every estimator
  dc <- .gen_cat(7, klev = 3L)
  expect_identical(mgm_fit(dc, native = FALSE)$weights,
                   mgm_fit(dc, native = FALSE, na_method = "listwise")$weights)
})

test_that("a missing k-level code is imputed with the modal level, never a mean", {
  m <- cbind(g = c(1.5, NA, 2.5, 3.5), c3 = c(2, 2, 0, NA), b = c(0, 1, NA, 1))
  out <- .na_prep_nodewise(m, "pairwise", categorical = c(FALSE, TRUE, TRUE))
  expect_equal(unname(out[4, "c3"]), 2)         # mode of {2, 2, 0}, not 4/3
  expect_equal(unname(out[3, "b"]), 1)          # 0/1 rule unchanged (mean >= .5)
  expect_equal(unname(out[2, "g"]), mean(c(1.5, 2.5, 3.5)))
  # without the flag the historical mean rule still applies
  expect_equal(unname(.na_prep_nodewise(m, "pairwise")[4, "c3"]), 4 / 3)
})

test_that("a constant or all-NA column is dropped, on both engines, without error", {
  d <- .gen_cat(8, klev = 3L)[, c("g1", "g2", "b1")]
  d$k_int <- 1L                                  # constant integer
  d$k_dbl <- 2.5                                 # constant double
  d$k_na  <- NA_real_                            # all missing
  fit <- mgm_fit(d)                              # no promotion warning, no error
  expect_equal(fit$nodes$label, c("g1", "g2", "b1"))
  expect_equal(unname(fit$types), c("g", "g", "c"))
  # an explicit `types` covering the original columns is subset alongside
  fit_t <- mgm_fit(d, types = c("g", "g", "c", "g", "g", "g"))
  expect_equal(fit_t$nodes$label, c("g1", "g2", "b1"))
  skip_if_not_installed("glmnet")
  fit_g <- mgm_fit(d, native = FALSE)
  expect_equal(fit_g$nodes$label, c("g1", "g2", "b1"))
})

test_that("the multi-level certificate is the multinomial (softmax) residual", {
  skip_if_not_installed("glmnet")
  # A 3-level node with a real dependence on g1, so its fit has active
  # coefficients: grading it per class as a logistic regression reported ~0.2
  # for a fit glmnet had converged; the multinomial residual is ~1e-5.
  set.seed(11); n <- 400
  f <- sample(c("lo", "mid", "hi"), n, TRUE)
  d <- data.frame(g1 = stats::rnorm(n) + c(lo = -1, mid = 0, hi = 1)[f],
                  g2 = stats::rnorm(n), b = stats::rbinom(n, 1, 0.5),
                  f = factor(f, levels = c("lo", "mid", "hi")))
  fit <- mgm_fit(d, native = FALSE)
  expect_gt(abs(fit$weights["f", "g1"]), 0.5)
  expect_lt(fit$kkt, 1e-4)
  expect_true(certificate(fit, tol = 1e-4)$certified)
})

test_that("an empty selected model certifies as optimal, not as a violation", {
  skip_if_not_installed("glmnet")
  # Independent binaries: EBIC selects the null model for every node. The
  # self-lambda heuristic has no active gradient to read the penalty from; the
  # empty solution is exactly optimal at lambda_max, so the residual is ~0.
  set.seed(2)
  bd <- as.data.frame(matrix(stats::rbinom(300 * 4, 1, 0.5), 300, 4))
  fit <- ising_fit(bd, native = FALSE)
  expect_equal(sum(fit$weights != 0), 0)
  expect_lt(fit$kkt, 1e-8)
})

test_that("types = 'c' on a detected-gaussian column keeps its true level count", {
  set.seed(4); n <- 200
  d <- data.frame(g = stats::rnorm(n), x = sample(c(1.5, 2.5, 3.5), n, TRUE),
                  b = stats::rbinom(n, 1, 0.5))
  # x is non-integer, so detection says gaussian; the caller overrides to "c".
  # It has 3 levels, so the base kernel must refuse it by name ...
  expect_error(mgm_fit(d, types = c("g", "c", "c")), "x")
  skip_if_not_installed("glmnet")
  # ... and the glmnet path must model it as multi-level (levels = 3) and
  # therefore refuse net_predict(), not silently predict from a zero row.
  fit <- mgm_fit(d, types = c("g", "c", "c"), native = FALSE)
  expect_equal(unname(fit$levels), c(1L, 3L, 2L))
  expect_error(net_predict(fit, data = d), "multi-level")
})

test_that("a high-cardinality character column is refused as an identifier", {
  skip_if_not_installed("glmnet")
  d <- .gen_cat(9, klev = 3L)
  d$id <- sprintf("id%03d", seq_len(nrow(d)))
  expect_error(mgm_fit(d, native = FALSE), "id")
  expect_error(mgm_fit(d, native = FALSE), "factor\\(\\)")
})

test_that("the promotion warning carries a condition class the resamplers can muffle", {
  d <- .gen_cat(10, klev = 3L)
  d$fc <- as.integer(d$fc)
  expect_warning(expect_error(mgm_fit(d), "native = FALSE"),
                 class = "psychnet_type_promotion")
})

test_that("dropping a constant column re-targets moderators and labels correctly", {
  d <- .gen_cat(12, klev = 3L)[, c("g1", "g2", "b1")]
  d <- data.frame(k = 1, d)                        # constant FIRST column
  # moderators is a positional index into the ORIGINAL columns: 3 = g2.
  fit <- mgm_fit(d, moderators = 3)
  ref <- mgm_fit(d[, -1], moderators = 2)
  expect_equal(condition(fit, value = 1)$weights, condition(ref, value = 1)$weights)
  # a moderator that is itself the constant column is refused by name
  expect_error(mgm_fit(d, moderators = 1), "fewer than two distinct")
  # labels supplied for every original column are subset alongside
  fit_l <- mgm_fit(d, labels = c("K", "A", "B", "C"))
  expect_equal(fit_l$nodes$label, c("A", "B", "C"))
})

test_that("an explicit factor() is the escape hatch for a > 10-level categorical", {
  skip_if_not_installed("glmnet")
  set.seed(13); n <- 600
  d <- data.frame(g1 = stats::rnorm(n), g2 = stats::rnorm(n),
                  many = factor(sample(letters[1:12], n, TRUE)))
  fit <- mgm_fit(d, native = FALSE)                # accepted, one 12-level node
  expect_equal(unname(fit$levels), c(1L, 1L, 12L))
  # ... whereas the same values as character, or as declared 'c' integers, are
  # refused
  dc <- d; dc$many <- as.character(dc$many)
  expect_error(mgm_fit(dc, native = FALSE), "identifiers")
  di <- d; di$many <- as.integer(di$many)
  expect_error(mgm_fit(di, types = c("g", "g", "c"), native = FALSE),
               "more than 10 distinct")
})

test_that("listwise deletion that leaves a one-level node drops it on both engines", {
  d <- .gen_cat(14, klev = 3L)
  # every row where fc != 2 has a missing g1: listwise leaves fc constant
  d$g1[d$fc != 2] <- NA
  fit_b <- mgm_fit(d[, c("g1", "g2", "b1")], na_method = "listwise")
  expect_equal(fit_b$nodes$label, c("g1", "g2", "b1"))
  skip_if_not_installed("glmnet")
  fit_g <- mgm_fit(d, native = FALSE, na_method = "listwise")
  expect_equal(fit_g$nodes$label, c("g1", "g2", "b1"))  # fc dropped, no error
  expect_equal(fit_g$n, sum(stats::complete.cases(d)))
})
