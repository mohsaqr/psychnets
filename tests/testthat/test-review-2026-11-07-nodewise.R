# Review 2026-11-07, nodewise estimators (ising_fit, ising_sampler, mgm_fit,
# the moderated MGM). Finding ids N1-N13 refer to the review; every error path
# is asserted by condition class, not by message.

.mod_data <- function(codes = 1:2, seed = 1, n = 400) {
  set.seed(seed)
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  g <- rep(1:2, each = n / 2)
  y <- x1 * (g == 2) + stats::rnorm(n)        # x1-y edge only in the 2nd group
  data.frame(x1 = x1, x2 = x2, y = y, mod = codes[g])
}

.common_factor_data <- function(seed = 1, n = 300, p = 8) {
  set.seed(seed)
  z <- matrix(stats::rnorm(n * p), n, p)
  as.data.frame(z + stats::rnorm(n))          # one shared component
}

.binary_data <- function(seed = 5, n = 600) {
  set.seed(seed)
  z <- matrix(stats::rnorm(n * 2), n, 2)
  x <- cbind(z[, 1], z[, 1], z[, 2], z[, 2]) + matrix(stats::rnorm(n * 4), n)
  b <- (x > 0) * 1
  colnames(b) <- paste0("V", 1:4)
  b
}

# --- N1: a categorical moderator's value is read in the user's coding --------

test_that("N1: condition() reads a 1/2-coded moderator in its own coding", {
  fit <- mgm_fit(.mod_data(1:2), types = c("g", "g", "g", "c"), moderators = 4)
  expect_identical(fit$moderator_levels, c("1", "2"))
  expect_equal(condition(fit, 1)$weights["x1", "y"], 0)
  expect_gt(condition(fit, 2)$weights["x1", "y"], 0.5)
})

test_that("N1: the conditioned network does not depend on the moderator coding", {
  fit01 <- mgm_fit(.mod_data(0:1), types = c("g", "g", "g", "c"), moderators = 4)
  fit12 <- mgm_fit(.mod_data(1:2), types = c("g", "g", "g", "c"), moderators = 4)
  fit59 <- mgm_fit(.mod_data(c(5, 9)), types = c("g", "g", "g", "c"),
                   moderators = 4)
  expect_identical(condition(fit01, 0)$weights, condition(fit12, 1)$weights)
  expect_identical(condition(fit01, 1)$weights, condition(fit12, 2)$weights)
  expect_identical(condition(fit01, 1)$weights, condition(fit59, 9)$weights)
  expect_identical(condition(fit01, 0)$weights, condition(fit59, 5)$weights)
})

test_that("N1: factor labels condition the glmnet moderated fit", {
  skip_if_not_installed("glmnet")
  d <- .mod_data(1:2)
  dl <- d
  dl$mod <- factor(c("ctrl", "treat")[d$mod])
  fit_l <- mgm_fit(dl, types = c("g", "g", "g", "c"), moderators = "mod",
                   native = FALSE)
  fit_n <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = 4,
                   native = FALSE)
  expect_identical(condition(fit_l, "treat")$weights,
                   condition(fit_n, 2)$weights)
  expect_identical(condition(fit_l, factor("ctrl"))$weights,
                   condition(fit_n, 1)$weights)
  expect_error(condition(fit_l, 1), class = "psychnet_bad_input")
})

test_that("N1: a value that is not a moderator level is refused, not read as the reference", {
  fit <- mgm_fit(.mod_data(1:2), types = c("g", "g", "g", "c"), moderators = 4)
  expect_error(condition(fit, 5), class = "psychnet_bad_input")
  expect_error(condition(fit, 0), class = "psychnet_bad_input")   # an internal code
  expect_error(condition(fit, "treat"), class = "psychnet_bad_input")
})

# --- N7: a continuous moderator's value is on its raw scale -------------------

test_that("N7: a continuous moderator is conditioned on its raw scale", {
  set.seed(2); n <- 400
  m <- stats::rnorm(n); x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  y <- 0.7 * x1 * m + stats::rnorm(n)
  d <- data.frame(x1 = x1, x2 = x2, y = y, m = m)
  d_shift <- d; d_shift$m <- 10 + 3 * d$m          # same moderator, new units
  fit <- mgm_fit(d, moderators = "m")
  fit_s <- mgm_fit(d_shift, moderators = "m")
  # invariance: the same person-level moderator value gives the same network
  for (v in c(-1, 0, 1.5)) {
    expect_equal(condition(fit, v)$weights, condition(fit_s, 10 + 3 * v)$weights,
                 tolerance = 1e-10)
  }
  # the raw mean is the no-moderation point (pure-interaction model)
  expect_equal(condition(fit_s, mean(d_shift$m))$weights["x1", "y"], 0,
               tolerance = 1e-8)
  expect_gt(condition(fit_s, 16)$weights["x1", "y"], 0.5)
  expect_equal(fit_s$moderator_center, mean(d_shift$m))
  expect_equal(fit_s$moderator_scale, stats::sd(d_shift$m))
})

test_that("N7: raw-scale condition() equals mgm::condition() at the z value", {
  skip_if(Sys.getenv("PSYCHNET_EQUIV_TESTS") == "")
  skip_if_not_installed("mgm")
  skip_if_not_installed("glmnet")
  set.seed(2); n <- 400
  m <- stats::rnorm(n, 10, 2); x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  d <- data.frame(x1 = x1, x2 = x2, y = x1 * (m - 10) / 2 + stats::rnorm(n),
                  m = m)
  fit <- mgm_fit(d, moderators = 4, native = FALSE)
  ref <- suppressWarnings(mgm::mgm(
    as.matrix(d), type = rep("g", 4), level = rep(1, 4), moderators = 4,
    lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND", threshold = "LW",
    overparameterize = FALSE, scale = TRUE, pbar = FALSE, signInfo = FALSE,
    warnings = FALSE))
  for (v in c(8, 10, 13)) {
    rr <- mgm::condition(ref, list("4" = (v - mean(m)) / stats::sd(m)))
    expect_lt(max(abs(unname(condition(fit, v)$weights) - rr$pairwise$wadj)),
              1e-10)
  }
})

# --- N9: condition() validates `value` ----------------------------------------

test_that("N9: condition() refuses a missing, vector, or non-numeric value", {
  fit <- mgm_fit(.mod_data(0:1), types = c("g", "g", "g", "c"), moderators = 4)
  expect_error(condition(fit, NA), class = "psychnet_bad_input")
  expect_error(condition(fit, c(0, 1)), class = "psychnet_bad_input")
  expect_error(condition(fit, numeric(0)), class = "psychnet_bad_input")
  set.seed(3)
  d <- data.frame(a = stats::rnorm(200), b = stats::rnorm(200),
                  c = stats::rnorm(200), m = stats::rnorm(200))
  fc <- mgm_fit(d, moderators = 4)
  expect_error(condition(fc, "1"), class = "psychnet_bad_input")
  expect_error(condition(fc, Inf), class = "psychnet_bad_input")
  expect_error(condition(list(), 1), class = "psychnet_bad_input")
})

# --- N10: as.data.frame() on a moderated fit ----------------------------------

test_that("N10: as.data.frame() returns main and moderation effects per edge", {
  fit <- mgm_fit(.mod_data(1:2), types = c("g", "g", "g", "c"), moderators = 4)
  tab <- as.data.frame(fit)
  expect_s3_class(tab, "data.frame")
  expect_named(tab, c("from", "to", "main_effect", "interaction", "moderated",
                      "involves_moderator"))
  xy <- subset(tab, from == "x1" & to == "y")
  expect_equal(nrow(xy), 1L)
  expect_true(xy$moderated)
  expect_gt(xy$interaction, 0)
  expect_true(all(is.na(subset(tab, involves_moderator)$interaction)))
  expect_true(all(tab$main_effect > 0 | tab$interaction > 0, na.rm = TRUE))
  # with a value it is the tidy edge list of the conditioned network
  expect_identical(as.data.frame(fit, value = 2),
                   as.data.frame(condition(fit, 2)))
})

# --- N2: observation weights are scale-free -----------------------------------

test_that("N2: mgm_fit is invariant to the scale of the weights", {
  d <- .mod_data(0:1)
  set.seed(4); w <- stats::runif(nrow(d), 0.2, 3)
  f1 <- mgm_fit(d, weights = w)
  expect_equal(mgm_fit(d, weights = 7 * w)$weights, f1$weights, tolerance = 1e-10)
  expect_equal(mgm_fit(d, weights = w / 1000)$weights, f1$weights,
               tolerance = 1e-10)
  # uniform weights of any size are the unweighted fit (rep(1/n, n) used to
  # select the empty network)
  f0 <- mgm_fit(d)
  expect_equal(mgm_fit(d, weights = rep(1 / nrow(d), nrow(d)))$weights,
               f0$weights, tolerance = 1e-10)
  expect_equal(mgm_fit(d, weights = rep(2, nrow(d)))$weights, f0$weights,
               tolerance = 1e-10)
  expect_gt(sum(f0$weights != 0), 0)
})

test_that("N2: ising_fit and ising_sampler are invariant to the scale of the weights", {
  b <- .binary_data()
  set.seed(6); w <- stats::runif(nrow(b), 0.5, 2)
  expect_equal(ising_fit(b, weights = 7 * w)$weights,
               ising_fit(b, weights = w)$weights, tolerance = 1e-10)
  expect_equal(ising_fit(b, weights = rep(1 / nrow(b), nrow(b)))$weights,
               ising_fit(b)$weights, tolerance = 1e-10)
  s1 <- ising_sampler(b, weights = w)
  s7 <- ising_sampler(b, weights = 7 * w)
  expect_equal(s7$weights, s1$weights, tolerance = 1e-8)
  expect_equal(s7$p_values, s1$p_values, tolerance = 1e-8)  # Wald SEs too
})

test_that("N2: invalid weights are refused by class", {
  b <- .binary_data()
  n <- nrow(b)
  bad <- list(rep(1, 3), c(Inf, rep(1, n - 1)), c(NA, rep(1, n - 1)),
              c(-1, rep(1, n - 1)), rep(0, n), rep("a", n))
  for (w in bad) {
    expect_error(ising_fit(b, weights = w), class = "psychnet_bad_input")
    expect_error(mgm_fit(as.data.frame(b), weights = w),
                 class = "psychnet_bad_input")
  }
})

# --- N3: a thresholded mgm network is not presented as certified ----------------

test_that("N3: when LW removes a coefficient, $kkt is NA and the fit residual is kept", {
  d <- .common_factor_data()
  fit <- mgm_fit(d)                                   # threshold = "LW"
  expect_true(is.na(fit$kkt))
  expect_lt(fit$fit_kkt, 1e-6)
  expect_identical(fit$certificate_target, "pre_threshold_fit")
  expect_true(is.na(certificate(fit)$certified))
  # unthresholded: certified, and kkt == fit_kkt
  full <- mgm_fit(d, threshold = "none")
  expect_identical(full$certificate_target, "returned_network")
  expect_identical(full$kkt, full$fit_kkt)
  expect_true(certificate(full)$certified)
  # a threshold that removed nothing leaves the certificate in place
  g <- golden_data_mixed()
  kept <- mgm_fit(g)
  expect_identical(kept$certificate_target, "returned_network")
  expect_true(certificate(kept)$certified)
})

test_that("N3: the glmnet mgm path follows the same certificate rule", {
  skip_if_not_installed("glmnet")
  d <- .common_factor_data()
  fit <- mgm_fit(d, native = FALSE)
  expect_true(is.na(fit$kkt))
  expect_true(is.finite(fit$fit_kkt))
  expect_identical(fit$certificate_target, "pre_threshold_fit")
  full <- mgm_fit(d, native = FALSE, threshold = "none")
  expect_identical(full$kkt, full$fit_kkt)
})

test_that("N3: the moderated fit reports a thresholded network as uncertified", {
  d <- .common_factor_data(p = 5)
  d$mod <- rep(0:1, length.out = nrow(d))
  fit <- mgm_fit(d, types = c(rep("g", 5), "c"), moderators = 6)
  full <- mgm_fit(d, types = c(rep("g", 5), "c"), moderators = 6,
                  threshold = "none")
  expect_identical(full$certificate_target, "returned_network")
  expect_lt(full$kkt, 1e-6)
  # on this design LW removes coefficients (verified), so it is uncertified
  expect_identical(fit$certificate_target, "pre_threshold_fit")
  expect_true(is.na(fit$kkt))
  expect_lt(fit$fit_kkt, 1e-6)
})

# --- N4: separation in the unregularized Ising model ---------------------------

test_that("N4: perfect separation warns, and is not certified", {
  set.seed(3)
  a <- stats::rbinom(200, 1, 0.5)
  X <- cbind(a, b = a, c = stats::rbinom(200, 1, 0.5),
             d = stats::rbinom(200, 1, 0.5))
  classes <- character(0)
  fit <- withCallingHandlers(ising_sampler(X), warning = function(w) {
    classes <<- c(classes, class(w))
    invokeRestart("muffleWarning")
  })
  expect_true("psychnet_separation" %in% classes)
  expect_true("psychnet_singular_information" %in% classes)
  expect_true(is.na(fit$kkt))
  expect_identical(fit$certificate_target, "nonconverged_fit")
  expect_true(is.na(certificate(fit)$certified))
})

test_that("N4: a well-posed unregularized Ising fit converges silently", {
  b <- .binary_data()
  expect_no_warning(fit <- ising_sampler(b))
  expect_identical(fit$certificate_target, "returned_network")
  expect_true(certificate(fit)$certified)
})

# --- N5: the gaussian nodewise lasso reaches its optimum -------------------------

test_that("N5: collinear gaussian nodes are solved to the optimum", {
  set.seed(1); n <- 300
  f <- stats::rnorm(n)
  d <- as.data.frame(vapply(1:10, function(i) f + stats::rnorm(n, sd = 0.15),
                            numeric(n)))
  expect_no_warning(fit <- mgm_fit(d, threshold = "none"))
  expect_lt(fit$kkt, 1e-9)                 # was 4.7e-05, certified FALSE
  expect_true(certificate(fit)$certified)
})

test_that("N5: the exact gaussian solver matches least squares at lambda = 0", {
  set.seed(9); n <- 120
  X <- .standardize(matrix(stats::rnorm(n * 4), n, 4))$X
  X[, 4] <- X[, 3] + 0.05 * X[, 4]               # strongly collinear pair
  X <- .standardize(X)$X
  y <- as.numeric(X %*% c(0.5, -0.2, 0.3, 0.1)) + stats::rnorm(n)
  fit <- .glm_lasso_fit(X, y, "gaussian", lambda = 0)
  ols <- stats::lm.fit(cbind(1, X), y)$coefficients
  expect_true(fit$converged)
  expect_equal(c(fit$b0, fit$beta), unname(ols), tolerance = 1e-10)
  # and the penalized optimum satisfies its own KKT conditions
  pen <- .glm_lasso_fit(X, y, "gaussian", lambda = 0.05)
  expect_lt(glm_lasso_kkt(X, y, pen$b0, pen$beta, 0.05), 1e-12)
})

test_that("N5: a nodewise fit that misses its optimum warns by class", {
  d <- golden_data_mixed()
  testthat::local_mocked_bindings(
    .gauss_lasso = function(X, z, w, lambda, b0, beta, gram = NULL, ...) {
      list(b0 = b0, beta = beta, converged = FALSE, iter = 1L)
    })
  expect_warning(mgm_fit(d), class = "psychnet_no_converge")
})

# --- N6: logical and two-level factor columns are binary nodes, never dropped -----

test_that("N6: logical and 2-level factor columns are coded 0/1, not dropped", {
  b <- as.data.frame(.binary_data())
  ref <- ising_fit(b)
  mixed <- b
  mixed$V1 <- b$V1 == 1
  mixed$V2 <- factor(c("no", "yes")[b$V2 + 1], levels = c("no", "yes"))
  fit <- ising_fit(mixed)
  expect_identical(fit$nodes$label, paste0("V", 1:4))
  expect_equal(fit$weights, ref$weights)
  # an all-logical data frame works like a logical matrix
  all_lgl <- as.data.frame(lapply(b, as.logical))
  expect_equal(ising_fit(all_lgl)$weights, ref$weights)
  expect_equal(ising_sampler(all_lgl)$weights, ising_sampler(b)$weights)
})

test_that("N6: a non-binary storage type is refused by name", {
  b <- as.data.frame(.binary_data())
  b3 <- b; b3$V1 <- factor(sample(c("x", "y", "z"), nrow(b), TRUE))
  expect_error(ising_fit(b3), "V1", class = "psychnet_bad_input")
  bc <- b; bc$V2 <- as.character(bc$V2)
  expect_error(ising_fit(bc), "V2", class = "psychnet_bad_input")
  expect_error(ising_sampler(bc), "V2", class = "psychnet_bad_input")
})

# --- N8 / N13: argument validation ---------------------------------------------

test_that("N8: moderators must be one whole column index or one column name", {
  d <- .mod_data(0:1)
  types <- c("g", "g", "g", "c")
  expect_error(mgm_fit(d, types = types, moderators = 2.5),
               class = "psychnet_bad_input")
  expect_error(mgm_fit(d, types = types, moderators = 0),
               class = "psychnet_bad_input")
  expect_error(mgm_fit(d, types = types, moderators = 9),
               class = "psychnet_bad_input")
  expect_error(mgm_fit(d, types = types, moderators = c(3, 4)),
               class = "psychnet_bad_input")
  expect_error(mgm_fit(d, types = types, moderators = "m"),
               class = "psychnet_bad_input")
  by_name <- mgm_fit(d, types = types, moderators = "mod")
  by_index <- mgm_fit(d, types = types, moderators = 4)
  expect_identical(condition(by_name, 1)$weights, condition(by_index, 1)$weights)
})

test_that("N13: nlambda, gamma and lambda_min_ratio are validated by class", {
  b <- .binary_data()
  d <- as.data.frame(b)
  for (f in list(function(...) ising_fit(b, ...), function(...) mgm_fit(d, ...))) {
    expect_error(f(nlambda = 2.5), class = "psychnet_bad_input")
    expect_error(f(nlambda = 1), class = "psychnet_bad_input")
    expect_error(f(gamma = -1), class = "psychnet_bad_input")
    expect_error(f(gamma = c(0.1, 0.2)), class = "psychnet_bad_input")
    expect_error(f(lambda_min_ratio = 1), class = "psychnet_bad_input")
    expect_error(f(lambda_min_ratio = 0), class = "psychnet_bad_input")
  }
  expect_error(ising_sampler(b, alpha = 1.5), class = "psychnet_bad_input")
  expect_error(mgm_fit(d, types = c("g", "x", "g", "g")),
               class = "psychnet_bad_input")
})

# --- N9: binary-input errors are classed and name the caller ----------------------

test_that("N9: 1/2-coded data is refused by class, naming the caller and columns", {
  b <- .binary_data() + 1
  expect_error(ising_fit(b), "ising_fit\\(\\).*V1", class = "psychnet_bad_input")
  expect_error(ising_sampler(b), "ising_sampler\\(\\)",
               class = "psychnet_bad_input")
})

test_that("N9: labels for every input column survive a dropped constant column", {
  b <- as.data.frame(.binary_data())
  b$k <- 1
  lab <- c("a", "b", "c", "d", "konst")
  quiet <- function(expr) withCallingHandlers(expr,
    psychnet_constant_variable = function(w) invokeRestart("muffleWarning"))
  expect_identical(quiet(ising_fit(b, labels = lab))$nodes$label, lab[1:4])
  expect_identical(quiet(ising_sampler(b, labels = lab))$nodes$label, lab[1:4])
  expect_error(quiet(ising_fit(b, labels = c("a", "b"))),
               class = "psychnet_bad_input")
})

# --- N11: a two-valued column is not a "type promotion" ---------------------------

test_that("N11: a 1/2-coded binary is not warned about; a 3-valued integer is", {
  set.seed(8); n <- 300
  d <- data.frame(g1 = stats::rnorm(n), g2 = stats::rnorm(n),
                  b12 = sample(1:2, n, TRUE))
  expect_no_warning(f12 <- mgm_fit(d), class = "psychnet_type_promotion")
  d01 <- d; d01$b12 <- d$b12 - 1L
  expect_equal(f12$weights, mgm_fit(d01)$weights)
  d3 <- d; d3$b12 <- sample(1:3, n, TRUE)
  expect_warning(.detect_types(d3), class = "psychnet_type_promotion")
})

test_that("glm_lasso_kkt() refuses mismatched inputs by class", {
  X <- matrix(stats::rnorm(20), 10, 2)
  expect_error(glm_lasso_kkt(X, 1:9, 0, c(0, 0), 0), class = "psychnet_bad_input")
  expect_error(glm_lasso_kkt(X, 1:10, 0, 0, 0), class = "psychnet_bad_input")
  expect_error(glm_lasso_kkt(X, 1:10, 0, c(0, 0), 0, weights = 1:3),
               class = "psychnet_bad_input")
})
