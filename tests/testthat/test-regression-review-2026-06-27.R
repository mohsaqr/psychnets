# Regression tests for the 2026-06-27 package-review fixes. Each test is written
# to FAIL on the pre-fix code and pass on the fix:
#  (A) ising_fit / mgm_fit must store the WEIGHTED centers/scales (the same ones
#      the nodewise fits standardize on), not the unweighted moments.
#  (B) cor_auto must not emit NaN for a degenerate ordinal/continuous pair that
#      collapses to a single level under pairwise NA handling.

# --- Helpers -----------------------------------------------------------------

# Weighted population center/scale, matching lasso_glm.R `.standardize(X, w)`.
.weighted_moments <- function(mat, w) {
  sw  <- sum(w)
  ctr <- colSums(w * mat) / sw
  matc <- sweep(mat, 2L, ctr, "-")
  scl <- sqrt(colSums(w * matc^2) / sw)
  scl[scl < 1e-12] <- 1
  list(center = ctr, scale = scl)
}

.gen_binary <- function(seed, n = 400, p = 4) {
  set.seed(seed)
  f1 <- stats::rnorm(n); f2 <- stats::rnorm(n)
  lin <- cbind(f1, f1, f2, f2)[, seq_len(p)] +
    matrix(stats::rnorm(n * p, sd = 0.5), n, p)
  b <- (lin > 0) * 1L
  colnames(b) <- paste0("V", seq_len(p))
  b
}

# --- (A) weighted standardization -------------------------------------------

test_that("ising_fit stores the weighted centers/scales under non-uniform weights", {
  b <- .gen_binary(11)
  set.seed(101)
  w <- stats::runif(nrow(b), 0.3, 2)             # genuinely non-uniform
  fit <- ising_fit(b, weights = w)
  expected <- .weighted_moments(b, w)

  expect_equal(unname(fit$nodewise$center), unname(expected$center), tolerance = 1e-10)
  expect_equal(unname(fit$nodewise$scale),  unname(expected$scale),  tolerance = 1e-10)
  # Guard against the pre-fix behavior: the UNWEIGHTED mean must NOT match, so a
  # silent revert to `.standardize(mat)` re-breaks this test.
  expect_false(isTRUE(all.equal(unname(fit$nodewise$center), unname(colMeans(b)))))
})

test_that("mgm_fit stores the weighted centers/scales under non-uniform weights", {
  b <- .gen_binary(12)                            # all-binary: composite == std
  set.seed(102)
  w <- stats::runif(nrow(b), 0.3, 2)
  fit <- mgm_fit(b, weights = w)
  expected <- .weighted_moments(b, w)

  expect_equal(unname(fit$nodewise$center), unname(expected$center), tolerance = 1e-10)
  expect_equal(unname(fit$nodewise$scale),  unname(expected$scale),  tolerance = 1e-10)
  expect_false(isTRUE(all.equal(unname(fit$nodewise$center), unname(colMeans(b)))))
})

test_that("weights = NULL is byte-identical to uniform weights (equivalence held)", {
  b <- .gen_binary(13)
  i0 <- ising_fit(b);            i1 <- ising_fit(b, weights = rep(1, nrow(b)))
  m0 <- mgm_fit(b);              m1 <- mgm_fit(b, weights = rep(1, nrow(b)))
  expect_equal(i0$weights, i1$weights)
  expect_equal(unname(i0$thresholds), unname(i1$thresholds))
  expect_equal(unname(i0$nodewise$center), unname(i1$nodewise$center))
  expect_equal(m0$weights, m1$weights)
  expect_equal(unname(m0$nodewise$center), unname(m1$nodewise$center))
})

# --- (B) cor_auto degenerate pair -------------------------------------------

test_that(".polyserial returns 0 for a single-level ordinal or constant partner", {
  set.seed(201)
  expect_equal(psychnets:::.polyserial(stats::rnorm(20), rep(2L, 20)), 0)  # 1 level
  expect_equal(psychnets:::.polyserial(rep(1, 20), sample(1:4, 20, TRUE)), 0)  # const cont
})

test_that("cor_auto yields a finite, symmetric matrix for a degenerate pair", {
  set.seed(202)
  n <- 120
  ord <- sample(1:4, n, replace = TRUE)          # ordinal (few integer levels)
  y   <- stats::rnorm(n)
  y[ord != 2] <- NA                              # y observed only where ord == 2
  x   <- stats::rnorm(n)                          # a well-behaved continuous anchor
  d   <- data.frame(x = x, ord = ord, y = y)

  R <- cor_auto(d)                               # pre-fix: NaN -> eigen() error
  expect_false(is.null(dimnames(R)))             # must return a named matrix
  expect_true(all(is.finite(R)))
  expect_equal(unname(R), unname(t(R)), tolerance = 1e-12)
  expect_equal(unname(diag(R)), rep(1, 3), tolerance = 1e-12)
  # the degenerate ord<->y pair carries no estimable association (allow a tiny
  # nudge from the nearest-PD projection)
  expect_lt(abs(R["ord", "y"]), 1e-8)
})
