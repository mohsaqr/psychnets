# Correctness is certified against the convex objective itself (KKT residual),
# so these tests need no reference solver.

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))
compound <- function(p, rho) { m <- matrix(rho, p, p); diag(m) <- 1; m }

test_that("ebic_glasso returns the certified global optimum", {
  for (S in list(ar1(6, 0.5), ar1(8, 0.6), compound(5, 0.4))) {
    fit <- ebic_glasso(cor_matrix = S, n = 250)
    expect_s3_class(fit, "psychnet")
    expect_lt(fit$kkt, 1e-7)
    expect_equal(glasso_kkt(fit$precision, S, fit$lambda), fit$kkt)
  }
})

test_that("glasso_kkt flags a non-optimal precision matrix", {
  S <- ar1(5, 0.5)
  fit <- ebic_glasso(cor_matrix = S, n = 250)
  bad <- fit$precision
  bad[1, 2] <- bad[2, 1] <- bad[1, 2] + 0.3
  expect_gt(glasso_kkt(bad, S, fit$lambda),
            glasso_kkt(fit$precision, S, fit$lambda))
})

test_that("ebic_glasso graph is a symmetric partial-correlation matrix", {
  S <- ar1(7, 0.5)
  fit <- ebic_glasso(cor_matrix = S, n = 300)
  g <- fit$graph
  expect_equal(g, t(g))
  expect_true(all(diag(g) == 0))
  expect_true(all(abs(g) <= 1 + 1e-8))
})

test_that("ebic_glasso runs from raw data and respects threshold", {
  set.seed(1)
  p <- 6
  X <- matrix(stats::rnorm(400 * p), 400, p) %*% chol(ar1(p, 0.5))
  colnames(X) <- paste0("V", seq_len(p))
  fit  <- ebic_glasso(as.data.frame(X), gamma = 0.5)
  fitT <- ebic_glasso(as.data.frame(X), gamma = 0.5, threshold = 0.05)
  expect_lt(fit$kkt, 1e-6)
  expect_true(all(abs(fitT$graph[fitT$graph != 0]) >= 0.05))
  expect_lte(fitT$n_edges, fit$n_edges)
})

test_that("node order is the input column order, never sorted", {
  # A documented promise, not an accident: consumers reindex defensively
  # (m[nodes, nodes]) precisely because it was only ever observed behaviour.
  nm <- c("zeta", "alpha", "mid", "Beta", "q1", "a10", "a2", "Z")
  S <- ar1(length(nm), 0.45)
  dimnames(S) <- list(nm, nm)
  fit <- ebic_glasso(cor_matrix = S, n = 300)
  expect_identical(colnames(fit$weights), nm)
  expect_identical(rownames(fit$weights), nm)
  expect_identical(colnames(fit$precision), nm)
  expect_identical(fit$nodes$label, nm)

  set.seed(4)
  X <- matrix(stats::rnorm(300 * length(nm)), 300, length(nm))
  colnames(X) <- nm
  expect_identical(colnames(ebic_glasso(as.data.frame(X))$weights), nm)
  expect_identical(colnames(cor_network(as.data.frame(X))$weights), nm)
  expect_identical(colnames(mgm_fit(as.data.frame(X))$weights), nm)
})

# --- resampling support (lambda_path / refit / penalize_diagonal) ------------

fixed_path <- function(S, nlambda = 60L) {
  lmax <- max(abs(S[upper.tri(S)]))
  exp(seq(log(lmax), log(lmax * 0.01), length.out = nlambda))
}

test_that("lambda_path overrides the computed grid and is echoed unchanged", {
  S <- ar1(6, 0.5)
  lam <- fixed_path(S)
  fit <- ebic_glasso(cor_matrix = S, n = 250, lambda_path = lam)
  expect_identical(fit$lambda_path, lam)
  expect_true(fit$lambda %in% lam)
  expect_lt(fit$kkt, 1e-7)
})

test_that("a supplied lambda_path keeps penalties comparable across resamples", {
  # The point of the argument: every resample is scored on the SAME grid, so a
  # selected lambda means the same thing across draws. Recomputing the path per
  # resample lets the grid drift with each sample's largest |correlation|.
  set.seed(21)
  p <- 6; n <- 200
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(ar1(p, 0.5))
  colnames(X) <- paste0("V", seq_len(p))
  lam <- fixed_path(stats::cor(X))
  sel <- vapply(seq_len(5), function(i) {
    Sb <- stats::cor(X[sample(n, replace = TRUE), ])
    ebic_glasso(cor_matrix = Sb, n = n, lambda_path = lam, refit = FALSE)$lambda
  }, numeric(1))
  expect_true(all(sel %in% lam))
})

test_that("a non-PD path fit is skipped, not raised", {
  # Resampling hits near-singular draws routinely; the caller must get a fit
  # back rather than an error it has to tryCatch away.
  set.seed(3)
  p <- 15; n <- 9                                  # n < p: S is singular
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))
  S <- stats::cor(X)
  lam <- exp(seq(log(0.9), log(1e-4), length.out = 40))
  sel <- .select_ebic(S, lam, n, 0.5, refit = FALSE)
  expect_gt(sum(!is.finite(sel$ebic_path)), 0)     # at least one skipped
  expect_true(is.finite(sel$lambda))               # selection still succeeded
})

test_that("refit = FALSE returns the scanned fit at the same lambda", {
  S <- ar1(8, 0.5)
  lam <- fixed_path(S)
  tight <- ebic_glasso(cor_matrix = S, n = 250, lambda_path = lam, refit = TRUE)
  scan  <- ebic_glasso(cor_matrix = S, n = 250, lambda_path = lam, refit = FALSE)
  expect_identical(tight$lambda, scan$lambda)
  # The refit is the strictly better answer -- that is why refit = FALSE is a
  # bit-compatibility switch and not a performance one.
  expect_lte(tight$kkt, scan$kkt)
  expect_lt(max(abs(tight$weights - scan$weights)), 1e-3)
})

test_that("refit = 'unregularized' keeps the support and switches certificate", {
  S <- ar1(7, 0.5)
  pen <- ebic_glasso(cor_matrix = S, n = 250)
  unr <- ebic_glasso(cor_matrix = S, n = 250, refit = "unregularized")
  expect_identical(unr$weights != 0, pen$weights != 0)
  expect_lt(unr$kkt, 1e-7)
  expect_equal(unr$kkt, ggm_support_kkt(unr$precision, S, unr$support))
})

test_that("penalize_diagonal is certified against the matching condition", {
  S <- ar1(8, 0.5)
  fit <- ebic_glasso(cor_matrix = S, n = 250, penalize_diagonal = TRUE)
  expect_lt(fit$kkt, 1e-7)
  # Grading the same fit with the wrong flag reports a spurious violation of
  # exactly rho -- which is why the flag exists on the certificate at all.
  expect_equal(glasso_kkt(fit$precision, S, fit$lambda,
                          penalize_diagonal = FALSE),
               fit$lambda, tolerance = 1e-3)
})

test_that("penalize_diagonal matches the glasso Fortran package", {
  skip_if_not_installed("glasso")
  S <- ar1(8, 0.5)
  ours <- .glasso_fit(S, rho = 0.15, penalize_diagonal = TRUE)$wi
  ref  <- glasso::glasso(S, 0.15, penalize.diagonal = TRUE, thr = 1e-10)$wi
  expect_lt(max(abs(ours - (ref + t(ref)) / 2)), 1e-6)
})

test_that("lambda_path and refit reject malformed input", {
  S <- ar1(5, 0.4)
  expect_error(ebic_glasso(cor_matrix = S, n = 200, lambda_path = c(-1, 0.2)),
               "strictly positive")
  expect_error(ebic_glasso(cor_matrix = S, n = 200, lambda_path = c(0.1, NA)),
               "strictly positive")
  expect_error(ebic_glasso(cor_matrix = S, n = 200, lambda_path = "x"),
               "numeric vector")
  expect_error(ebic_glasso(cor_matrix = S, n = 200, lambda_path = numeric(0)),
               "numeric vector")
  expect_error(ebic_glasso(cor_matrix = S, n = 200, refit = "nope"),
               "must be TRUE, FALSE")
})
