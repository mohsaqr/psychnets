# Gated external-equivalence suite.
#
# These tests cross-check psychnet against the established reference packages
# (qgraph, IsingFit, mgm) at independent-solver precision. They are the
# complement to the self-certification tests: certification proves each fit
# solves its own objective; this proves the objective is the one the field's
# reference computes.
#
# The suite is GATED: it runs only when PSYCHNET_EQUIV_TESTS is set to a
# non-empty value AND the relevant reference package is installed. By default
# (and on CRAN) every test skips, so `R CMD check` stays self-contained and
# dependency-free. Run it with:
#
#   PSYCHNET_EQUIV_TESTS=true Rscript -e 'testthat::test_file(
#     "tests/testthat/test-equivalence.R")'
#
# Equivalence convention. Independent convex solvers do not byte-match (they
# differ in stopping rule and in EBIC tie-breaking between adjacent path
# points), so we assert agreement of the zero/nonzero edge pattern (`struct`)
# and a small maximum absolute weight difference, not equality. For Ising and
# mgm the comparison is on edge MAGNITUDE: mgm and IsingFit leave a
# binary-binary edge unsigned (mgm stores `NA` in its sign matrix), whereas
# psychnet keeps a signed coefficient. The magnitudes are the shared quantity;
# psychnet's sign is an addition, not a disagreement.

skip_equiv <- function(pkg) {
  if (!nzchar(Sys.getenv("PSYCHNET_EQUIV_TESTS")))
    skip("PSYCHNET_EQUIV_TESTS not set; gated external-equivalence suite skipped.")
  skip_if_not_installed(pkg)
}

# Off-diagonal agreement of two weighted adjacency matrices, on magnitude.
off_compare <- function(A, E, eps = 1e-6) {
  A <- as.matrix(A); E <- as.matrix(E)
  ut <- upper.tri(A)
  a <- abs(A[ut]); e <- abs(E[ut])
  list(struct = mean((a > eps) == (e > eps)),
       max_abs = max(abs(a - e)),
       n_a = sum(a > eps), n_e = sum(e > eps))
}

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

test_that("EBICglasso matches qgraph::EBICglasso on a clear chain", {
  skip_equiv("qgraph")
  set.seed(1)
  p <- 8; n <- 1000
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(ar1(p, 0.5))
  S <- stats::cor(X)

  pn <- ebic_glasso(cor_matrix = S, n = n, gamma = 0.5)
  qg <- suppressWarnings(
    qgraph::EBICglasso(S, n = n, gamma = 0.5, returnAllResults = FALSE))

  cmp <- off_compare(pn$graph, qg)
  expect_equal(cmp$struct, 1)                 # identical edge set
  expect_lt(cmp$max_abs, 0.02)                # weights agree to qgraph's tolerance
  # and psychnet is at the certified optimum
  expect_lt(pn$kkt, 1e-7)
})

test_that("engine = 'glasso' is byte-identical to the glasso package", {
  skip_equiv("glasso")
  set.seed(10)
  p <- 8; n <- 500
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(ar1(p, 0.5))
  S <- stats::cor(X)

  fast <- ebic_glasso(cor_matrix = S, n = n, engine = "glasso")
  g <- glasso::glasso(S, fast$lambda, penalize.diagonal = FALSE)
  expect_equal(fast$precision, (g$wi + t(g$wi)) / 2, ignore_attr = TRUE)
  expect_identical(fast$engine, "glasso")

  # same edge set as the certified base engine; the certificate reflects the
  # solver: base reaches ~1e-9, the glasso engine sits at its own tolerance.
  base <- ebic_glasso(cor_matrix = S, n = n)
  expect_equal(off_compare(base$weights, fast$weights)$struct, 1)
  expect_lt(base$kkt, 1e-7)
  expect_gt(fast$kkt, base$kkt)
})

test_that("ising_fit matches IsingFit::IsingFit on a binary chain", {
  skip_equiv("IsingFit")
  set.seed(2)
  p <- 6; n <- 2000
  f <- matrix(stats::rnorm(n * p), n, p) %*% chol(ar1(p, 0.5))
  X <- (f > 0) * 1L
  colnames(X) <- paste0("V", seq_len(p))

  pn <- ising_fit(X, gamma = 0.25, rule = "AND")
  isf <- suppressWarnings(IsingFit::IsingFit(X, gamma = 0.25, AND = TRUE,
                                             plot = FALSE, progressbar = FALSE))

  cmp <- off_compare(pn$graph, isf$weiadj)
  expect_gte(cmp$struct, 0.9)                 # at most one borderline edge differs
  expect_lt(cmp$max_abs, 0.25)                # logit-scale weights, generous tol

  # raw-scale node thresholds agree with IsingFit's (the scale-fix regression):
  # both are logit intercepts on the 0/1 scale, so they share sign and order.
  expect_equal(sign(pn$thresholds), sign(isf$thresholds))
  expect_gt(stats::cor(pn$thresholds, isf$thresholds), 0.95)
})

test_that("mgm_fit matches mgm::mgm magnitudes on mixed and binary data", {
  skip_equiv("mgm")
  set.seed(3)
  n <- 2000
  lat <- matrix(stats::rnorm(n * 5), n, 5) %*% chol(ar1(5, 0.5))
  # two gaussian, three binary
  D <- data.frame(g1 = lat[, 1], g2 = lat[, 2],
                  b1 = (lat[, 3] > 0) * 1, b2 = (lat[, 4] > 0) * 1,
                  b3 = (lat[, 5] > 0) * 1)
  types <- c("g", "g", "c", "c", "c")
  levels <- c(1, 1, 2, 2, 2)

  pn <- mgm_fit(D, gamma = 0.25)
  mg <- mgm::mgm(as.matrix(D), type = types, level = levels,
                 lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND",
                 pbar = FALSE, signInfo = FALSE)

  cmp <- off_compare(pn$graph, mg$pairwise$wadj)
  expect_lt(cmp$max_abs, 0.1)                 # magnitudes agree closely
  # Clear-signal edges (|w| > 0.1) agree on presence exactly; weak ~0.05 edges
  # can flip at the LW/lambda boundary, which psychnet selects on an independent
  # base-R path (documented; not a byte-match on borderline edges).
  s_pn <- abs(pn$graph) > 0.1; s_mg <- abs(mg$pairwise$wadj) > 0.1
  ut <- upper.tri(s_pn)
  expect_equal(mean(s_pn[ut] == s_mg[ut]), 1)

  # gaussian-gaussian edge (positions 1,2): psychnet and mgm agree on SIGN too,
  # because a gaussian-gaussian edge has a defined sign (mgm's wadj carries no
  # dimnames, so index by position).
  gg_pn <- pn$graph[1, 2]
  gg_mg <- mg$pairwise$wadj[1, 2] * mg$pairwise$signs[1, 2]
  if (abs(gg_pn) > 1e-6 && !is.na(gg_mg) && abs(gg_mg) > 1e-6)
    expect_equal(sign(gg_pn), sign(gg_mg))
})

test_that("mgm gaussian edges match mgm on NON-unit-variance columns", {
  skip_equiv("mgm")
  set.seed(7)
  n <- 2000
  lat <- matrix(stats::rnorm(n * 4), n, 4) %*% chol(ar1(4, 0.5))
  # gaussian columns deliberately off the unit-variance scale (SD 1.8, 0.6) plus
  # an arbitrary mean: mgm standardizes internally and so must mgm_fit().
  D <- data.frame(g1 = lat[, 1] * 1.8 + 3, g2 = lat[, 2] * 0.6 - 1,
                  b1 = (lat[, 3] > 0) * 1, b2 = (lat[, 4] > 0) * 1)
  pn <- mgm_fit(D)
  mg <- mgm::mgm(as.matrix(D), type = c("g", "g", "c", "c"),
                 level = c(1, 1, 2, 2), lambdaSel = "EBIC", lambdaGam = 0.25,
                 ruleReg = "AND", pbar = FALSE, signInfo = FALSE, scale = TRUE)$pairwise$wadj
  # the gaussian-gaussian edge now matches to solver precision (the scale fix);
  # without it psychnet was ~28% too large on this data.
  if (abs(pn$graph[1, 2]) > 1e-6 && mg[1, 2] > 1e-6)
    expect_equal(abs(pn$graph[1, 2]), mg[1, 2], tolerance = 0.02)
})

test_that("mgm gaussian-binary edges use mgm's categorical scale", {
  skip_equiv("mgm")
  set.seed(11)
  n <- 3000
  f <- stats::rnorm(n)
  D <- data.frame(
    g1 = f * 1.4 + stats::rnorm(n) * 0.5 + 2,
    b1 = (f + stats::rnorm(n) * 0.4 > 0) * 1,
    g2 = stats::rnorm(n),
    b2 = (stats::rnorm(n) > 0) * 1
  )

  pn <- mgm_fit(D)
  mg <- mgm::mgm(as.matrix(D), type = c("g", "c", "g", "c"),
                 level = c(1, 2, 1, 2), lambdaSel = "EBIC",
                 lambdaGam = 0.25, ruleReg = "AND", pbar = FALSE,
                 signInfo = FALSE, scale = TRUE)$pairwise$wadj

  expect_equal(abs(pn$graph[1, 2]), mg[1, 2], tolerance = 0.01)
  expect_lt(pn$kkt, 1e-6)
})
