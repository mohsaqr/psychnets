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

test_that("native = FALSE is byte-identical to the glasso package", {
  skip_equiv("glasso")
  set.seed(10)
  p <- 8; n <- 500
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(ar1(p, 0.5))
  S <- stats::cor(X)

  fast <- ebic_glasso(cor_matrix = S, n = n, native = FALSE)
  g <- glasso::glasso(S, fast$lambda, penalize.diagonal = FALSE)
  expect_equal(fast$precision, (g$wi + t(g$wi)) / 2, ignore_attr = TRUE)
  expect_false(fast$native)

  # same edge set as the certified base engine; the certificate reflects the
  # solver: base reaches ~1e-9, the glasso engine sits at its own tolerance.
  base <- ebic_glasso(cor_matrix = S, n = n)
  expect_equal(off_compare(base$weights, fast$weights)$struct, 1)
  expect_lt(base$kkt, 1e-7)
  expect_gt(fast$kkt, base$kkt)
})

test_that("ising_fit(native=FALSE) byte-matches IsingFit::IsingFit", {
  skip_equiv("IsingFit")
  skip_if_not_installed("glmnet")
  set.seed(2)
  p <- 6; n <- 2000
  f <- matrix(stats::rnorm(n * p), n, p) %*% chol(ar1(p, 0.5))
  X <- (f > 0) * 1L
  colnames(X) <- paste0("V", seq_len(p))

  pn <- ising_fit(X, gamma = 0.25, rule = "AND", native = FALSE)
  isf <- suppressWarnings(IsingFit::IsingFit(X, gamma = 0.25, AND = TRUE,
                                             plot = FALSE, progressbar = FALSE))

  # The glmnet engine reproduces IsingFit's exact lambda path + EBIC selection,
  # so the symmetric weights and the raw-scale node thresholds byte-match.
  cmp <- off_compare(pn$graph, isf$weiadj)
  expect_equal(cmp$struct, 1)                 # identical edge set
  expect_lt(cmp$max_abs, 1e-7)                # byte-match (~0 to solver noise)
  expect_equal(unname(pn$thresholds), unname(isf$thresholds), tolerance = 1e-7)

  # the OR rule also byte-matches IsingFit's OR symmetrization
  pn_or <- ising_fit(X, gamma = 0.25, rule = "OR", native = FALSE)
  isf_or <- suppressWarnings(IsingFit::IsingFit(X, gamma = 0.25, AND = FALSE,
                                                plot = FALSE, progressbar = FALSE))
  expect_lt(off_compare(pn_or$graph, isf_or$weiadj)$max_abs, 1e-7)
})

test_that("mgm_fit(native=FALSE) byte-matches mgm::mgm magnitudes (mixed/binary)", {
  skip_equiv("mgm")
  skip_if_not_installed("glmnet")
  set.seed(3)
  n <- 2000
  lat <- matrix(stats::rnorm(n * 5), n, 5) %*% chol(ar1(5, 0.5))
  # two gaussian, three binary
  D <- data.frame(g1 = lat[, 1], g2 = lat[, 2],
                  b1 = (lat[, 3] > 0) * 1, b2 = (lat[, 4] > 0) * 1,
                  b3 = (lat[, 5] > 0) * 1)
  types <- c("g", "g", "c", "c", "c")
  levels <- c(1, 1, 2, 2, 2)

  pn <- mgm_fit(D, gamma = 0.25, native = FALSE)
  mg <- mgm::mgm(as.matrix(D), type = types, level = levels,
                 lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND",
                 pbar = FALSE, signInfo = FALSE)

  # edge magnitudes byte-match mgm's wadj, and the support is identical
  cmp <- off_compare(pn$graph, mg$pairwise$wadj)
  expect_lt(cmp$max_abs, 1e-6)
  ut <- upper.tri(pn$graph)
  expect_equal((abs(pn$graph) > 1e-9)[ut], (abs(mg$pairwise$wadj) > 1e-9)[ut])

  # gaussian-gaussian edge (positions 1,2) agrees on SIGN too (mgm's wadj
  # carries no dimnames, so index by position).
  gg_pn <- pn$graph[1, 2]
  gg_mg <- mg$pairwise$wadj[1, 2] * mg$pairwise$signs[1, 2]
  if (abs(gg_pn) > 1e-6 && !is.na(gg_mg) && abs(gg_mg) > 1e-6)
    expect_equal(sign(gg_pn), sign(gg_mg))
})

test_that("mgm(native=FALSE) byte-matches mgm on NON-unit-variance columns", {
  skip_equiv("mgm")
  skip_if_not_installed("glmnet")
  set.seed(7)
  n <- 2000
  lat <- matrix(stats::rnorm(n * 4), n, 4) %*% chol(ar1(4, 0.5))
  # gaussian columns deliberately off the unit-variance scale (SD 1.8, 0.6) plus
  # an arbitrary mean: mgm standardizes internally and so must mgm_fit().
  D <- data.frame(g1 = lat[, 1] * 1.8 + 3, g2 = lat[, 2] * 0.6 - 1,
                  b1 = (lat[, 3] > 0) * 1, b2 = (lat[, 4] > 0) * 1)
  pn <- mgm_fit(D, native = FALSE)
  mg <- mgm::mgm(as.matrix(D), type = c("g", "g", "c", "c"),
                 level = c(1, 1, 2, 2), lambdaSel = "EBIC", lambdaGam = 0.25,
                 ruleReg = "AND", pbar = FALSE, signInfo = FALSE, scale = TRUE)$pairwise$wadj
  # the gaussian-gaussian edge byte-matches to solver precision
  if (abs(pn$graph[1, 2]) > 1e-6 && mg[1, 2] > 1e-6)
    expect_equal(abs(pn$graph[1, 2]), mg[1, 2], tolerance = 1e-6)
})

test_that("mgm(native=FALSE) byte-matches mgm gaussian-binary edges", {
  skip_equiv("mgm")
  skip_if_not_installed("glmnet")
  set.seed(11)
  n <- 3000
  f <- stats::rnorm(n)
  D <- data.frame(
    g1 = f * 1.4 + stats::rnorm(n) * 0.5 + 2,
    b1 = (f + stats::rnorm(n) * 0.4 > 0) * 1,
    g2 = stats::rnorm(n),
    b2 = (stats::rnorm(n) > 0) * 1
  )

  pn <- mgm_fit(D, native = FALSE)
  mg <- mgm::mgm(as.matrix(D), type = c("g", "c", "g", "c"),
                 level = c(1, 2, 1, 2), lambdaSel = "EBIC",
                 lambdaGam = 0.25, ruleReg = "AND", pbar = FALSE,
                 signInfo = FALSE, scale = TRUE)$pairwise$wadj

  expect_equal(abs(pn$graph[1, 2]), mg[1, 2], tolerance = 1e-6)
})

test_that("moderated mgm_fit + condition() matches mgm::mgm(moderators) + condition()", {
  skip_equiv("mgm")
  skip_if_not_installed("glmnet")
  gen <- function(seed, n, p_g, p_c = 0L) {
    set.seed(seed)
    Z <- matrix(stats::rnorm(n * p_g), n, p_g)
    Z <- Z + 0.3 * matrix(rowSums(Z), n, p_g)
    G <- sample(0:1, n, replace = TRUE)
    Z[, p_g] <- Z[, p_g] + ifelse(G == 1, 0.6 * Z[, 1], -0.1 * Z[, 1])
    dat <- as.data.frame(Z)
    if (p_c > 0L)
      for (k in seq_len(p_c))
        dat[[paste0("C", k)]] <- sample(seq_len(sample(2:4, 1)), n, replace = TRUE)
    dat$G <- G
    type <- c(rep("g", p_g), rep("c", p_c), "c")
    level <- c(rep(1L, p_g),
               if (p_c > 0L)
                 vapply(dat[, (p_g + 1):(p_g + p_c), drop = FALSE],
                        function(x) length(unique(x)), integer(1)) else integer(0),
               2L)
    colnames(dat) <- paste0("V", seq_len(ncol(dat)))
    list(dat = dat, type = type, level = level, mod = ncol(dat))
  }
  cfgs <- list(c(200, 4, 0), c(250, 3, 1), c(300, 3, 2))
  for (cf in cfgs) {
    d <- gen(cf[1] * 7L, cf[1], cf[2], cf[3])
    # This gate checks the optional glmnet-backed reference engine for exact
    # parity with mgm.  The default base engine is independently optimized and
    # has its own structure/magnitude/KKT comparison in test-mgm.R; it also
    # intentionally rejects categorical nodes with more than two levels.
    fit <- mgm_fit(d$dat, types = d$type, moderators = d$mod, gamma = 0.25,
                   rule = "AND", threshold = "LW", native = FALSE)
    ref <- suppressWarnings(suppressMessages(mgm::mgm(
      as.matrix(d$dat), type = d$type, level = d$level, moderators = d$mod,
      lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND", threshold = "LW",
      overparameterize = FALSE, scale = TRUE, pbar = FALSE, signInfo = FALSE,
      warnings = FALSE)))
    for (v in c(0, 1)) {
      ours <- unname(condition(fit, v)$weights)
      rr <- unname(mgm::condition(
        ref, values = stats::setNames(list(v), as.character(d$mod)))$pairwise$wadj)
      expect_lt(max(abs(ours - rr)), 1e-10)
    }
  }
})

# ---- batch: graph metrics ported from networktools / qgraph / cocor ----------

test_that("net_bridge matches networktools::bridge (undirected)", {
  skip_if(Sys.getenv("PSYCHNET_EQUIV_TESTS") == "")
  skip_if_not_installed("networktools")
  set.seed(7); p <- 8
  S <- cov2cor(crossprod(matrix(stats::rnorm(p * p), p)) + diag(p))
  fit <- ebic_glasso(cor_matrix = S, n = 500); W <- fit$weights
  comm <- c(1, 1, 1, 2, 2, 2, 3, 3)
  ob <- net_bridge(fit, communities = comm)
  rb <- networktools::bridge(W, communities = comm)
  expect_lt(max(abs(ob$bridge_strength    - rb[["Bridge Strength"]])), 1e-8)
  expect_lt(max(abs(ob$bridge_betweenness - rb[["Bridge Betweenness"]])), 1e-8)
  expect_lt(max(abs(ob$bridge_closeness   - rb[["Bridge Closeness"]])), 1e-8)
  expect_lt(max(abs(ob$bridge_ei1 - rb[["Bridge Expected Influence (1-step)"]])), 1e-8)
  expect_lt(max(abs(ob$bridge_ei2 - rb[["Bridge Expected Influence (2-step)"]])), 1e-8)
})

test_that("2-step expected influence matches networktools::expectedInf", {
  skip_if(Sys.getenv("PSYCHNET_EQUIV_TESTS") == "")
  skip_if_not_installed("networktools")
  set.seed(7); p <- 8
  S <- cov2cor(crossprod(matrix(stats::rnorm(p * p), p)) + diag(p))
  fit <- ebic_glasso(cor_matrix = S, n = 500)
  ours <- net_centralities(fit, measures = "expected_influence_2step")$expected_influence_2step
  ref  <- as.numeric(networktools::expectedInf(fit$weights, step = 2)$step2)
  expect_lt(max(abs(ours - ref)), 1e-8)
})

test_that("net_clustering matches qgraph::clustcoef_auto", {
  skip_if(Sys.getenv("PSYCHNET_EQUIV_TESTS") == "")
  skip_if_not_installed("qgraph")
  set.seed(7); p <- 8
  S <- cov2cor(crossprod(matrix(stats::rnorm(p * p), p)) + diag(p))
  fit <- ebic_glasso(cor_matrix = S, n = 500)
  oc <- net_clustering(fit); rc <- qgraph::clustcoef_auto(fit$weights)
  for (col in intersect(names(oc), names(rc)))
    expect_lt(max(abs(oc[[col]] - rc[[col]]), na.rm = TRUE), 1e-8)
})

test_that("Hittner2003 test matches cocor", {
  skip_if(Sys.getenv("PSYCHNET_EQUIV_TESTS") == "")
  skip_if_not_installed("cocor")
  set.seed(3)
  for (. in 1:5) {
    r <- stats::runif(3, -0.7, 0.7); n <- sample(80:400, 1)
    ours <- .psn_hittner2003(r[1], r[2], r[3], n)
    ref  <- suppressWarnings(cocor::cocor.dep.groups.overlap(
      r[1], r[2], r[3], n, test = "hittner2003")@hittner2003$p.value)
    expect_lt(abs(ours - ref), 1e-10)
  }
})

test_that("redundancy matches networktools::goldbricker", {
  skip_if(Sys.getenv("PSYCHNET_EQUIV_TESTS") == "")
  skip_if_not_installed("networktools")
  og <- redundancy(SRL_Claude, cor_method = "auto")
  rg <- tryCatch(
    suppressWarnings(networktools::goldbricker(SRL_Claude, progressbar = FALSE)),
    error = function(e) {
      skip(paste("networktools::goldbricker unavailable with the installed lavaan:",
                 conditionMessage(e)))
    })
  expect_lt(max(abs(attr(og, "proportion_matrix") - rg$proportion_matrix),
                na.rm = TRUE), 1e-10)
})

test_that("net_aggregate PCA composite matches networktools::net_reduce", {
  skip_if(Sys.getenv("PSYCHNET_EQUIV_TESTS") == "")
  skip_if_not_installed("networktools")
  set.seed(1)
  d <- as.data.frame(matrix(stats::rnorm(200 * 4), 200, 4)); names(d) <- paste0("V", 1:4)
  ours <- net_aggregate(d, c(1, 1, 2, 2), method = "pca")[["1"]]
  red  <- networktools::net_reduce(d, badpairs = list(c("V1", "V2")), method = "PCA")
  expect_equal(abs(stats::cor(ours, red[["PCA.V1.V2"]])), 1, tolerance = 1e-8)
})

test_that("net_edge_betweenness matches igraph on an unweighted graph", {
  skip_if(Sys.getenv("PSYCHNET_EQUIV_TESTS") == "")
  skip_if_not_installed("igraph")
  set.seed(2); A <- matrix(0, 7, 7); on <- sample(which(upper.tri(A)), 9)
  A[on] <- 1; A <- A + t(A); colnames(A) <- rownames(A) <- paste0("n", 1:7)
  ours <- net_edge_betweenness(A)
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  el <- igraph::as_edgelist(g)
  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  refm  <- stats::setNames(igraph::edge_betweenness(g), key(el[, 1], el[, 2]))
  oursm <- stats::setNames(ours$edge_betweenness, key(ours$from, ours$to))
  common <- intersect(names(refm), names(oursm))
  expect_lt(max(abs(oursm[common] - refm[common])), 1e-8)
})
