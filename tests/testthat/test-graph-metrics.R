# Self-contained tests for the qgraph/networktools graph-metric ports
# (bridge, 2-step EI, clustering, smallworldness, redundancy). Reference-package
# equivalence lives in the gated test-equivalence.R.

mk_fit <- function(p = 6, seed = 1) {
  set.seed(seed)
  S <- cov2cor(crossprod(matrix(stats::rnorm(p * p), p)) + diag(p))
  ebic_glasso(cor_matrix = S, n = 400)
}

test_that("2-step expected influence equals ei1 + W %*% ei1", {
  fit <- mk_fit()
  W <- fit$weights; diag(W) <- 0
  ei1 <- rowSums(W)
  manual <- as.numeric(ei1 + W %*% ei1)
  ours <- net_centralities(fit, measures = "expected_influence_2step")$expected_influence_2step
  expect_equal(ours, manual)
})

test_that("net_bridge returns a tidy per-node frame; strength is cross-community", {
  fit <- mk_fit(6)
  comm <- c(1, 1, 1, 2, 2, 2)
  b <- net_bridge(fit, communities = comm)
  expect_s3_class(b, "psychnet_bridge")
  expect_s3_class(b, "data.frame")
  expect_named(b, c("node", "community", "bridge_strength", "bridge_betweenness",
                    "bridge_closeness", "bridge_ei1", "bridge_ei2"))
  expect_equal(nrow(b), 6L)
  W <- fit$weights; diag(W) <- 0
  cross <- outer(comm, comm, "!=")
  expect_equal(b$bridge_strength, unname(rowSums(abs(W) * cross)))
  expect_equal(b$bridge_ei1, unname(rowSums(W * cross)))
})

test_that("net_bridge accepts named communities and errors on <2 groups", {
  fit <- mk_fit(6)
  named <- stats::setNames(c(1, 1, 1, 2, 2, 2), fit$labels)
  expect_equal(net_bridge(fit, named)$bridge_strength,
               net_bridge(fit, c(1, 1, 1, 2, 2, 2))$bridge_strength)
  expect_error(net_bridge(fit, rep(1, 6)), "at least two communities")
})

test_that("net_clustering is tidy and correct on a known unweighted triangle", {
  fit <- mk_fit()
  cc <- net_clustering(fit)
  expect_s3_class(cc, "psychnet_clustering")
  expect_true(all(c("clustWS", "clustZhang", "clustOnnela", "clustBarrat") %in% names(cc)))
  tri <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
  colnames(tri) <- c("a", "b", "c")
  ct <- net_clustering(tri)
  expect_equal(ct$clustWS, c(1, 1, 1))
  expect_equal(ct$clustBarrat, c(1, 1, 1))
})

test_that("net_smallworld returns a sane one-row summary", {
  fit <- mk_fit(8)
  sw <- net_smallworld(fit, n_rand = 30, seed = 1)
  expect_s3_class(sw, "psychnet_smallworld")
  expect_equal(nrow(sw), 1L)
  expect_named(sw, c("smallworldness", "transitivity", "aspl",
                     "transitivity_rand", "aspl_rand"))
  expect_true(sw$transitivity >= 0 && sw$transitivity <= 1)
  expect_true(is.finite(sw$smallworldness) && sw$smallworldness > 0)
})

test_that("Hittner2003 test is symmetric, in [0,1], and 1 when correlations equal", {
  expect_equal(.psn_hittner2003(0.3, 0.3, 0.5, 200), 1)
  p1 <- .psn_hittner2003(0.6, 0.2, 0.4, 150)
  p2 <- .psn_hittner2003(0.2, 0.6, 0.4, 150)
  expect_equal(p1, p2)
  expect_true(p1 >= 0 && p1 <= 1)
  expect_true(is.na(.psn_hittner2003(1, 0.3, 0.4, 100)))
})

test_that("redundancy returns a tidy pairs frame with a symmetric proportion matrix", {
  r <- redundancy(SRL_Claude, cor_method = "pearson")
  expect_s3_class(r, "psychnet_redundancy")
  expect_named(r, c("item1", "item2", "proportion", "correlation"))
  pm <- attr(r, "proportion_matrix")
  expect_equal(pm, t(pm))
  expect_true(all(is.na(diag(pm))))
})

test_that("net_edge_betweenness returns a tidy per-edge frame", {
  fit <- mk_fit(6)
  eb <- net_edge_betweenness(fit)
  expect_s3_class(eb, "data.frame")
  expect_named(eb, c("from", "to", "edge_betweenness"))
  W <- fit$weights; diag(W) <- 0
  expect_equal(nrow(eb), sum(W != 0 & upper.tri(W)))   # one row per undirected edge
  expect_true(all(eb$edge_betweenness >= 0))
})

test_that("net_edge_betweenness: middle edge of a path is the most central", {
  # path a-b-c-d-e: edge b-c and c-d carry the most shortest paths
  A <- matrix(0, 5, 5); for (i in 1:4) A[i, i + 1] <- A[i + 1, i] <- 1
  colnames(A) <- letters[1:5]
  eb <- net_edge_betweenness(A)
  mid <- eb$edge_betweenness[eb$from == "b" & eb$to == "c"]
  end <- eb$edge_betweenness[eb$from == "a" & eb$to == "b"]
  expect_gt(mid, end)
})
