# Fixes to the literature (2026-11-22): a planarity test that holds for any graph in tmfg_certificate(), and
# glmnet's default penalty path for the native nodewise solvers of ising_fit() and mgm_fit().

graph_of <- \(p, edges) {
  adjacency <- matrix(FALSE, p, p)
  adjacency[edges] <- TRUE
  adjacency[edges[, 2:1, drop = FALSE]] <- TRUE
  adjacency
}
complete_graph <- \(p) !diag(p) > 0
k33 <- graph_of(6L, as.matrix(expand.grid(1:3, 4:6)))
petersen <- graph_of(10L, rbind(cbind(1:5, c(2:5, 1L)), cbind(1:5, 6:10), cbind(6:10, c(8L, 9L, 10L, 6L, 7L))))
cube <- graph_of(8L, rbind(cbind(1:4, c(2:4, 1L)), cbind(5:8, c(6:8, 5L)), cbind(1:4, 5:8)))
octahedron <- complete_graph(6L) & !graph_of(6L, rbind(c(1L, 2L), c(3L, 4L), c(5L, 6L)))
drop_edge <- \(adjacency, a, b) { adjacency[a, b] <- adjacency[b, a] <- FALSE; adjacency }
disjoint <- \(a, b) {
  out <- matrix(FALSE, nrow(a) + nrow(b), nrow(a) + nrow(b))
  out[seq_len(nrow(a)), seq_len(nrow(a))] <- a
  out[nrow(a) + seq_len(nrow(b)), nrow(a) + seq_len(nrow(b))] <- b
  out
}
as_network_of <- \(adjacency) {
  p <- ncol(adjacency)
  fit <- tmfg_network(cor_matrix = 0.5^abs(outer(seq_len(p), seq_len(p), "-")))
  dimnames(adjacency) <- dimnames(fit$adjacency)
  fit$adjacency <- adjacency
  fit
}

test_that("the planarity test decides Kuratowski's graphs and planar solids correctly", {
  expect_false(.is_planar(complete_graph(5L)))
  expect_false(.is_planar(k33))
  expect_false(.is_planar(petersen))
  expect_false(.is_planar(disjoint(complete_graph(4L), k33)))   # non-planar in the second component only
  expect_true(.is_planar(drop_edge(complete_graph(5L), 1L, 2L)))
  expect_true(.is_planar(drop_edge(k33, 1L, 4L)))
  expect_true(.is_planar(cube))
  expect_true(.is_planar(octahedron))
  expect_true(.is_planar(disjoint(complete_graph(4L), complete_graph(4L))))
})

test_that("tmfg_certificate() counts K3,3 as violating edge count, chordality and planarity", {
  # K3,3: 9 edges, not 3(6 - 2) = 12; the 4-cycle 1-4-2-5 has no chord; it is a Kuratowski graph. The old
  # stacked-triangulation check, valid only for chordal graphs, counted it planar and returned 2.
  expect_identical(tmfg_certificate(as_network_of(k33)), 3L)
  expect_identical(tmfg_certificate(as_network_of(complete_graph(5L))), 2L)   # 10 edges, not 9; non-planar
})

test_that("a TMFG is planar, and any added edge makes it non-planar", {
  set.seed(7)
  fit <- tmfg_network(cor_matrix = stats::cor(matrix(stats::rnorm(80 * 12), 80, 12)))
  adjacency <- fit$adjacency
  expect_true(.is_planar(adjacency))
  expect_identical(tmfg_certificate(fit), 0L)
  gaps <- which(!adjacency & upper.tri(adjacency), arr.ind = TRUE)
  # A TMFG is a maximal planar graph (3(p - 2) edges), so every added edge breaks planarity.
  added <- vapply(seq_len(nrow(gaps)), \(k) {
    a <- adjacency
    a[gaps[k, 1L], gaps[k, 2L]] <- a[gaps[k, 2L], gaps[k, 1L]] <- TRUE
    .is_planar(a)
  }, logical(1))
  expect_false(any(added))
})

test_that("planarity does not depend on the order of the nodes", {
  set.seed(11)
  graphs <- lapply(seq_len(40), \(k) {
    p <- sample(6:14, 1L)
    a <- matrix(stats::runif(p * p) < 0.35, p, p)
    a[lower.tri(a, diag = TRUE)] <- FALSE
    a | t(a)
  })
  same <- vapply(graphs, \(a) {
    perm <- sample(nrow(a))
    identical(.is_planar(a), .is_planar(a[perm, perm]))
  }, logical(1))
  expect_true(all(same))
})

test_that("glmnet's default lambda.min.ratio and path termination are reproduced", {
  # glmnet(): lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-4)
  expect_identical(.glmnet_default_ratio(10L, 20L), 0.01)
  expect_identical(.glmnet_default_ratio(20L, 20L), 1e-4)
  expect_identical(.glmnet_default_ratio(2632L, 5L), 1e-4)
  # glmnet.control(): no stop before the 5th penalty; then a binomial path stops when the explained deviance
  # grows by less than 1e-5, a gaussian path when it grows by less than 1e-5 of itself, either above 0.999.
  expect_false(.glmnet_path_stops(4L, 100L, "binomial", 0.5, 0.5))
  expect_true(.glmnet_path_stops(5L, 100L, "binomial", 0.500005, 0.5))
  expect_false(.glmnet_path_stops(5L, 100L, "binomial", 0.50002, 0.5))
  expect_true(.glmnet_path_stops(5L, 100L, "gaussian", 0.500004, 0.5))      # (0.500004 - 0.5) / 0.500004 < 1e-5
  expect_false(.glmnet_path_stops(5L, 100L, "gaussian", 0.500006, 0.5))     # 1.2e-5 > 1e-5
  expect_true(.glmnet_path_stops(6L, 100L, "binomial", 0.9995, 0.99))
  expect_true(.glmnet_path_stops(2L, 2L, "gaussian", 0.3, 0.3))             # mnlam is capped at nlambda
})

test_that("the native nodewise path is glmnet's default path", {
  skip_if_not_installed("glmnet")
  set.seed(5)
  X <- matrix(stats::rbinom(600 * 4, 1, 0.5), 600, 4)
  y <- stats::rbinom(600, 1, stats::plogis(-0.3 + X %*% c(0.8, -0.5, 0, 0.3)))
  fit <- glmnet::glmnet(X, y, family = "binomial")
  ours <- .nodewise_ebic(X, y, "binomial", gamma = 0.25, nlambda = 100L)
  # the selected penalty is one of glmnet's, chosen on the same grid
  expect_equal(min(abs(fit$lambda - ours$lambda)), 0, tolerance = 1e-10)
})

test_that("native ising_fit() and mgm_fit() agree with their glmnet engine on the default path", {
  skip_if_not_installed("glmnet")
  b <- dichotomize(SRL_GPT)
  native <- ising_fit(b)
  engine <- ising_fit(b, native = FALSE)
  expect_equal(native$weights, engine$weights, tolerance = 1e-4)
  expect_identical(ising_fit(b, lambda_min_ratio = 1e-4)$weights, native$weights)
  set.seed(3)
  latent <- matrix(stats::rnorm(2000L * 5L), 2000L, 5L) %*% chol(0.5^abs(outer(1:5, 1:5, "-")))
  d <- data.frame(g1 = latent[, 1L], g2 = latent[, 2L], b1 = (latent[, 3L] > 0) * 1,
                  b2 = (latent[, 4L] > 0) * 1, b3 = (latent[, 5L] > 0) * 1)
  m_native <- mgm_fit(d, gamma = 0.25)
  m_engine <- mgm_fit(d, gamma = 0.25, native = FALSE)
  expect_equal(abs(m_native$weights), abs(m_engine$weights), tolerance = 1e-4)
  # The old grid (down to 0.01 of lambda_max) missed the g2-b2 edge that mgm finds (0.0497).
  old <- mgm_fit(d, gamma = 0.25, lambda_min_ratio = 0.01)
  expect_identical(old$weights["g2", "b2"], 0)
  expect_gt(abs(m_native$weights["g2", "b2"]), 0.04)
})

test_that("lambda_min_ratio is still validated by class", {
  b <- dichotomize(SRL_GPT)
  expect_error(ising_fit(b, lambda_min_ratio = 0), class = "psychnet_bad_input")
  expect_error(mgm_fit(data.frame(a = stats::rnorm(50), b = stats::rnorm(50), c = stats::rnorm(50)),
                       lambda_min_ratio = 2), class = "psychnet_bad_input")
})
