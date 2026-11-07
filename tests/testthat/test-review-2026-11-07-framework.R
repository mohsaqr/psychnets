# Regression tests for the 2026-11-07 review of the framework verbs:
# centrality (W2, W3, W11), bridge (W1, W4, W7), aggregate / redundancy (W4,
# W5, W6), predictability (W10) and the base-graphics plots (W8, W9, N13).

rv_fit <- function(p = 6, seed = 1) {
  set.seed(seed)
  S <- cov2cor(crossprod(matrix(stats::rnorm(p * p), p)) + diag(p))
  ebic_glasso(cor_matrix = S, n = 400)
}

rv_data <- function(seed = 1, n = 200, p = 5) {
  set.seed(seed)
  d <- as.data.frame(matrix(stats::rnorm(n * p), n, p) %*%
                       chol(0.4^abs(outer(seq_len(p), seq_len(p), "-"))))
  names(d) <- paste0("V", seq_len(p))
  d
}

# Another agent makes .as_numeric_matrix() warn `psychnet_constant_variable`
# when it drops a constant column; that warning is expected here and is not
# what these tests are about, so exactly that class is muffled.
quiet_constant <- function(expr) {
  withCallingHandlers(expr, psychnet_constant_variable = function(w)
    invokeRestart("muffleWarning"))
}

# Draw to a throwaway PDF and return the graphics parameters before / after.
draw_par <- function(expr) {
  grDevices::pdf(file.path(tempdir(), "psn-review-plot.pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  keys <- c("mar", "oma", "mfrow", "mfcol", "fig", "cex", "xpd", "mgp", "las")
  before <- graphics::par(keys)
  force(expr)
  after <- graphics::par(keys)
  list(before = before, after = after)
}

# ---- W1: bridge normalize matches networktools' divisors ---------------------

test_that("W1: bridge normalize uses the node count p in the 2-step divisor", {
  fit <- rv_fit(8)
  comm <- c(1, 1, 1, 2, 2, 2, 3, 3)
  raw  <- net_bridge(fit, comm)
  nrm  <- net_bridge(fit, comm, normalize = TRUE)
  p <- 8
  q <- vapply(comm, function(cc) sum(comm != cc), numeric(1))
  expect_equal(nrm$bridge_strength, raw$bridge_strength / q)
  expect_equal(nrm$bridge_ei1, raw$bridge_ei1 / q)
  expect_equal(nrm$bridge_betweenness, raw$bridge_betweenness / ((p - 1) * q))
  expect_equal(nrm$bridge_ei2, raw$bridge_ei2 / (q + (p - 1) * (q - q / (p - 1))))
  expect_equal(nrm$bridge_closeness, raw$bridge_closeness)
  expect_error(net_bridge(fit, comm, normalize = NA), class = "psychnet_bad_input")
})

# ---- W2: bare asymmetric matrices are directed --------------------------------

test_that("W2: net_edge_betweenness treats an asymmetric bare matrix as directed", {
  rl <- relimp_network(SRL_Claude)
  from_obj <- net_edge_betweenness(rl)
  from_mat <- net_edge_betweenness(rl$weights)
  expect_identical(nrow(from_mat), nrow(from_obj))
  expect_equal(from_mat, from_obj)
  # a symmetric bare matrix stays undirected: one row per unordered edge
  W <- rv_fit()$weights
  eb <- net_edge_betweenness(W)
  expect_identical(nrow(eb), sum(abs(W) > 0 & upper.tri(W)))
  expect_error(net_edge_betweenness(W, labels = letters[1:3]),
               class = "psychnet_bad_input")
})

# ---- W3: centrality_fn names are honoured -------------------------------------

test_that("W3: a named centrality_fn vector is aligned by node label", {
  fit <- rv_fit()
  labs <- fit$nodes$label
  fn_rev <- function(W) {
    v <- rowSums(abs(W))
    list(mine = rev(v))                     # named, reversed order
  }
  fn_plain <- function(W) list(mine = unname(rowSums(abs(W))))
  a <- net_centralities(fit, measures = c("strength", "mine"),
                        centrality_fn = fn_rev)
  b <- net_centralities(fit, measures = c("strength", "mine"),
                        centrality_fn = fn_plain)
  expect_equal(a$mine, a$strength)          # aligned, not reversed
  expect_equal(a, b)
  # invariance: any permutation of a named vector gives the same table
  perm <- function(W) { v <- rowSums(abs(W)); list(mine = v[c(3, 1, 6, 2, 5, 4)]) }
  expect_equal(net_centralities(fit, measures = "mine", centrality_fn = perm)$mine,
               a$strength)
})

test_that("W3: bad centrality_fn results raise psychnet_bad_input", {
  fit <- rv_fit()
  short <- function(W) list(mine = 1:3)
  wrong <- function(W) list(mine = stats::setNames(rowSums(W), paste0("x", seq_len(nrow(W)))))
  chr   <- function(W) list(mine = rep("a", nrow(W)))
  expect_error(net_centralities(fit, measures = "mine", centrality_fn = short),
               class = "psychnet_bad_input")
  expect_error(net_centralities(fit, measures = "mine", centrality_fn = wrong),
               class = "psychnet_bad_input")
  expect_error(net_centralities(fit, measures = "mine", centrality_fn = chr),
               class = "psychnet_bad_input")
  expect_error(net_centralities(fit, measures = "mine", centrality_fn = "x"),
               class = "psychnet_bad_input")
})

# ---- W11: documented scaling --------------------------------------------------

test_that("W11: edge betweenness is a path count; node betweenness is normalized", {
  A <- matrix(0, 5, 5)
  A[cbind(1:4, 2:5)] <- 1; A <- A + t(A)
  colnames(A) <- letters[1:5]
  expect_equal(net_edge_betweenness(A)$edge_betweenness, c(4, 6, 6, 4))
  # star: the hub lies on every leaf-leaf geodesic -> normalized betweenness 1
  star <- matrix(0, 5, 5); star[1, -1] <- star[-1, 1] <- 1
  colnames(star) <- letters[1:5]
  cen <- net_centralities(star, measures = c("betweenness", "closeness"))
  expect_equal(cen$betweenness, c(1, 0, 0, 0, 0))
  expect_equal(cen$closeness[1], 1)         # 4 reachable / distance sum 4
})

test_that("W11: betweenness/closeness equal igraph's normalized versions", {
  skip_if_not_installed("igraph")
  W <- rv_fit(7, seed = 7)$weights
  cen <- net_centralities(W, measures = c("betweenness", "closeness"))
  G <- igraph::graph_from_adjacency_matrix(abs(W), mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  igraph::E(G)$weight <- 1 / igraph::E(G)$weight
  expect_equal(cen$betweenness, unname(igraph::betweenness(G, normalized = TRUE)))
  expect_equal(cen$closeness, unname(igraph::closeness(G, normalized = TRUE)))
})

# ---- W4: community specification ----------------------------------------------

test_that("W4: bridge communities as a list of label vectors equal the vector form", {
  fit <- rv_fit()
  labs <- fit$nodes$label
  lst <- list(a = labs[1:3], b = labs[4:6])
  expect_equal(net_bridge(fit, lst)$bridge_strength,
               net_bridge(fit, c(1, 1, 1, 2, 2, 2))$bridge_strength)
  expect_identical(net_bridge(fit, lst)$community, rep(c("a", "b"), each = 3))
})

test_that("W4: duplicated, unknown and missing community labels are errors", {
  fit <- rv_fit()
  labs <- fit$nodes$label
  expect_error(net_bridge(fit, list(a = labs[1:4], b = labs[4:6])),
               "more than one community", class = "psychnet_bad_input")
  expect_error(net_bridge(fit, list(a = c(labs[1:3], "zz"), b = labs[4:6])),
               "not in the network", class = "psychnet_bad_input")
  expect_error(net_bridge(fit, list(a = labs[1:2], b = labs[4:6])),
               "no community", class = "psychnet_bad_input")
  # the old documented "list keyed by node label" form now says why it fails
  keyed <- stats::setNames(as.list(c(1, 1, 1, 2, 2, 2)), labs)
  expect_error(net_bridge(fit, keyed), class = "psychnet_bad_input")
  named <- stats::setNames(c(1, 1, 1, 2, 2, 2, 3), c(labs, "zz"))
  expect_error(net_bridge(fit, named), "not in the network",
               class = "psychnet_bad_input")
  expect_error(net_bridge(fit, c(1, 2)), "one entry per node",
               class = "psychnet_bad_input")
})

test_that("W4: bridge centrality is permutation-equivariant", {
  fit <- rv_fit(7, seed = 3)
  W <- fit$weights
  comm <- c(1, 1, 2, 2, 2, 3, 3)
  ref <- net_bridge(W, comm, normalize = TRUE)
  o <- c(5, 2, 7, 1, 4, 6, 3)
  prm <- net_bridge(W[o, o], comm[o], normalize = TRUE)
  ref_o <- ref[o, , drop = FALSE]; rownames(ref_o) <- NULL
  expect_equal(as.data.frame(prm), as.data.frame(ref_o))
})

test_that("W4: aggregate accepts the list form and rejects bad labels", {
  d <- rv_data()
  vec <- net_aggregate(d, c(1, 1, 2, 2, 2))
  lst <- net_aggregate(d, list(`1` = c("V1", "V2"), `2` = c("V3", "V4", "V5")))
  expect_equal(lst, vec)
  expect_error(net_aggregate(d, list(a = c("V1", "V2"), b = c("V2", "V3", "V4", "V5"))),
               class = "psychnet_bad_input")
  expect_error(net_aggregate(d, list(a = c("V1", "V9"), b = c("V3", "V4", "V5"))),
               class = "psychnet_bad_input")
})

# ---- W5: dropped rows / columns ------------------------------------------------

test_that("W5: a constant column is matched before it is dropped", {
  d <- rv_data()
  d6 <- d; d6$K <- 1
  # the bare 5-vector no longer silently shifts onto the wrong columns
  expect_error(quiet_constant(net_aggregate(d6, c(1, 1, 2, 2, 2))),
               "one entry per item", class = "psychnet_bad_input")
  # a community on the dropped column is named in the error
  expect_error(quiet_constant(net_aggregate(d6, c(1, 1, 2, 2, 2, 3))),
               "K", class = "psychnet_dropped_column")
  named <- stats::setNames(c(1, 1, 2, 2, 2, 3), names(d6))
  expect_error(quiet_constant(net_aggregate(d6, named)),
               class = "psychnet_dropped_column")
  # NA for the dropped column (or leaving it out of a named vector) works
  ok <- quiet_constant(net_aggregate(d6, c(1, 1, 2, 2, 2, NA)))
  expect_equal(ok, net_aggregate(d, c(1, 1, 2, 2, 2)))
  ok2 <- quiet_constant(net_aggregate(d6, stats::setNames(c(1, 1, 2, 2, 2), names(d))))
  expect_equal(ok2, ok)
})

test_that("W5: listwise row dropping is announced with its count", {
  d <- rv_data()
  d[c(3, 10, 11), "V2"] <- NA
  expect_warning(r <- net_aggregate(d, c(1, 1, 2, 2, 2)),
                 "dropped 3 of 200", class = "psychnet_rows_dropped")
  expect_identical(nrow(r), 197L)
  expect_warning(redundancy(d, cor_method = "pearson"), "dropped 3 of 200",
                 class = "psychnet_rows_dropped")
  # complete data: no warning
  expect_no_warning(net_aggregate(rv_data(), c(1, 1, 2, 2, 2)))
})

test_that("W5: redundancy validates its thresholds", {
  expect_error(redundancy(SRL_Claude, p = 2), class = "psychnet_bad_input")
  expect_error(redundancy(SRL_Claude, threshold = -1), class = "psychnet_bad_input")
  expect_error(redundancy(SRL_Claude, cor_min = "a"), class = "psychnet_bad_input")
})

# ---- W6: factor scores do not fail silently -------------------------------------

test_that("W6: factor method falls back to PCA on < 3 items without a warning", {
  d <- rv_data()
  expect_no_warning(f <- net_aggregate(d, c(1, 1, 2, 2, 2), method = "factor"))
  p <- net_aggregate(d, c(1, 1, 2, 2, 2), method = "pca")
  expect_equal(f[["1"]], p[["1"]])
})

test_that("W6: a factanal failure is a classed warning naming the community", {
  set.seed(2)
  n <- 120
  a <- stats::rnorm(n); b <- stats::rnorm(n)
  d <- data.frame(A1 = a, A2 = b, A3 = a + b,            # singular: factanal fails
                  B1 = stats::rnorm(n), B2 = stats::rnorm(n), B3 = stats::rnorm(n))
  comm <- c("bad", "bad", "bad", "ok", "ok", "ok")
  expect_warning(f <- net_aggregate(d, comm, method = "factor"),
                 "community 'bad'", class = "psychnet_fit_failed")
  p <- net_aggregate(d, comm, method = "pca")
  expect_equal(f[["bad"]], p[["bad"]])
})

test_that("minor: numeric communities are ordered numerically", {
  d <- rv_data()
  r <- net_aggregate(d, c(10, 10, 2, 2, 1))
  expect_named(r, c("1", "2", "10"))
  r2 <- net_aggregate(d, c("z", "z", "a", "a", "a"))
  expect_named(r2, c("a", "z"))                 # character: alphabetical, as before
  r3 <- net_aggregate(d, factor(c("z", "z", "a", "a", "a"), levels = c("z", "a")))
  expect_named(r3, c("z", "a"))
})

# ---- W10: predictability errors are not swallowed --------------------------------

test_that("W10: predictability is attached only where defined, errors propagate", {
  expect_null(cor_network(SRL_Claude)$nodes$predictability)
  expect_false(is.null(ebic_glasso(SRL_Claude)$nodes$predictability))
  net <- ebic_glasso(SRL_Claude)
  net$nodes$predictability <- NULL
  net$precision <- matrix("x", 5, 5)            # defined case, broken input
  expect_error(.attach_predictability(net))
  net$precision <- NULL                          # undefined case: no attempt
  expect_null(.attach_predictability(net)$nodes$predictability)
  # nodewise model without data: undefined at construction, skipped
  b <- dichotomize(SRL_GPT)
  isg <- ising_fit(b)
  expect_false(.psn_has_predictability(isg, NULL))
  expect_true(.psn_has_predictability(isg, b))
})

# ---- W7 / W8 / W9 / N13: plots ----------------------------------------------------

okabe_plus <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7", "#999999", "#000000",
                "#D33F6A", "#4A6FE3")

test_that("W8: the shared palette is Okabe-Ito (plus the diverging pair)", {
  expect_true(all(unlist(.psn_pal) %in% okabe_plus))
  expect_true(all(.psn_series(12) %in% okabe_plus))
  key <- .psn_group_key(9)
  expect_identical(anyDuplicated(paste(key$col, key$pch)), 0L)
  expect_identical(length(unique(key$pch[1:6])), 6L)
  expect_identical(length(unique(key$lty[1:6])), 6L)
})

test_that("W9: centrality plot validates measures with classed conditions", {
  cen <- net_centralities(rv_fit())
  expect_error(draw_par(plot(cen, measures = "nope")), "unknown",
               class = "psychnet_bad_input")
  cen$empty <- NA_real_
  expect_warning(draw_par(plot(cen, measures = c("strength", "empty"))),
                 class = "psychnet_plot_skipped")
  expect_error(draw_par(plot(cen, measures = "empty")),
               class = "psychnet_nothing_to_plot")
  # a partially missing measure still draws, in both types
  cen$part <- c(NA, 1, 2, 3, 4, 5)
  expect_no_warning(draw_par(plot(cen, measures = c("strength", "part"))))
  expect_no_warning(draw_par(plot(cen, measures = "part", type = "line", scale = "z")))
})

test_that("W7/W9: bridge plot draws per-community styling and validates measures", {
  br <- net_bridge(rv_fit(), c(1, 1, 2, 2, 3, 3))
  pp <- expect_no_warning(draw_par(plot(br)))
  expect_identical(pp$after, pp$before)
  expect_error(draw_par(plot(br, measures = "community")),
               class = "psychnet_bad_input")
  expect_error(draw_par(plot(br, measures = "nope")), class = "psychnet_bad_input")
})

test_that("W8/N13: every plot restores par() and draws without warnings", {
  cen <- net_centralities(rv_fit(), measures = c("strength", "expected_influence",
                                                 "betweenness"))
  restores <- function(expr) {
    pp <- expect_no_warning(draw_par(expr))
    expect_identical(pp$after, pp$before)
  }
  restores(plot(cen, type = "bar"))
  restores(plot(cen, type = "line"))
  stab <- structure(list(
    table = data.frame(measure = rep(c("strength", "expected_influence",
                                       "betweenness"), each = 3),
                       drop_prop = rep(c(0.1, 0.5, 0.9), 3),
                       mean_cor = c(0.99, 0.9, 0.6, 0.98, 0.85, 0.5,
                                    0.9, 0.6, 0.2),
                       sd_cor = 0.05),
    cs = c(strength = 0.75, expected_influence = 0.5, betweenness = 0.1),
    threshold = 0.7), class = "psychnet_stability")
  restores(plot(stab))
  skip_on_cran()
  set.seed(1)
  bs <- net_boot(as.matrix(rv_data()), n_boot = 40, cores = 1)
  invisible(lapply(c("edges", "centrality", "edge_diff", "centrality_diff"),
                   function(ty) restores(plot(bs, type = ty))))
  expect_error(draw_par(plot(bs, type = "centrality", measure = "nope")),
               class = "psychnet_bad_input")
  dt <- difference_test(bs, type = "strength")
  restores(plot(dt, style = "box"))
  restores(plot(dt, style = "forest"))
  d <- rv_data(seed = 4, n = 120)
  set.seed(2)
  nct <- net_compare(d[1:60, ], d[61:120, ], iter = 20)
  invisible(lapply(c("strength", "structure", "edges"),
                   function(ty) restores(plot(nct, type = ty))))
})

test_that("W10: unusable column names fall back to positional prediction", {
  # cbind(a, <expr>) leaves an empty column name; it used to make
  # .nodewise_predict() index by "" and the error was swallowed at fit time.
  set.seed(1)
  a <- stats::rbinom(100, 1, 0.5)
  b <- cbind(a, ifelse(stats::runif(100) < 0.9, a, 1 - a))
  fit <- ising_fit(b)
  expect_length(fit$nodes$predictability, 2L)
  expect_true(all(is.finite(fit$nodes$predictability)))
  expect_equal(net_predict(fit, data = b)$predictability,
               net_predict(fit, data = unname(b))$predictability)
})

test_that("W8: significance is carried by shape and line type, not colour only", {
  sty <- .psn_sig_style(c(TRUE, FALSE, NA), 3)
  expect_identical(sty$pch, c(19, 1, 1))
  expect_identical(sty$lty, c(1, 2, 2))
  expect_false(identical(sty$col[1], sty$col[2]))
  expect_identical(.psn_sig_style(NULL, 2)$pch, c(19, 19))
})
