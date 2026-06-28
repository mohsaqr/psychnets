# certificate() and node_predictability() are the tidy accessors the vignettes
# use in place of $kkt / method-specific certificate functions / $r2.

test_that("certificate() reports a tidy, certified residual for every method", {
  set.seed(1)
  X <- matrix(stats::rnorm(300 * 6), 300, 6)
  colnames(X) <- paste0("V", 1:6)
  for (m in c("glasso", "ggm", "tmfg", "logo", "relimp")) {
    net <- psychnet(X, method = m)
    cert <- certificate(net)
    expect_s3_class(cert, "data.frame")
    expect_identical(nrow(cert), 1L)
    expect_named(cert, c("method", "certificate", "kind", "certified"))
    expect_identical(cert$method, m)
    expect_true(cert$certified)
  }
})

test_that("certificate() reports NA / 'none' for cor and pcor", {
  S <- 0.4^abs(outer(1:5, 1:5, "-"))
  cert <- certificate(cor_network(cor_matrix = S))
  expect_true(is.na(cert$certificate))
  expect_identical(cert$kind, "none")
  expect_true(cert$certified)
})

test_that("node_predictability() returns a node-order vector clamped to [0,1]", {
  S <- 0.5^abs(outer(1:6, 1:6, "-"))
  net <- ebic_glasso(cor_matrix = S, n = 300)
  v <- node_predictability(net)
  expect_length(v, 6L)
  expect_true(all(v >= 0 & v <= 1))
  expect_identical(names(v), net$nodes$label)
})

test_that("net_predict() and node_predictability() work for relimp via stored r2", {
  S <- 0.4^abs(outer(1:5, 1:5, "-"))
  ri <- relimp_network(cor_matrix = S)
  pr <- net_predict(ri)
  expect_identical(nrow(pr), 5L)
  expect_true(all(pr$metric == "R2"))
  expect_length(node_predictability(ri), 5L)
})

test_that("predictability is stored on the node table at fit time", {
  S <- 0.5^abs(outer(1:6, 1:6, "-"))
  g <- ebic_glasso(cor_matrix = S, n = 300)
  expect_true("predictability" %in% names(g$nodes))
  expect_true(all(g$nodes$predictability >= 0 & g$nodes$predictability <= 1))
  # glasso is the one network flagged to draw its ring by default
  expect_true(isTRUE(g$meta$predictability_default))
  expect_false(isTRUE(logo_network(cor_matrix = S, n = 300)$meta$predictability_default))
})

test_that("nodewise predictability is stored without re-passing data", {
  set.seed(3)
  b <- matrix(rbinom(300 * 5, 1, 0.5), 300, 5); colnames(b) <- paste0("B", 1:5)
  fit <- ising_fit(b)
  expect_true("predictability" %in% names(fit$nodes))
  expect_length(fit$nodes$predictability, 5L)
})

test_that("psychnet() accepts both our short names and the qgraph aliases", {
  X <- matrix(stats::rnorm(200 * 5), 200, 5)
  pairs <- list(c("glasso", "EBICglasso"), c("ggm", "ggmModSelect"),
                c("tmfg", "TMFG"), c("logo", "LoGo"))
  for (p in pairs) {
    short <- psychnet(X, p[1])$method
    alias <- psychnet(X, p[2])$method
    expect_identical(short, p[1])      # canonical stored name is the short one
    expect_identical(alias, p[1])      # alias resolves to the same canonical
  }
})
