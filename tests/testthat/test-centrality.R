test_that("centrality returns a tidy one-row-per-node data frame", {
  S <- 0.4^abs(outer(1:6, 1:6, "-"))
  fit <- ebic_glasso(cor_matrix = S, n = 300)
  ct <- net_centralities(fit)
  expect_s3_class(ct, "data.frame")
  expect_named(ct, c("node", "strength", "expected_influence"))
  expect_equal(nrow(ct), 6)
  expect_equal(ct$node, fit$labels)
})

test_that("strength and expected_influence match direct row sums", {
  S <- 0.5^abs(outer(1:5, 1:5, "-"))
  fit <- ebic_glasso(cor_matrix = S, n = 300)
  g <- fit$graph
  ct <- net_centralities(fit)
  expect_equal(ct$strength, unname(rowSums(abs(g))))
  expect_equal(ct$expected_influence, unname(rowSums(g)))
})

test_that("centrality accepts a bare matrix", {
  m <- matrix(c(0, 0.3, -0.2, 0.3, 0, 0.1, -0.2, 0.1, 0), 3, 3)
  colnames(m) <- c("a", "b", "c")
  ct <- net_centralities(m)
  expect_equal(ct$node, c("a", "b", "c"))
  expect_equal(ct$strength, c(0.5, 0.4, 0.3))
})

test_that("closeness of an isolated node is 0, as documented", {
  # A path A - B - C plus an isolated node D: D reaches no other node.
  w <- matrix(0, 4, 4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  w["A", "B"] <- w["B", "A"] <- 0.5
  w["B", "C"] <- w["C", "B"] <- 0.25
  closeness <- net_centralities(w, measures = "closeness")
  expect_identical(closeness$closeness[closeness$node == "D"], 0)
  # B reaches A (distance 2) and C (distance 4): 2 / 6.
  expect_equal(closeness$closeness[closeness$node == "B"], 2 / 6, tolerance = 1e-12)
})
