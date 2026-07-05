# Regression tests for the 2026-06-30 code-review fixes.

test_that("#1 group net_boot/net_stability forward the group's labels", {
  set.seed(1)
  g <- data.frame(matrix(stats::rnorm(160 * 4), 160, 4), grp = rep(c("A", "B"), 80))
  gp <- psychnet(g, group = "grp", method = "glasso", labels = c("W", "X", "Y", "Z"))
  bb <- net_boot(gp, n_boot = 20, cores = 1)
  expect_identical(bb[["A"]]$node_labels, c("W", "X", "Y", "Z"))
  # edge labels in the bootstrap table use the forwarded node labels, not V1..V4
  expect_true(all(grepl("^[WXYZ]--[WXYZ]$", bb[["A"]]$edges$from |> paste0("--", bb[["A"]]$edges$to))))
  expect_s3_class(net_stability(gp, drop_prop = c(0.3, 0.5), iter = 10),
                  "psychnet_stability_group")
})

test_that("#2 tmfg_network accepts and ignores a stray n", {
  S <- stats::cor(matrix(stats::rnorm(200 * 5), 200, 5))
  expect_s3_class(tmfg_network(cor_matrix = S, n = 200), "psychnet")
})

test_that("#3 net_centralities/net_predict reject multilevel groups clearly", {
  set.seed(1)
  ev <- do.call(rbind, lapply(1:30, function(a) do.call(rbind, lapply(1:3,
    function(s) data.frame(actor = a, session = s,
      grp = if (a <= 15) "A" else "B",
      action = sample(letters[1:4], 5, TRUE))))))
  mlg <- psychnet(ev, actor = "actor", session = "session", action = "action",
                  group = "grp", standardize = FALSE)
  expect_error(net_centralities(mlg), "multilevel")
  expect_error(net_predict(mlg), "multilevel")
})

test_that("#4 net_aggregate(canonical) survives a rank-deficient/wide block", {
  set.seed(1)
  d <- as.data.frame(matrix(stats::rnorm(8 * 12), 8, 12)); names(d) <- paste0("v", 1:12)
  agg <- net_aggregate(d, communities = c(rep(1, 2), rep(2, 10)), method = "canonical")
  expect_s3_class(agg, "psychnet")
  expect_true(all(is.finite(agg$weights)))
})

test_that("#5 .psn_refit pads a dropped (constant) column back to the full node set", {
  set.seed(1)
  m <- matrix(stats::rnorm(60 * 4), 60, 4); colnames(m) <- paste0("V", 1:4)
  m[, 3] <- 1                                            # constant -> psychnet drops it
  r <- psychnets:::.psn_refit(m, "glasso", colnames(m))
  expect_equal(dim(r$weights), c(4L, 4L))               # re-expanded, not skipped
  expect_equal(unname(r$weights[3, ]), rep(0, 4))       # dropped node is isolated
})

test_that("#6 group net_compare ignores a stray non-numeric column", {
  set.seed(1)
  g <- data.frame(matrix(stats::rnorm(160 * 4), 160, 4),
                  note = sample(letters, 160, TRUE), grp = rep(c("A", "B"), 80))
  gp <- psychnet(g, group = "grp", method = "pcor")
  cmp <- net_compare(gp, c("A", "B"), iter = 30)
  expect_s3_class(cmp, "psychnet_nct")
})

test_that("#7 redundancy print advertises r (not |r|), matching the raw filter", {
  rd <- redundancy(SRL_Claude)
  expect_true(any(grepl("r > ", capture.output(print(rd)), fixed = TRUE)))
  expect_false(any(grepl("|r|", capture.output(print(rd)), fixed = TRUE)))
})

test_that("#8 net_smallworld returns NA (not Inf/NaN) when reference has no triangles", {
  A <- matrix(0, 5, 5); A[1, 2] <- A[2, 1] <- 1; A[3, 4] <- A[4, 3] <- 1
  colnames(A) <- paste0("n", 1:5)
  expect_warning(sw <- net_smallworld(A, n_rand = 20, seed = 1), "undefined")
  expect_true(is.na(sw$smallworldness))
})

test_that("#9 difference box matrix orders items by observed value", {
  set.seed(1)
  bs <- net_boot(matrix(stats::rnorm(120 * 4), 120, 4), n_boot = 30, cores = 1)
  pdf(file = tempfile(fileext = ".pdf")); on.exit(dev.off())
  expect_no_error(plot(difference_test(bs, type = "strength")))
})
