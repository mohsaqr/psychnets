# cor_auto() and cor_network(cor_method = "auto"): polychoric correlations of ordinal items. References:
# psych::polychoric() and qgraph::cor_auto() (which calls lavaan::lavCor()), on the ten bfi items (six
# ordered answers each, 2,632 rows) and on four simulated items cut from a latent AR(1) normal into four
# categories.

items <- utils::read.csv(reference_data("data/bfi10.csv"))

# Two-step estimators of the same bivariate likelihood (thresholds from the margins, then rho). The tolerance
# is that of the reference's one-dimensional optimizer (tightened 2026-11-23 from the inherited 1e-3):
#   psych::polychoric() maximizes over rho with stats::optimize() at its default tol,
#     .Machine$double.eps^0.25 = 1.22e-4, so it may stop up to that far from the optimum;
#   qgraph::cor_auto() calls lavaan::lavCor(), whose optimizer converges far inside 1e-6 (largest difference
#     1.3e-7 in the 2026-11-23 audit, on sparse five-category items).
psych_tolerance <- .Machine$double.eps^0.25
lavaan_tolerance <- 1e-6

simulated_ordinal <- function() {
  set.seed(5)
  latent <- matrix(stats::rnorm(1000L * 4L), 1000L, 4L) %*% chol(0.5^abs(outer(seq_len(4L), seq_len(4L), "-")))
  out <- apply(latent, 2L, \(column) as.integer(cut(column, c(-Inf, -1, 0, 1, Inf))))
  colnames(out) <- paste0("V", seq_len(4L))
  out
}

test_that("polychoric correlations equal psych::polychoric()", {
  skip_without_reference("psych")
  expect_reference(ref_edges(cor_auto(items), "r"), ref_edges(psych::polychoric(items)$rho, "r"),
                   reference = "psych", covers = "cor_auto()", quantity = "polychoric correlation, bfi items",
                   tolerance = psych_tolerance, key = c("from", "to"),
                   note = "cor_auto() returns a bare matrix (API gap); its upper triangle is compared.")
  expect_reference(as.data.frame(cor_network(data = items, cor_method = "auto"), include_zero = TRUE),
                   ref_edges(psych::polychoric(items)$rho), reference = "psych", covers = "cor_network()",
                   quantity = "cor_network(cor_method = \"auto\"), bfi items", tolerance = psych_tolerance,
                   key = c("from", "to"))
  simulated <- simulated_ordinal()
  expect_reference(ref_edges(cor_auto(simulated), "r"), ref_edges(psych::polychoric(simulated)$rho, "r"),
                   reference = "psych", covers = "cor_auto()", quantity = "polychoric correlation, four cut normals",
                   tolerance = psych_tolerance, key = c("from", "to"))
})

test_that("polychoric correlations equal qgraph::cor_auto() (lavaan)", {
  skip_without_reference(c("qgraph", "lavaan"))
  expect_reference(ref_edges(cor_auto(items), "r"), ref_edges(qgraph::cor_auto(items, verbose = FALSE), "r"),
                   reference = "qgraph", covers = "cor_auto()", quantity = "polychoric correlation, bfi items",
                   tolerance = lavaan_tolerance, key = c("from", "to"))
  simulated <- simulated_ordinal()
  expect_reference(ref_edges(cor_auto(simulated), "r"), ref_edges(qgraph::cor_auto(simulated, verbose = FALSE), "r"),
                   reference = "qgraph", covers = "cor_auto()", quantity = "polychoric correlation, four cut normals",
                   tolerance = lavaan_tolerance, key = c("from", "to"))
})

test_that("continuous columns: cor_auto() is the Pearson correlation", {
  expect_reference(ref_edges(cor_auto(SRL_Claude), "r"), ref_edges(stats::cor(SRL_Claude), "r"), reference = "stats",
                   covers = "cor_auto()", quantity = "Pearson correlation of continuous scales, SRL_Claude",
                   tolerance = 1e-12, key = c("from", "to"))
})
