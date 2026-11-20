# huge_network(): the nonparanormal graphical model. The transformed correlation matrix of each `npn`
# option against huge::huge.npn(), then the network against qgraph::EBICglasso() run on huge's transformed
# correlations (the model huge::huge() fits with the nonparanormal transform, selected by the EBIC).

items <- utils::read.csv(reference_data("data/bfi10.csv"))
n_rows <- nrow(items)

huge_correlation <- function(method) {
  transformed <- huge::huge.npn(as.matrix(items), npn.func = method, verbose = FALSE)
  if (identical(method, "skeptic")) transformed else stats::cor(transformed)
}

test_that("the network on each nonparanormal transform equals qgraph::EBICglasso() on huge.npn()", {
  skip_without_reference(c("huge", "qgraph"))
  networks <- lapply(c(shrinkage = "shrinkage", truncation = "truncation", skeptic = "skeptic"), \(method) {
    S <- huge_correlation(method)
    qgraph::EBICglasso(S, n = n_rows, gamma = 0.5, refit = FALSE, threshold = FALSE, thr = 1e-10, maxit = 1e5,
                       verbose = FALSE)
  })
  # Two independent graphical lasso solvers of one convex problem (glasso run to thr = 1e-10); the transforms
  # are closed-form. 1e-6, as for ebic_glasso().
  expect_reference(as.data.frame(huge_network(data = items, npn = "shrinkage"), include_zero = TRUE),
                   ref_edges(networks$shrinkage), reference = "huge", covers = "huge_network()",
                   quantity = "partial correlation, shrinkage transform", tolerance = 1e-6, key = c("from", "to"),
                   note = "Reference: huge::huge.npn(npn.func = \"shrinkage\"), then qgraph::EBICglasso(gamma = 0.5).")
  expect_reference(as.data.frame(huge_network(data = items, npn = "truncation"), include_zero = TRUE),
                   ref_edges(networks$truncation), reference = "huge", covers = "huge_network()",
                   quantity = "partial correlation, truncation transform", tolerance = 1e-6, key = c("from", "to"))
  expect_reference(as.data.frame(huge_network(data = items, npn = "skeptic"), include_zero = TRUE),
                   ref_edges(networks$skeptic), reference = "huge", covers = "huge_network()",
                   quantity = "partial correlation, skeptic", tolerance = 1e-6, key = c("from", "to"))
  expect_reference(as.data.frame(psychnet(items, method = "huge"), include_zero = TRUE), ref_edges(networks$shrinkage),
                   reference = "huge", covers = "psychnet()", quantity = "partial correlation, psychnet(method = \"huge\")",
                   tolerance = 1e-6, key = c("from", "to"))
})

test_that("the transformed correlation matrices equal huge::huge.npn()", {
  skip_without_reference("huge")
  # Closed-form: normal scores of ranks, or 2 sin(pi rho_S / 6). The skeptic matrix of these items is
  # positive definite, so no projection is applied on either side.
  expect_reference(ref_edges(huge_network(data = items, npn = "shrinkage")$cor_matrix, "r"),
                   ref_edges(huge_correlation("shrinkage"), "r"), reference = "huge",
                   covers = c("huge_network()", "$.psychnet()"), quantity = "correlation after the shrinkage transform",
                   tolerance = 1e-12, key = c("from", "to"))
  expect_reference(ref_edges(huge_network(data = items, npn = "truncation")$cor_matrix, "r"),
                   ref_edges(huge_correlation("truncation"), "r"), reference = "huge", covers = "huge_network()",
                   quantity = "correlation after the truncation transform", tolerance = 1e-12, key = c("from", "to"))
  expect_reference(ref_edges(huge_network(data = items, npn = "skeptic")$cor_matrix, "r"),
                   ref_edges(huge_correlation("skeptic"), "r"), reference = "huge", covers = "huge_network()",
                   quantity = "skeptic correlation", tolerance = 1e-12, key = c("from", "to"))
})
