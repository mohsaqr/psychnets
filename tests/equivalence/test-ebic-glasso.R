# ebic_glasso(), glasso_kkt() and certificate(): the EBIC graphical lasso network of the ten bfi items
# (Pearson correlations, 2,632 rows) against qgraph::EBICglasso() and glasso::glasso(); the KKT residual
# against the optimality conditions written out in base R. psychnets' own solver (native = TRUE, the default)
# is independent of glasso; the native = FALSE hand-over to glasso is in test-ebic-glasso-native-false.R.

items <- utils::read.csv(reference_data("data/bfi10.csv"))
correlations <- stats::cor(items)
n_rows <- nrow(items)

#' The KKT residual of the graphical lasso at precision `theta` (Friedman, Hastie and Tibshirani 2008):
#' W = theta^-1 must satisfy W_ii = S_ii, W_ij - S_ij = rho sign(theta_ij) where theta_ij != 0, and
#' |W_ij - S_ij| <= rho where theta_ij = 0. Returns the largest violation.
kkt_residual <- function(theta, S, rho, active_tol = 1e-8) {
  W <- solve(theta)
  off <- upper.tri(theta)
  gap <- (W - S)[off]
  active <- abs(theta[off]) > active_tol
  max(abs(diag(W) - diag(S)),
      abs(gap[active] - rho * sign(theta[off][active])),
      pmax(abs(gap[!active]) - rho, 0))
}

test_that("the EBIC graphical lasso network equals qgraph::EBICglasso() with a tight glasso threshold", {
  skip_without_reference("qgraph")
  fit <- ebic_glasso(data = items, gamma = 0.5, nlambda = 100L, lambda_min_ratio = 0.01, threshold = 0,
                     penalize_diagonal = FALSE, refit = TRUE)
  # qgraph's glasso stops at thr = 1e-4 by default; thr = 1e-10 (passed through EBICglasso's dots) lets the
  # two solvers be read against each other. refit = FALSE in qgraph keeps the regularized network, which is
  # psychnets' refit = TRUE (the penalized optimum re-solved to machine precision).
  reference <- qgraph::EBICglasso(correlations, n = n_rows, gamma = 0.5, nlambda = 100, lambda.min.ratio = 0.01,
                                  penalize.diagonal = FALSE, refit = FALSE, threshold = FALSE,
                                  returnAllResults = TRUE, thr = 1e-10, maxit = 1e5, verbose = FALSE)
  # Two independent solvers of one strictly convex problem: psychnets certifies a KKT residual below 1e-7
  # and glasso runs to thr = 1e-10, so partial correlations agree far below 1e-6.
  expect_reference(as.data.frame(fit, include_zero = TRUE), ref_edges(reference$optnet), reference = "qgraph",
                   covers = c("ebic_glasso()", "as.data.frame.psychnet()"), quantity = "partial correlation",
                   tolerance = 1e-6, key = c("from", "to"))
  # The selected penalty: both compute the same grid in closed form and pick the EBIC minimum.
  expect_reference(fit$lambda, reference$lambda[which.min(reference$ebic)], reference = "qgraph",
                   covers = c("ebic_glasso()", "$.psychnet()"), quantity = "selected penalty", tolerance = 1e-10,
                   measure = "relative")
  expect_reference(as.data.frame(psychnet(items, method = "glasso"), include_zero = TRUE),
                   ref_edges(reference$optnet), reference = "qgraph", covers = "psychnet()",
                   quantity = "partial correlation, psychnet(method = \"glasso\")", tolerance = 1e-6,
                   key = c("from", "to"))
})

test_that("a sparse chain: the EBIC graphical lasso equals qgraph::EBICglasso(), zeros included", {
  skip_without_reference("qgraph")
  set.seed(1)
  chain <- 0.5^abs(outer(seq_len(8L), seq_len(8L), "-"))
  simulated <- matrix(stats::rnorm(1000L * 8L), 1000L, 8L) %*% chol(chain)
  colnames(simulated) <- paste0("V", seq_len(8L))
  S <- stats::cor(simulated)
  reference <- qgraph::EBICglasso(S, n = 1000L, gamma = 0.5, refit = FALSE, threshold = FALSE,
                                  returnAllResults = TRUE, thr = 1e-10, maxit = 1e5, verbose = FALSE)
  # The same convex problem as above; 1e-6 for the same reason.
  expect_reference(as.data.frame(ebic_glasso(cor_matrix = S, n = 1000L), include_zero = TRUE),
                   ref_edges(reference$optnet), reference = "qgraph", covers = "ebic_glasso()",
                   quantity = "partial correlation, AR(1) chain, n = 1000", tolerance = 1e-6, key = c("from", "to"))
})

test_that("the precision matrix at a fixed penalty equals glasso::glasso()", {
  skip_without_reference("glasso")
  chain <- 0.5^abs(outer(seq_len(8L), seq_len(8L), "-"))
  dimnames(chain) <- list(paste0("V", seq_len(8L)), paste0("V", seq_len(8L)))
  # A single-penalty path fits the graphical lasso at that penalty.
  fit <- ebic_glasso(cor_matrix = chain, n = 250, lambda_path = 0.15, penalize_diagonal = TRUE)
  reference <- glasso::glasso(chain, rho = 0.15, penalize.diagonal = TRUE, thr = 1e-10, maxit = 1e5)$wi
  dimnames(reference) <- dimnames(chain)
  # glasso's precision is symmetric to its threshold only: the mean of the two triangles is compared.
  expect_reference(fit$precision, (reference + t(reference)) / 2, reference = "glasso",
                   covers = c("ebic_glasso()", "$.psychnet()"), quantity = "precision, rho = 0.15, diagonal penalized",
                   tolerance = 1e-6)
  fit_items <- ebic_glasso(cor_matrix = correlations, n = n_rows, lambda_path = 0.05)
  reference_items <- glasso::glasso(correlations, rho = 0.05, penalize.diagonal = FALSE, thr = 1e-10, maxit = 1e5)$wi
  dimnames(reference_items) <- dimnames(correlations)
  expect_reference(fit_items$precision, (reference_items + t(reference_items)) / 2, reference = "glasso",
                   covers = "ebic_glasso()", quantity = "precision, bfi items, rho = 0.05", tolerance = 1e-6)
})

test_that("glasso_kkt() and certificate() equal the KKT conditions written out", {
  skip_without_reference("glasso")
  # The residual of glasso's own solution at its default threshold (1e-4): a clearly non-zero value.
  loose <- glasso::glasso(correlations, rho = 0.05, penalize.diagonal = FALSE)$wi
  loose <- (loose + t(loose)) / 2
  # Both sides solve the same 10 x 10 matrix and take a maximum; relative 1e-8 covers the inverse.
  expect_reference(glasso_kkt(loose, correlations, rho = 0.05), kkt_residual(loose, correlations, 0.05),
                   reference = "independent: the graphical lasso KKT conditions at glasso's solution, base R",
                   covers = "glasso_kkt()", quantity = "KKT residual at glasso::glasso(rho = 0.05)", tolerance = 1e-8,
                   measure = "relative")
  fit <- ebic_glasso(data = items)
  expect_reference(certificate(fit)$certificate, kkt_residual(fit$precision, correlations, fit$lambda),
                   reference = "independent: the graphical lasso KKT conditions at glasso's solution, base R",
                   covers = "certificate()", quantity = "certificate of ebic_glasso()",
                   tolerance = 1e-12, note = "Both values are near 1e-10; absolute 1e-12 compares them.")
})
