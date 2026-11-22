# cor_network(), pcor_network() and psychnet(method = "cor" / "pcor"): Pearson and Spearman correlations and
# full-order partial correlations of the ten bfi items (2,632 complete rows), with and without the
# significance rule `alpha`. References: stats::cor() and the inverse of the correlation matrix (the
# definition), qgraph(graph = "cor" / "pcor"), corpcor::cor2pcor() and ppcor::pcor() (partial correlations
# and their p-values), stats::cor.test() (correlation p-values).

items <- utils::read.csv(reference_data("data/bfi10.csv"))
srl <- SRL_Claude

# Closed-form quantities (a correlation; one inverse of a 10 x 10 matrix): only the last bits may differ.
exact <- 1e-12

test_that("Pearson and Spearman correlation networks equal stats::cor()", {
  expect_reference(as.data.frame(cor_network(data = items, cor_method = "pearson"), include_zero = TRUE),
                   ref_edges(stats::cor(items)), reference = "stats",
                   covers = c("cor_network()", "as.data.frame.psychnet()"), quantity = "Pearson correlation",
                   tolerance = exact, key = c("from", "to"))
  expect_reference(as.data.frame(cor_network(data = items, cor_method = "spearman"), include_zero = TRUE),
                   ref_edges(stats::cor(items, method = "spearman")), reference = "stats",
                   covers = "cor_network()", quantity = "Spearman correlation", tolerance = exact,
                   key = c("from", "to"))
  expect_reference(as.data.frame(psychnet(srl, method = "cor"), include_zero = TRUE), ref_edges(stats::cor(srl)),
                   reference = "stats", covers = c("psychnet()", "as.data.frame.psychnet()"),
                   quantity = "Pearson correlation, psychnet(method = \"cor\"), SRL_Claude", tolerance = exact,
                   key = c("from", "to"))
  expect_reference(summary(cor_network(data = srl)), ref_edges(stats::cor(srl)), reference = "stats",
                   covers = "summary.psychnet()", quantity = "Pearson correlation, summary(), SRL_Claude",
                   tolerance = exact, key = c("from", "to"))
})

test_that("partial correlations equal the negated, rescaled inverse of the correlation matrix", {
  # pcor_ij = -P_ij / sqrt(P_ii P_jj), P the inverse correlation matrix
  partial <- -stats::cov2cor(solve(stats::cor(items)))
  expect_reference(as.data.frame(pcor_network(data = items), include_zero = TRUE), ref_edges(partial),
                   reference = "independent: -cov2cor(solve(cor(x))), base R", covers = "pcor_network()",
                   quantity = "partial correlation", tolerance = exact, key = c("from", "to"))
  srl_partial <- -stats::cov2cor(solve(stats::cor(srl)))
  expect_reference(as.data.frame(psychnet(srl, method = "pcor"), include_zero = TRUE), ref_edges(srl_partial),
                   reference = "independent: -cov2cor(solve(cor(x))), base R", covers = "psychnet()",
                   quantity = "partial correlation, psychnet(method = \"pcor\"), SRL_Claude", tolerance = exact,
                   key = c("from", "to"))
})

test_that("correlation and partial correlation networks equal qgraph(graph = \"cor\" / \"pcor\")", {
  skip_without_reference("qgraph")
  correlations <- qgraph::getWmat(qgraph::qgraph(stats::cor(items), graph = "cor", DoNotPlot = TRUE))
  partials <- qgraph::getWmat(qgraph::qgraph(stats::cor(items), graph = "pcor", DoNotPlot = TRUE))
  # qgraph's weights matrix carries the node labels as dimnames
  expect_reference(as.data.frame(cor_network(data = items), include_zero = TRUE), ref_edges(correlations),
                   reference = "qgraph", covers = "cor_network()", quantity = "Pearson correlation",
                   tolerance = exact, key = c("from", "to"))
  # qgraph computes the partial correlations through corpcor's pseudo-inverse: 1e-10 allows for that
  # different route to the same inverse.
  expect_reference(as.data.frame(pcor_network(data = items), include_zero = TRUE), ref_edges(partials),
                   reference = "qgraph", covers = "pcor_network()", quantity = "partial correlation",
                   tolerance = 1e-10, key = c("from", "to"))
})

test_that("partial correlations and their p-values equal corpcor and ppcor", {
  skip_without_reference(c("corpcor", "ppcor"))
  from_corpcor <- corpcor::cor2pcor(stats::cor(items))
  dimnames(from_corpcor) <- list(names(items), names(items))
  expect_reference(as.data.frame(pcor_network(data = items), include_zero = TRUE), ref_edges(from_corpcor),
                   reference = "corpcor", covers = "pcor_network()", quantity = "partial correlation",
                   tolerance = 1e-10, key = c("from", "to"))
  test <- ppcor::pcor(items)
  dimnames(test$estimate) <- dimnames(test$p.value) <- list(names(items), names(items))
  fit <- pcor_network(data = items, alpha = 0.05)
  expect_reference(as.data.frame(pcor_network(data = items), include_zero = TRUE), ref_edges(test$estimate),
                   reference = "ppcor", covers = "pcor_network()", quantity = "partial correlation",
                   tolerance = 1e-10, key = c("from", "to"))
  # The t test of a partial correlation on n - 2 - (p - 2) degrees of freedom. p-values run down to 1e-50:
  # relative, 1e-8 for the tail of pt() through two routes.
  expect_reference(ref_edges(fit$p_values, "p"), ref_edges(test$p.value, "p"), reference = "ppcor",
                   covers = c("pcor_network()", "$.psychnet()"), quantity = "partial correlation p-value",
                   tolerance = 1e-8, measure = "relative", key = c("from", "to"),
                   note = "psychnets keeps the p-values as a matrix reached with `$p_values`: no accessor (API gap).")
  significant <- ifelse(test$p.value < 0.05, test$estimate, 0)
  expect_reference(as.data.frame(fit, include_zero = TRUE), ref_edges(significant), reference = "ppcor",
                   covers = "pcor_network()", quantity = "partial correlation kept at alpha = 0.05",
                   tolerance = 1e-10, key = c("from", "to"))
})

test_that("correlation p-values and the alpha rule equal stats::cor.test()", {
  pairs <- ref_edges(stats::cor(srl), "r")
  pairs$p <- vapply(seq_len(nrow(pairs)), \(k) stats::cor.test(srl[[pairs$from[[k]]]], srl[[pairs$to[[k]]]])$p.value,
                    numeric(1))
  fit <- cor_network(data = srl, alpha = 0.05, adjust = "holm")
  expect_reference(ref_edges(fit$p_values, "p"), pairs[c("from", "to", "p")], reference = "stats",
                   covers = c("cor_network()", "$.psychnet()"), quantity = "correlation p-value, SRL_Claude",
                   tolerance = 1e-8, measure = "relative", key = c("from", "to"),
                   note = "stats::cor.test(): t = r sqrt((n - 2) / (1 - r^2)); relative, the p-values reach 1e-100.")
  kept <- data.frame(from = pairs$from, to = pairs$to,
                     weight = ifelse(stats::p.adjust(pairs$p, method = "holm") < 0.05, pairs$r, 0))
  expect_reference(as.data.frame(fit, include_zero = TRUE), kept, reference = "stats", covers = "cor_network()",
                   quantity = "correlation kept at alpha = 0.05 after Holm's adjustment, SRL_Claude",
                   tolerance = exact, key = c("from", "to"))
})

# Parameter recovery: data drawn by simulab from a known AR(1) correlation matrix (rho = 0.4, five variables,
# n = 5000). The sampling SE of a correlation is at most 1 / sqrt(n) (at rho = 0): the bound, fixed before the
# run, is 4 / sqrt(5000) = 0.057.
test_that("psychnet(method = 'cor') recovers the correlation matrix simulab drew the data from", {
  invisible(lapply(c(1L, 2L, 3L), \(seed) {
    drawn <- simulab::simulate_correlation(n = 5000L, means = rep(0, 5L), rho = 0.4, structure = "ar1",
                                           variable_names = paste0("v", seq_len(5L)), seed = seed)
    truth <- as.data.frame(drawn, what = "correlation")
    names_ <- paste0("v", seq_len(5L))
    truth_matrix <- matrix(0, 5L, 5L, dimnames = list(names_, names_))
    truth_matrix[cbind(match(truth$row, names_), match(truth$column, names_))] <- truth$correlation
    expect_reference(as.data.frame(psychnet(subset(as.data.frame(drawn), select = -id), method = "cor"), include_zero = TRUE),
                     ref_edges(truth_matrix),
                     reference = "independent: recovery of the parameters of simulab simulate_correlation()",
                     covers = "psychnet()", quantity = sprintf("correlation, AR(1) rho 0.4, n 5000, seed %d", seed),
                     tolerance = 4 / sqrt(5000), key = c("from", "to"))
  }))
})
