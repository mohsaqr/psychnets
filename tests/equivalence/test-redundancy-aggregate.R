# redundancy() and net_aggregate(): the goldbricker procedure (Hittner et al. 2003 test of dependent
# overlapping correlations) against networktools::goldbricker(); composite scores and community networks
# against networktools::net_reduce() and their definitions.

#' Six continuous items: two pairs of near-duplicates and two distinct items
redundant_items <- function() {
  set.seed(11)
  f <- matrix(stats::rnorm(500L * 3L), 500L, 3L) %*% chol(0.4^abs(outer(seq_len(3L), seq_len(3L), "-")))
  data.frame(a1 = f[, 1L] + stats::rnorm(500L, sd = 0.3), a2 = f[, 1L] + stats::rnorm(500L, sd = 0.3),
             b1 = f[, 2L] + stats::rnorm(500L, sd = 0.3), b2 = f[, 2L] + stats::rnorm(500L, sd = 0.3),
             c1 = f[, 3L] + stats::rnorm(500L, sd = 0.5), d1 = stats::rnorm(500L))
}

test_that("the proportion of significantly different correlations equals networktools::goldbricker()", {
  skip_without_reference("networktools")
  items <- redundant_items()
  reference <- networktools::goldbricker(items, p = 0.05, method = "hittner2003", threshold = 0.25, corMin = 0.5,
                                         progressbar = FALSE)
  ours <- redundancy(items, p = 0.05, threshold = 0.25, cor_min = 0.5, cor_method = "pearson")
  # A share of correlations whose Hittner test is below p: equal unless a p-value sits at 0.05 within rounding.
  expect_reference(attr(ours, "proportion_matrix"), reference$proportion_matrix, reference = "networktools",
                   covers = "redundancy()", quantity = "proportion of correlations that differ", tolerance = 1e-10,
                   note = "The matrix is reached with attr(, \"proportion_matrix\"): no accessor (API gap).")
})

test_that("the PCA composite equals networktools::net_reduce() and stats::prcomp()", {
  skip_without_reference("networktools")
  items <- SRL_Claude
  composite <- net_aggregate(items, communities = c(1, 1, 2, 2, 2), method = "pca")
  reduced <- networktools::net_reduce(items, badpairs = list(c("CSU", "IV")), method = "PCA")
  # Both are the first principal component of the two standardized items: scores identical up to sign and the
  # scale of the component; their correlation is the comparable number.
  expect_reference(stats::cor(composite[["1"]], reduced[["PCA.CSU.IV"]])^2, 1, reference = "networktools",
                   covers = "net_aggregate()", quantity = "squared correlation of the CSU-IV composite with networktools'",
                   tolerance = 1e-10)
  pc <- stats::prcomp(scale(items[c("CSU", "IV")]))$x[, 1L]
  expect_reference(stats::cor(composite[["1"]], pc)^2, 1, reference = "independent: first principal component, stats::prcomp()",
                   covers = "net_aggregate()", quantity = "squared correlation of the CSU-IV composite with prcomp()'s",
                   tolerance = 1e-10)
})

test_that("mean and sum composites equal their definitions on standardized items", {
  items <- SRL_Claude
  standardized <- scale(items)
  means <- data.frame(`1` = rowMeans(standardized[, c("CSU", "IV")]), `2` = rowMeans(standardized[, c("SE", "SR", "TA")]),
                      check.names = FALSE)
  definition <- "independent: row means and sums of the standardized items of each community, base R"
  expect_reference(as.matrix(net_aggregate(items, communities = c(1, 1, 2, 2, 2), method = "mean")), as.matrix(means), reference = definition,
                   covers = "net_aggregate()", quantity = "mean composite", tolerance = 1e-12)
  sums <- data.frame(`1` = rowSums(standardized[, c("CSU", "IV")]), `2` = rowSums(standardized[, c("SE", "SR", "TA")]),
                     check.names = FALSE)
  expect_reference(as.matrix(net_aggregate(items, communities = c(1, 1, 2, 2, 2), method = "sum")), as.matrix(sums), reference = definition,
                   covers = "net_aggregate()", quantity = "sum composite", tolerance = 1e-12)
})

test_that("the RV coefficient between communities equals Escoufier's definition", {
  items <- scale(SRL_Claude, scale = FALSE)
  first <- items[, c("CSU", "IV")]
  second <- items[, c("SE", "SR", "TA")]
  # RV(X, Y) = tr(XX'YY') / sqrt(tr(XX'XX') tr(YY'YY')) on column-centred blocks (Escoufier 1973)
  cross <- \(a, b) sum(crossprod(a, b)^2)
  rv <- cross(first, second) / sqrt(cross(first, first) * cross(second, second))
  expect_reference(as.data.frame(net_aggregate(SRL_Claude, communities = c(1, 1, 2, 2, 2), method = "rv")),
                   data.frame(from = "1", to = "2", weight = rv),
                   reference = "independent: Escoufier's RV coefficient, base R", covers = "net_aggregate()",
                   quantity = "RV coefficient between the two communities", tolerance = 1e-12, key = c("from", "to"),
                   note = "On centred, unstandardized items: the RV of net_aggregate() is the same with scale = TRUE or FALSE.")
})
