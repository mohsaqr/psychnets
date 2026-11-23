# net_boot(), net_stability(), net_casedrop_reliability(), net_split_reliability() and their as.data.frame()
# methods: the bootstrap and case-dropping procedures of Epskamp, Borsboom and Fried (2018) on a partial
# correlation network (method = "pcor": both sides estimate the same network of a resample, so only the draws
# differ). psychnets draws its resamples inside the functions, so the evidence is statistical: distribution
# summaries against bootnet 1.9.1 (bootnet(type = "nonparametric") and bootnet(type = "case"), each with its own
# seeded draws), and for split halves, which bootnet does not offer, against split halves written in base R.
#
# Tolerances, fixed before the run from the reference's own spread, 4 standard errors of the difference of two
# Monte Carlo estimates (psychnets' B1 draws against the reference's B2), per cell (2026-11-23: each edge, node
# or subset size has its own standard error, from its own sd; before, the largest sd of the table was used for
# every row):
#   a mean: 4 sd sqrt(1/B1 + 1/B2);
#   a 2.5% or 97.5% percentile: 4 sd sqrt(0.025 * 0.975) / dnorm(qnorm(0.975)) sqrt(1/B1 + 1/B2), the
#     large-sample SE of a sample quantile of a normal distribution (the bootstrap distributions of partial
#     correlations and strengths are close to normal at n = 300);
#   a share of subsets above 0.7: 4 sqrt(q (1 - q) (1/B1 + 1/B2)), q the reference's share shrunk by Laplace's
#     rule (x + 1) / (B2 + 2), so that a share of 1 still has a positive standard error.
# The observed (sample) values are deterministic and are compared separately, at rounding precision.
# The CS-coefficient (Epskamp, Borsboom and Fried 2018) is compared with bootnet::corStability() on a grid of
# drop proportions; it is a step function of the shares, so the comparison is meaningful only where no
# reference share lies within 4 standard errors of the 0.95 certainty, which the test asserts first.
# A probe of 10 seeds of psychnets against the bootnet references passed 10 of 10 for each bootstrap and
# case-dropping summary before the full run; the split-half comparison was not probed beyond the three seeds here.

data <- local({
  set.seed(61)
  x <- matrix(stats::rnorm(300 * 6), 300, 6) %*% chol(0.4^abs(outer(1:6, 1:6, "-")))
  colnames(x) <- paste0("V", 1:6)
  x
})
quantile_se <- sqrt(0.025 * 0.975) / stats::dnorm(stats::qnorm(0.975))

#' One Monte Carlo comparison per row: each row of `actual` against the same row of `expected` at 4 times
#' that row's standard error (`se`, one per row of `expected`), the row's key in the quantity
expect_per_row <- function(actual, expected, se, key, quantity, ...) {
  invisible(lapply(seq_len(nrow(expected)), \(i) {
    row <- expected[i, , drop = FALSE]
    label <- paste(unlist(row[key]), collapse = " ")
    expect_reference(merge(actual, row[key], by = key), row, quantity = sprintf("%s, %s", quantity, label),
                     tolerance = 4 * se[[i]], key = key, ...)
  }))
}

#' A bootnet bootstrap of the pcor network, seeded and restoring the caller's stream
bootnet_run <- function(type, n_boots, seed, ...) {
  network <- bootnet::estimateNetwork(data, default = "pcor", corMethod = "cor", verbose = FALSE)
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()),
          add = TRUE, after = FALSE)
  set.seed(seed)
  bootnet::bootnet(network, nBoots = n_boots, type = type, statistics = c("edge", "strength"), verbose = FALSE, ...)
}

#' bootnet's resamples as one value per draw (columns) and edge or node (rows), with its sample values
bootnet_draws <- function(boot, statistic) {
  table <- as.data.frame(boot$bootTable)
  table <- table[table$type == statistic, ]
  key <- if (statistic == "edge") paste(table$node1, table$node2) else table$node1
  values <- tapply(table$value, list(key, table$name), sum)
  sample <- as.data.frame(boot$sampleTable)
  sample <- sample[sample$type == statistic, ]
  sample_key <- if (statistic == "edge") paste(sample$node1, sample$node2) else sample$node1
  sizes <- tapply(table$nPerson, table$name, unique)
  list(values = values[sample_key, , drop = FALSE], sample = stats::setNames(sample$value, sample_key),
       sizes = sizes[colnames(values)])
}

test_that("net_boot() edge means and percentile intervals agree with bootnet's nonparametric bootstrap", {
  skip_without_reference("bootnet")
  b1 <- 1000L
  b2 <- 2000L
  reference <- bootnet_run("nonparametric", b2, seed = 1)
  edges <- bootnet_draws(reference, "edge")
  strength <- bootnet_draws(reference, "strength")
  edge_sd <- apply(edges$values, 1L, stats::sd)
  strength_sd <- apply(strength$values, 1L, stats::sd)
  scale <- sqrt(1 / b1 + 1 / b2)
  expected_edges <- data.frame(
    from = sub(" .*", "", rownames(edges$values)), to = sub(".* ", "", rownames(edges$values)),
    mean = rowMeans(edges$values),
    lower = apply(edges$values, 1L, stats::quantile, 0.025, names = FALSE),
    upper = apply(edges$values, 1L, stats::quantile, 0.975, names = FALSE))
  expected_strength <- data.frame(node = names(strength$sample),
                                  strength_lower = apply(strength$values, 1L, stats::quantile, 0.025, names = FALSE),
                                  strength_upper = apply(strength$values, 1L, stats::quantile, 0.975, names = FALSE))
  covers <- c("net_boot()", "as.data.frame.psychnet_bootstrap()")
  lapply(1:3, \(seed) {
    set.seed(seed)
    boot <- net_boot(data, method = "pcor", n_boot = b1, cores = 1)
    table <- as.data.frame(boot)
    expect_per_row(table, expected_edges[c("from", "to", "mean")], se = edge_sd * scale, key = c("from", "to"),
                   quantity = sprintf("Monte Carlo summary: bootstrap mean of the edge, seed %d", seed),
                   reference = "bootnet", covers = covers)
    expect_per_row(table, expected_edges[c("from", "to", "lower", "upper")], se = edge_sd * quantile_se * scale,
                   key = c("from", "to"),
                   quantity = sprintf("Monte Carlo summary: 95%% percentile interval of the edge, seed %d", seed),
                   reference = "bootnet", covers = covers)
    centrality <- as.data.frame(boot, what = "centrality")
    expect_per_row(centrality, expected_strength, se = strength_sd * quantile_se * scale, key = "node",
                   quantity = sprintf("Monte Carlo summary: 95%% percentile interval of strength, seed %d", seed),
                   reference = "bootnet", covers = covers)
    # The observed strength is the sample network's, not a resampled value: rounding precision.
    expect_reference(centrality, data.frame(node = names(strength$sample), strength = strength$sample),
                     reference = "bootnet", covers = covers,
                     quantity = sprintf("observed strength of the sample network, seed %d", seed), tolerance = 1e-12,
                     key = "node")
  })
})

#' Mean correlation with the sample values, and share above 0.7, per subset size, of bootnet's case-drop draws
casedrop_summary <- function(draws) {
  correlation <- apply(draws$values, 2L, \(v) stats::cor(v, draws$sample))
  by_size <- split(correlation, draws$sizes)
  list(sizes = as.integer(names(by_size)), mean = vapply(by_size, mean, numeric(1)),
       sd = vapply(by_size, stats::sd, numeric(1)), above = vapply(by_size, \(r) mean(r > 0.7), numeric(1)),
       count = lengths(by_size), mean_abs_dev = vapply(split(colMeans(abs(draws$values - draws$sample)), draws$sizes),
                                                      mean, numeric(1)))
}

test_that("net_stability() and net_casedrop_reliability() agree with bootnet's case-dropping bootstrap", {
  skip_without_reference("bootnet")
  b1 <- 200L
  reference <- bootnet_run("case", 2000L, seed = 2)
  strength <- casedrop_summary(bootnet_draws(reference, "strength"))
  edges <- casedrop_summary(bootnet_draws(reference, "edge"))
  drop_prop <- 1 - strength$sizes / nrow(data)
  scale <- sqrt(1 / b1 + 1 / min(strength$count))
  covers <- c("net_stability()", "as.data.frame.psychnet_stability()")
  shrunk <- (strength$above * strength$count + 1) / (strength$count + 2)
  share_scale <- sqrt(1 / b1 + 1 / strength$count)
  lapply(1:2, \(seed) {
    set.seed(seed)
    stability <- net_stability(data, method = "pcor", measures = "strength", drop_prop = drop_prop, iter = b1)
    table <- transform(as.data.frame(stability), size = as.integer(round(nrow(data) * (1 - drop_prop))))
    expect_per_row(table, data.frame(size = strength$sizes, mean_cor = strength$mean),
                   se = strength$sd * sqrt(1 / b1 + 1 / strength$count), key = "size",
                   quantity = sprintf("Monte Carlo summary: mean strength correlation, seed %d, subset size", seed),
                   reference = "bootnet", covers = covers)
    expect_per_row(table, data.frame(size = strength$sizes, prop_above = strength$above),
                   se = sqrt(shrunk * (1 - shrunk)) * share_scale, key = "size",
                   quantity = sprintf("Monte Carlo summary: share of strength correlations above 0.7, seed %d, subset size", seed),
                   reference = "bootnet", covers = covers)
    set.seed(seed + 10L)
    casedrop <- as.data.frame(net_casedrop_reliability(data, method = "pcor", drop_prop = drop_prop, iter = b1))
    rows <- transform(subset(casedrop, metric == "correlation"), size = as.integer(round(nrow(data) * (1 - drop_prop))))
    expect_per_row(rows, data.frame(size = edges$sizes, mean = edges$mean),
                   se = edges$sd * sqrt(1 / b1 + 1 / edges$count), key = "size",
                   quantity = sprintf("Monte Carlo summary: mean edge correlation, seed %d, subset size", seed),
                   reference = "bootnet", covers = "net_casedrop_reliability()")
  })
})

test_that("net_stability()'s CS-coefficient equals bootnet::corStability()", {
  skip_without_reference("bootnet")
  b1 <- 200L
  # Drop proportions 0.3, 0.6 and 0.9 (bootnet's caseMin, caseMax, caseN): far enough from the CS boundary of
  # this network that the step function is decided by the data, not by the Monte Carlo draws.
  reference <- bootnet_run("case", 2000L, seed = 4, caseMin = 0.3, caseMax = 0.9, caseN = 3)
  strength <- casedrop_summary(bootnet_draws(reference, "strength"))
  drop_prop <- 1 - strength$sizes / nrow(data)
  shrunk <- (strength$above * strength$count + 1) / (strength$count + 2)
  # The precondition of a meaningful comparison: every reference share is more than 4 standard errors of the
  # difference of two Monte Carlo shares away from the certainty 0.95.
  expect_true(all(abs(strength$above - 0.95) > 4 * sqrt(shrunk * (1 - shrunk) * (1 / b1 + 1 / strength$count))))
  expected <- data.frame(measure = "strength",
                         cs = unname(bootnet::corStability(reference, cor = 0.7, statistics = "strength",
                                                           verbose = FALSE)[["strength"]]))
  lapply(1:3, \(seed) {
    set.seed(seed)
    stability <- net_stability(data, method = "pcor", measures = "strength", drop_prop = drop_prop, iter = b1)
    # Both are one of the drop proportions: equal or not.
    expect_reference(as.data.frame(stability, what = "cs"), expected, reference = "bootnet",
                     covers = c("net_stability()", "as.data.frame.psychnet_stability()"),
                     quantity = sprintf("Monte Carlo summary: CS-coefficient of strength, seed %d", seed),
                     tolerance = 1e-12, key = "measure")
  })
})

test_that("net_split_reliability() agrees with split halves written in base R", {
  b1 <- 300L
  b2 <- 2000L
  pcor <- \(x) {
    precision <- solve(stats::cor(x))
    p <- -stats::cov2cor(precision)
    p[upper.tri(p)]
  }
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv())
  set.seed(3)
  splits <- t(vapply(seq_len(b2), \(i) {
    half <- sample.int(nrow(data), nrow(data) / 2)
    d <- abs(pcor(data[half, ]) - pcor(data[-half, ]))
    c(mean_abs_dev = mean(d), median_abs_dev = stats::median(d),
      correlation = stats::cor(pcor(data[half, ]), pcor(data[-half, ])), max_abs_dev = max(d))
  }, numeric(4)))
  if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
  lapply(1:3, \(seed) {
    set.seed(seed)
    table <- as.data.frame(net_split_reliability(data, method = "pcor", iter = b1))
    lapply(colnames(splits), \(metric) {
      expect_reference(table$mean[table$metric == metric], mean(splits[, metric]),
                       reference = "independent: split halves of the partial correlation network, base R",
                       covers = "net_split_reliability()",
                       quantity = sprintf("Monte Carlo summary: mean %s over splits, seed %d", metric, seed),
                       tolerance = 4 * stats::sd(splits[, metric]) * sqrt(1 / b1 + 1 / b2))
    })
  })
})
