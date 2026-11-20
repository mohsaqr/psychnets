# net_boot(), net_stability(), net_casedrop_reliability(), net_split_reliability() and their as.data.frame()
# methods: the bootstrap and case-dropping procedures of Epskamp, Borsboom and Fried (2018) on a partial
# correlation network (method = "pcor": both sides estimate the same network of a resample, so only the draws
# differ). psychnets draws its resamples inside the functions, so the evidence is statistical: distribution
# summaries against bootnet 1.9.1 (bootnet(type = "nonparametric") and bootnet(type = "case"), each with its own
# seeded draws), and for split halves, which bootnet does not offer, against split halves written in base R.
#
# Tolerances, fixed before the run from the reference's own spread, 4 standard errors of the difference of two
# Monte Carlo estimates (psychnets' B1 draws against the reference's B2):
#   a mean: 4 max sd sqrt(1/B1 + 1/B2);
#   a 2.5% or 97.5% percentile: 4 max sd sqrt(0.025 * 0.975) / dnorm(qnorm(0.975)) sqrt(1/B1 + 1/B2), the
#     large-sample SE of a sample quantile of a normal distribution (the bootstrap distributions of partial
#     correlations and strengths are close to normal at n = 300);
#   a share of subsets above 0.7: 4 sqrt(0.25 (1/B1 + 1/B2)).
# A probe of 10 seeds of psychnets against the bootnet references passed 10 of 10 for each bootstrap and
# case-dropping summary before the full run; the split-half comparison was not probed beyond the three seeds here.

data <- local({
  set.seed(61)
  x <- matrix(stats::rnorm(300 * 6), 300, 6) %*% chol(0.4^abs(outer(1:6, 1:6, "-")))
  colnames(x) <- paste0("V", 1:6)
  x
})
quantile_se <- sqrt(0.025 * 0.975) / stats::dnorm(stats::qnorm(0.975))

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
  lapply(1:3, \(seed) {
    set.seed(seed)
    boot <- net_boot(data, method = "pcor", n_boot = b1, cores = 1)
    table <- as.data.frame(boot)
    expect_reference(table, expected_edges[c("from", "to", "mean")], reference = "bootnet",
                     covers = c("net_boot()", "as.data.frame.psychnet_bootstrap()"),
                     quantity = sprintf("Monte Carlo summary: bootstrap mean of the edge, seed %d", seed),
                     tolerance = 4 * max(edge_sd) * scale, key = c("from", "to"))
    expect_reference(table, expected_edges[c("from", "to", "lower", "upper")], reference = "bootnet",
                     covers = c("net_boot()", "as.data.frame.psychnet_bootstrap()"),
                     quantity = sprintf("Monte Carlo summary: 95%% percentile interval of the edge, seed %d", seed),
                     tolerance = 4 * max(edge_sd) * quantile_se * scale, key = c("from", "to"))
    centrality <- as.data.frame(boot, what = "centrality")
    expect_reference(centrality, data.frame(node = names(strength$sample), strength = strength$sample,
                                            strength_lower = apply(strength$values, 1L, stats::quantile, 0.025, names = FALSE),
                                            strength_upper = apply(strength$values, 1L, stats::quantile, 0.975, names = FALSE)),
                     reference = "bootnet", covers = c("net_boot()", "as.data.frame.psychnet_bootstrap()"),
                     quantity = sprintf("Monte Carlo summary: strength and its 95%% percentile interval, seed %d", seed),
                     tolerance = 4 * max(strength_sd) * quantile_se * scale, key = "node")
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
  lapply(1:2, \(seed) {
    set.seed(seed)
    stability <- net_stability(data, method = "pcor", measures = "strength", drop_prop = drop_prop, iter = b1)
    table <- as.data.frame(stability)
    expect_reference(stats::setNames(table$mean_cor, strength$sizes), stats::setNames(strength$mean, strength$sizes),
                     reference = "bootnet", covers = c("net_stability()", "as.data.frame.psychnet_stability()"),
                     quantity = sprintf("Monte Carlo summary: mean strength correlation per subset size, seed %d", seed),
                     tolerance = 4 * max(strength$sd) * scale)
    expect_reference(stats::setNames(table$prop_above, strength$sizes), stats::setNames(strength$above, strength$sizes),
                     reference = "bootnet", covers = c("net_stability()", "as.data.frame.psychnet_stability()"),
                     quantity = sprintf("Monte Carlo summary: share of strength correlations above 0.7, seed %d", seed),
                     tolerance = 4 * sqrt(0.25) * scale)
    set.seed(seed + 10L)
    casedrop <- as.data.frame(net_casedrop_reliability(data, method = "pcor", drop_prop = drop_prop, iter = b1))
    rows <- casedrop[casedrop$metric == "correlation", ]
    expect_reference(stats::setNames(rows$mean, edges$sizes), stats::setNames(edges$mean, edges$sizes),
                     reference = "bootnet", covers = "net_casedrop_reliability()",
                     quantity = sprintf("Monte Carlo summary: mean edge correlation per subset size, seed %d", seed),
                     tolerance = 4 * max(edges$sd) * scale)
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
