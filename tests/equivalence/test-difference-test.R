# difference_test(): the bootstrapped difference test of edge weights and of node strengths (Epskamp,
# Borsboom and Fried 2018), written out from its definition on the resamples a net_boot() object stores:
# for every pair of items, the percentile interval of the difference of their bootstrapped values and the
# two-sided bootstrap p-value 2 min(P(d > 0), P(d < 0)), with ties at 0 split evenly. net_boot() draws its
# resamples inside the function and is not compared here (see verify.json); this checks what
# difference_test() computes from them.

set.seed(2)
boot <- net_boot(SRL_Claude, method = "pcor", n_boot = 200L, measures = c("strength", "expected_influence"))

by_definition <- function(draws, labels, observed, level = 0.95, adjust = "none") {
  alpha <- (1 - level) / 2
  pairs <- which(upper.tri(diag(length(labels))), arr.ind = TRUE)
  rows <- lapply(seq_len(nrow(pairs)), \(k) {
    d <- draws[, pairs[k, 1L]] - draws[, pairs[k, 2L]]
    ties <- mean(d == 0)
    data.frame(item1 = labels[pairs[k, 1L]], item2 = labels[pairs[k, 2L]],
               obs_diff = observed[pairs[k, 1L]] - observed[pairs[k, 2L]],
               lower = stats::quantile(d, alpha, names = FALSE), upper = stats::quantile(d, 1 - alpha, names = FALSE),
               p_value = min(1, 2 * min(mean(d > 0) + ties / 2, mean(d < 0) + ties / 2)))
  })
  out <- do.call(rbind, rows)
  out$p_value <- stats::p.adjust(out$p_value, method = adjust)
  out
}

definition <- "independent: percentile interval and two-sided bootstrap p-value of the difference, base R"

test_that("edge differences equal the bootstrapped difference test written out", {
  expected <- by_definition(boot$edge_boot, boot$edge_labels, as.data.frame(boot)$observed)
  expect_reference(difference_test(boot, type = "edge"), expected, reference = definition, covers = "difference_test()",
                   quantity = "edge difference, 95% interval and p-value", tolerance = 1e-12, key = c("item1", "item2"),
                   note = "The resamples are read with `$edge_boot` and `$edge_labels`: no accessor (API gap).")
  holm <- by_definition(boot$edge_boot, boot$edge_labels, as.data.frame(boot)$observed, level = 0.9, adjust = "holm")
  expect_reference(difference_test(boot, type = "edge", ci = 0.9, p_adjust = "holm"), holm, reference = definition,
                   covers = "difference_test()", quantity = "edge difference, 90% interval, Holm-adjusted p-value",
                   tolerance = 1e-12, key = c("item1", "item2"))
})

test_that("strength differences equal the bootstrapped difference test written out", {
  strength <- as.data.frame(boot, what = "centrality")$strength
  expected <- by_definition(boot$centrality_boot$strength, boot$node_labels, strength)
  expect_reference(difference_test(boot, type = "strength"), expected, reference = definition,
                   covers = "difference_test()", quantity = "strength difference, 95% interval and p-value",
                   tolerance = 1e-12, key = c("item1", "item2"))
})
