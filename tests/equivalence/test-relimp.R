# relimp_network() and lmg_certificate(): the directed relative-importance network (LMG shares, Lindeman,
# Merenda and Gold 1980) of six bfi items against relaimpo::calc.relimp(type = "lmg"), each item regressed
# on the other five; the shares' sum against the R-squared of stats::lm().

items <- utils::read.csv(reference_data("data/bfi10.csv"))[c(paste0("A", seq_len(5L)), "C1")]
outcomes <- names(items)

#' One row per (predictor, outcome): relaimpo's LMG share of the outcome's R-squared, raw or normalized
relaimpo_table <- function(normalized) {
  do.call(rbind, lapply(outcomes, \(outcome) {
    fit <- relaimpo::calc.relimp(stats::lm(stats::reformulate(setdiff(outcomes, outcome), outcome), data = items),
                                 type = "lmg", rela = normalized)
    data.frame(from = names(fit@lmg), to = outcome, weight = unname(fit@lmg))
  }))
}

test_that("LMG shares equal relaimpo::calc.relimp()", {
  skip_without_reference("relaimpo")
  # Both average R-squared increments over all 120 orders from the same correlation matrix: closed form.
  expect_reference(as.data.frame(relimp_network(data = items), include_zero = TRUE), relaimpo_table(FALSE),
                   reference = "relaimpo", covers = "relimp_network()", quantity = "LMG share of R-squared",
                   tolerance = 1e-12, key = c("from", "to"))
  expect_reference(as.data.frame(relimp_network(data = items, normalized = TRUE), include_zero = TRUE),
                   relaimpo_table(TRUE), reference = "relaimpo", covers = "relimp_network()",
                   quantity = "LMG share, normalized to sum to 1", tolerance = 1e-12, key = c("from", "to"))
  expect_reference(as.data.frame(psychnet(items, method = "relimp"), include_zero = TRUE), relaimpo_table(FALSE),
                   reference = "relaimpo", covers = "psychnet()", quantity = "LMG share, psychnet(method = \"relimp\")",
                   tolerance = 1e-12, key = c("from", "to"))
})

test_that("the R-squared and lmg_certificate() equal stats::lm() and the shares' sums", {
  fit <- relimp_network(data = items)
  r_squared <- vapply(outcomes, \(outcome)
    summary(stats::lm(stats::reformulate(setdiff(outcomes, outcome), outcome), data = items))$r.squared, numeric(1))
  expect_reference(fit$r2, r_squared, reference = "stats", covers = c("relimp_network()", "$.psychnet()"),
                   quantity = "R-squared of each item on the others", tolerance = 1e-12)
  shares <- as.data.frame(fit, include_zero = TRUE)
  summed <- vapply(outcomes, \(outcome) sum(shares$weight[shares$to == outcome]), numeric(1))
  # The certificate is the largest |sum of an outcome's shares - its R-squared|: a value near 1e-16.
  expect_reference(lmg_certificate(fit), max(abs(summed - r_squared)),
                   reference = "independent: the LMG shares summed per outcome against stats::lm()'s R-squared",
                   covers = "lmg_certificate()", quantity = "LMG certificate", tolerance = 1e-12)
  expect_reference(certificate(fit)$certificate, max(abs(summed - r_squared)),
                   reference = "independent: the LMG shares summed per outcome against stats::lm()'s R-squared",
                   covers = "certificate()", quantity = "certificate of relimp_network()", tolerance = 1e-12)
})
