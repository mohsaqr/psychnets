# relimp_network() and lmg_certificate(): the directed relative-importance network (LMG shares, Lindeman,
# Merenda and Gold 1980) of six bfi items against relaimpo::calc.relimp(type = "lmg"), each item regressed
# on the other five; the R-squared against stats::lm(); lmg_certificate() on a share table that is not the LMG
# decomposition (relaimpo's first-entry shares) against relaimpo's own gap.

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

test_that("the R-squared equals stats::lm(); at the fit the certificates are rounding residues", {
  fit <- relimp_network(data = items)
  r_squared <- vapply(outcomes, \(outcome)
    summary(stats::lm(stats::reformulate(setdiff(outcomes, outcome), outcome), data = items))$r.squared, numeric(1))
  expect_reference(fit$r2, r_squared, reference = "stats", covers = c("relimp_network()", "$.psychnet()"),
                   quantity = "R-squared of each item on the others", tolerance = 1e-12)
  shares <- as.data.frame(fit, include_zero = TRUE)
  summed <- vapply(outcomes, \(outcome) sum(shares$weight[shares$to == outcome]), numeric(1))
  # At the fit both certificates are rounding residues of an exact identity (about 1e-16 on each side): a
  # function returning any value below a tolerance would agree, so nothing is judged there. The certificate
  # is compared below, where it is clearly not zero.
  expect_not_comparable(lmg_certificate(fit), max(abs(summed - r_squared)),
                        reference = "independent: the LMG shares summed per outcome against stats::lm()'s R-squared",
                        covers = "lmg_certificate()", quantity = "LMG certificate at the fit (rounding residue)",
                        reason = "Both values are rounding residues near 1e-16 of an identity that holds exactly; their difference cannot tell a right certificate from a constant 0.")
  expect_not_comparable(certificate(fit)$certificate, max(abs(summed - r_squared)),
                        reference = "independent: the LMG shares summed per outcome against stats::lm()'s R-squared",
                        covers = "certificate()", quantity = "certificate of relimp_network() at the fit (rounding residue)",
                        reason = "Both values are rounding residues near 1e-16 of an identity that holds exactly; certificate() is compared where it is not zero in test-ebic-glasso.R, test-ggm-modselect.R and test-tmfg-logo.R.")
})

test_that("lmg_certificate() of a share table that is not the LMG decomposition equals relaimpo's gap", {
  skip_without_reference("relaimpo")
  # The first-entry shares (each predictor's squared correlation with the outcome, relaimpo type = "first")
  # do not sum to the R-squared when predictors correlate, so the certificate is far from 0 (about 0.33 here).
  # relaimpo computes both the shares and the R-squared; the certificate is the largest gap over outcomes.
  first <- lapply(outcomes, \(outcome) relaimpo::calc.relimp(
    stats::lm(stats::reformulate(setdiff(outcomes, outcome), outcome), data = items), type = "first", rela = FALSE))
  gap <- max(vapply(first, \(fit) abs(sum(fit@first) - fit@R2), numeric(1)))
  shares <- vapply(first, \(fit) replace(stats::setNames(numeric(length(outcomes)), outcomes), names(fit@first), fit@first),
                   numeric(length(outcomes)))
  dimnames(shares) <- list(outcomes, outcomes)
  # lmg_certificate() has no entry point for a share table (API gap): the fitted network's shares are replaced.
  fit <- relimp_network(data = items)
  fit$weights <- shares
  # Sums of six numbers of order 0.1: relative 1e-10.
  expect_reference(lmg_certificate(fit), gap, reference = "relaimpo", covers = "lmg_certificate()",
                   quantity = "certificate of relaimpo's first-entry shares against its R-squared", tolerance = 1e-10,
                   measure = "relative",
                   note = "Shares replaced through `$weights<-`: lmg_certificate() takes only a fitted network (API gap).")
})
