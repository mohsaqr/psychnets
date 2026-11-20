# ebic_glasso(native = FALSE): psychnets calls glasso::glasso() itself at the EBIC-selected penalty, so the
# agreement checks the hand-over, not a second algorithm, and is declared internal (it counts as no
# independent evidence). Kept in its own file because one case gives one reference one authority.

items <- utils::read.csv(reference_data("data/bfi10.csv"))
correlations <- stats::cor(items)

test_that("native = FALSE hands the fit to glasso::glasso() (internal: psychnets calls glasso)", {
  skip_without_reference("glasso")
  fit <- ebic_glasso(data = items, native = FALSE)
  reference <- glasso::glasso(correlations, rho = fit$lambda, penalize.diagonal = FALSE)$wi
  dimnames(reference) <- dimnames(correlations)
  expect_reference(fit$precision, (reference + t(reference)) / 2, reference = "glasso", authority = "internal",
                   covers = "ebic_glasso()", quantity = "precision, native = FALSE, at the selected penalty",
                   tolerance = 1e-12, note = "With native = FALSE psychnets calls glasso::glasso(); this checks the hand-over.")
})
