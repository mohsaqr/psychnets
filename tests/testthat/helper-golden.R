# Fixed-seed data generators shared by test-golden.R and the fixture generator
# (tests/testthat/fixtures/make-golden.R). Kept in one place so the frozen
# fixtures and the tests that check them can never drift apart.

golden_data_mixed <- function(n = 400) {
  set.seed(20260807)
  f <- stats::rnorm(n)
  data.frame(
    g1 = f + stats::rnorm(n, sd = 0.6),
    g2 = f + stats::rnorm(n, sd = 0.6),
    g3 = stats::rnorm(n),
    b1 = ((f + stats::rnorm(n, sd = 0.6)) > 0) * 1L,
    b2 = ((f + stats::rnorm(n, sd = 0.9)) > 0) * 1L
  )
}

golden_data_binary <- function(n = 400) {
  d <- golden_data_mixed(n)
  data.frame(
    b1 = (d$g1 > 0) * 1L, b2 = (d$g2 > 0) * 1L,
    b3 = (d$g3 > 0) * 1L, b4 = d$b1, b5 = d$b2
  )
}

golden_cor <- function() 0.45^abs(outer(1:8, 1:8, "-"))
