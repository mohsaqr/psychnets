# The comparisons with reference packages that used to live here (qgraph, glasso, IsingFit, mgm,
# networktools, igraph, psychonetrics) are declared equivalence tests in tests/equivalence/ (run by the
# `verify` job). What stays is the check of an internal function, which the declared tests cannot call:
# the Hittner et al. (2003) test that redundancy() applies to each pair of correlations, against cocor.
# Gated: runs only when PSYCHNET_EQUIV_TESTS is set and cocor is installed.

test_that("Hittner2003 test matches cocor", {
  skip_if(Sys.getenv("PSYCHNET_EQUIV_TESTS") == "")
  skip_if_not_installed("cocor")
  set.seed(3)
  invisible(lapply(seq_len(5L), function(.) {
    r <- stats::runif(3, -0.7, 0.7); n <- sample(80:400, 1)
    ours <- .psn_hittner2003(r[1], r[2], r[3], n)
    ref  <- suppressWarnings(cocor::cocor.dep.groups.overlap(
      r[1], r[2], r[3], n, test = "hittner2003")@hittner2003$p.value)
    expect_lt(abs(ours - ref), 1e-10)
  }))
})
