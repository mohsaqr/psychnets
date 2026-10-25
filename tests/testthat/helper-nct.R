# Shared data generator and configurations for the network comparison test.
# fixtures/make-nct.R freezes net_compare() output and the
# NetworkComparisonTest::NCT() reference on these, and test-nct.R checks
# against them, so the generator and the checks can never drift apart.
nct_groups <- function(seed, n = 120, p = 5) {
  set.seed(seed)
  a <- matrix(stats::rnorm(n * p), n, p)
  a <- a + 0.4 * matrix(rowSums(a), n, p)
  b <- matrix(stats::rnorm(n * p), n, p)
  b <- b + 0.15 * matrix(rowSums(b), n, p)
  b[, 2] <- b[, 2] + 0.5 * b[, 1]
  colnames(a) <- colnames(b) <- paste0("V", seq_len(p))
  list(a = a, b = b)
}

nct_configs <- data.frame(
  seed     = c(1, 2, 3, 4),
  paired   = c(FALSE, TRUE, FALSE, FALSE),
  abs      = c(TRUE, TRUE, FALSE, TRUE),
  weighted = c(TRUE, TRUE, TRUE, FALSE)
)
nct_iter <- 15L

# The components of a net_compare() result that carry numbers: M (network
# structure), S (global strength) and E (edges), in NCT notation.
nct_numbers <- function(fit) {
  list(nw1 = fit$nw1, nw2 = fit$nw2,
       structure = fit$structure[c("observed", "perm", "p_value")],
       strength = fit$strength[c("observed", "perm", "p_value")],
       E = fit$E[c("observed", "perm", "p_value")])
}
