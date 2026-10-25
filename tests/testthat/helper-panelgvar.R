# Fixed-seed panel-data generator shared by test-panelgvar.R, the golden /
# equivalence tests, and tests/testthat/fixtures/make-panelgvar.R, so the frozen
# fixtures and the tests that check them can never drift apart.
#
# Stationary random-intercept VAR(1): y_it = mu + b_i + w_it,
# w_t = B w_{t-1} + z_t, with w_1 drawn from the stationary distribution.
# Returns long data (id, wave, V1..Vp). The RNG state is restored on exit.
# `trait_scale` multiplies each variable's stable trait AFTER all draws (so the
# random stream, and every fixture built from it, is unchanged); 0 removes a
# variable's between-person variance.
panel_sim <- function(n = 400, p = 3, waves = 4, seed = 1, trait_scale = 1) {
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  b <- diag(0.4, p)
  b[2, 1] <- 0.25
  b[3, 2] <- -0.2
  sz <- diag(p)
  sz[1, 2] <- sz[2, 1] <- 0.3
  sb <- diag(0.6, p)
  sb[1, 3] <- sb[3, 1] <- 0.25
  s0 <- matrix(solve(diag(p * p) - b %x% b, as.vector(sz)), p, p)
  trait <- matrix(stats::rnorm(n * p), n) %*% chol(sb)
  w <- matrix(stats::rnorm(n * p), n) %*% chol(s0)
  states <- Reduce(function(prev, t)
    prev %*% t(b) + matrix(stats::rnorm(n * p), n) %*% chol(sz),
    seq_len(waves - 1L), w, accumulate = TRUE)
  trait <- sweep(trait, 2L, rep_len(trait_scale, p), "*")
  rows <- lapply(seq_len(waves), function(t) {
    y <- sweep(states[[t]] + trait, 2L, seq_len(p), "+")
    colnames(y) <- paste0("V", seq_len(p))
    data.frame(id = seq_len(n), wave = t, y)
  })
  do.call(rbind, rows)
}

# Wide form of panel_sim() plus the psychonetrics-style design matrix.
panel_sim_wide <- function(...) {
  long <- panel_sim(...)
  vars <- setdiff(names(long), c("id", "wave"))
  waves <- sort(unique(long$wave))
  wide <- Reduce(function(acc, t) {
    part <- long[long$wave == t, c("id", vars)]
    names(part) <- c("id", paste0(vars, "_w", t))
    merge(acc, part, by = "id")
  }, waves[-1L], {
    part <- long[long$wave == waves[1L], c("id", vars)]
    names(part) <- c("id", paste0(vars, "_w", waves[1L]))
    part
  })
  design <- matrix(paste0(rep(vars, length(waves)), "_w",
                          rep(waves, each = length(vars))),
                   nrow = length(vars), dimnames = list(vars, NULL))
  list(data = wide, design = design)
}
