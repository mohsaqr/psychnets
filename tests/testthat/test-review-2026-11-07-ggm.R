# Regression tests for the 2026-11-07 review of the GGM / structural / front-door
# estimators (findings G1-G10, S1-S15). Each block names the finding it pins.

sim_gauss <- function(n = 300, p = 5, seed = 1) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p) %*% chol(0.4^abs(outer(seq_len(p), seq_len(p), "-")))
  colnames(x) <- paste0("V", seq_len(p))
  x
}

# ---- G1: sample size for significance tests -----------------------------------

test_that("G1: an invalid n is refused, never stored", {
  S <- stats::cor(sim_gauss())
  expect_error(cor_network(cor_matrix = S, n = -3, alpha = 0.05),
               class = "psychnet_bad_input")
  expect_error(cor_network(cor_matrix = S, n = 2, alpha = 0.05),
               class = "psychnet_bad_input")
  expect_error(cor_network(cor_matrix = S, n = 50.5, alpha = 0.05),
               class = "psychnet_bad_input")
  expect_error(cor_network(cor_matrix = S, n = -5), class = "psychnet_bad_input")
  expect_error(pcor_network(cor_matrix = S, n = NA, alpha = 0.05),
               class = "psychnet_bad_input")
  expect_error(ebic_glasso(cor_matrix = S, n = -5), class = "psychnet_bad_input")
  expect_error(huge_network(cor_matrix = S, n = 0), class = "psychnet_bad_input")
  expect_error(ggm_modselect(cor_matrix = S, n = "a"), class = "psychnet_bad_input")
  expect_error(logo_network(cor_matrix = S, n = -1), class = "psychnet_bad_input")
})

test_that("G1: df <= 0 is an error, not an all-p = 1 empty network", {
  S <- stats::cor(sim_gauss(p = 5))
  # full-order partials of 5 variables need n > 2 + 3 = 5
  expect_error(pcor_network(cor_matrix = S, n = 5, alpha = 0.05),
               class = "psychnet_degenerate_sample")
  expect_error(pcor_network(cor_matrix = S, n = 5, alpha = 0.05),
               class = "psychnet_bad_input")
  # the smallest valid n still tests (df = 1)
  fit <- pcor_network(cor_matrix = S, n = 6, alpha = 0.5)
  expect_true(all(is.finite(fit$p_values)))
  # pairwise data with a pair co-observed only twice
  x <- sim_gauss(n = 30, p = 3)
  x[3:30, 1] <- NA
  expect_error(suppressWarnings(cor_network(x, alpha = 0.05)),
               class = "psychnet_degenerate_sample")
})

# ---- G2: singular correlation in pcor_network -----------------------------------

test_that("G2: n <= p partial correlations warn psychnet_singular", {
  set.seed(1)
  x <- matrix(stats::rnorm(40), 5, 8)
  expect_warning(fit <- pcor_network(x), class = "psychnet_singular")
  expect_s3_class(fit, "psychnet")                    # result still returned
  # an exactly duplicated column is singular too
  y <- sim_gauss()
  expect_warning(pcor_network(cbind(y, dup = y[, 1])), class = "psychnet_singular")
  # a well-conditioned matrix is silent and unchanged
  S <- stats::cor(y)
  expect_no_warning(fit2 <- pcor_network(cor_matrix = S))
  expect_equal(fit2$precision, solve(S), tolerance = 1e-12,
               ignore_attr = TRUE)
})

# ---- G3: ordered factors, dropped non-numeric columns ---------------------------

test_that("G3: ordered factors enter as ordinal codes, identical to integer input", {
  set.seed(2)
  z <- sim_gauss(n = 300, p = 4)
  lik <- apply(z, 2, function(col) as.integer(cut(col, 5)))
  df_int <- as.data.frame(lik)
  df_ord <- df_int
  df_ord$V2 <- factor(df_ord$V2, levels = 1:5, ordered = TRUE)
  expect_identical(cor_auto(df_ord), cor_auto(df_int))
  expect_identical(ebic_glasso(df_ord, cor_method = "auto")$weights,
                   ebic_glasso(df_int, cor_method = "auto")$weights)
  expect_identical(colnames(cor_auto(df_ord)), colnames(df_int))
})

test_that("G3: unordered factor / character / logical columns drop with a classed message", {
  x <- as.data.frame(sim_gauss(p = 4))
  x$grp <- factor(rep(c("a", "b", "c"), 100))
  x$id <- as.character(seq_len(300))
  x$flag <- rep(c(TRUE, FALSE), 150)
  expect_message(fit <- cor_network(x), class = "psychnet_dropped_columns")
  expect_identical(fit$nodes$label, paste0("V", 1:4))
  msg <- tryCatch(cor_network(x), message = function(m) conditionMessage(m))
  expect_match(msg, "grp")
  expect_match(msg, "id")
  expect_match(msg, "flag")
  # an all-numeric frame is silent
  expect_silent(cor_network(as.data.frame(sim_gauss(p = 4))))
})

# ---- G4: constant columns + labels --------------------------------------------

test_that("G4: constant / all-NA columns warn psychnet_constant_variable by name", {
  x <- cbind(sim_gauss(p = 4), konst = 1, empty = NA_real_)
  w <- tryCatch(cor_network(x), warning = function(w) w)
  expect_s3_class(w, "psychnet_constant_variable")
  expect_match(conditionMessage(w), "konst")
  expect_match(conditionMessage(w), "empty")
  expect_warning(ebic_glasso(x), class = "psychnet_constant_variable")
})

test_that("G4: labels of the wrong length give the classed hint in every verb", {
  x <- cbind(sim_gauss(p = 5), konst = 1)
  labs <- paste0("L", 1:6)
  verbs <- list(cor_network, pcor_network, ebic_glasso, huge_network,
                ggm_modselect, logo_network, tmfg_network, relimp_network)
  for (f in verbs) {
    expect_error(suppressWarnings(f(x, labels = labs)),
                 class = "psychnet_bad_input")
    expect_error(suppressWarnings(f(x, labels = labs)), "may have been dropped")
  }
  # cor_matrix input checks labels too
  S <- stats::cor(sim_gauss(p = 5))
  expect_error(ebic_glasso(cor_matrix = S, n = 300, labels = letters[1:4]),
               class = "psychnet_bad_input")
  expect_error(tmfg_network(cor_matrix = S, labels = rep("a", 5)),
               class = "psychnet_bad_input")
  # correct labels still work
  expect_identical(ebic_glasso(cor_matrix = S, n = 300,
                               labels = letters[1:5])$nodes$label, letters[1:5])
})

test_that("G4: .psn_with_quiet_columns mutes the column notices and restores state", {
  x <- cbind(sim_gauss(p = 4), konst = 1)
  expect_no_warning(.psn_with_quiet_columns(.as_numeric_matrix(x)))
  expect_false(.psn_state$quiet_columns)
  expect_warning(.as_numeric_matrix(x), class = "psychnet_constant_variable")
  # restored on error too
  try(.psn_with_quiet_columns(stop("boom")), silent = TRUE)
  expect_false(.psn_state$quiet_columns)
})

test_that("G3/G4: numeric input gives byte-identical matrices", {
  x <- sim_gauss()
  expect_identical(.as_numeric_matrix(x), {
    m <- x; storage.mode(m) <- "double"; m
  })
  expect_identical(.as_numeric_matrix(as.data.frame(x)), as.matrix(as.data.frame(x)))
})

# ---- G5: ggm_modselect argument checks ------------------------------------------

test_that("G5: ggm_modselect validates gamma, nlambda, lambda_min_ratio, stepwise", {
  x <- sim_gauss()
  expect_error(ggm_modselect(x, gamma = -1), class = "psychnet_bad_input")
  expect_error(ggm_modselect(x, gamma = NA), class = "psychnet_bad_input")
  expect_error(ggm_modselect(x, nlambda = 1), class = "psychnet_bad_input")
  expect_error(ggm_modselect(x, nlambda = 10.5), class = "psychnet_bad_input")
  expect_error(ggm_modselect(x, lambda_min_ratio = 0), class = "psychnet_bad_input")
  expect_error(ggm_modselect(x, lambda_min_ratio = 1), class = "psychnet_bad_input")
  expect_error(ggm_modselect(x, stepwise = NA), class = "psychnet_bad_input")
  expect_error(ggm_modselect(x, stepwise = "yes"), class = "psychnet_bad_input")
})

# ---- G6 / S1 / S12: threshold and alpha ---------------------------------------

test_that("G6: threshold must be one finite number >= 0 in every verb", {
  x <- sim_gauss()
  verbs <- list(cor_network, pcor_network, ebic_glasso, huge_network,
                ggm_modselect, logo_network)
  bad <- list(TRUE, "LW", "a", NA, c(0.1, 0.2), -0.1, Inf)
  for (f in verbs) for (b in bad) {
    expect_error(f(x, threshold = b), class = "psychnet_bad_input")
  }
  expect_error(ebic_glasso(x, threshold = TRUE), "qgraph")
  expect_error(psychnet(x, "logo", threshold = TRUE), class = "psychnet_bad_input")
  expect_error(psychnet(x, "glasso", threshold = "LW"), class = "psychnet_bad_input")
  # a valid numeric threshold still thresholds
  fit <- cor_network(x, threshold = 0.2)
  expect_true(all(abs(as.data.frame(fit)$weight) >= 0.2))
})

test_that("G6: alpha is NULL or in (0, 1)", {
  x <- sim_gauss()
  for (a in list(NA, 5, 0, 1, "0.05", c(0.01, 0.05)))
    expect_error(cor_network(x, alpha = a), class = "psychnet_bad_input")
  expect_error(pcor_network(x, alpha = NA), class = "psychnet_bad_input")
  expect_s3_class(cor_network(x, alpha = 0.05), "psychnet")
})

# ---- G8: certificate inputs ---------------------------------------------------

test_that("G8: glasso_kkt refuses a negative / missing rho", {
  S <- 0.5^abs(outer(1:5, 1:5, "-"))
  expect_error(glasso_kkt(solve(S), S, rho = -1), class = "psychnet_bad_input")
  expect_error(glasso_kkt(solve(S), S, rho = NA), class = "psychnet_bad_input")
  expect_error(glasso_kkt(solve(S), S, rho = c(0, 1)), class = "psychnet_bad_input")
  expect_lt(glasso_kkt(solve(S), S, rho = 0), 1e-10)
})

test_that("G8: ggm_support_kkt symmetrizes a triangular support", {
  S <- 0.4^abs(outer(1:6, 1:6, "-"))
  chain <- abs(row(S) - col(S)) <= 1
  up <- chain & upper.tri(chain)
  lo <- chain & lower.tri(chain)
  ref <- ggm_support_kkt(solve(S), S, support = chain)
  expect_identical(ggm_support_kkt(solve(S), S, support = up), ref)
  expect_identical(ggm_support_kkt(solve(S), S, support = lo), ref)
  expect_lt(ref, 1e-10)
  # the full inverse is NOT a chain-restricted MLE of a denser S
  S2 <- stats::cor(sim_gauss(p = 6))
  expect_gt(ggm_support_kkt(solve(S2), S2, support = lo), 1e-3)
})

test_that("G8: ggm_support_kkt refuses malformed supports", {
  S <- 0.4^abs(outer(1:4, 1:4, "-"))
  sup <- abs(row(S) - col(S)) <= 1
  expect_error(ggm_support_kkt(solve(S), S, support = sup * 1),
               class = "psychnet_bad_input")
  expect_error(ggm_support_kkt(solve(S), S, support = sup[1:3, 1:3]),
               class = "psychnet_bad_input")
  sup_na <- sup; sup_na[1, 2] <- NA
  expect_error(ggm_support_kkt(solve(S), S, support = sup_na),
               class = "psychnet_bad_input")
})

test_that("G8: certificate() validates tol", {
  fit <- ebic_glasso(sim_gauss())
  for (t in list("a", NA, -1, c(1e-6, 1e-3)))
    expect_error(certificate(fit, tol = t), class = "psychnet_bad_input")
  expect_error(certificate(list()), class = "psychnet_bad_input")
  expect_true(certificate(fit, tol = 1e-6)$certified)
})

# ---- G9: cor_auto ordinal_max_levels --------------------------------------------

test_that("G9: cor_auto validates ordinal_max_levels", {
  x <- sim_gauss()
  for (v in list("a", NA, 1, 2.5, c(3, 4)))
    expect_error(cor_auto(x, ordinal_max_levels = v), class = "psychnet_bad_input")
  expect_true(is.matrix(cor_auto(x, ordinal_max_levels = 2)))
})

# ---- G10: ggm_modselect failed fits ---------------------------------------------

test_that("G10: a failed candidate fit is reported as psychnet_fit_failed", {
  x <- sim_gauss()
  real_fit <- .ggm_fit_support
  calls <- 0L
  local_mocked_bindings(.ggm_fit_support = function(S, support, ...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("Lapack routine dgesv: system is exactly singular")
    real_fit(S, support, ...)
  })
  w <- tryCatch(ggm_modselect(x), psychnet_fit_failed = function(w) w)
  expect_s3_class(w, "psychnet_fit_failed")
  expect_match(conditionMessage(w), "1 candidate")
  expect_match(conditionMessage(w), "singular")
})

test_that("G10: all candidate fits failing is a classed error", {
  x <- sim_gauss()
  local_mocked_bindings(.ggm_fit_support = function(S, support, ...)
    stop("system is exactly singular"))
  expect_error(ggm_modselect(x, stepwise = FALSE), class = "psychnet_no_converge")
})

test_that("G10: an uncertified final fit warns psychnet_no_converge", {
  x <- sim_gauss()
  local_mocked_bindings(ggm_support_kkt = function(...) 0.5)
  expect_warning(ggm_modselect(x), class = "psychnet_no_converge")
})

test_that("G10: a normal ggm_modselect fit is warning-free and certified", {
  x <- sim_gauss()
  expect_no_warning(fit <- ggm_modselect(x))
  expect_lt(fit$kkt, 1e-6)
})

# ---- S2 / S13: dichotomize -------------------------------------------------------

test_that("S2: a median split that leaves one group warns psychnet_constant_variable", {
  x <- cbind(floor = c(0, 0, 0, 0, 0, 1, 2), ok = 1:7)
  expect_warning(b <- dichotomize(x, method = "median"),
                 class = "psychnet_constant_variable")
  expect_true(all(b[, "floor"] == 1L))
  msg <- tryCatch(dichotomize(x), warning = function(w) conditionMessage(w))
  expect_match(msg, "floor")
  expect_false(grepl("ok", sub("^.*: ", "", sub("\\. A constant.*$", "", msg))))
  # mean and rank split it
  expect_no_warning(bm <- dichotomize(x, method = "mean"))
  expect_setequal(unique(bm[, "floor"]), c(0L, 1L))
  expect_no_warning(br <- dichotomize(x, method = "rank"))
  expect_setequal(unique(br[, "floor"]), c(0L, 1L))
})

test_that("S13: dichotomize keeps NA as NA and validates input with a class", {
  x <- cbind(a = c(1, 2, NA, 4, 5), b = c(5, NA, 3, 2, 1))
  for (m in c("median", "mean", "rank")) {
    b <- dichotomize(x, method = m)
    expect_identical(is.na(b), is.na(x))
    expect_true(all(b[!is.na(b)] %in% c(0L, 1L)))
  }
  expect_error(dichotomize(matrix(letters[1:6], 3)), class = "psychnet_bad_input")
  expect_error(dichotomize(cbind(a = c(1, 1, 1), b = 1:3)),
               class = "psychnet_bad_input")
})

# ---- S3 / S11: relimp_network ------------------------------------------------------

test_that("S3: relimp_network accepts and ignores n", {
  S <- 0.4^abs(outer(1:5, 1:5, "-"))
  a <- relimp_network(cor_matrix = S, n = 300)
  b <- relimp_network(cor_matrix = S)
  expect_identical(a$weights, b$weights)
  expect_s3_class(psychnet(sim_gauss(), "relimp"), "psychnet")
})

test_that("S11: relimp_network validates max_nodes", {
  S <- 0.4^abs(outer(1:5, 1:5, "-"))
  for (v in list(NA, "a", 1, 3.5))
    expect_error(relimp_network(cor_matrix = S, max_nodes = v),
                 class = "psychnet_bad_input")
  expect_error(relimp_network(cor_matrix = S, max_nodes = 4),
               class = "psychnet_bad_input")      # p = 5 > 4
})

# ---- S4: crosswalk refit row ----------------------------------------------------

test_that("S4: the EBICglasso refit row maps qgraph TRUE to 'unregularized'", {
  cw <- net_crosswalk("EBICglasso")
  row <- subset(cw, ref_arg == "refit")
  expect_equal(nrow(row), 1L)
  expect_match(row$note, "unregularized", fixed = TRUE)
  expect_match(row$note, "qgraph refit = TRUE", fixed = TRUE)
})

# ---- S5: lmg_certificate ---------------------------------------------------------

test_that("S5: lmg certificate anchors the full-model R^2 to the closed form", {
  S <- 0.4^abs(outer(1:5, 1:5, "-"))
  fit <- relimp_network(cor_matrix = S)
  expect_lt(lmg_certificate(fit), 1e-12)
  expect_equal(unname(fit$r2), 1 - 1 / diag(solve(S)), tolerance = 1e-12)
})

test_that("S5 mutation: a subset solver off by a factor fails the certificate", {
  S <- 0.4^abs(outer(1:5, 1:5, "-"))
  real_node <- .lmg_node
  local_mocked_bindings(.lmg_node = function(S, j) {
    res <- real_node(S, j)        # every subset R^2 halved => shares halved too
    res$lmg <- 0.5 * res$lmg
    res$r2_full <- 0.5 * res$r2_full
    res
  })
  bad <- relimp_network(cor_matrix = S)
  # Shapley efficiency still holds exactly (it telescopes) ...
  expect_lt(max(abs(colSums(bad$weights) - bad$r2)), 1e-12)
  # ... but the closed-form anchor exposes the wrong solver
  expect_gt(lmg_certificate(bad), 0.05)
  expect_gt(bad$kkt, 0.05)
  expect_false(certificate(bad)$certified)
})

# ---- S7: LoGo near-duplicate column ------------------------------------------------

test_that("S7: LoGo refuses a near-singular clique block with psychnet_singular", {
  x <- sim_gauss(n = 300, p = 5)
  set.seed(3)
  x_dup <- cbind(x, near = x[, 1] + stats::rnorm(300, sd = 1e-7))
  expect_error(logo_network(x_dup), class = "psychnet_singular")
  expect_error(logo_network(cbind(x, dup = x[, 1])), class = "psychnet_singular")
  fit <- logo_network(x)
  expect_lt(fit$kkt, 1e-10)
})

# ---- S8 / S9: the front door -----------------------------------------------------

test_that("S8: an explicit threshold / gamma for a method that ignores it warns", {
  x <- sim_gauss()
  expect_warning(psychnet(x, "tmfg", threshold = 0.3),
                 class = "psychnet_ignored_argument")
  expect_warning(psychnet(x, "relimp", gamma = 0.5),
                 class = "psychnet_ignored_argument")
  expect_warning(psychnet(x, "cor", gamma = 0.5),
                 class = "psychnet_ignored_argument")
  b <- dichotomize(x)
  expect_warning(psychnet(b, "ising", threshold = 0.1),
                 class = "psychnet_ignored_argument")
  # defaults and accepted arguments are silent
  expect_no_warning(psychnet(x, "tmfg"))
  expect_no_warning(psychnet(x, "tmfg", threshold = 0))
  expect_no_warning(psychnet(x, "glasso", threshold = 0.05, gamma = 0.25))
  expect_no_warning(psychnet(b, "ising", gamma = 0.1))
})

test_that("S8: group mode reports an ignored argument once, not per group", {
  d <- data.frame(g = rep(c("A", "B", "C"), each = 100), sim_gauss())
  n_warn <- 0L
  withCallingHandlers(
    psychnet(d, method = "tmfg", group = "g", threshold = 0.3),
    psychnet_ignored_argument = function(w) {
      n_warn <<- n_warn + 1L
      invokeRestart("muffleWarning")
    })
  expect_identical(n_warn, 1L)
  expect_no_warning(psychnet(d, method = "tmfg", group = "g"))
})

test_that("S8: resampling verbs with default arguments stay warning-free", {
  skip_on_cran()
  x <- sim_gauss(n = 200, p = 5)
  expect_no_warning(net_boot(x, method = "tmfg", n_boot = 5, cores = 1))
  expect_no_warning(net_boot(x, method = "relimp", n_boot = 5, cores = 1))
})

test_that("S9: several methods at once are an error; the untouched default is fine", {
  x <- sim_gauss()
  expect_error(psychnet(x, method = c("tmfg", "logo")), class = "psychnet_bad_input")
  expect_error(psychnet(x, method = character()), class = "psychnet_bad_input")
  expect_identical(psychnet(x)$method, "glasso")
})

# ---- S14 / S15: TMFG zero edges and certificate labels ------------------------

test_that("S14: a retained TMFG edge of correlation 0 is in $adjacency, not $edges", {
  S <- diag(6)
  S[1:3, 1:3] <- 0.5; S[4:6, 4:6] <- 0.4; diag(S) <- 1
  fit <- tmfg_network(cor_matrix = S)
  expect_equal(sum(fit$adjacency[upper.tri(fit$adjacency)]), 3 * (6 - 2))
  expect_lt(nrow(as.data.frame(fit)), 3 * (6 - 2))
  expect_identical(tmfg_certificate(fit), 0L)
})

test_that("S15: print labels structural certificates as structural, not KKT", {
  S <- 0.4^abs(outer(1:6, 1:6, "-"))
  out_r <- capture.output(print(relimp_network(cor_matrix = S)))
  out_t <- capture.output(print(tmfg_network(cor_matrix = S)))
  out_g <- capture.output(print(ebic_glasso(cor_matrix = S, n = 250)))
  expect_false(any(grepl("KKT", out_r)))
  expect_true(any(grepl("structural certificate", out_r)))
  expect_false(any(grepl("KKT", out_t)))
  expect_true(any(grepl("violated TMFG invariants): 0", out_t, fixed = TRUE)))
  expect_true(any(grepl("optimality (KKT residual)", out_g, fixed = TRUE)))
  expect_identical(certificate(relimp_network(cor_matrix = S))$kind, "structural")
})

# ---- invariants ------------------------------------------------------------------

test_that("invariant: column permutation permutes the validated networks", {
  x <- sim_gauss(p = 5)
  perm <- c(3, 1, 5, 2, 4)
  a <- pcor_network(x)$weights
  b <- pcor_network(x[, perm])$weights
  expect_equal(b, a[perm, perm], tolerance = 1e-12)
  ra <- relimp_network(x)
  rb <- relimp_network(x[, perm])
  expect_equal(rb$weights, ra$weights[perm, perm], tolerance = 1e-12)
  expect_lt(lmg_certificate(rb), 1e-12)
})
