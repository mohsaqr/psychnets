# Panel graphical VAR (panel-GVAR): the stationary random-intercept VAR(1) for
# few-wave, many-person panel data (Epskamp 2020). Each person's scores are a
# stable trait (between-person part) plus a within-person deviation that follows
# a stationary lag-1 VAR. The stacked (wave x variable) covariance is therefore
# block Toeplitz:
#
#   Sigma[t, s] = B^(t - s) Sigma0 + Sigma_between      (t >= s)
#   vec(Sigma0) = (I - B %x% B)^-1 vec(Sigma_zeta)      (stationary within cov)
#
# with Sigma_zeta (contemporaneous innovations) and Sigma_between both
# parameterized as Gaussian graphical models, Sigma = D (I - Omega)^-1 D, so the
# off-diagonal of Omega is directly the partial-correlation network. Means are
# stationary (one mean per variable, shared by every wave).
#
# Estimation is maximum likelihood on the sample moments (covariance with the
# ML denominator n), the same fit function as psychonetrics::panelgvar():
#
#   F = tr(S K) + (m - mu)' K (m - mu) - log det K,     K = Sigma^-1
#
# Certification: the ML optimum is a stationary point, so the score (gradient of
# F) on the free parameters is ~0. The optimizer uses a reverse-mode gradient;
# the certificate recomputes the score with a forward-mode Jacobian -- an
# independent code path -- and stores the worst absolute value as `kkt`.

# ---- parameter layout --------------------------------------------------------

# The full parameter vector, in psychonetrics' order: mu, omega_zeta_within
# (lower triangle, column-major), delta_zeta_within, beta (column-major,
# beta[i, j] = effect of j at t-1 on i at t), omega_zeta_between,
# delta_zeta_between.
#' @noRd
.pgv_layout <- function(p) {
  m <- p * (p - 1L) / 2L
  sizes <- c(mu = p, omega_zeta_within = m, delta_zeta_within = p,
             beta = p * p, omega_zeta_between = m, delta_zeta_between = p)
  ends <- cumsum(sizes)
  idx <- Map(function(s, e) seq.int(to = e, length.out = s), sizes, ends)
  lower <- which(lower.tri(diag(p)), arr.ind = TRUE)
  diag_rc <- cbind(seq_len(p), seq_len(p))
  full_rc <- which(matrix(TRUE, p, p), arr.ind = TRUE)
  rc <- rbind(cbind(seq_len(p), 1L), lower, diag_rc, full_rc, lower, diag_rc)
  table <- data.frame(matrix = rep(names(sizes), sizes),
                      row = as.integer(rc[, 1L]), col = as.integer(rc[, 2L]),
                      stringsAsFactors = FALSE)
  list(p = p, idx = idx, n_par = sum(sizes), lower = lower, table = table)
}

#' @noRd
.pgv_unpack <- function(theta, lay) {
  p <- lay$p
  ix <- lay$idx
  sym <- function(v) {
    out <- matrix(0, p, p)
    out[lay$lower] <- v
    out + t(out)
  }
  list(mu = theta[ix$mu],
       omega_w = sym(theta[ix$omega_zeta_within]),
       delta_w = theta[ix$delta_zeta_within],
       beta = matrix(theta[ix$beta], p, p),
       omega_b = sym(theta[ix$omega_zeta_between]),
       delta_b = theta[ix$delta_zeta_between])
}

# ---- implied moments ---------------------------------------------------------

# GGM covariance D (I - Omega)^-1 D; NULL when I - Omega is not positive
# definite (an infeasible point the optimizer must back away from).
#' @noRd
.pgv_ggm_cov <- function(omega, delta) {
  a <- diag(length(delta)) - omega
  if (min(eigen(a, symmetric = TRUE, only.values = TRUE)$values) <= 1e-10)
    return(NULL)
  w <- chol2inv(chol(a))
  list(sigma = outer(delta, delta) * w, w = w)
}

# Stack lag blocks into the block-Toeplitz (waves*p) square: block (i, j) is
# lags[[i - j + 1]] below the diagonal and its transpose above, plus `between`
# on every block.
#' @noRd
.pgv_assemble <- function(lags, between, waves) {
  block <- function(i, j) {
    if (i >= j) lags[[i - j + 1L]] + between else t(lags[[j - i + 1L]]) + between
  }
  do.call(rbind, lapply(seq_len(waves), function(i)
    do.call(cbind, lapply(seq_len(waves), function(j) block(i, j)))))
}

# Implied moments at `par`; NULL at an infeasible point (non-stationary B, or a
# non-positive-definite innovation / between / stacked covariance).
#' @noRd
.pgv_moments <- function(par, waves) {
  p <- length(par$mu)
  b <- par$beta
  if (any(par$delta_w <= 0) || any(par$delta_b <= 0)) return(NULL)
  if (max(Mod(eigen(b, only.values = TRUE)$values)) >= 1 - 1e-10) return(NULL)
  zw <- .pgv_ggm_cov(par$omega_w, par$delta_w)
  zb <- .pgv_ggm_cov(par$omega_b, par$delta_b)
  if (is.null(zw) || is.null(zb)) return(NULL)
  lyap <- diag(p * p) - b %x% b
  s0 <- matrix(solve(lyap, as.vector(zw$sigma)), p, p)
  s0 <- (s0 + t(s0)) / 2
  lags <- Reduce(function(a, k) b %*% a, seq_len(waves - 1L), s0,
                 accumulate = TRUE)
  sigma <- .pgv_assemble(lags, zb$sigma, waves)
  sigma <- (sigma + t(sigma)) / 2
  list(sigma = sigma, mu = rep(par$mu, waves), lags = lags, s0 = s0,
       zw = zw, zb = zb, lyap = lyap)
}

# Everything the objective, gradient, and information need at one parameter
# vector; NULL at an infeasible point.
#' @noRd
.pgv_state <- function(theta, ctx) {
  par <- .pgv_unpack(theta, ctx$lay)
  mo <- .pgv_moments(par, ctx$waves)
  if (is.null(mo)) return(NULL)
  # Relative, not absolute, floor: a covariance conditioned beyond 1e12 makes K
  # overflow and the gradient NaN, so it is treated as outside the domain.
  ev <- eigen(mo$sigma, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(ev)) return(NULL)
  ch <- chol(mo$sigma)
  k <- chol2inv(ch)
  r <- ctx$means - mo$mu
  f <- sum(ctx$cov * k) + sum(r * (k %*% r)) + 2 * sum(log(diag(ch)))
  c(mo, list(par = par, k = k, r = r, f = f))
}

# ---- objective and reverse-mode gradient -------------------------------------

#' @noRd
.pgv_full <- function(theta_free, ctx) {
  theta <- ctx$theta
  theta[ctx$free] <- theta_free
  theta
}

# dF/dSigma (as a symmetric matrix G with dF = tr(G dSigma) - 2 (K r)' dmu).
#' @noRd
.pgv_dsigma <- function(st, ctx) {
  kr <- st$k %*% st$r
  st$k - st$k %*% ctx$cov %*% st$k - tcrossprod(kr)
}

# Gradient of a GGM covariance D W D (W = (I - Omega)^-1) w.r.t. its omega
# (lower triangle) and delta, given R with dF = tr(R dSigma).
#' @noRd
.pgv_ggm_grad <- function(r_mat, delta, w, lower) {
  z <- w %*% (outer(delta, delta) * r_mat) %*% w
  list(omega = 2 * z[lower], delta = 2 * diag(w %*% (delta * r_mat)))
}

# Reverse-mode gradient of F w.r.t. the FULL parameter vector.
#' @noRd
.pgv_gradient_full <- function(st, ctx) {
  p <- ctx$lay$p
  waves <- ctx$waves
  b <- st$par$beta
  g <- .pgv_dsigma(st, ctx)
  rows <- function(i) (i - 1L) * p + seq_len(p)
  blk <- function(i, j) g[rows(i), rows(j), drop = FALSE]
  g_mu <- -2 * rowSums(matrix(st$k %*% st$r, p, waves))
  # Between part: every block carries Sigma_between, so its gradient is the sum
  # of all blocks of G (symmetric).
  pairs <- expand.grid(i = seq_len(waves), j = seq_len(waves))
  g_between <- Reduce(`+`, Map(blk, pairs$i, pairs$j))
  # Lag-k block A_k = B^k Sigma0 appears below the diagonal (and transposed
  # above), so dF = tr(C_0 dA_0) + 2 sum_k tr(C_k dA_k), C_k = sum_j G[j, j+k].
  # Store the gradient w.r.t. each A_k, i.e. t(multiplier * C_k).
  ga <- lapply(seq_len(waves) - 1L, function(k) {
    ck <- Reduce(`+`, lapply(seq_len(waves - k), function(j) blk(j, j + k)))
    t(if (k == 0L) ck else 2 * ck)
  })
  g_beta <- matrix(0, p, p)
  # Back-propagate through A_k = B A_{k-1}. Each step needs the adjoint the
  # previous (higher-lag) step produced -- a genuinely sequential recursion.
  for (k in rev(seq_len(waves - 1L))) {
    g_beta <- g_beta + ga[[k + 1L]] %*% t(st$lags[[k]])
    ga[[k]] <- ga[[k]] + t(b) %*% ga[[k + 1L]]
  }
  g_s0 <- (ga[[1L]] + t(ga[[1L]])) / 2
  # Adjoint of the Lyapunov map Sigma0 = Sigma_zeta + B Sigma0 B':
  # Q - B' Q B = dF/dSigma0.
  q <- matrix(solve(t(st$lyap), as.vector(g_s0)), p, p)
  q <- (q + t(q)) / 2
  g_beta <- g_beta + 2 * q %*% b %*% st$s0
  gw <- .pgv_ggm_grad(q, st$par$delta_w, st$zw$w, ctx$lay$lower)
  gb <- .pgv_ggm_grad(g_between, st$par$delta_b, st$zb$w, ctx$lay$lower)
  c(g_mu, gw$omega, gw$delta, as.vector(g_beta), gb$omega, gb$delta)
}

# ---- forward-mode Jacobian (information + certificate) -----------------------

# Derivative of every implied moment w.r.t. each full parameter in `which`:
# a list of (waves*p)^2 x q `dsigma` (vec per column) and (waves*p) x q `dmu`.
#' @noRd
.pgv_jacobian <- function(st, ctx, which) {
  lay <- ctx$lay
  p <- lay$p
  waves <- ctx$waves
  b <- st$par$beta
  ggm_d <- function(domega, ddelta, delta, w) {
    t1 <- (ddelta * w) * rep(delta, each = p)
    t1 + t(t1) + outer(delta, delta) * (w %*% domega %*% w)
  }
  one <- function(a) {
    e <- numeric(lay$n_par)
    e[a] <- 1
    d <- .pgv_unpack(e, lay)
    dzw <- ggm_d(d$omega_w, d$delta_w, st$par$delta_w, st$zw$w)
    dzb <- ggm_d(d$omega_b, d$delta_b, st$par$delta_b, st$zb$w)
    rhs <- dzw + d$beta %*% st$s0 %*% t(b) + b %*% st$s0 %*% t(d$beta)
    ds0 <- matrix(solve(st$lyap, as.vector(rhs)), p, p)
    dlags <- Reduce(function(da, k) d$beta %*% st$lags[[k]] + b %*% da,
                    seq_len(waves - 1L), ds0, accumulate = TRUE)
    list(dsigma = as.vector(.pgv_assemble(dlags, dzb, waves)),
         dmu = rep(d$mu, waves))
  }
  parts <- lapply(which, one)
  list(dsigma = vapply(parts, `[[`, numeric((waves * p)^2), "dsigma"),
       dmu = vapply(parts, `[[`, numeric(waves * p), "dmu"))
}

# Expected Fisher information per observation on the parameters in `which`:
# J_mu' K J_mu + 1/2 tr(K dSigma_a K dSigma_b).
#' @noRd
.pgv_information <- function(st, jac) {
  k <- st$k
  dim_s <- nrow(k)
  kdk <- apply(jac$dsigma, 2L, function(v)
    as.vector(k %*% matrix(v, dim_s, dim_s) %*% k))
  kdk <- matrix(kdk, ncol = ncol(jac$dsigma))
  info <- crossprod(jac$dmu, k %*% jac$dmu) + 0.5 * crossprod(jac$dsigma, kdk)
  (info + t(info)) / 2
}

# Score of F by forward mode (the certificate's independent route).
#' @noRd
.pgv_score <- function(st, ctx, jac) {
  g <- .pgv_dsigma(st, ctx)
  as.vector(crossprod(jac$dsigma, as.vector(g)) -
              2 * crossprod(jac$dmu, st$k %*% st$r))
}

# ---- starting values ---------------------------------------------------------

#' @noRd
.pgv_spectral_shift <- function(x, floor = 1e-3) {
  x <- (x + t(x)) / 2
  ev <- min(eigen(x, symmetric = TRUE, only.values = TRUE)$values)
  if (ev < floor) x <- x + diag(floor - ev, nrow(x))
  x
}

#' @noRd
.pgv_ggm_from_cov <- function(sigma) {
  k <- solve(sigma)
  omega <- -stats::cov2cor(k)
  diag(omega) <- 0
  list(omega = omega, delta = 1 / sqrt(diag(k)))
}

# Moment-based starts in the spirit of psychonetrics' "version2": the averaged
# same-wave block as the stationary covariance, the first-last cross-wave block
# as the between part, and a lag-1 regression for beta (shrunk to be stationary).
# Starts affect only the optimizer path, not the optimum.
#' @noRd
.pgv_start <- function(ctx) {
  p <- ctx$lay$p
  waves <- ctx$waves
  rows <- function(i) (i - 1L) * p + seq_len(p)
  blk <- function(i, j) ctx$cov[rows(i), rows(j), drop = FALSE]
  stationary <- .pgv_spectral_shift(
    Reduce(`+`, lapply(seq_len(waves), function(t) blk(t, t))) / waves)
  lag1 <- Reduce(`+`, lapply(seq_len(waves - 1L), function(t) blk(t + 1L, t))) /
    (waves - 1L)
  far <- blk(waves, 1L)
  between <- .pgv_spectral_shift(pmin(far, t(far)))
  within <- .pgv_spectral_shift(stationary - between)
  beta <- (lag1 - between) %*% solve(within)
  shrink <- function(b) {
    rho <- max(Mod(eigen(b, only.values = TRUE)$values))
    if (rho > 0.9) b * 0.9 / rho else b
  }
  beta <- shrink(pmin(pmax(shrink(beta), -0.5), 0.5))
  innov <- .pgv_spectral_shift(within - beta %*% within %*% t(beta))
  zw <- .pgv_ggm_from_cov(innov)
  zb <- .pgv_ggm_from_cov(between)
  lay <- ctx$lay
  c(rowMeans(matrix(ctx$means, p, waves)), zw$omega[lay$lower], zw$delta,
    as.vector(beta), zb$omega[lay$lower], zb$delta)
}

# ---- optimizer ---------------------------------------------------------------

# Minimize F over the free parameters from one start: nlminb with the analytic
# gradient, then Fisher-scoring polish to drive the score to ~machine precision
# (the same scan-loose / refit-tight split the glasso path uses). Returns NULL
# when the start, or the point nlminb ends at, is outside the model's domain.
#' @noRd
.pgv_optimize <- function(ctx, start, max_iter = 2000L, polish_steps = 25L,
                          tol = 1e-9) {
  free <- which(ctx$free)
  start[!ctx$free] <- 0
  if (is.null(.pgv_state(start, ctx))) return(NULL)
  ctx$theta <- start
  objective <- function(x) {
    st <- .pgv_state(.pgv_full(x, ctx), ctx)
    if (is.null(st)) Inf else st$f
  }
  gradient <- function(x) {
    st <- .pgv_state(.pgv_full(x, ctx), ctx)
    if (is.null(st)) return(rep(NaN, length(x)))
    .pgv_gradient_full(st, ctx)[free]
  }
  is_delta <- ctx$lay$table$matrix[free] %in%
    c("delta_zeta_within", "delta_zeta_between")
  lower <- ifelse(is_delta, 1e-8, -Inf)
  opt <- stats::nlminb(start[free], objective, gradient, lower = lower,
                       control = list(eval.max = 4L * max_iter,
                                      iter.max = max_iter, rel.tol = 1e-14))
  theta <- .pgv_full(opt$par, ctx)
  st <- .pgv_state(theta, ctx)
  if (is.null(st)) return(NULL)
  # Fisher-scoring polish. E[Hessian of F] = 2 * per-observation information;
  # each step depends on the previous iterate, hence the loop.
  score <- .pgv_gradient_full(st, ctx)[free]
  steps <- 0L
  for (it in seq_len(polish_steps)) {
    if (max(abs(score)) < tol * 1e-3) break
    jac <- .pgv_jacobian(st, ctx, free)
    # A boundary solution (a near-singular between/innovation covariance) makes
    # the information near-singular; the pseudo-inverse step stays finite and
    # still reduces F along the identified directions.
    step <- as.vector(MASS_ginv(2 * .pgv_information(st, jac), tol = 1e-12) %*% score)
    accepted <- FALSE
    # Step halving: stop at the first fraction that does not increase F.
    for (half in 0:20) {
      cand <- theta
      cand[free] <- theta[free] - step / 2^half
      cst <- .pgv_state(cand, ctx)
      if (!is.null(cst) && cst$f <= st$f + 1e-13 * abs(st$f)) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) break
    cscore <- .pgv_gradient_full(cst, ctx)[free]
    if (max(abs(cscore)) >= max(abs(score)) && cst$f >= st$f) break
    theta <- cand
    st <- cst
    score <- cscore
    steps <- steps + 1L
  }
  list(theta = theta, state = st, iterations = opt$iterations,
       polish_steps = steps, message = opt$message)
}

# Run every start and keep the lowest objective. `starts` is a named list of
# full parameter vectors; the record of each start's objective is returned so
# the agreement of the starts (a check on local optima) can be reported.
#' @noRd
.pgv_best_fit <- function(ctx, starts) {
  fits <- lapply(starts, function(st) .pgv_optimize(ctx, st))
  f <- vapply(fits, function(x) if (is.null(x)) Inf else x$state$f, numeric(1))
  if (!any(is.finite(f))) {
    stop(errorCondition(
      "panel_gvar(): no start reached a feasible optimum (non-stationary temporal matrix or non-positive-definite covariance at every start).",
      class = c("psychnet_no_converge", "psychnet_error"), call = NULL))
  }
  best <- fits[[which.min(f)]]
  best$starts <- data.frame(start = names(starts), objective = f,
                            row.names = NULL, stringsAsFactors = FALSE)
  best
}

# A neutral start on standardized data: no temporal or network structure, half
# of each variable's variance within and half between persons.
#' @noRd
.pgv_simple_start <- function(ctx) {
  lay <- ctx$lay
  p <- lay$p
  c(rowMeans(matrix(ctx$means, p, ctx$waves)), numeric(nrow(lay$lower)),
    rep(sqrt(0.5 * (1 - 0.1^2)), p), as.vector(diag(0.1, p)),
    numeric(nrow(lay$lower)), rep(sqrt(0.5), p))
}

# Per-parameter factor mapping the standardized-scale parameters to the data
# scale (the model is scale-equivariant): mu and delta scale with s, beta[i, j]
# with s_i / s_j, the partial correlations are invariant.
#' @noRd
.pgv_scale_factor <- function(lay, s) {
  tab <- lay$table
  out <- rep(1, nrow(tab))
  by_s <- tab$matrix %in% c("mu", "delta_zeta_within", "delta_zeta_between")
  out[by_s] <- s[tab$row[by_s]]
  is_beta <- tab$matrix == "beta"
  out[is_beta] <- s[tab$row[is_beta]] / s[tab$col[is_beta]]
  out
}

# ---- step-up search -----------------------------------------------------------

# Modification indices (univariate score tests) for fixed parameters, exactly as
# psychonetrics::addMIs(): with g the gradient of F and I the per-observation
# expected information over the free parameters plus the candidates,
# MI_j = n g_j^2 / (4 V_jj), V = I_cc - I_cf I_ff^-1 I_fc (the candidates'
# information with the free parameters partialled out), and the expected
# parameter change EPC_j = -g_j / (2 V_jj) (the Newton step). A V_jj below
# machine tolerance (on psychonetrics' 4 n V scale) gives MI = NA.
#' @noRd
.pgv_modification <- function(st, ctx, candidates) {
  free <- which(ctx$free & !.pgv_at_bound(st, ctx))
  jac <- .pgv_jacobian(st, ctx, c(free, candidates))
  info <- .pgv_information(st, jac)
  fi <- seq_along(free)
  ci <- length(free) + seq_along(candidates)
  i_ff <- info[fi, fi, drop = FALSE]
  i_ff_inv <- if (min(eigen(i_ff, symmetric = TRUE, only.values = TRUE)$values) > 0)
    chol2inv(chol(i_ff)) else MASS_ginv(i_ff)
  i_cf <- info[ci, fi, drop = FALSE]
  v <- diag(info[ci, ci, drop = FALSE]) - rowSums((i_cf %*% i_ff_inv) * i_cf)
  g <- .pgv_gradient_full(st, ctx)[candidates]
  identified <- 4 * ctx$n * v >= sqrt(.Machine$double.eps)
  data.frame(index = candidates,
             mi = ifelse(identified, ctx$n * g^2 / (4 * v), NA_real_),
             epc = ifelse(identified, -g / (2 * v), NA_real_))
}

# psychonetrics::stepup() (non-greedy): free the fixed network parameter with
# the largest modification index while that index exceeds the chi-square
# critical value at `alpha`; keep it only if the refit improves the model by a
# chi-square difference test at `alpha` and, with criterion = "bic", lowers BIC
# (one extra parameter: BIC falls iff the LR statistic exceeds log n). The first
# rejected candidate ends the search. Each step depends on the model the
# previous one left, hence the loop.
#' @noRd
.pgv_stepup <- function(fit, ctx, sctx, network_rows, alpha, criterion) {
  added <- integer(0)
  # Scores, information and the LR statistic are invariant to the variables'
  # units, so the search runs on the standardized scale (data-scale moments can
  # be numerically singular when units differ by orders of magnitude).
  repeat {
    st <- .pgv_state(fit$theta, sctx)
    candidates <- which(network_rows & !sctx$free)
    if (!length(candidates)) break
    mod <- .pgv_modification(st, sctx, candidates)
    mi <- ifelse(is.na(mod$mi), 0, mod$mi)
    if (max(mi) <= stats::qchisq(alpha, 1, lower.tail = FALSE)) break
    best <- which.max(mi)
    j <- candidates[best]
    trial_ctx <- ctx
    trial_ctx$free[j] <- TRUE
    trial_sctx <- sctx
    trial_sctx$free <- trial_ctx$free
    epc_start <- fit$theta
    epc_start[j] <- mod$epc[best]
    trial <- .pgv_best_fit(trial_sctx, list(previous = fit$theta,
                                            epc = epc_start))
    lr <- ctx$n * (st$f - .pgv_state(trial$theta, sctx)$f)
    if (!(stats::pchisq(lr, 1, lower.tail = FALSE) < alpha)) break
    if (identical(criterion, "bic") && lr < log(ctx$n)) break
    fit <- trial
    ctx <- trial_ctx
    sctx <- trial_sctx
    added <- c(added, j)
  }
  list(fit = fit, ctx = ctx, sctx = sctx, added = added)
}

# ---- data preparation --------------------------------------------------------

#' @noRd
.pgv_bad_input <- function(msg) {
  stop(errorCondition(msg, class = c("psychnet_bad_input", "psychnet_error"),
                      call = NULL))
}

# Resolve `data` + (`vars`, `id`, `wave`) into a wide numeric matrix whose
# columns are ordered wave-major (wave 1 variables, wave 2 variables, ...),
# plus the variable labels and wave levels.
#' @noRd
.pgv_wide <- function(data, vars, id, wave) {
  if (!is.data.frame(data) && !is.matrix(data))
    .pgv_bad_input("`data` must be a data.frame or matrix.")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (is.matrix(vars)) {
    if (!is.character(vars))
      .pgv_bad_input("A wide `vars` design must be a character matrix of column names (rows = variables, columns = waves).")
    miss <- setdiff(as.vector(vars), names(data))
    if (length(miss))
      .pgv_bad_input(sprintf("Columns not found in `data`: %s.",
                             paste(miss, collapse = ", ")))
    labels <- rownames(vars) %||% sub("[._]?[wW]?[0-9]+$", "", vars[, 1L])
    if (anyDuplicated(labels) || any(!nzchar(labels)))
      labels <- paste0("V", seq_len(nrow(vars)))
    waves <- colnames(vars) %||% as.character(seq_len(ncol(vars)))
    cols <- as.vector(vars)
    wide <- data[, cols, drop = FALSE]
  } else {
    if (!is.character(vars) || !length(vars))
      .pgv_bad_input("`vars` must be a character vector of variable names (long data) or a character design matrix (wide data).")
    if (is.null(id) || is.null(wave))
      .pgv_bad_input("Long data needs `id` and `wave` column names; for wide data pass `vars` as a design matrix (rows = variables, columns = waves).")
    need <- c(vars, id, wave)
    miss <- setdiff(need, names(data))
    if (length(miss))
      .pgv_bad_input(sprintf("Columns not found in `data`: %s.",
                             paste(miss, collapse = ", ")))
    if (anyNA(data[[id]]))
      .pgv_bad_input("The person identifier column `id` must not contain missing values; unknown persons cannot be matched across waves.")
    key <- paste(data[[id]], data[[wave]], sep = "\r")
    if (anyDuplicated(key))
      .pgv_bad_input("Each `id` must have at most one row per `wave`.")
    wave_levels <- if (is.factor(data[[wave]])) levels(droplevels(data[[wave]]))
                   else as.character(sort(unique(data[[wave]])))
    ids <- unique(data[[id]])
    per_wave <- lapply(wave_levels, function(w) {
      sub <- data[as.character(data[[wave]]) == w, , drop = FALSE]
      sub[match(ids, sub[[id]]), vars, drop = FALSE]
    })
    wide <- do.call(cbind, per_wave)
    labels <- vars
    waves <- wave_levels
  }
  num <- vapply(wide, is.numeric, logical(1))
  if (!all(num))
    .pgv_bad_input("Every panel variable must be numeric.")
  list(wide = as.matrix(wide), labels = labels, waves = waves)
}

# Standardize the wide matrix (columns wave-major: wave 1 variables, wave 2
# variables, ...) exactly as psychonetrics' panelvar(standardize =): "z" uses
# one mean and SD per variable pooled over its waves, "z_per_wave" z-scores
# every variable-by-wave column separately; both use the available values (SD
# with denominator n - 1) before listwise deletion, as psychonetrics does.
#' @noRd
.pgv_standardize <- function(wide, p, waves, standardize) {
  if (identical(standardize, "none")) return(wide)
  z <- function(x, m, s) (x - m) / s
  if (identical(standardize, "z_per_wave")) {
    sds <- apply(wide, 2L, stats::sd, na.rm = TRUE)
    if (any(!is.finite(sds) | sds <= 0))
      .pgv_bad_input("standardize = \"z_per_wave\" needs every variable to vary within every wave.")
    return(sweep(sweep(wide, 2L, colMeans(wide, na.rm = TRUE)), 2L, sds, "/"))
  }
  var_of <- rep(seq_len(p), waves)
  pooled_mean <- vapply(seq_len(p), function(v)
    mean(wide[, var_of == v], na.rm = TRUE), numeric(1))
  pooled_sd <- vapply(seq_len(p), function(v)
    stats::sd(as.vector(wide[, var_of == v]), na.rm = TRUE), numeric(1))
  if (any(!is.finite(pooled_sd) | pooled_sd <= 0))
    .pgv_bad_input("standardize = \"z\" needs every variable to vary.")
  sweep(sweep(wide, 2L, pooled_mean[var_of]), 2L, pooled_sd[var_of], "/")
}

# ---- fit indices -------------------------------------------------------------

#' @noRd
# `f`, `f_base` are data-scale objectives; `sctx` the standardized context and
# `log_jacobian` the data-minus-standardized objective offset.
.pgv_fit_indices <- function(f, npar, f_base, npar_base, sctx, log_jacobian) {
  n <- sctx$n
  k <- nrow(sctx$cov)
  nobs <- k + k * (k + 1L) / 2L
  f_sat <- as.numeric(determinant(sctx$cov, logarithm = TRUE)$modulus) + k +
    log_jacobian
  logl <- -n / 2 * (f + k * log(2 * pi))
  chisq <- n * (f - f_sat)
  df <- nobs - npar
  b_chisq <- n * (f_base - f_sat)
  b_df <- nobs - npar_base
  excess <- max(chisq - df, 0)
  cfi_denom <- max(b_chisq - b_df, chisq - df, 0)
  data.frame(
    logl = logl, npar = npar, df = df, chisq = chisq,
    pvalue = stats::pchisq(chisq, df, lower.tail = FALSE),
    baseline_chisq = b_chisq, baseline_df = b_df,
    cfi = if (cfi_denom == 0) 1 else 1 - excess / cfi_denom,
    tli = (b_chisq / b_df - chisq / df) / (b_chisq / b_df - 1),
    rmsea = sqrt(excess / (df * n)),
    aic = -2 * logl + 2 * npar,
    bic = -2 * logl + npar * log(n),
    row.names = NULL)
}

# ---- the verb ----------------------------------------------------------------

#' Panel graphical VAR network (panel-GVAR)
#'
#' Estimates the three networks of the panel graphical vector-autoregressive
#' model (Epskamp 2020) from few-wave, many-person panel data: a directed
#' **temporal** network (within-person lag-1 effects), an undirected
#' **contemporaneous** network (partial correlations of the within-person
#' innovations), and an undirected **between**-person network (partial
#' correlations of the stable person means). The within-person process is a
#' stationary VAR(1) around a random intercept; means are stationary across
#' waves. This is the model fitted by `psychonetrics::panelgvar()`, estimated by
#' maximum likelihood in base R.
#'
#' By default the full (saturated) model is returned: every edge with its
#' standard error, p-value, 95% Wald interval and a `significant` flag at
#' `alpha`, drawn solid (significant) or dashed (not) by `plot()`.
#'
#' `select = "ebic"` instead minimizes F + lambda * sum|theta| over the
#' temporal (autoregressive diagonal included), contemporaneous and between
#' parameters, for 50 log-spaced lambdas from the smallest penalty that empties
#' every network, keeping the lambda with the lowest EBIC = -2 logL + k log n
#' + 4 k gamma log V (k free parameters, V wave x variable columns). The solver
#' (accelerated proximal gradient, then an active-set Newton polish) stays
#' inside the model's domain and is certified by the lasso KKT conditions; the
#' selected lambda is re-solved from three starts. The likelihood is not
#' convex, so the certificate proves a stationary point, and `print()` reports
#' how many starts agree.
#'
#' With `select = "prune"` the saturated model is fitted, every temporal,
#' contemporaneous, and between edge whose Wald p-value exceeds `alpha` is fixed
#' to zero in a single pass, and the model is refitted -- the
#' `psychonetrics::prune()` default workflow (`recursive = TRUE` repeats until
#' every remaining edge is significant). Then, with `stepup = TRUE` (default),
#' the `psychonetrics::stepup()` search adds pruned edges back one at a time
#' while modification indices, chi-square difference tests and BIC support
#' them, so `select = "prune"` equals psychonetrics'
#' `panelgvar() |> runmodel() |> prune() |> stepup()`.
#'
#' @param data A data.frame (or matrix). Either **long** (one row per person per
#'   wave; give `vars`, `id`, `wave`) or **wide** (one row per person; give
#'   `vars` as a design matrix).
#' @param vars Long data: character vector of the variable columns. Wide data: a
#'   character matrix of column names, one row per variable and one column per
#'   wave (the `psychonetrics` `vars` design); its row names become the labels.
#' @param id,wave Long data only: names of the person-identifier and wave
#'   columns. Person identifiers must be non-missing. Waves are ordered by
#'   `sort()` (or factor level order).
#' @param select Edge selection. `"none"` (default): every edge is estimated
#'   and reported with its standard error and p-value -- significance is shown,
#'   not imposed (the `psychonetrics` default). `"prune"`: remove edges with
#'   p > `alpha` and refit, then step up (see `stepup`) -- psychonetrics'
#'   `prune() |> stepup()`. `"ebic"`: penalized maximum likelihood (lasso on
#'   every temporal, contemporaneous and between parameter) with the penalty
#'   chosen by EBIC along a path -- psychonetrics'
#'   `panelgvar(estimator = "PML")`, the panel counterpart of [ebic_glasso()].
#' @param gamma,nlambda,lambda_min_ratio For `select = "ebic"`: the EBIC
#'   hyperparameter (default 0.5), the number of penalty values on the path
#'   (50; a finite integer >= 2) and the smallest as a fraction of the largest (1e-4) --
#'   psychonetrics' `find_penalized_lambda()` defaults.
#' @param refit For `select = "ebic"`: refit the selected network by
#'   unpenalized maximum likelihood (standard errors and p-values then come
#'   from that refit). Default `FALSE`: the penalized estimates are returned,
#'   and standard errors, p-values and `significant` are `NA` (they do not
#'   apply to lasso estimates).
#' @param alpha Significance level for the `significant` flag, the solid/dashed
#'   edges of `plot()`, and `select = "prune"`. Default `0.01`
#'   (`psychonetrics::prune()`'s default).
#' @param adjust Multiplicity correction applied to the edge p-values before
#'   pruning, one of [stats::p.adjust.methods]. Default `"none"`
#'   (`psychonetrics::prune()`'s default).
#' @param recursive Logical; repeat pruning until no edge has `p > alpha`.
#'   Default `FALSE`.
#' @param stepup Logical; after pruning, run the `psychonetrics::stepup()`
#'   search: re-add the pruned edge with the largest modification index while
#'   it is significant at `alpha` and the refit passes a chi-square difference
#'   test at `alpha` (and lowers BIC, with `criterion = "bic"`). Default `TRUE`,
#'   because one-pass pruning can remove parameters that substitute for each
#'   other (e.g. an autoregression and correlated cross-lagged effects), each
#'   non-significant alone but not jointly. `FALSE` reproduces
#'   `psychonetrics::prune()` alone. Used only with `select = "prune"`.
#' @param criterion Step-up acceptance criterion besides the chi-square test:
#'   `"bic"` (default, as `psychonetrics::stepup()`) or `"none"`.
#' @param standardize Rescale the data before fitting, as psychonetrics'
#'   `panelvar(standardize =)`: `"none"` (default), `"z"` (one mean and SD per
#'   variable, pooled over waves), or `"z_per_wave"` (z-score every variable
#'   separately within every wave). `"z_per_wave"` removes wave-to-wave drift
#'   in means and variances -- useful when the stationary-means assumption
#'   fails, at the price of no longer modelling that drift -- and makes the
#'   edges standardized. Standardization uses the available values (SD with
#'   denominator n - 1) before listwise deletion.
#' @param labels Optional node labels (default: `vars` for long data, the design
#'   row names for wide data).
#' @param na_method Missing-data handling. Only `"listwise"` is supported: a
#'   person missing any wave/variable is dropped (with a warning giving the
#'   count). A pairwise covariance is not a likelihood for this model and full-
#'   information ML is not yet implemented, so `"pairwise"` is refused.
#' @return An object of class `c("psychnet_panelgvar", "netobject_group")`: a
#'   named list of three [psychnet] networks, `temporal` (directed; edge
#'   `from -> to` weight `beta[to, from]`, autoregressive effects kept in
#'   `$beta`), `contemporaneous`, and `between` (undirected partial
#'   correlations; the innovation / between covariance in `$sigma`). Being a
#'   `netobject_group`, it plots as a grid with `cograph::splot()`. Read it with
#'   `as.data.frame()`: `what = "edges"` (one row per edge with standard error,
#'   p-value and 95% Wald interval), `what = "parameters"` (every model
#'   parameter), `what = "fit"` (one row of fit indices: chi-square, CFI,
#'   TLI, RMSEA, AIC, BIC), `what = "variance"` (one row per variable: the
#'   between-person share of its variance), or `what = "diagnostics"` (one
#'   row: the score
#'   certificate -- largest absolute gradient of the ML fit function on the
#'   free parameters, or the lasso KKT residual for `select = "ebic"` --
#'   boundary flag, SE availability, and how many of the three optimizer
#'   starts reached the best optimum). Fitting, inference and the certificate
#'   run on standardized variables (pooled SD per variable), so they are
#'   unit-free and the certificate does not depend on the variables' units;
#'   the model is scale-equivariant, so estimates and standard errors map back
#'   exactly and are returned on the data scale.
#'
#'   Errors are classed: `psychnet_bad_input` for unusable input (fewer than
#'   three waves, too few persons, non-numeric or missing columns, a requested
#'   `"pairwise"`), `psychnet_no_converge` if the optimizer ends at an
#'   infeasible point. Warnings are classed too: `psychnet_no_converge` (score
#'   above `1e-6` at an interior solution), `psychnet_boundary` (the final
#'   contemporaneous or between covariance is near-singular, i.e. partial
#'   correlations near +/-1 -- common with three waves; those estimates are
#'   unreliable, and the score, which need not vanish there, is quoted in the
#'   message),
#'   `psychnet_singular_information` (some parameters are not identified at
#'   the returned estimate: their standard errors are NA, the others exact),
#'   `psychnet_approximate_se` (the fit being pruned had a singular information;
#'   the parameters it does not identify were treated as having infinite
#'   standard errors, so pruning removes them), and
#'   `psychnet_listwise_drop` (persons removed for missing values).
#' @references Epskamp, S. (2020). Psychometric network models from time-series
#'   and panel data. *Psychometrika*, 85(1), 206-231.
#'   \doi{10.1007/s11336-020-09697-3}
#' @examples
#' # srl_panel: 100 students x 4 weeks, simulated from a panel-GVAR in which
#' # effort in one week drives efficacy in the next.
#' fit <- panel_gvar(srl_panel,
#'                   vars = c("efficacy", "planning", "effort", "enjoyment",
#'                            "anxiety"),
#'                   id = "student", wave = "week")
#' fit
#' as.data.frame(fit)
#' as.data.frame(fit, what = "fit")
#' as.data.frame(fit, what = "variance")
#' @export
panel_gvar <- function(data, vars, id = NULL, wave = NULL,
                       select = c("none", "prune", "ebic"),
                       alpha = 0.01, adjust = "none", recursive = FALSE,
                       stepup = TRUE, criterion = c("bic", "none"),
                       gamma = 0.5, nlambda = 50L, lambda_min_ratio = 1e-4,
                       refit = FALSE,
                       standardize = c("none", "z", "z_per_wave"),
                       labels = NULL, na_method = "listwise") {
  select <- match.arg(select)
  standardize <- match.arg(standardize)
  criterion <- match.arg(criterion)
  prune <- identical(select, "prune")
  if (!is.numeric(nlambda) || is.complex(nlambda) || length(nlambda) != 1L ||
      !is.finite(nlambda) || nlambda < 2 || nlambda > .Machine$integer.max ||
      nlambda != floor(nlambda))
    .pgv_bad_input("`nlambda` must be a single finite integer >= 2 and within the integer range.")
  stopifnot(
    "`refit` must be TRUE or FALSE" = isTRUE(refit) || isFALSE(refit),
    "`gamma` must be a single non-negative number" =
      is.numeric(gamma) && length(gamma) == 1L && is.finite(gamma) && gamma >= 0,
    "`lambda_min_ratio` must be a single number in (0, 1)" =
      is.numeric(lambda_min_ratio) && length(lambda_min_ratio) == 1L &&
      lambda_min_ratio > 0 && lambda_min_ratio < 1,
    "`stepup` must be TRUE or FALSE" = isTRUE(stepup) || isFALSE(stepup),
    "`recursive` must be TRUE or FALSE" = isTRUE(recursive) || isFALSE(recursive),
    "`alpha` must be a single number in (0, 1)" =
      is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
      alpha > 0 && alpha < 1
  )
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  if (!identical(na_method, "listwise"))
    .pgv_bad_input("panel_gvar() supports na_method = \"listwise\" only: a pairwise covariance is not a likelihood for the panel model, and full-information ML is not yet implemented.")
  wd <- .pgv_wide(data, vars, id, wave)
  p <- length(wd$labels)
  waves <- length(wd$waves)
  if (p < 2L) .pgv_bad_input("panel_gvar() needs at least two variables.")
  if (waves < 3L)
    .pgv_bad_input(sprintf("panel_gvar() needs at least three waves to separate within- from between-person variance; got %d.", waves))
  if (!is.null(labels)) {
    if (length(labels) != p)
      .pgv_bad_input(sprintf("`labels` must have length %d.", p))
    wd$labels <- as.character(labels)
  }
  wd$wide <- .pgv_standardize(wd$wide, p, waves, standardize)
  complete <- stats::complete.cases(wd$wide)
  if (!all(complete)) {
    warning(warningCondition(
      sprintf("panel_gvar(): dropped %d of %d persons with missing values (listwise).",
              sum(!complete), length(complete)),
      class = "psychnet_listwise_drop"))
  }
  y <- wd$wide[complete, , drop = FALSE]
  n <- nrow(y)
  if (n <= ncol(y))
    .pgv_bad_input(sprintf("panel_gvar() needs more complete persons (%d) than wave x variable columns (%d).", n, ncol(y)))
  means <- colMeans(y)
  cov_ml <- crossprod(sweep(y, 2L, means)) / n
  # Judged on the correlation scale, so the variables' units cannot trigger it.
  if (any(diag(cov_ml) <= 0) ||
      min(eigen(stats::cov2cor(cov_ml), symmetric = TRUE,
                only.values = TRUE)$values) <= 1e-10)
    .pgv_bad_input("The sample covariance of the stacked waves is singular (a constant or collinear variable).")
  lay <- .pgv_layout(p)
  network_rows <- lay$table$matrix %in%
    c("beta", "omega_zeta_within", "omega_zeta_between")
  # `ctx` is the data scale (inference, certificate, fit indices). The optimizer
  # works on standardized variables (pooled SD over waves) in `sctx`: the model
  # is scale-equivariant, and a common scale keeps raw 0-100 data from making
  # the parameters differ by orders of magnitude.
  ctx <- list(lay = lay, waves = waves, cov = cov_ml, means = means, n = n,
              free = rep(TRUE, lay$n_par))
  s <- sqrt(rowMeans(matrix(diag(cov_ml), p, waves)))
  s_full <- rep(s, waves)
  sctx <- ctx
  sctx$cov <- cov_ml / outer(s_full, s_full)
  sctx$means <- means / s_full
  to_data <- .pgv_scale_factor(lay, s)

  # The baseline (no temporal or network structure) is nested in every model
  # here, so it serves as a start and as a floor: no fit may end above it.
  base_ctx <- sctx
  base_ctx$free <- !network_rows
  base_fit <- .pgv_best_fit(base_ctx, list(simple = .pgv_simple_start(sctx)))
  fit <- .pgv_best_fit(sctx, list(moments = .pgv_start(sctx),
                                  simple = .pgv_simple_start(sctx),
                                  baseline = base_fit$theta))
  log_jacobian <- 2 * waves * sum(log(s))   # F(data) = F(standardized) + this
  starts <- fit$starts
  starts$objective <- starts$objective + log_jacobian
  starts$best <- abs(starts$objective - min(starts$objective)) <=
    1e-6 * max(1, abs(min(starts$objective)))

  pruned_any <- FALSE
  pruned_unidentified <- FALSE
  repeat {
    if (!prune) break
    inf <- .pgv_inference(.pgv_state(fit$theta, sctx), sctx,
                          approximate = TRUE)
    pruned_unidentified <- pruned_unidentified || any(inf$unidentified)
    test <- which(network_rows & ctx$free)
    p_adj <- stats::p.adjust(inf$p[test], method = adjust)
    drop <- test[!is.na(p_adj) & p_adj > alpha]
    if (!length(drop)) break
    ctx$free[drop] <- FALSE
    sctx$free <- ctx$free
    theta <- fit$theta
    theta[drop] <- 0
    pruned_any <- TRUE
    fit <- .pgv_best_fit(sctx, list(previous = theta, baseline = base_fit$theta))
    if (!recursive) break
  }
  added <- integer(0)
  if (prune && stepup && pruned_any) {
    su <- .pgv_stepup(fit, ctx, sctx, network_rows, alpha, criterion)
    fit <- su$fit
    ctx <- su$ctx
    sctx <- su$sctx
    added <- su$added
  }
  ebic <- NULL
  if (identical(select, "ebic")) {
    ebic <- .pgv_ebic_path(sctx, ctx, to_data, network_rows, base_fit$theta,
                           list(saturated = fit$theta,
                                moments = .pgv_start(sctx)),
                           gamma, as.integer(nlambda), lambda_min_ratio,
                           log_jacobian)
    starts <- ebic$starts
    zero <- ebic$pen & ebic$theta == 0
    ctx$free <- ctx$free & !zero
    sctx$free <- ctx$free
    fit <- list(theta = ebic$theta)
    if (refit) {
      fit <- .pgv_best_fit(sctx, list(selected = ebic$theta,
                                      baseline = base_fit$theta))
    }
  }

  # Every post-fit computation runs on the standardized scale: scores,
  # information, p-values, the certificate and boundary checks are unit-free
  # there, and data-scale moments can be numerically singular when the
  # variables' units differ by orders of magnitude. Estimates and standard
  # errors map back exactly (multiply by `to_data`).
  sctx$free <- ctx$free
  st <- .pgv_state(fit$theta, sctx)
  base_st <- .pgv_state(base_fit$theta, sctx)
  if (st$f > base_st$f + 1e-8 * abs(base_st$f) &&
      !(identical(select, "ebic") && !refit)) {
    stop(errorCondition(
      sprintf("panel_gvar(): the fitted model (F = %.6g) is worse than its nested baseline (F = %.6g); the optimization failed.",
              st$f + log_jacobian, base_st$f + log_jacobian),
      class = c("psychnet_no_converge", "psychnet_error"), call = NULL))
  }
  inf <- .pgv_inference(st, sctx)
  if (pruned_unidentified) {
    warning(warningCondition(
      "panel_gvar(): the Fisher information at a fit being pruned was singular; parameters the data do not identify were treated as non-significant (infinite standard error) and pruned.",
      class = "psychnet_approximate_se"))
  }
  if (any(inf$unidentified) && !(identical(select, "ebic") && !refit)) {
    tab <- lay$table
    who <- sprintf("%s[%s, %s]", tab$matrix, wd$labels[tab$row],
                   wd$labels[tab$col])[inf$unidentified]
    warning(warningCondition(
      sprintf("panel_gvar(): %d parameter(s) are not identified at the returned estimate (singular Fisher information): %s. Their standard errors and p-values are NA; all other standard errors are exact.",
              length(who), paste(who, collapse = ", ")),
      class = "psychnet_singular_information"))
  }
  if (identical(select, "ebic") && !refit) {
    # Lasso estimates: Wald standard errors and p-values do not apply, and the
    # certificate is the lasso KKT residual (penalty weights map the
    # data-scale penalty to this scale).
    inf$se[] <- NA_real_
    inf$p[] <- NA_real_
    inf$kkt <- .pgv_pen_kkt(.pgv_gradient_full(st, sctx), inf$est,
                            sctx$free | ebic$pen, ebic$pen, to_data,
                            ebic$lambda, inf$bound)
  }
  boundary <- .pgv_boundary(st, sctx, wd$labels, inf$kkt)
  # At a near-singular GGM block the optimum can sit on the edge of the
  # parameter space, where the score need not vanish -- so a large residual
  # there is reported inside the boundary warning, not as non-convergence.
  if (length(boundary)) {
    warning(warningCondition(
      sprintf("panel_gvar(): boundary solution -- %s. Estimates and standard errors involving it are unreliable (KKT residual %.2e; the optimum lies on the edge of the parameter space). More waves or persons, or pruning, usually resolve this.",
              paste(boundary, collapse = "; "), inf$kkt),
      class = "psychnet_boundary"))
  } else if (inf$kkt > 1e-6) {
    warning(warningCondition(
      sprintf("panel_gvar(): the ML score is %.2e at the returned estimate (expected < 1e-6); the fit may not have converged.",
              inf$kkt),
      class = "psychnet_no_converge"))
  }
  fit_tab <- .pgv_fit_indices(st$f + log_jacobian, sum(ctx$free),
                              base_st$f + log_jacobian, sum(base_ctx$free),
                              sctx, log_jacobian)
  inf$est <- inf$est * to_data
  inf$se <- inf$se * to_data
  .pgv_object(st, inf, ctx, s, wd$labels, wd$waves, fit_tab, boundary,
              list(select = select, standardize = standardize,
                   prune = prune, alpha = alpha,
                   adjust = adjust, recursive = recursive, pruned = pruned_any,
                   stepup = prune && stepup, criterion = criterion,
                   added = length(added),
                   lambda = ebic$lambda, gamma = gamma, refit = refit,
                   path = ebic$path), starts)
}

# Rank of an information matrix judged on its correlation scale,
# D^-1/2 h D^-1/2 (D = diag(h)): a near-zero variance makes some derivatives
# enormous, and a relative eigenvalue cut-off on the raw matrix then declares
# almost every direction null (31 of 40 on real data, means included); on the
# correlation scale the deficiency is the genuine one (there: the 10
# between-person parameters, 4 null directions, an eigenvalue gap of 1e-9).
#' @noRd
.pgv_scaled_eigen <- function(h) {
  h <- (h + t(h)) / 2
  d <- diag(h)
  scale <- ifelse(d > 0, 1 / sqrt(d), 0)
  ev <- eigen(h * outer(scale, scale), symmetric = TRUE)
  list(values = ev$values, vectors = ev$vectors, scale = scale, d = d)
}

# Standard errors from a singular information matrix `h` (n * I), on its
# correlation scale (see .pgv_scaled_eigen). A parameter with a non-negligible
# loading on the null space (eigenvalues below 1e-10 of the largest) is
# unidentified and gets Inf. Every other coordinate is an identifiable
# function, whose variance is e' G e for any generalized inverse G -- here
# D^-1/2 (scaled)^+ D^-1/2. (Dropping the unidentified coordinates and
# inverting the rest would also drop identifiable combinations of them, e.g.
# b2 + b3 when only b2 - b3 is unidentified, and understate the SEs.)
#' @noRd
.pgv_identified_se <- function(h, tol = 1e-10, load_tol = 1e-4) {
  ev <- .pgv_scaled_eigen(h)
  null <- ev$values < tol * max(ev$values)
  loading <- rowSums(ev$vectors[, null, drop = FALSE]^2)
  range <- ev$vectors[, !null, drop = FALSE]
  variance <- ev$scale^2 *
    rowSums(range^2 * rep(1 / ev$values[!null], each = nrow(h)))
  ifelse(loading > load_tol | ev$d <= 0, Inf, sqrt(variance))
}

# Free variance (delta) parameters ON the optimizer's lower bound (1e-8 on the
# standardized scale every post-fit computation uses): a variance estimated as
# zero (a Heywood case). Only these are active constraints. `ctx` must be the
# standardized context.
#' @noRd
.pgv_at_bound <- function(st, ctx, bound = 1e-8) {
  tab <- ctx$lay$table
  is_delta <- tab$matrix %in% c("delta_zeta_within", "delta_zeta_between")
  est <- .pgv_theta_from_state(st, ctx$lay)
  ctx$free & is_delta & est <= bound * (1 + 1e-6)
}

# Diagnostic only (never used to relax the certificate): free variance scales
# below `tol` on the standardized scale, i.e. a variance near zero.
#' @noRd
.pgv_near_zero <- function(st, ctx, tol = 0.01) {
  tab <- ctx$lay$table
  is_delta <- tab$matrix %in% c("delta_zeta_within", "delta_zeta_between")
  ctx$free & is_delta & .pgv_theta_from_state(st, ctx$lay) < tol
}

# Standard errors (inverse expected Fisher information / n), two-sided Wald
# p-values, and the score certificate. A variance on its lower bound is an
# active constraint: it is left out of the information (its SE is NA), and the
# certificate is the KKT condition -- score ~0 on interior parameters, score
# >= 0 (cannot improve by moving into the domain) on bound ones.
#' @noRd
.pgv_inference <- function(st, ctx, approximate = FALSE) {
  bound <- .pgv_at_bound(st, ctx)
  interior <- which(ctx$free & !bound)
  jac <- .pgv_jacobian(st, ctx, which(ctx$free))
  keep <- which(ctx$free) %in% interior
  info <- .pgv_information(st, list(dsigma = jac$dsigma[, keep, drop = FALSE],
                                    dmu = jac$dmu[, keep, drop = FALSE]))
  se <- rep(NA_real_, ctx$lay$n_par)
  unidentified <- rep(FALSE, ctx$lay$n_par)
  scaled <- .pgv_scaled_eigen(info)
  if (all(scaled$d > 0) &&
      min(scaled$values) > 1e-10 * max(scaled$values)) {
    se[interior] <- sqrt(diag(chol2inv(chol(ctx$n * info))))
  } else {
    # Singular information: parameters loading on its null space are not
    # identified (infinite SE); every other parameter keeps its exact SE, the
    # variance of an identifiable function. For pruning (`approximate`) the
    # unidentified ones keep Inf, so they are removed -- the outcome
    # psychonetrics reaches through huge SEs; for reporting they become NA.
    se[interior] <- .pgv_identified_se(ctx$n * info)
    unidentified <- !is.finite(se)
    if (!approximate) se[unidentified] <- NA_real_
  }
  est <- .pgv_theta_from_state(st, ctx$lay)
  p <- rep(NA_real_, ctx$lay$n_par)
  p[interior] <- 2 * stats::pnorm(abs(est[interior]) / se[interior],
                                  lower.tail = FALSE)
  score <- .pgv_score(st, ctx, jac)
  kkt <- max(abs(score[keep]), pmax(0, -score[!keep]), 0)
  list(est = est, se = se, p = p, kkt = kkt, bound = bound,
       unidentified = unidentified & ctx$free)
}

# Describe every part of the solution on (or heading for) a boundary of the
# parameter space. On it: a near-singular GGM block (smallest eigenvalue of
# I - Omega below 1e-3, i.e. partial correlations near +/-1) or a variance
# estimated at zero. Heading for it, reported only when the KKT residual says
# the fit has not converged: eigenvalue below 0.05 or a variance scale below
# 0.01 on the standardized scale -- the GGM parameterization reaches a rank-deficient
# covariance only in the limit, so the optimizer crawls along that ridge.
#' @noRd
.pgv_boundary <- function(st, ctx, labels, kkt = 0, tol = 1e-3,
                          near_tol = 0.05, near_sd = 0.01) {
  par <- st$par
  lam <- c(contemporaneous = min(eigen(diag(nrow(par$omega_w)) - par$omega_w,
                                       symmetric = TRUE, only.values = TRUE)$values),
           "between-person" = min(eigen(diag(nrow(par$omega_b)) - par$omega_b,
                                        symmetric = TRUE, only.values = TRUE)$values))
  block_name <- c(contemporaneous = "contemporaneous (innovation)",
                  "between-person" = "between-person")
  stalled <- kkt > 1e-6
  blocks <- c(sprintf("%s covariance near-singular", block_name[lam < tol]),
              if (stalled) sprintf("%s covariance approaching singularity",
                                   block_name[lam >= tol & lam < near_tol]))
  tab <- ctx$lay$table
  bound <- .pgv_at_bound(st, ctx)
  near <- if (stalled) .pgv_near_zero(st, ctx, tol = near_sd) & !bound
          else rep(FALSE, length(bound))
  what <- c(delta_zeta_within = "innovation", delta_zeta_between = "between-person")
  c(blocks,
    sprintf("%s variance of %s at zero", what[tab$matrix[bound]],
            labels[tab$row[bound]]),
    sprintf("%s variance of %s approaching zero", what[tab$matrix[near]],
            labels[tab$row[near]]))
}

#' @noRd
.pgv_theta_from_state <- function(st, lay) {
  par <- st$par
  c(par$mu, par$omega_w[lay$lower], par$delta_w, as.vector(par$beta),
    par$omega_b[lay$lower], par$delta_b)
}

# ---- result object -----------------------------------------------------------

# The result object. `st` is the standardized-scale state and `s` the
# per-variable scales: beta[i, j] maps by s_i / s_j, covariances by s_i s_j,
# partial correlations are unchanged; the variance split is computed here, on
# the standardized scale, and stored.
#' @noRd
.pgv_object <- function(st, inf, ctx, s, labels, waves, fit_tab, boundary,
                        settings, starts) {
  lay <- ctx$lay
  par <- st$par
  z <- stats::qnorm(0.975)
  tab <- lay$table
  params <- data.frame(
    matrix = tab$matrix, row = labels[tab$row],
    col = ifelse(tab$matrix == "mu", NA_character_, labels[tab$col]),
    estimate = inf$est, se = inf$se, p = inf$p,
    ci_lower = inf$est - z * inf$se, ci_upper = inf$est + z * inf$se,
    fixed = !ctx$free, stringsAsFactors = FALSE, row.names = NULL)
  ss <- outer(s, s)
  sz <- st$zw$sigma * ss
  sb <- st$zb$sigma * ss
  dimnames(sz) <- dimnames(sb) <- list(labels, labels)
  beta <- par$beta * outer(s, 1 / s)
  dimnames(beta) <- list(labels, labels)
  between <- diag(st$zb$sigma) * s^2
  within <- diag(st$s0) * s^2
  variance <- data.frame(variable = labels, between = between, within = within,
                         share_between = between / (between + within),
                         row.names = NULL, stringsAsFactors = FALSE)
  nets <- list(
    temporal = .new_psychnet(t(beta), labels, "panel_gvar", directed = TRUE,
                             n_obs = ctx$n, extra = list(beta = beta)),
    contemporaneous = .new_psychnet(par$omega_w, labels, "panel_gvar",
                                    directed = FALSE, n_obs = ctx$n,
                                    extra = list(sigma = sz)),
    between = .new_psychnet(par$omega_b, labels, "panel_gvar",
                            directed = FALSE, n_obs = ctx$n,
                            extra = list(sigma = sb)))
  structure(nets,
            parameters = params, fit = fit_tab, kkt = inf$kkt, n = ctx$n,
            starts = starts, boundary = boundary, variance = variance,
            waves = waves, settings = settings, group_col = "network",
            class = c("psychnet_panelgvar", "netobject_group"))
}

# Whether p-values (and so `significant`) apply: not for lasso estimates
# (select = "ebic" without refit).
#' @noRd
.pgv_sig_applies <- function(x) {
  st <- attr(x, "settings")
  !(identical(st$select, "ebic") && !isTRUE(st$refit))
}

# Tidy edge rows (one per estimated edge) from the parameter table.
#' @noRd
.pgv_edges <- function(params, include_zero, alpha, applies = TRUE) {
  net_of <- c(beta = "temporal", omega_zeta_within = "contemporaneous",
              omega_zeta_between = "between")
  e <- params[params$matrix %in% names(net_of), , drop = FALSE]
  if (!include_zero) e <- e[!e$fixed, , drop = FALSE]
  # beta[row, col] is the effect of col (t-1) on row (t); the omegas are stored
  # lower-triangular (row > col). Either way the edge runs col -> row.
  out <- data.frame(
    network = unname(net_of[e$matrix]),
    from = e$col, to = e$row,
    weight = e$estimate, se = e$se, p = e$p,
    ci_lower = e$ci_lower, ci_upper = e$ci_upper,
    significant = if (applies) !is.na(e$p) & e$p < alpha
                  else rep(NA, nrow(e)),
    stringsAsFactors = FALSE, row.names = NULL)
  out$network <- factor(out$network,
                        levels = c("temporal", "contemporaneous", "between"))
  out <- out[order(out$network, out$from, out$to), , drop = FALSE]
  out$network <- as.character(out$network)
  rownames(out) <- NULL
  out
}

#' Read a panel graphical VAR fit
#'
#' Methods for the object returned by [panel_gvar()].
#'
#' @param x,object A `psychnet_panelgvar` object.
#' @param what For `as.data.frame()`: `"edges"` (default; one row per estimated
#'   edge across the three networks, with standard error, p-value and 95% Wald
#'   interval -- temporal rows include the autoregressive `from == to`
#'   effects), `"parameters"` (every model parameter, including means, scaling
#'   `delta`s, and fixed/pruned entries), `"fit"` (one row of fit indices), or
#'   `"diagnostics"` (one row: persons `n`, `waves`, the ML `score`
#'   certificate, `boundary` / `boundary_part`, `se_available`, and
#'   `starts_agree` of `starts` optimizer starts reaching the best saturated
#'   optimum), or `"variance"` (one row per variable: the model-implied
#'   `between`-person variance, the stationary `within`-person variance, and
#'   `share_between`, their ratio to the total -- an intraclass correlation).
#' @param include_zero For `what = "edges"`: also return pruned (fixed-zero)
#'   edges. Default `FALSE`.
#' @param row.names,optional Ignored (S3 consistency).
#' @param ... Passed to `cograph::splot()` by `plot()`; named arguments
#'   override the plotting defaults, with unique argument abbreviations
#'   accepted. Otherwise ignored.
#' @return `as.data.frame()`: a data.frame as described under `what`. The edge
#'   table has columns `network`, `from`, `to`, `weight`, `se`, `p`,
#'   `ci_lower`, `ci_upper`, `significant` (p < the fit's `alpha`); the
#'   parameter table has `matrix`, `row`, `col`,
#'   `estimate`, `se`, `p`, `ci_lower`, `ci_upper`, `fixed`; the diagnostics
#'   row has `n`, `waves`, `score`, `boundary`, `boundary_part`,
#'   `se_available`, `starts_agree`, `starts`; the variance table has
#'   `variable`, `between`, `within`, `share_between`. `print()` and
#'   `plot()` return `x` invisibly. `summary()` returns a data.frame with one
#'   row per network (`network`, `directed`, `nodes`, `edges`, `significant`,
#'   `mean_abs_weight`).
#' @examples
#' # srl_panel: 100 students x 4 weeks, simulated from a panel-GVAR in which
#' # effort in one week drives efficacy in the next.
#' fit <- panel_gvar(srl_panel,
#'                   vars = c("efficacy", "planning", "effort", "enjoyment",
#'                            "anxiety"),
#'                   id = "student", wave = "week")
#' summary(fit)
#' as.data.frame(fit, what = "parameters")
#' as.data.frame(fit, what = "diagnostics")
#' if (requireNamespace("cograph", quietly = TRUE)) plot(fit)
#' @name psychnet_panelgvar
#' @export
as.data.frame.psychnet_panelgvar <- function(x, row.names = NULL,
                                             optional = FALSE, ...,
                                             what = c("edges", "parameters",
                                                      "fit", "diagnostics",
                                                      "variance"),
                                             include_zero = FALSE) {
  what <- match.arg(what)
  switch(what,
    edges = .pgv_edges(attr(x, "parameters"), include_zero,
                       attr(x, "settings")$alpha, .pgv_sig_applies(x)),
    parameters = attr(x, "parameters"),
    fit = attr(x, "fit"),
    diagnostics = .pgv_diagnostics(x),
    variance = .pgv_variance(x))
}

# Per-variable split of the model-implied variance into the stable
# between-person part (Sigma_between) and the stationary within-person part
# (Sigma0, the variance of the VAR(1) process), computed at fit time on the
# standardized scale and rescaled.
#' @noRd
.pgv_variance <- function(x) attr(x, "variance")

# One row of fit diagnostics: sample, certificate, boundary, SE availability,
# and the agreement of the optimizer starts.
#' @noRd
.pgv_diagnostics <- function(x) {
  params <- attr(x, "parameters")
  starts <- attr(x, "starts")
  boundary <- attr(x, "boundary")
  data.frame(n = attr(x, "n"), waves = length(attr(x, "waves")),
             score = attr(x, "kkt"), boundary = length(boundary) > 0L,
             boundary_part = if (length(boundary))
               paste(boundary, collapse = "; ") else NA_character_,
             se_available = !anyNA(params$se[!params$fixed]),
             starts_agree = sum(starts$best), starts = nrow(starts),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' @rdname psychnet_panelgvar
#' @export
print.psychnet_panelgvar <- function(x, ...) {
  ft <- attr(x, "fit")
  st <- attr(x, "settings")
  cat("<psychnet_panelgvar> panel graphical VAR (temporal + contemporaneous + between)\n")
  cat(sprintf("  persons: %d   waves: %d   nodes: %d\n", attr(x, "n"),
              length(attr(x, "waves")), nrow(x$temporal$nodes)))
  edges <- .pgv_edges(attr(x, "parameters"), include_zero = FALSE,
                     attr(x, "settings")$alpha, .pgv_sig_applies(x))
  applies <- .pgv_sig_applies(x)
  tally <- function(rows) if (applies)
    sprintf("%d/%d", sum(edges$significant[rows]), sum(rows))
    else as.character(sum(rows))
  cross <- edges$from != edges$to
  cat(sprintf("  edges (%s): temporal %s (+%s autoregressive)   contemporaneous %s   between %s\n",
              if (applies) sprintf("significant/estimated, p < %g", st$alpha)
              else "EBIC-selected",
              tally(edges$network == "temporal" & cross),
              tally(edges$network == "temporal" & !cross),
              tally(edges$network == "contemporaneous"),
              tally(edges$network == "between")))
  if (!identical(st$standardize %||% "none", "none"))
    cat(sprintf("  standardize: %s\n", st$standardize))
  cat(sprintf("  selection: %s\n", switch(st$select %||% "none",
    prune = sprintf("prune (alpha = %g, adjust = %s%s%s)", st$alpha, st$adjust,
                    if (isTRUE(st$recursive)) ", recursive" else "",
                    if (isTRUE(st$stepup))
                      sprintf("; stepup (%s) re-added %d", st$criterion, st$added)
                    else ""),
    ebic = sprintf("EBIC lasso (gamma = %g, lambda = %.4g%s)", st$gamma,
                   st$lambda, if (isTRUE(st$refit)) ", refitted" else ""),
    "none (full model)")))
  cat(sprintf("  fit: chisq(%d) = %.2f, p = %.3g   CFI = %.3f   TLI = %.3f   RMSEA = %.3f\n",
              as.integer(ft$df), ft$chisq, ft$pvalue, ft$cfi, ft$tli, ft$rmsea))
  is_ebic <- identical(st$select, "ebic")
  cat(sprintf("  optimality (%s): %.2e\n",
              if (is_ebic && !isTRUE(st$refit)) "lasso KKT residual" else "ML score",
              attr(x, "kkt")))
  starts <- attr(x, "starts")
  cat(sprintf("  starts: %d of %d reached the best %s optimum\n",
              sum(starts$best), nrow(starts),
              if (is_ebic) "penalized" else "saturated"))
  invisible(x)
}

#' @rdname psychnet_panelgvar
#' @export
summary.psychnet_panelgvar <- function(object, ...) {
  edges <- .pgv_edges(attr(object, "parameters"), include_zero = FALSE,
                     attr(object, "settings")$alpha, .pgv_sig_applies(object))
  edges <- edges[edges$from != edges$to, , drop = FALSE]
  nets <- c("temporal", "contemporaneous", "between")
  data.frame(
    network = nets,
    directed = c(TRUE, FALSE, FALSE),
    nodes = nrow(object$temporal$nodes),
    edges = vapply(nets, function(nm) sum(edges$network == nm), integer(1)),
    significant = vapply(nets, function(nm)
      if (.pgv_sig_applies(object))
        sum(edges$significant[edges$network == nm]) else NA_integer_,
      integer(1)),
    mean_abs_weight = vapply(nets, function(nm) {
      w <- edges$weight[edges$network == nm]
      if (length(w)) mean(abs(w)) else 0
    }, numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE)
}

# The temporal network as drawn: the stored network has a zero diagonal (the
# psychnet contract), so the autoregressive effects are added back as self-loop
# edges -- cograph draws a loop only from an `edges` row with from == to.
#' @noRd
.pgv_plot_temporal <- function(net) {
  ar <- diag(net$beta)
  keep <- abs(ar) > 1e-12
  net$weights <- t(net$beta)
  dimnames(net$weights) <- list(net$nodes$label, net$nodes$label)
  net$edges <- rbind(net$edges,
                     data.frame(from = which(keep), to = which(keep),
                                weight = ar[keep]))
  net
}

# Significance of each edge of a network as drawn, in the order of its `edges`
# rows (integer node indices; undirected edges stored from < to, which is the
# from/to orientation of the tidy edge table).
#' @noRd
.pgv_edge_significance <- function(net, name, edges) {
  labels <- net$nodes$label
  tab <- edges[edges$network == name, , drop = FALSE]
  # Match ordered integer pairs, so spaces or separators in labels cannot
  # collapse two different edges onto the same key.
  pair <- function(from, to) (from - 1) * length(labels) + to
  hit <- match(pair(net$edges$from, net$edges$to),
               pair(match(tab$from, labels), match(tab$to, labels)))
  !is.na(hit) & tab$significant[hit]
}

#' @rdname psychnet_panelgvar
#' @param edge_labels For `plot()`: print each edge's weight on it. Default
#'   `TRUE`.
#' @details By default, `plot()` draws the three networks side by side on one shared circle
#'   layout (so a node sits in the same place in every panel), with the
#'   temporal autoregressive effects as self-loops. Line type shows
#'   significance: solid edges have p < `alpha` (the fit's `alpha`), dashed
#'   edges do not (an edge without a p-value counts as not significant). Colour
#'   shows sign -- blue positive, vermillion negative (Okabe-Ito) -- and the
#'   signed weight is printed on each edge, so sign never rests on colour
#'   alone. A network with nothing estimated shows "no edges". It needs the
#'   `cograph` package. If `edge_style` is overridden, the legend no longer
#'   assigns significance to solid or dashed lines.
#' @export
plot.psychnet_panelgvar <- function(x, edge_labels = TRUE, ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop("Network plotting needs the 'cograph' package. Install it, or read ",
         "the networks with as.data.frame().", call. = FALSE)
  }
  alpha <- attr(x, "settings")$alpha
  applies <- .pgv_sig_applies(x)
  edges <- .pgv_edges(attr(x, "parameters"), include_zero = FALSE, alpha,
                      applies)
  nets <- list(temporal = .pgv_plot_temporal(x$temporal),
               contemporaneous = x$contemporaneous, between = x$between)
  p <- nrow(x$temporal$nodes)
  angle <- pi / 2 - 2 * pi * (seq_len(p) - 1) / p
  coords <- cbind(cos(angle), sin(angle))
  op <- graphics::par(mfrow = c(1L, 3L), oma = c(2, 0, 0, 0))
  on.exit(graphics::par(op), add = TRUE)
  dots <- list(...)
  # Resolve abbreviations before combining with defaults: otherwise supplying
  # node_fill explicitly would leave a caller's node_f unmatched in splot().
  dot_names <- names(dots) %||% rep("", length(dots))
  formal_names <- setdiff(names(formals(cograph::splot)), "...")
  named <- nzchar(dot_names)
  for (i in which(named)) {
    if (dot_names[i] %in% formal_names) next
    matches <- formal_names[startsWith(formal_names, dot_names[i])]
    if (length(matches) > 1L)
      .pgv_bad_input(sprintf("Plot argument `%s` is ambiguous; use a full argument name.", dot_names[i]))
    if (length(matches) == 1L) dot_names[i] <- matches
  }
  if (anyDuplicated(dot_names[named]))
    .pgv_bad_input("Plot arguments must not name the same option more than once.")
  names(dots) <- dot_names
  Map(function(net, name) {
    n_edges <- nrow(net$edges)
    significant <- if (applies) .pgv_edge_significance(net, name, edges)
                   else rep(TRUE, n_edges)
    defaults <- list(layout = coords,
                   title = if (applies)
                     sprintf("%s: %d of %d significant", name,
                             sum(significant), n_edges)
                   else sprintf("%s: %d edges (EBIC)", name, n_edges),
                   layout_margin = 0.3,
                   edge_labels = edge_labels && n_edges > 0L,
                   edge_label_size = 1,
                   edge_positive_color = "#0072B2",
                   edge_negative_color = "#D55E00",
                   edge_style = ifelse(significant, 1L, 2L),
                   node_fill = "#56B4E9")
    defaults[names(dots)[named]] <- dots[named]
    do.call(cograph::splot, c(list(net), defaults, dots[!named]))
    # The circle layout is centred, so the middle of the panel is free.
    if (!n_edges) {
      usr <- graphics::par("usr")
      graphics::text(mean(usr[1:2]), mean(usr[3:4]), "no edges",
                     col = "#999999")
    }
  }, nets, names(nets))
  sign_key <- sprintf("%s: positive   %s: negative",
    paste(dots$edge_positive_color %||% "blue", collapse = ", "),
    paste(dots$edge_negative_color %||% "vermillion", collapse = ", "))
  key <- if (applies && !"edge_style" %in% dot_names)
    sprintf("solid: p < %g   dashed: not significant   %s", alpha, sign_key)
  else if (applies) sprintf("custom line styles   %s", sign_key)
  else sprintf("EBIC-selected edges (lambda = %.3g; penalized estimates)   %s",
               attr(x, "settings")$lambda, sign_key)
  graphics::mtext(key, side = 1L, outer = TRUE, line = 0.5, cex = 0.8,
                  col = "#555555")
  invisible(x)
}
