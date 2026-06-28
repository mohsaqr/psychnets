# Clean-room graphical lasso (Friedman, Hastie & Tibshirani 2008, Biostatistics)
# and EBIC model selection (Foygel & Drton 2010). Pure base R, no compiled
# dependency. The objective
#     min_{Theta > 0}  -log det Theta + tr(S Theta) + rho * sum_{i != j}|Theta_ij|
# is strictly convex, so its minimiser is unique; glasso_kkt() grades any
# candidate against that unique optimum with no reference solver.
#
# The covariance (W) block-coordinate descent solves, for each column, an
# ordinary lasso by Gauss-Seidel soft-thresholding, then reconstructs the
# precision matrix from the converged W and the lasso coefficients. The outer
# column sweep and inner coordinate descent are Gauss-Seidel: each update
# consumes the latest values, a genuine sequential dependency (the inner dot
# products are vectorised).

# --- soft-threshold ----------------------------------------------------------
.soft <- function(z, g) sign(z) * pmax(abs(z) - g, 0)

# --- single-column lasso (coordinate descent) --------------------------------
# Pathwise coordinate descent with covariance/residual updates (Friedman,
# Hastie & Tibshirani 2010): keep a running r = W11 %*% beta so each coordinate
# reads r[k] in O(1), and apply an O(p) rank-1 update to r only when a
# coordinate actually moves -- in a sparse network most stay zero, so the inner
# sweeps are far cheaper than recomputing a full dot product per coordinate. The
# update is algebraically identical to the naive sweep, so it converges to the
# same unique optimum of the strictly convex objective.
.glasso_lasso_column <- function(W11, s12, beta, rho, max_inner, tol_inner, dgg) {
  pp  <- length(s12)
  r   <- as.numeric(W11 %*% beta)           # running W11 %*% beta
  for (inner in seq_len(max_inner)) {
    max_diff <- 0
    for (k in seq_len(pp)) {
      wkk <- dgg[k]
      if (wkk < 1e-12) {
        new_k <- 0
      } else {
        partial <- s12[k] - (r[k] - wkk * beta[k])
        thr <- abs(partial) - rho           # soft-threshold, inlined (scalar)
        new_k <- if (thr > 0) (if (partial > 0) thr else -thr) / wkk else 0
      }
      delta <- new_k - beta[k]
      if (delta != 0) {
        r <- r + W11[, k] * delta           # rank-1 residual update
        beta[k] <- new_k
        ad <- abs(delta); if (ad > max_diff) max_diff <- ad
      }
    }
    if (max_diff < tol_inner) break
  }
  beta
}

# --- single graphical-lasso fit at fixed penalty -----------------------------
.glasso_fit <- function(S, rho,
                        max_outer = 1e4, tol_outer = 1e-8,
                        max_inner = 1e4, tol_inner = 1e-10,
                        w_init = NULL, beta_init = NULL) {
  p <- ncol(S)
  W <- if (is.null(w_init)) S else w_init
  dS <- diag(S)
  diag(W) <- dS                                        # penalize.diagonal = FALSE
  Beta <- if (is.null(beta_init)) matrix(0, p, p) else beta_init

  # diag(W) is held fixed at diag(S) (only off-diagonals update), so each
  # column's W11[k, k] is just dS[idx] -- precompute it instead of diag(W11).
  for (outer in seq_len(max_outer)) {
    max_diff <- 0
    for (j in seq_len(p)) {
      idx  <- seq_len(p)[-j]
      W11  <- W[idx, idx, drop = FALSE]
      s12  <- S[idx, j]
      beta <- .glasso_lasso_column(W11, s12, Beta[j, idx], rho,
                                   max_inner, tol_inner, dS[idx])
      w12  <- as.numeric(W11 %*% beta)
      d <- max(abs(w12 - W[idx, j]))
      if (d > max_diff) max_diff <- d
      W[idx, j] <- w12
      W[j, idx] <- w12
      Beta[j, idx] <- beta
    }
    if (max_diff < tol_outer) break
  }

  Theta <- matrix(0, p, p)
  for (j in seq_len(p)) {
    idx   <- seq_len(p)[-j]
    beta  <- Beta[j, idx]
    denom <- W[j, j] - sum(W[idx, j] * beta)
    tjj   <- if (abs(denom) > 1e-12) 1 / denom else 1e6
    Theta[j, j]   <- tjj
    Theta[idx, j] <- -beta * tjj
  }
  Theta <- (Theta + t(Theta)) / 2
  dimnames(Theta) <- dimnames(W) <- dimnames(S)
  list(wi = Theta, w = W, beta = Beta)
}

# --- log-spaced lambda path --------------------------------------------------
.compute_lambda_path <- function(S, nlambda, lambda_min_ratio) {
  lambda_max <- max(abs(S[upper.tri(S)]))
  if (lambda_max <= 0) {
    stop("All off-diagonal correlations are zero; nothing to regularise.",
         call. = FALSE)
  }
  exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
          length.out = nlambda))
}

# Empty-graph fit. When no off-diagonal association exists the optimum is the
# diagonal precision diag(1/S_ii) at lambda 0 (FHT 2008 KKT: W = S, every
# inactive residual 0 <= rho). Returns the same shape as .select_ebic().
.empty_glasso <- function(S, n) {
  wi <- diag(1 / diag(S), nrow = ncol(S))
  dimnames(wi) <- dimnames(S)
  loglik <- (n / 2) * (sum(log(1 / diag(S))) - ncol(S))   # tr(S wi) = p
  ebic <- -2 * loglik                                      # npar = 0
  list(wi = wi, lambda = 0, ebic = ebic, ebic_path = ebic)
}

# Validate the solve engine. "base" is the pure-R FHT-2008 solver (default,
# dependency-free, self-certified); "glasso" delegates each fixed-penalty solve
# to the glasso Fortran package (opt-in, Suggests) for speed and byte-identical
# parity with qgraph/bootnet, at the cost of glasso's looser convergence.
#' @noRd
.check_engine <- function(engine) {
  engine <- match.arg(engine, c("base", "glasso"))
  if (engine == "glasso" && !requireNamespace("glasso", quietly = TRUE)) {
    stop("engine = \"glasso\" needs the 'glasso' package; install it or use ",
         "engine = \"base\".", call. = FALSE)
  }
  engine
}

# One fixed-penalty graphical-lasso solve, dispatched by engine. Returns
# list(wi, w, beta) (beta is NULL for the glasso engine). The base path is
# unchanged, so engine = "base" is byte-identical to the previous behaviour.
#' @noRd
.glasso_solve <- function(S, rho, engine, tol_outer, tol_inner, warm = NULL) {
  if (engine == "base") {
    return(.glasso_fit(S, rho, tol_outer = tol_outer, tol_inner = tol_inner,
                       w_init = warm$w, beta_init = warm$beta))
  }
  g <- glasso::glasso(S, rho, penalize.diagonal = FALSE, thr = tol_outer)
  wi <- (g$wi + t(g$wi)) / 2
  dimnames(wi) <- dimnames(S)
  list(wi = wi, w = g$w, beta = NULL)
}

# --- two-tier EBIC selection: fast scan, then tight refit at the winner ------
.select_ebic <- function(S, lambda_path, n, gamma, engine = "base",
                         scan_tol = 1e-4, refit_tol = 1e-8) {
  p <- ncol(S)
  ebic_vals <- numeric(length(lambda_path))
  w_prev <- NULL; beta_prev <- NULL
  best_idx <- 1L; best_ebic <- Inf; have_best <- FALSE

  for (i in seq_along(lambda_path)) {
    fit <- tryCatch(
      .glasso_solve(S, lambda_path[i], engine, scan_tol, scan_tol,
                    warm = list(w = w_prev, beta = beta_prev)),
      error = function(e) NULL
    )
    if (is.null(fit)) { ebic_vals[i] <- Inf; next }
    w_prev <- fit$w; beta_prev <- fit$beta

    ld <- determinant(fit$wi, logarithm = TRUE)
    if (ld$sign <= 0) { ebic_vals[i] <- Inf; next }
    loglik <- (n / 2) * (as.numeric(ld$modulus) - sum(diag(S %*% fit$wi)))
    npar   <- sum(abs(fit$wi[upper.tri(fit$wi)]) > 1e-10)
    ebic_vals[i] <- -2 * loglik + npar * log(n) + 4 * npar * gamma * log(p)

    if (ebic_vals[i] < best_ebic) {
      best_ebic <- ebic_vals[i]; best_idx <- i; have_best <- TRUE
    }
  }
  if (!have_best) stop("All glasso fits failed; check the input data.",
                       call. = FALSE)

  # base refits to the certified optimum; glasso stays at its own tolerance so
  # the result is exactly what the glasso package returns.
  refit_outer <- if (engine == "glasso") scan_tol else refit_tol
  refit_inner <- if (engine == "glasso") scan_tol else refit_tol * 1e-2
  best_wi <- .glasso_solve(S, lambda_path[best_idx], engine,
                           refit_outer, refit_inner)$wi
  dimnames(best_wi) <- dimnames(S)
  list(wi = best_wi, lambda = lambda_path[best_idx], ebic = best_ebic,
       ebic_path = ebic_vals)
}

# --- precision -> partial correlation ----------------------------------------
.precision_to_pcor <- function(wi) {
  x <- -stats::cov2cor(wi)
  diag(x) <- 0
  (x + t(x)) / 2
}

#' Graphical-lasso stationarity (KKT) residual
#'
#' A dependency-free correctness certificate for a fitted Gaussian graphical
#' model. For the convex objective
#' \deqn{\min_{\Theta \succ 0} -\log\det\Theta + \mathrm{tr}(S\Theta) +
#'       \rho \sum_{i \neq j} |\Theta_{ij}|}
#' (off-diagonal penalty), let \eqn{W = \Theta^{-1}}. The subgradient
#' optimality conditions are \eqn{W_{ii} = S_{ii}};
#' \eqn{W_{ij} - S_{ij} = \rho\,\mathrm{sign}(\Theta_{ij})} where
#' \eqn{\Theta_{ij} \neq 0}; and \eqn{|W_{ij} - S_{ij}| \le \rho} otherwise. By
#' strict convexity, a precision matrix with zero violation is the unique global
#' optimum, so a near-zero return certifies correctness independently of any
#' reference solver.
#'
#' @param theta Precision matrix to test.
#' @param cor_matrix Correlation / covariance the model was fit to.
#' @param rho Scalar penalty.
#' @param active_tol Magnitude above which an off-diagonal entry is "active".
#' @return Maximum absolute stationarity violation (scalar); 0 = exact optimum.
#' @examples
#' S <- 0.5^abs(outer(1:5, 1:5, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 200)
#' glasso_kkt(fit$precision, S, fit$lambda)
#' @export
glasso_kkt <- function(theta, cor_matrix, rho, active_tol = 1e-8) {
  theta <- as.matrix(theta); cor_matrix <- as.matrix(cor_matrix)
  stopifnot(nrow(theta) == ncol(theta), all(dim(theta) == dim(cor_matrix)),
            all(is.finite(theta)), all(is.finite(cor_matrix)),
            is.numeric(rho), length(rho) == 1L)
  # Feasibility is part of the optimality conditions: the objective is defined
  # only on symmetric Theta > 0, so an asymmetric, indefinite, or numerically
  # non-invertible matrix is infeasible, not optimal.
  if (max(abs(theta - t(theta))) > 1e-8) return(Inf)
  if (min(eigen(theta, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    return(Inf)
  W <- tryCatch(solve(theta), error = function(e) NULL)
  if (is.null(W)) return(Inf)
  diag_v <- max(abs(diag(W) - diag(cor_matrix)))
  off <- upper.tri(theta)
  r  <- (W - cor_matrix)[off]
  th <- theta[off]
  active <- abs(th) > active_tol
  v_a <- if (any(active))  max(abs(r[active] - rho * sign(th[active]))) else 0
  v_i <- if (any(!active)) max(pmax(abs(r[!active]) - rho, 0)) else 0
  max(diag_v, v_a, v_i)
}

#' Constrained Gaussian-MRF (graph-restricted MLE) stationarity residual
#'
#' Certificate for an *unregularized* Gaussian graphical model whose precision
#' is constrained to a fixed graph (the estimator behind [ggm_modselect()] and
#' [logo_network()]). The maximum-likelihood / maximum-entropy conditions for a
#' Gaussian Markov random field on a graph \eqn{G} are exact:
#' \eqn{W_{ij} = S_{ij}} for every \eqn{(i,j)} on the graph and on the diagonal
#' (\eqn{W = \Theta^{-1}}), and \eqn{\Theta_{ij} = 0} for every \eqn{(i,j)} not
#' on the graph. A near-zero return certifies the constrained optimum with no
#' reference solver.
#'
#' @param theta Precision matrix to test.
#' @param cor_matrix Correlation / covariance the model was fit to.
#' @param support Logical p x p matrix; `TRUE` where an edge is allowed.
#' @param active_tol Magnitude above which an off-support entry counts as a
#'   nonzero violation.
#' @return Maximum absolute stationarity violation (scalar); 0 = exact optimum.
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' fit <- ggm_modselect(cor_matrix = S, n = 250)
#' ggm_support_kkt(fit$precision, S, fit$support)
#' @export
ggm_support_kkt <- function(theta, cor_matrix, support, active_tol = 1e-8) {
  theta <- as.matrix(theta); cor_matrix <- as.matrix(cor_matrix)
  stopifnot(nrow(theta) == ncol(theta), all(dim(theta) == dim(cor_matrix)),
            all(dim(support) == dim(theta)), is.logical(support),
            all(is.finite(theta)), all(is.finite(cor_matrix)))
  if (max(abs(theta - t(theta))) > 1e-8) return(Inf)
  if (min(eigen(theta, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    return(Inf)
  W <- tryCatch(solve(theta), error = function(e) NULL)
  if (is.null(W)) return(Inf)
  diag_v <- max(abs(diag(W) - diag(cor_matrix)))
  off <- upper.tri(theta)
  sup <- support[off]
  r   <- (W - cor_matrix)[off]
  th  <- theta[off]
  v_match <- if (any(sup))  max(abs(r[sup])) else 0       # W = S on the graph
  v_zero  <- if (any(!sup)) max(abs(th[!sup])) else 0     # Theta = 0 off-graph
  max(diag_v, v_match, v_zero)
}

#' EBIC-regularized Gaussian graphical model (graphical lasso)
#'
#' Selects an L1 penalty by the extended BIC (Foygel & Drton 2010) over a
#' log-spaced path, then refits the chosen penalty to machine precision so the
#' returned network is the certified global optimum of the convex objective.
#' Equivalent in purpose to `qgraph::EBICglasso()` / `bootnet`'s `"EBICglasso"`
#' default, but pure base R and self-certified (see [glasso_kkt()]).
#'
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied.
#' @param cor_matrix Optional correlation matrix; if given, `n` is required and
#'   `data` is ignored.
#' @param n Sample size (required when `cor_matrix` is supplied).
#' @param gamma EBIC hyperparameter. Default 0.5.
#' @param nlambda Number of penalties on the path. Default 100.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest.
#'   Default 0.01.
#' @param threshold Partial correlations with absolute value below this are set
#'   to zero. Default 0.
#' @param cor_method Correlation used when `data` is supplied: `"pearson"`
#'   (default), `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial
#'   for ordinal items, the `qgraph::cor_auto` / `bootnet` default). See
#'   [cor_auto()].
#' @param na_method Missing-data handling when `data` is supplied: `"pairwise"`
#'   (default, pairwise-complete correlations + nearest-PD projection) or
#'   `"listwise"` (drop incomplete rows). Identical for complete data.
#' @param engine Solver: `"base"` (default) is the pure-R, dependency-free,
#'   self-certified solver; `"glasso"` delegates each fixed-penalty solve to the
#'   `glasso` Fortran package (in `Suggests`) for speed and byte-identical
#'   `glasso`/`qgraph` output, at its looser convergence (the reported `$kkt`
#'   then shows glasso's tolerance rather than ~1e-11).
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$weights` is the partial-correlation matrix,
#'   with `$precision`, `$lambda`, `$gamma`, `$cor_matrix`, `$ebic`, `$engine`,
#'   and `$kkt` (the stationarity residual of the returned network).
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 250)
#' fit
#' as.data.frame(fit)
#' @export
ebic_glasso <- function(data = NULL, cor_matrix = NULL, n = NULL,
                        gamma = 0.5, nlambda = 100L, lambda_min_ratio = 0.01,
                        threshold = 0,
                        cor_method = c("pearson", "spearman", "kendall", "auto"),
                        na_method = c("pairwise", "listwise"),
                        engine = c("base", "glasso"), labels = NULL) {
  na_method <- match.arg(na_method)
  cor_method <- match.arg(cor_method)
  engine <- .check_engine(engine)
  if (is.null(cor_matrix)) {
    ci <- .cor_input(data, method = cor_method, na_method = na_method)
    S <- ci$S; n <- ci$n
    if (is.null(labels)) labels <- ci$labels
  } else {
    S <- .check_cor_matrix(cor_matrix)
    if (is.null(n)) stop("`n` is required when `cor_matrix` is supplied.",
                         call. = FALSE)
    stopifnot(is.numeric(n), length(n) == 1L, is.finite(n), n > 0)
    if (is.null(labels)) {
      labels <- colnames(S)
      if (is.null(labels)) labels <- paste0("V", seq_len(ncol(S)))
    }
  }
  stopifnot(is.numeric(gamma), length(gamma) == 1L, gamma >= 0,
            nlambda >= 2L, lambda_min_ratio > 0, lambda_min_ratio < 1)

  dimnames(S) <- list(labels, labels)
  if (max(abs(S[upper.tri(S)])) <= 1e-12) {     # no association: empty optimum
    sel <- .empty_glasso(S, n); lambda_path <- 0
  } else {
    lambda_path <- .compute_lambda_path(S, nlambda, lambda_min_ratio)
    sel <- .select_ebic(S, lambda_path, n, gamma, engine = engine)
  }

  pcor <- .precision_to_pcor(sel$wi)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)

  .new_psychnet(
    graph = pcor, labels = labels, method = "EBICglasso",
    directed = FALSE, n_obs = n,
    extra = list(
      precision = sel$wi, lambda = sel$lambda, gamma = gamma,
      cor_matrix = S, ebic = sel$ebic, ebic_path = sel$ebic_path,
      lambda_path = lambda_path, na_method = na_method, engine = engine,
      kkt = glasso_kkt(sel$wi, S, sel$lambda)
    )
  )
}
