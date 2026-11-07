# Penalized maximum likelihood (lasso) for the panel graphical VAR, with the
# penalty chosen by EBIC -- the panel counterpart of ebic_glasso(), and the
# model psychonetrics fits with panelgvar(estimator = "PML").
#
# Objective (data scale): F(theta) + lambda * sum_{j in P} |theta_j|, with F the
# ML fit function of panelgvar.R and P every temporal (beta, autoregressive
# diagonal included), contemporaneous and between partial-correlation parameter
# -- psychonetrics' default penalized set. The optimizer works on standardized
# variables, where the data-scale penalty becomes a weighted lasso with weight
# to_data[j] (beta[i, j] scales by s_i / s_j; partial correlations are
# scale-free), so the two problems have the same solution.
#
# Solver: FISTA (accelerated proximal gradient) with backtracking that rejects
# any step leaving the model's domain, then an active-set Newton polish on the
# selected support. Certificate: the lasso KKT conditions -- on a non-zero
# penalized parameter the gradient equals -lambda * w * sign; on a zero one
# |gradient| <= lambda * w; unpenalized parameters have zero gradient (or a
# non-negative one on a variance's lower bound).
#
# The panel likelihood is not convex, so the KKT residual certifies a
# stationary point, not the global optimum; the selected lambda is re-solved
# from several starts and their agreement reported, as for the ML fit.

# Penalized objective at a state.
#' @noRd
.pgv_pen_objective <- function(st, theta, pen, w, lambda) {
  st$f + lambda * sum(w[pen] * abs(theta[pen]))
}

# Lasso KKT residual over the free parameters (see the file header).
#' @noRd
.pgv_pen_kkt <- function(g, theta, free, pen, w, lambda, bound) {
  free_pen <- free & pen
  nonzero <- free_pen & theta != 0
  zero <- free_pen & theta == 0
  smooth <- free & !pen & !bound
  max(c(abs(g[smooth]),
        pmax(0, -g[free & bound]),
        abs(g[nonzero] + lambda * w[nonzero] * sign(theta[nonzero])),
        pmax(0, abs(g[zero]) - lambda * w[zero]), 0))
}

# Proximal map: soft-threshold the penalized coordinates (threshold
# step * lambda * w), keep variance scales above zero.
#' @noRd
.pgv_prox <- function(x, step, pen, w, lambda, is_delta) {
  x[pen] <- sign(x[pen]) * pmax(abs(x[pen]) - step * lambda * w[pen], 0)
  x[is_delta] <- pmax(x[is_delta], 1e-8)
  x
}

# One penalized fit at fixed lambda from `theta0` (standardized scale): FISTA
# with backtracking and adaptive restart until the KKT residual is below `tol`
# or `max_iter` iterations. Returns the estimate, its state and its residual.
#' @noRd
.pgv_pen_fista <- function(ctx, theta0, pen, w, lambda, tol = 1e-7,
                           max_iter = 5000L) {
  is_delta <- ctx$lay$table$matrix %in%
    c("delta_zeta_within", "delta_zeta_between")
  fixed <- !ctx$free
  theta0[fixed] <- 0
  x <- .pgv_prox(theta0, 0, pen, w, lambda, is_delta)
  st_x <- .pgv_state(x, ctx)
  if (is.null(st_x)) return(NULL)
  obj_x <- .pgv_pen_objective(st_x, x, pen, w, lambda)
  y <- x
  st_y <- st_x
  momentum <- 1
  step <- 1
  bound <- rep(FALSE, length(x))
  kkt <- Inf
  # Proximal gradient iterations are sequential by construction.
  for (it in seq_len(max_iter)) {
    g <- .pgv_gradient_full(st_y, ctx)
    g[fixed] <- 0
    f_y <- st_y$f
    accepted <- FALSE
    for (bt in seq_len(60L)) {
      z <- .pgv_prox(y - step * g, step, pen, w, lambda, is_delta)
      z[fixed] <- 0
      st_z <- .pgv_state(z, ctx)
      if (!is.null(st_z)) {
        d <- z - y
        if (st_z$f <= f_y + sum(g * d) + sum(d * d) / (2 * step) + 1e-14) {
          accepted <- TRUE
          break
        }
      }
      step <- step / 2
    }
    if (!accepted) break
    obj_z <- .pgv_pen_objective(st_z, z, pen, w, lambda)
    # Adaptive restart (O'Donoghue & Candes 2015): drop momentum whenever the
    # penalized objective goes up.
    if (obj_z > obj_x) {
      momentum <- 1
      y <- x
      st_y <- st_x
      next
    }
    next_momentum <- (1 + sqrt(1 + 4 * momentum^2)) / 2
    y_new <- z + ((momentum - 1) / next_momentum) * (z - x)
    y_new[fixed] <- 0
    y_new <- .pgv_prox(y_new, 0, pen, w, lambda, is_delta)
    st_y_new <- .pgv_state(y_new, ctx)
    if (is.null(st_y_new)) {
      y_new <- z
      st_y_new <- st_z
      next_momentum <- 1
    }
    x <- z
    st_x <- st_z
    obj_x <- obj_z
    y <- y_new
    st_y <- st_y_new
    momentum <- next_momentum
    step <- step * 1.5
    bound <- is_delta & x <= 1e-8 * 1.0001
    kkt <- .pgv_pen_kkt(.pgv_gradient_full(st_x, ctx), x, ctx$free, pen, w,
                        lambda, bound)
    if (kkt < tol) break
  }
  bound <- is_delta & x <= 1e-8 * 1.0001
  kkt <- .pgv_pen_kkt(.pgv_gradient_full(st_x, ctx), x, ctx$free, pen, w,
                      lambda, bound)
  list(theta = x, state = st_x, objective = obj_x, kkt = kkt)
}

# Active-set Newton polish at fixed lambda: on the non-zero penalized and the
# interior unpenalized parameters the penalized objective is smooth, so take
# Fisher-scoring steps (E[Hessian of F] = 2 I) on it, rejecting steps that flip
# a sign or raise the objective; then hand back to FISTA if the zero set is no
# longer optimal. Returns the polished fit.
#' @noRd
.pgv_pen_polish <- function(ctx, fit, pen, w, lambda, tol = 1e-10,
                            max_steps = 30L) {
  is_delta <- ctx$lay$table$matrix %in%
    c("delta_zeta_within", "delta_zeta_between")
  theta <- fit$theta
  st <- fit$state
  obj <- fit$objective
  # Newton steps are sequential by construction.
  for (it in seq_len(max_steps)) {
    bound <- is_delta & theta <= 1e-8 * 1.0001
    active <- which(ctx$free & !bound & (!pen | theta != 0))
    g <- .pgv_gradient_full(st, ctx)
    h <- g[active] + ifelse(pen[active], lambda * w[active] * sign(theta[active]), 0)
    kkt <- .pgv_pen_kkt(g, theta, ctx$free, pen, w, lambda, bound)
    if (kkt < tol) break
    jac <- .pgv_jacobian(st, ctx, active)
    step <- as.vector(MASS_ginv(2 * .pgv_information(st, jac), tol = 1e-12) %*% h)
    accepted <- FALSE
    for (half in 0:30) {
      cand <- theta
      cand[active] <- theta[active] - step / 2^half
      flips <- pen[active] & sign(cand[active]) != sign(theta[active])
      if (any(flips)) next
      cst <- .pgv_state(cand, ctx)
      if (is.null(cst)) next
      cobj <- .pgv_pen_objective(cst, cand, pen, w, lambda)
      if (cobj <= obj + 1e-13 * abs(obj)) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) break
    theta <- cand
    st <- cst
    obj <- cobj
  }
  bound <- is_delta & theta <= 1e-8 * 1.0001
  kkt <- .pgv_pen_kkt(.pgv_gradient_full(st, ctx), theta, ctx$free, pen, w,
                      lambda, bound)
  list(theta = theta, state = st, objective = obj, kkt = kkt)
}

# Penalized fit at one lambda: a short FISTA run identifies the support (the
# zeros of the lasso solution), the active-set Newton polish then converges on
# it quadratically; if the polished point violates a zero's KKT condition the
# support is wrong, so FISTA runs again from there. FISTA alone crawls near a
# boundary (a variance heading to zero): run to its iteration cap it left KKT
# residuals up to ~0.5 on real data.
#' @noRd
.pgv_pen_fit <- function(ctx, theta0, pen, w, lambda, tol = 1e-9,
                         fista_iter = 60L) {
  fit <- .pgv_pen_fista(ctx, theta0, pen, w, lambda, tol = tol,
                        max_iter = fista_iter)
  if (is.null(fit)) return(NULL)
  # Each pass either certifies the support or moves it, hence the loop. A pass
  # that no longer lowers the objective means the fit is crawling along a
  # boundary ridge (no interior stationary point at this lambda): stop there and
  # let the caller see the residual.
  for (pass in seq_len(6L)) {
    before <- fit$objective
    fit <- .pgv_pen_polish(ctx, fit, pen, w, lambda, tol = tol)
    if (fit$kkt < tol) break
    again <- .pgv_pen_fista(ctx, fit$theta, pen, w, lambda, tol = tol,
                            max_iter = fista_iter)
    if (!is.null(again) && again$objective <= fit$objective) fit <- again
    if (before - fit$objective <= 1e-10 * max(1, abs(before))) break
  }
  fit
}

# The EBIC lambda path (psychonetrics::find_penalized_lambda()): lambda_max is
# the largest weighted gradient of a penalized parameter at the fitted baseline
# (all penalized parameters zero), the grid is `nlambda` log-spaced values down
# to lambda_max * lambda_min_ratio with warm starts, and
# EBIC = -2 logL + k log n + 4 k gamma log V, with k the free unpenalized
# parameters plus the penalized ones with |theta| >= 1e-5 on the data scale and
# V the number of observed (wave x variable) columns. The selected lambda is
# re-solved, polished, from the path solution and two further starts.
#' @noRd
.pgv_ebic_path <- function(sctx, ctx, to_data, network_rows, base_theta,
                           extra_starts, gamma, nlambda, lambda_min_ratio,
                           log_jacobian, beta_min = 1e-5) {
  pen <- network_rows & sctx$free
  w <- to_data
  base_state <- .pgv_state(base_theta, sctx)
  g0 <- .pgv_gradient_full(base_state, sctx)
  lambda_max <- max(abs(g0[pen]) / w[pen])
  grid <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
                  length.out = nlambda))
  n <- ctx$n
  k_obs <- nrow(ctx$cov)
  n_unpen <- sum(sctx$free & !pen)
  ebic_of <- function(fit) {
    f_data <- fit$state$f + log_jacobian
    logl <- -n / 2 * (f_data + k_obs * log(2 * pi))
    nonzero <- sum(abs(fit$theta[pen] * to_data[pen]) >= beta_min)
    k <- n_unpen + nonzero
    list(ebic = -2 * logl + k * log(n) + 4 * k * gamma * log(k_obs),
         edges = nonzero)
  }
  path <- vector("list", length(grid))
  current <- base_theta
  # Warm starts: each lambda starts from the previous solution, hence the loop.
  for (i in seq_along(grid)) {
    fit <- .pgv_pen_fit(sctx, current, pen, w, grid[i], tol = 1e-8)
    if (is.null(fit)) next
    crit <- ebic_of(fit)
    path[[i]] <- list(lambda = grid[i], ebic = crit$ebic, edges = crit$edges,
                      kkt = fit$kkt, theta = fit$theta)
    current <- fit$theta
  }
  ok <- !vapply(path, is.null, logical(1))
  if (!any(ok)) {
    stop(errorCondition(
      "panel_gvar(): the penalized fit failed at every lambda on the path.",
      class = c("psychnet_no_converge", "psychnet_error"), call = NULL))
  }
  path <- path[ok]
  # Solve one lambda to certification from the path solution and the extra
  # starts, keep the best, apply beta_min, and score it with EBIC.
  solve_final <- function(i) {
    starts <- c(list(path = path[[i]]$theta), extra_starts)
    # The path solution gets the full budget (it is the certified fit); the
    # other starts only need to reach their own optimum to test for competing
    # ones.
    fits <- Map(function(s0, nm)
      .pgv_pen_fit(sctx, s0, pen, w, path[[i]]$lambda, tol = 1e-10,
                   fista_iter = if (identical(nm, "path")) 3000L else 400L),
      starts, names(starts))
    obj <- vapply(fits, function(f) if (is.null(f)) Inf else f$objective,
                  numeric(1))
    best <- fits[[which.min(obj)]]
    # beta_min: a penalized estimate below 1e-5 on the data scale is zero.
    tiny <- pen & abs(best$theta * to_data) < beta_min
    best$theta[tiny] <- 0
    best$state <- .pgv_state(best$theta, sctx)
    # Thresholding changes the point being scored. Report its actual residual,
    # rather than the certificate of the pre-threshold optimizer solution.
    # With very different variable units a tiny data-scale beta can still
    # matter on the standardized scale, or even make the rounded model invalid.
    if (is.null(best$state)) {
      best$kkt <- Inf
      best$objective <- Inf
      crit <- list(ebic = Inf, edges = sum(pen & best$theta != 0))
    } else {
      best$kkt <- .pgv_pen_kkt(.pgv_gradient_full(best$state, sctx),
                               best$theta, sctx$free, pen, w,
                               path[[i]]$lambda, .pgv_at_bound(best$state, sctx))
      best$objective <- .pgv_pen_objective(best$state, best$theta, pen, w,
                                           path[[i]]$lambda)
      crit <- ebic_of(best)
    }
    list(fit = best, ebic = crit$ebic, edges = crit$edges,
         starts = .pgv_starts_table(
           data.frame(start = names(starts), objective = obj,
                      stringsAsFactors = FALSE), log_jacobian))
  }
  ebic <- vapply(path, `[[`, numeric(1), "ebic")
  edges <- vapply(path, `[[`, numeric(1), "edges")
  kkt <- vapply(path, `[[`, numeric(1), "kkt")
  final <- rep(FALSE, length(path))
  solved <- vector("list", length(path))
  # Selection is settled only when the EBIC minimum belongs to a finalized
  # candidate: a final solve can find a different (e.g. denser)
  # optimum than the short path fit, change that lambda's EBIC, and hand the
  # minimum to another candidate, which is then solved in turn. Sequential by
  # construction; at most `length(path)` passes.
  repeat {
    if (!any(is.finite(ebic))) {
      stop(errorCondition(
        "panel_gvar(): no feasible penalized model remained after coefficient thresholding.",
        class = c("psychnet_no_converge", "psychnet_error"), call = NULL))
    }
    i <- which.min(ebic)
    if (final[i]) break
    solved[[i]] <- solve_final(i)
    ebic[i] <- solved[[i]]$ebic
    edges[i] <- solved[[i]]$edges
    kkt[i] <- solved[[i]]$fit$kkt
    final[i] <- TRUE
  }
  chosen <- solved[[i]]
  table <- data.frame(
    lambda = vapply(path, `[[`, numeric(1), "lambda"),
    ebic = ebic, edges = edges, kkt = kkt, final = final,
    row.names = NULL)
  list(theta = chosen$fit$theta, lambda = path[[i]]$lambda,
       kkt = chosen$fit$kkt, path = table, pen = pen,
       starts = chosen$starts)
}
