# Automatic correlation detection (polychoric / polyserial / Pearson), the
# base-R counterpart of qgraph::cor_auto. Ordinal items are treated as coarse
# discretizations of latent normals: ordinal-ordinal pairs get a polychoric
# correlation, ordinal-continuous pairs a polyserial correlation, and
# continuous-continuous pairs the ordinary Pearson correlation. This is the
# correlation bootnet/qgraph use by default for Likert data, and it changes the
# estimated network relative to Pearson.
#
# The one primitive base R lacks is the bivariate-normal CDF (stats has only the
# univariate pnorm). We build it by Gauss-Legendre quadrature of the bivariate
# density along the correlation path (Drezner & Wesolowsky 1990):
#   Phi2(h, k; rho) = Phi(h) Phi(k)
#                     + (1/2pi) \int_0^rho (1-t^2)^{-1/2}
#                         exp(-(h^2 - 2 h k t + k^2) / (2(1-t^2))) dt.

# Gauss-Legendre nodes/weights on [-1, 1] via the Golub-Welsch eigenproblem.
#' @noRd
.gauss_legendre <- function(m) {
  i <- seq_len(m - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, m, m)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(x = e$values[ord], w = 2 * (e$vectors[1L, ord])^2)
}

# Standard bivariate-normal CDF P(Z1 <= h, Z2 <= k), correlation rho, with the
# thresholds allowed to be +/-Inf (so it doubles as the rectangle primitive).
#' @noRd
.pbivnorm <- function(h, k, rho, gl) {
  if (h == -Inf || k == -Inf) return(0)
  if (h == Inf && k == Inf) return(1)
  if (h == Inf) return(stats::pnorm(k))
  if (k == Inf) return(stats::pnorm(h))
  if (abs(rho) < 1e-12) return(stats::pnorm(h) * stats::pnorm(k))
  t <- (rho / 2) * (gl$x + 1)                       # map [-1,1] -> [0, rho]
  dens <- exp(-(h * h - 2 * h * k * t + k * k) / (2 * (1 - t * t))) /
    sqrt(1 - t * t)
  stats::pnorm(h) * stats::pnorm(k) +
    (rho / 2) * sum(gl$w * dens) / (2 * pi)
}

# Two-step polychoric correlation of two integer-coded ordinal vectors: fix the
# thresholds at the marginal normal quantiles, then maximise the contingency
# table likelihood over rho (Olsson 1979).
#' @noRd
.polychoric <- function(x, y, gl) {
  tab <- table(x, y)
  n <- sum(tab)
  ax <- c(-Inf, stats::qnorm(cumsum(rowSums(tab)) / n))
  ay <- c(-Inf, stats::qnorm(cumsum(colSums(tab)) / n))
  ax[length(ax)] <- Inf
  ay[length(ay)] <- Inf
  La <- length(ax); Lb <- length(ay)
  tabm <- matrix(tab, nrow(tab), ncol(tab))
  negll <- function(rho) {
    G <- matrix(0, La, Lb)
    for (a in seq_len(La)) for (b in seq_len(Lb)) {
      G[a, b] <- .pbivnorm(ax[a], ay[b], rho, gl)
    }
    P <- G[-1L, -1L] - G[-La, -1L] - G[-1L, -Lb] + G[-La, -Lb]
    -sum(tabm * log(pmax(P, 1e-12)))
  }
  stats::optimize(negll, c(-0.999, 0.999), tol = 1e-6)$minimum
}

# Two-step polyserial correlation of a continuous vector and an integer-coded
# ordinal vector (Olsson, Drasgow & Dorans 1982).
#' @noRd
.polyserial <- function(cont, ord) {
  # A pair carrying no estimable association (the ordinal collapsed to one level,
  # or the continuous partner has zero variance) has correlation 0; without this
  # guard `sum(dnorm(tau))` is 0 and `rho` becomes NaN, which later breaks the
  # nearest-PD projection's eigen decomposition.
  if (length(unique(ord)) < 2L || stats::sd(cont) == 0) return(0)
  n <- length(ord)
  ny <- as.numeric(table(ord))
  tau <- stats::qnorm(cumsum(ny[-length(ny)]) / n)
  r <- stats::cor(cont, as.numeric(ord))
  sy <- stats::sd(as.numeric(ord))
  rho <- r * sy / sum(stats::dnorm(tau))
  max(min(rho, 0.999), -0.999)
}

# A variable is treated as ordinal if it is integer-valued with few levels.
#' @noRd
.is_ordinal <- function(v, max_levels) {
  v <- v[!is.na(v)]
  length(unique(v)) <= max_levels && all(abs(v - round(v)) < 1e-8)
}

# Pairwise auto-correlation matrix, projected to the nearest PD matrix.
#' @noRd
.cor_auto_matrix <- function(mat, max_levels) {
  p <- ncol(mat)
  gl <- .gauss_legendre(60L)
  ord <- vapply(seq_len(p), function(j) .is_ordinal(mat[, j], max_levels),
                logical(1))
  R <- diag(p)
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    ok <- !is.na(mat[, i]) & !is.na(mat[, j])
    a <- mat[ok, i]; b <- mat[ok, j]
    r <- if (ord[i] && ord[j]) .polychoric(a, b, gl)
         else if (ord[i]) .polyserial(b, a)
         else if (ord[j]) .polyserial(a, b)
         else stats::cor(a, b)
    if (!is.finite(r)) r <- 0   # zero-variance / single-level pair -> no association
    R[i, j] <- R[j, i] <- max(min(r, 1), -1)
  }
  dimnames(R) <- list(colnames(mat), colnames(mat))
  .nearest_pd_cor(R)
}

#' Automatic correlation matrix (polychoric / polyserial / Pearson)
#'
#' Detects ordinal variables (integer-valued with at most `ordinal_max_levels`
#' levels) and returns the correlation matrix using a polychoric correlation for
#' ordinal-ordinal pairs, a polyserial correlation for ordinal-continuous pairs,
#' and Pearson otherwise, projected to the nearest positive-definite matrix. The
#' base-R counterpart of `qgraph::cor_auto()`; this is the correlation
#' `bootnet`/`qgraph` use by default for Likert data.
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param ordinal_max_levels Maximum distinct values for a variable to count as
#'   ordinal. Default 7.
#' @param na_method `"pairwise"` (default) or `"listwise"`.
#' @return A correlation matrix with the variable names as dimnames.
#' @examples
#' set.seed(1)
#' z <- matrix(stats::rnorm(300 * 4), 300, 4) %*% chol(0.5^abs(outer(1:4, 1:4, "-")))
#' x <- apply(z, 2, function(col) as.integer(cut(col, 5)))   # 5-level Likert
#' cor_auto(x)
#' @export
cor_auto <- function(data, ordinal_max_levels = 7L,
                     na_method = c("pairwise", "listwise")) {
  na_method <- match.arg(na_method)
  mat <- .as_numeric_matrix(data, drop_na = (na_method == "listwise"))
  .cor_auto_matrix(mat, ordinal_max_levels)
}
