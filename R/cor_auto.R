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
  if (n < 2L || nrow(tab) < 2L || ncol(tab) < 2L) return(0)
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

# Polyserial correlation of a continuous vector and an integer-coded ordinal
# vector by the Olsson, Drasgow & Dorans (1982) ad-hoc (moment) estimator:
# rho = r_{x,y} * sd(y) / sum(dnorm(tau)), thresholds tau at the marginal normal
# quantiles. It is closed-form, NOT the two-step maximum-likelihood estimator
# qgraph::cor_auto (via lavaan) uses; the two differ by up to ~0.006 in
# simulations (documented in cor_auto()'s details).
#' @noRd
.polyserial <- function(cont, ord) {
  # A pair carrying no estimable association (the ordinal collapsed to one level,
  # or the continuous partner has zero variance) has correlation 0; without this
  # guard `sum(dnorm(tau))` is 0 and `rho` becomes NaN, which later breaks the
  # nearest-PD projection's eigen decomposition.
  if (length(unique(ord)) < 2L || stats::sd(cont) == 0) return(0)
  # The moment denominator assumes consecutive ordered category indices.
  # Arbitrary numeric labels (e.g. 10, 20, 30) carry no spacing information.
  ord <- match(ord, sort(unique(ord)))
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
    if (length(a) < 2L) {
      .psn_stop(sprintf(
        "Automatic correlation is undefined: variables '%s' and '%s' have fewer than two joint observations.",
        colnames(mat)[i], colnames(mat)[j]), call. = FALSE)
    }
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

#' Compute correlations suited to ordinal and continuous variables
#'
#' @description
#' `cor_auto()` computes a correlation matrix suited to a mix of ordinal and
#' continuous variables: polychoric correlations between ordinal items,
#' polyserial correlations between an ordinal item and a continuous variable,
#' and Pearson correlations otherwise. The polychoric correlation assumes that
#' each ordinal item is a coarse categorization of a normally distributed
#' latent variable and estimates the correlation between the latent variables
#' (Olsson, 1979); the polyserial correlation does the same for an ordinal
#' item paired with a continuous variable (Olsson et al., 1982).
#'
#' Questionnaire items are usually recorded on a few ordered categories, such
#' as a five-point Likert scale, and treating these codes as continuous and
#' computing Pearson correlations understates the association between the
#' underlying traits. `cor_auto()` decides for each variable whether it is
#' ordinal (integer values with at most `ordinal_max_levels` distinct values)
#' and computes a polychoric correlation for each pair of ordinal variables, a
#' polyserial correlation for each ordinal-continuous pair and a Pearson
#' correlation otherwise.
#'
#' It returns the resulting matrix adjusted to the nearest positive-definite
#' correlation matrix, so that it can be used by the network estimators.
#'
#' Network estimators for ordinal data commonly start from this correlation
#' matrix. The estimators use the same correlations when given
#' `cor_method = "auto"`, for example [ebic_glasso()], [ggm_modselect()] and
#' [cor_network()]. The function corresponds to `qgraph::cor_auto()`.
#'
#' @details
#' The polychoric correlation is the two-step estimator of Olsson (1979): the
#' thresholds are fixed at the normal quantiles of the cumulative marginal
#' proportions, and the likelihood of the contingency table is then maximized
#' over the correlation. The polyserial correlation is the ad hoc estimator of
#' Olsson, Drasgow and Dorans (1982),
#' \eqn{r_{xy}\,\sigma_y / \sum_k \phi(\tau_k)}, rather than the two-step
#' maximum-likelihood estimator that `qgraph::cor_auto()` obtains through
#' `lavaan`. Both are consistent, but on finite samples they differ slightly
#' (by up to about 0.006 in simulations), so a matrix that mixes ordinal and
#' continuous variables does not reproduce `qgraph::cor_auto()` exactly.
#'
#' An ordered factor in a data frame enters as its integer level codes and is
#' ordinal whenever it has at most `ordinal_max_levels` levels. Unordered
#' factors, character and logical columns are dropped with a message, and
#' constant columns with a warning. A pair with no estimable association,
#' such as an item observed at a single level, has correlation 0.
#'
#' @param data A numeric data frame or matrix with one row per observation.
#' @param ordinal_max_levels Largest number of distinct values for an
#'   integer-valued variable to be treated as ordinal: a whole number >= 2.
#'   Default `7`.
#' @param na_method Handling of missing values: `"pairwise"` (default) uses
#'   the observations available for each pair; `"listwise"` drops incomplete
#'   rows first.
#' @return A symmetric correlation matrix with the variable names as row and
#'   column names, one row and column per variable kept.
#' @section Conditions:
#'   `psychnet_bad_input` when `ordinal_max_levels` or `na_method` is invalid,
#'   `data` has no numeric columns or fewer than two usable variables, or a
#'   pair of variables has fewer than two joint observations;
#'   `psychnet_constant_variable` (warning) when a constant or all-missing
#'   column is dropped. A `psychnet_dropped_columns` message names
#'   non-numeric columns that were dropped.
#' @seealso [ebic_glasso()], [cor_network()] and the other estimators, which
#'   accept `cor_method = "auto"`.
#' @references Olsson, U. (1979). Maximum likelihood estimation of the polychoric
#'   correlation coefficient. *Psychometrika*, 44(4), 443-460.
#'   \doi{10.1007/BF02296207}
#'
#'   Olsson, U., Drasgow, F., & Dorans, N. J. (1982). The polyserial correlation
#'   coefficient. *Psychometrika*, 47(3), 337-347. \doi{10.1007/BF02294164}
#' @examples
#' set.seed(1)
#' z <- matrix(stats::rnorm(300 * 4), 300, 4) %*% chol(0.5^abs(outer(1:4, 1:4, "-")))
#' x <- apply(z, 2, function(col) as.integer(cut(col, 5)))   # 5-level Likert
#' cor_auto(x)
#' @export
cor_auto <- function(data, ordinal_max_levels = 7L,
                     na_method = c("pairwise", "listwise")) {
  na_method <- .psn_match_arg(na_method)
  .psn_check_number(ordinal_max_levels, "ordinal_max_levels", lower = 2,
                    whole = TRUE)
  mat <- .as_numeric_matrix(data, drop_na = (na_method == "listwise"))
  .cor_auto_matrix(mat, ordinal_max_levels)
}
