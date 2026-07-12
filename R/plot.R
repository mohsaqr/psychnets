# Base-R visualisations for psychnet result objects.
#
# Every plot here is drawn with the base `graphics` engine (no ggplot2, no grid,
# no new hard dependency) so the package keeps its base-R-only contract. The
# chart grammar mirrors the bootnet / Nestimate vocabulary: sorted edge-weight
# confidence intervals, centrality CIs, case-dropping stability curves, the
# bootstrapped difference "significance box" matrix, and the network-comparison
# permutation nulls.

# ---- shared palette + small drawing helpers --------------------------------

.psn_pal <- list(
  sig    = "#2C6E8A",  # significant / emphasis
  nonsig = "#9AA5AD",  # non-significant / muted
  pos    = "#C0392B",  # positive difference (row > col)
  neg    = "#2C6E8A",  # negative difference (col > row)
  ref    = "#B5B5B5",  # reference lines
  fill   = "#2A6FBB"   # gradient high end
)

# A clean, well-separated qualitative palette for series (one per measure).
.psn_series <- function(k) {
  base <- c("#2C6E8A", "#C0392B", "#E1A52D", "#3A8C5F", "#7D5BA6", "#11707F",
            "#B5651D")
  if (k <= length(base)) base[seq_len(k)] else grDevices::hcl.colors(k, "Dynamic")
}

# Stars for a p-value, bootnet-style.
.psn_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", ""))))
}

# Horizontal point-range panel: one row per item, segment = interval, dot =
# point estimate, optional significance colouring, optional reference line.
.psn_pointrange <- function(value, lower, upper, labels, significant = NULL,
                            main = "", xlab = "", ref = 0, sort = TRUE) {
  n <- length(value)
  ord <- if (isTRUE(sort)) order(value) else seq_len(n)
  value <- value[ord]; lower <- lower[ord]; upper <- upper[ord]
  labels <- labels[ord]
  cols <- rep(.psn_pal$nonsig, n)
  if (!is.null(significant)) cols[significant[ord] %in% TRUE] <- .psn_pal$sig
  xr <- range(c(lower, upper, value, ref), na.rm = TRUE)
  y <- seq_len(n)
  graphics::plot(NA, xlim = xr, ylim = c(0.5, n + 0.5), yaxt = "n",
                 xlab = xlab, ylab = "", main = main, bty = "n")
  if (length(ref)) graphics::abline(v = ref, col = .psn_pal$ref, lty = 2)
  graphics::segments(lower, y, upper, y, col = cols, lwd = 2)
  graphics::points(value, y, pch = 19, col = cols,
                   cex = min(1, 0.9 * 30 / max(n, 1)) + 0.4)
  cx <- if (n > 40) 0.55 else if (n > 20) 0.7 else 0.85
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE,
                 cex.axis = cx)
  invisible(NULL)
}

# Wider left margin so long item labels fit; widest label drives the width.
.psn_left_margin <- function(labels) {
  w <- max(graphics::strwidth(labels, units = "inches", cex = 0.8))
  min(max(4, w / 0.2 * 4 + 1), 14)
}

# The bootstrapped "significance box" matrix (bootnet differenceTest grammar):
# items on both axes ordered by their observed value, the diagonal carries the
# observed value shaded by magnitude, an off-diagonal cell is filled when the
# pair differs significantly (red if the row item is larger, blue if smaller)
# and left faint otherwise, with significance stars.
.psn_box_matrix <- function(items, observed, df, main = "", show_p = TRUE) {
  # Order axes by observed value (high to low), as documented; df is keyed by
  # item name so the index map below follows the reordering automatically.
  o <- order(observed, decreasing = TRUE)
  items <- items[o]; observed <- observed[o]
  k <- length(items)
  idx <- stats::setNames(seq_len(k), items)
  D <- matrix(NA_real_, k, k); P <- matrix(NA_real_, k, k)
  S <- matrix(FALSE, k, k)
  ra <- idx[df$item1]; rb <- idx[df$item2]
  for (r in seq_len(nrow(df))) {
    a <- ra[r]; b <- rb[r]
    if (is.na(a) || is.na(b)) next
    D[a, b] <- df$obs_diff[r]; D[b, a] <- -df$obs_diff[r]
    P[a, b] <- P[b, a] <- df$p_value[r]
    S[a, b] <- S[b, a] <- isTRUE(df$significant[r])
  }
  dmax <- max(abs(D), na.rm = TRUE); if (!is.finite(dmax) || dmax == 0) dmax <- 1
  omax <- max(abs(observed), na.rm = TRUE); if (!is.finite(omax) || omax == 0) omax <- 1

  graphics::plot(NA, xlim = c(0.5, k + 0.5), ylim = c(0.5, k + 0.5),
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = main,
                 bty = "n", asp = 1)
  cx <- if (k > 25) 0.5 else if (k > 12) 0.65 else 0.8
  graphics::axis(1, at = seq_len(k), labels = items, las = 2, tick = FALSE,
                 cex.axis = cx)
  graphics::axis(2, at = rev(seq_len(k)), labels = items, las = 1,
                 tick = FALSE, cex.axis = cx)
  for (ci in seq_len(k)) {
    for (ri in seq_len(k)) {
      yy <- k - ri + 1
      if (ci == ri) {                                   # diagonal: observed
        shade <- grDevices::adjustcolor(.psn_pal$fill,
                                        alpha.f = 0.15 + 0.7 * abs(observed[ri]) / omax)
        graphics::rect(ci - 0.5, yy - 0.5, ci + 0.5, yy + 0.5,
                       col = shade, border = "white")
        graphics::text(ci, yy, formatC(observed[ri], format = "f", digits = 2),
                       cex = cx * 0.9, col = "#1A1A1A")
        next
      }
      d <- D[ri, ci]
      if (is.na(d)) next
      base_col <- if (d > 0) .psn_pal$pos else .psn_pal$neg
      a <- if (isTRUE(S[ri, ci])) 0.25 + 0.65 * abs(d) / dmax else 0.08
      graphics::rect(ci - 0.5, yy - 0.5, ci + 0.5, yy + 0.5,
                     col = grDevices::adjustcolor(base_col, alpha.f = a),
                     border = "white")
      if (show_p && isTRUE(S[ri, ci])) {
        st <- .psn_stars(P[ri, ci])
        if (nzchar(st))
          graphics::text(ci, yy, st, cex = cx, col = "white", font = 2)
      }
    }
  }
  invisible(NULL)
}

# ---- centrality ------------------------------------------------------------

#' Plot node centralities
#'
#' Draws the centrality table returned by [net_centralities()]. `type = "bar"`
#' gives one sorted horizontal lollipop panel per measure; `type = "line"` gives
#' the qgraph/bootnet centrality plot — one faceted panel per measure, nodes on a
#' shared vertical axis, a line-and-marker series within each panel.
#'
#' @param x A `psychnet_centrality` data frame from [net_centralities()].
#' @param type `"bar"` (default) for one sorted lollipop panel per measure, or
#'   `"line"` for the faceted qgraph-style centrality plot.
#' @param scale For `type = "line"`, the per-measure transform: `"raw"`
#'   (default — each faceted panel keeps its own axis, so no rescaling is
#'   needed), `"z"` (z-score per measure, centred at 0), or `"relative"` (each
#'   measure min-max scaled to \[0, 1\]). `type = "bar"` always shows raw values.
#' @param measures Which measure columns to draw. Default: all of them.
#' @param ... Unused.
#' @return `x`, invisibly. Called for the plot it draws.
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 300)
#' plot(net_centralities(fit))
#' plot(net_centralities(fit), type = "line")
#' @export
plot.psychnet_centrality <- function(x, type = c("bar", "line"),
                                     scale = c("raw", "z", "relative"),
                                     measures = NULL, ...) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  ms <- if (is.null(measures)) setdiff(names(x), "node") else measures
  nodes <- x$node
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))

  if (type == "line") {
    # One faceted panel per measure (qgraph centralityPlot layout): nodes share
    # a single vertical order across panels, each panel keeps its OWN x-axis, so
    # the measures never overlap and no cross-measure rescaling is forced.
    val <- vapply(ms, function(m) {
      v <- x[[m]]
      switch(scale,
        raw = v,
        z = { s <- stats::sd(v); if (is.na(s) || s == 0) v - mean(v) else (v - mean(v)) / s },
        relative = { rng <- range(v); d <- diff(rng); if (d == 0) v * 0 else (v - rng[1]) / d })
    }, numeric(length(nodes)))
    val <- matrix(val, nrow = length(nodes), dimnames = list(nodes, ms))
    # Stable node order shared by every panel: overall importance (mean z-score).
    zord <- vapply(ms, function(m) {
      v <- x[[m]]; s <- stats::sd(v)
      if (is.na(s) || s == 0) v - mean(v) else (v - mean(v)) / s
    }, numeric(length(nodes)))
    ord <- order(rowMeans(matrix(zord, nrow = length(nodes))))
    val <- val[ord, , drop = FALSE]; lab <- nodes[ord]
    np <- length(lab); y <- seq_len(np)
    cols <- .psn_series(length(ms))
    xlab <- switch(scale, raw = "centrality",
                   z = "standardised centrality (z)",
                   relative = "relative centrality (0-1)")
    lm <- .psn_left_margin(lab)
    graphics::par(mfrow = c(1, length(ms)), oma = c(3.5, lm, 1.5, 1))
    for (j in seq_along(ms)) {
      graphics::par(mar = c(0.6, 0.4, 2.2, 1.2))
      vv <- val[, j]
      graphics::plot(NA, xlim = grDevices::extendrange(vv, f = 0.08),
                     ylim = c(0.5, np + 0.5), yaxt = "n", bty = "n",
                     xlab = "", ylab = "", main = ms[j],
                     col.main = cols[j], font.main = 2, cex.main = 0.95)
      graphics::segments(graphics::par("usr")[1], y, graphics::par("usr")[2], y,
                         col = "#EDEDED", lwd = 1)
      if (scale == "z") graphics::abline(v = 0, col = .psn_pal$ref, lty = 2)
      graphics::lines(vv, y, col = cols[j], lwd = 2)
      graphics::points(vv, y, col = cols[j], pch = 19,
                       cex = if (np > 20) 0.8 else 1.2)
      if (j == 1L)
        graphics::axis(2, at = y, labels = lab, las = 1, tick = FALSE,
                       cex.axis = if (np > 20) 0.65 else 0.9, xpd = NA)
    }
    graphics::mtext(xlab, side = 1, outer = TRUE, line = 1.6, cex = 0.9)
    return(invisible(x))
  }

  graphics::par(mfrow = c(1, length(ms)),
                mar = c(4, .psn_left_margin(nodes), 3, 1), oma = c(0, 0, 0, 0))
  for (m in ms) {
    .psn_pointrange(x[[m]], x[[m]], x[[m]], nodes, main = m, xlab = "",
                    ref = numeric(0), sort = TRUE)
  }
  invisible(x)
}

# ---- bootstrap -------------------------------------------------------------

#' Plot a network bootstrap
#'
#' Visualises a [net_boot()] result. `type = "edges"` (default) draws the
#' bootstrapped edge-weight confidence intervals sorted by the observed weight
#' (bootnet's edge-accuracy plot); `type = "centrality"` draws the bootstrapped
#' centrality intervals, one sorted panel per measure; `type = "edge_diff"` and
#' `type = "centrality_diff"` draw the bootstrapped difference "significance box"
#' matrix for edges or for one centrality; `type = "predictability"` draws the
#' node predictability intervals (only when [net_boot()] was run with
#' `predictability = TRUE`).
#'
#' @param x A `psychnet_bootstrap` object from [net_boot()].
#' @param type One of `"edges"`, `"centrality"`, `"edge_diff"`,
#'   `"centrality_diff"`, `"predictability"`.
#' @param measure For `"centrality"`/`"centrality_diff"`, which measure(s) to
#'   draw. Default: all bootstrapped measures (`"centrality"`) or the first
#'   (`"centrality_diff"`).
#' @param ... Unused.
#' @return `x`, invisibly. Called for the plot it draws.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' bs <- net_boot(x, n_boot = 50, cores = 1)   # n_boot >= 1000 for real use
#' plot(bs)                       # edge-weight CIs
#' plot(bs, type = "centrality")  # centrality CIs
#' @export
plot.psychnet_bootstrap <- function(x, type = c("edges", "centrality",
                                                 "edge_diff", "centrality_diff",
                                                 "predictability"),
                                    measure = NULL, ...) {
  type <- match.arg(type)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))

  if (type == "edges") {
    ed <- x$edges
    graphics::par(mar = c(4, .psn_left_margin(x$edge_labels), 3, 1))
    .psn_pointrange(ed$observed, ed$lower, ed$upper, x$edge_labels,
                    significant = ed$significant,
                    main = sprintf("Edge weights (%.0f%% bootstrap CI)", 100 * x$ci),
                    xlab = "edge weight")
    return(invisible(x))
  }

  if (type == "centrality") {
    ms <- if (is.null(measure)) x$measures else measure
    graphics::par(mfrow = c(1, length(ms)),
                  mar = c(4, .psn_left_margin(x$node_labels), 3, 1))
    for (m in ms) {
      .psn_pointrange(x$centrality[[m]],
                      x$centrality[[paste0(m, "_lower")]],
                      x$centrality[[paste0(m, "_upper")]],
                      x$node_labels, main = m, xlab = "", ref = numeric(0))
    }
    return(invisible(x))
  }

  if (type == "predictability") {
    if (is.null(x$predictability))
      stop("No predictability stored. Re-run net_boot(predictability = TRUE).",
           call. = FALSE)
    pr <- x$predictability
    graphics::par(mar = c(4, .psn_left_margin(pr$node), 3, 1))
    .psn_pointrange(pr$value, pr$lower, pr$upper, pr$node,
                    main = "Predictability (R^2)", xlab = "R^2",
                    ref = numeric(0))
    return(invisible(x))
  }

  # difference box matrices, computed from the retained draws.
  if (type == "edge_diff") {
    df <- difference_test(x, type = "edge")
    graphics::par(mar = c(7, 7, 3, 1))
    .psn_box_matrix(x$edge_labels, x$edges$observed, df,
                    main = "Edge-weight differences")
    return(invisible(x))
  }
  # centrality_diff
  m <- if (is.null(measure)) x$measures[1] else measure[1]
  df <- difference_test(x, type = m)
  graphics::par(mar = c(7, 7, 3, 1))
  .psn_box_matrix(x$node_labels, x$centrality[[m]], df,
                  main = sprintf("%s differences", m))
  invisible(x)
}

#' Plot a bootstrapped difference test
#'
#' Draws the bootnet-style "significance box" matrix for the pairwise difference
#' test returned by [difference_test()]: items on both axes ordered by their
#' observed value, the diagonal showing the observed value, and each off-diagonal
#' cell filled when that pair differs significantly (red when the row item is
#' larger, blue when smaller). `style = "forest"` instead draws a forest plot:
#' one row per pair, the bootstrapped difference as a point with its confidence
#' interval, a reference line at zero, and significant pairs (interval excluding
#' zero) emphasised.
#'
#' @param x A `psychnet_difference` data frame from [difference_test()].
#' @param style `"box"` (default) for the significance-box matrix, or
#'   `"forest"` for a forest plot of the pairwise differences with their CIs.
#' @param ... Unused.
#' @return `x`, invisibly. Called for the plot it draws.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' bs <- net_boot(x, n_boot = 50, cores = 1)   # n_boot >= 1000 for real use
#' plot(difference_test(bs, type = "strength"))                   # box matrix
#' plot(difference_test(bs, type = "strength"), style = "forest") # forest plot
#' @export
plot.psychnet_difference <- function(x, style = c("box", "forest"), ...) {
  style <- match.arg(style)
  obs <- attr(x, "observed")
  ty  <- attr(x, "diff_type")
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  if (style == "forest") {
    labels <- paste(x$item1, x$item2, sep = " vs ")
    graphics::par(mar = c(4, .psn_left_margin(labels), 3, 1))
    .psn_pointrange(x$obs_diff, x$lower, x$upper, labels,
                    significant = x$significant,
                    main = sprintf("%s pairwise differences",
                                   if (is.null(ty)) "" else ty),
                    xlab = "difference (item1 - item2)", ref = 0)
    return(invisible(x))
  }
  graphics::par(mar = c(7, 7, 3, 1))
  .psn_box_matrix(names(obs), unname(obs), x,
                  main = sprintf("%s differences", if (is.null(ty)) "" else ty))
  invisible(x)
}

# ---- stability -------------------------------------------------------------

#' Plot centrality stability (case-dropping)
#'
#' Draws the case-dropping stability curves from [net_stability()]: mean rank
#' correlation with the full-sample centrality against the proportion of cases
#' dropped, one line per measure, with a +/- 1 SD band, the acceptance threshold,
#' and the CS-coefficient annotated in the legend.
#'
#' @param x A `psychnet_stability` object from [net_stability()].
#' @param ... Unused.
#' @return `x`, invisibly. Called for the plot it draws.
#' @examples
#' set.seed(1)
#' d <- matrix(stats::rnorm(200 * 5), 200, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' s <- net_stability(d, drop_prop = c(0.3, 0.6), iter = 10)
#' plot(s)
#' @export
plot.psychnet_stability <- function(x, ...) {
  tab <- x$table
  ms <- unique(tab$measure)
  cols <- grDevices::hcl.colors(max(length(ms), 2), "Dynamic")[seq_along(ms)]
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mar = c(4.5, 4.5, 3, 1))
  graphics::plot(NA, xlim = range(tab$drop_prop), ylim = c(0, 1),
                 xlab = "proportion of cases dropped",
                 ylab = "mean correlation with full sample",
                 main = "Centrality stability", bty = "n")
  graphics::abline(h = x$threshold, col = .psn_pal$ref, lty = 2)
  for (i in seq_along(ms)) {
    d <- tab[tab$measure == ms[i], , drop = FALSE]
    d <- d[order(d$drop_prop), , drop = FALSE]
    if (!is.null(d$sd_cor) && any(is.finite(d$sd_cor))) {
      lo <- pmax(0, d$mean_cor - d$sd_cor); hi <- pmin(1, d$mean_cor + d$sd_cor)
      graphics::polygon(c(d$drop_prop, rev(d$drop_prop)), c(lo, rev(hi)),
                        col = grDevices::adjustcolor(cols[i], alpha.f = 0.18),
                        border = NA)
    }
    graphics::lines(d$drop_prop, d$mean_cor, col = cols[i], lwd = 2)
    graphics::points(d$drop_prop, d$mean_cor, col = cols[i], pch = 19, cex = 0.8)
  }
  leg <- sprintf("%s (CS = %.2f)", ms, x$cs[ms])
  graphics::legend("bottomleft", legend = leg, col = cols, lwd = 2, pch = 19,
                   bty = "n", cex = 0.8)
  invisible(x)
}

# ---- network comparison test ----------------------------------------------

#' Plot a Network Comparison Test
#'
#' Visualises a [net_compare()] result. `type = "strength"` (default) and
#' `type = "structure"` draw the permutation null distribution for the global
#' strength invariance (M) and the maximum edge-difference (S) statistics, with
#' the observed value and p-value marked; `type = "edges"` draws the observed
#' per-edge absolute differences, coloured by whether each edge differs
#' significantly.
#'
#' @param x A `psychnet_nct` object from [net_compare()].
#' @param type One of `"strength"`, `"structure"`, `"edges"`.
#' @param alpha Significance level for colouring per-edge differences. Default
#'   `0.05`.
#' @param ... Unused.
#' @return `x`, invisibly. Called for the plot it draws.
#' @examples
#' set.seed(1)
#' mk <- function(s) { set.seed(s)
#'   matrix(stats::rnorm(120 * 4), 120, 4) %*% chol(0.3^abs(outer(1:4, 1:4, "-"))) }
#' cmp <- net_compare(mk(1), mk(2), iter = 25)
#' plot(cmp)                  # global strength permutation null
#' plot(cmp, type = "edges")  # per-edge differences
#' @export
plot.psychnet_nct <- function(x, type = c("strength", "structure", "edges"),
                              alpha = 0.05, ...) {
  type <- match.arg(type)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))

  if (type %in% c("strength", "structure")) {
    stat <- if (type == "strength") x$M else x$S
    lab  <- if (type == "strength") "global strength difference (M)"
            else "maximum edge difference (S)"
    graphics::par(mar = c(4.5, 4.5, 3, 1))
    graphics::hist(stat$perm, breaks = 30, col = .psn_pal$nonsig,
                   border = "white", xlab = lab,
                   xlim = range(c(stat$perm, stat$observed)),
                   main = sprintf("NCT permutation null (p = %.3f)", stat$p_value))
    graphics::abline(v = stat$observed, col = .psn_pal$pos, lwd = 3)
    graphics::legend("topright", legend = "observed", col = .psn_pal$pos,
                     lwd = 3, bty = "n", cex = 0.85)
    return(invisible(x))
  }

  en <- x$E$edge_names
  labels <- paste(en$from, en$to, sep = "--")
  sig <- x$E$p_value < alpha
  cols <- ifelse(sig, .psn_pal$sig, .psn_pal$nonsig)
  graphics::par(mar = c(4, .psn_left_margin(labels), 3, 1))
  ord <- order(x$E$observed)
  graphics::barplot(x$E$observed[ord], names.arg = labels[ord], horiz = TRUE,
                    las = 1, col = cols[ord], border = NA,
                    xlab = "|edge difference|",
                    main = sprintf("Per-edge differences (%d sig.)", sum(sig)),
                    cex.names = if (length(labels) > 20) 0.6 else 0.8)
  invisible(x)
}

# ---- network ---------------------------------------------------------------

#' Plot a psychnet network
#'
#' Renders the estimated network with `cograph::splot()` (a Suggested package);
#' `psychnet` objects inherit from `cograph_network` for exactly this purpose.
#' For the bootstrap / centrality / difference / stability diagnostics use the
#' dedicated `plot()` methods for those result objects instead.
#'
#' @param x A `psychnet` object.
#' @param ... Passed to [cograph::splot()].
#' @return The value of [cograph::splot()], invisibly.
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 300)
#' if (requireNamespace("cograph", quietly = TRUE)) {
#'   plot(fit)
#' }
#' @export
plot.psychnet <- function(x, ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop("Network plotting needs the 'cograph' package. ",
         "Install it, or use net_centralities()/net_boot() with their plot ",
         "methods for diagnostics.", call. = FALSE)
  }
  invisible(cograph::splot(x, ...))
}
