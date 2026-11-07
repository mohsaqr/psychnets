# Base-R visualisations for psychnet result objects.
#
# Every plot here is drawn with the base `graphics` engine (no ggplot2, no grid,
# no new hard dependency) so the package keeps its base-R-only contract. The
# chart grammar mirrors the bootnet / Nestimate vocabulary: sorted edge-weight
# confidence intervals, centrality CIs, case-dropping stability curves, the
# bootstrapped difference "significance box" matrix, and the network-comparison
# permutation nulls.

# ---- shared palette + small drawing helpers --------------------------------

# Okabe-Ito (colour-blind safe) is the only qualitative palette; the diverging
# pair (low = red, high = blue, white midpoint) codes a signed difference. No
# distinction is carried by colour alone: significance is also drawn as filled
# vs open points, series also differ in line type and point shape, and the
# sign of a difference cell is also printed as a +/- glyph.
.psn_okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                    "#D55E00", "#CC79A7", "#999999", "#000000")

.psn_pal <- list(
  sig    = "#0072B2",  # significant / emphasis (Okabe-Ito blue)
  nonsig = "#999999",  # non-significant / muted (Okabe-Ito grey)
  pos    = "#4A6FE3",  # positive difference (row > col): diverging high
  neg    = "#D33F6A",  # negative difference (col > row): diverging low
  ref    = "#999999",  # reference lines
  fill   = "#0072B2",  # magnitude shading of observed values
  mark   = "#D55E00"   # an observed statistic against its null
)

# Series colours (one per measure / group): Okabe-Ito, ordered for contrast on
# white (yellow, the weakest on white, last but one); recycled past nine, so a
# caller always pairs it with .psn_group_key()'s shapes and line types.
.psn_series <- function(k) {
  base <- c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9",
            "#000000", "#F0E442", "#999999")
  rep_len(base, k)
}

# Colour + point shape + line type per group, so a group is identifiable
# without colour (shapes and line types cycle on different periods).
.psn_group_key <- function(k) {
  list(col = .psn_series(k),
       pch = rep_len(c(19, 17, 15, 18, 8, 1, 2, 0, 5, 6), k),
       lty = rep_len(c(1, 2, 3, 4, 5, 6), k))
}

# Validate the measures a plot method is asked to draw: each must be one of
# `available` (numeric columns of `x`). Columns with no finite value are
# skipped with a classed warning; if none is left the call errors.
#' @noRd
.psn_plot_measures <- function(x, measures, available, fn) {
  if (!is.character(measures) || !length(measures) || anyNA(measures))
    .psn_stop("`measures` must be a character vector of column names to plot.")
  measures <- unique(measures)
  bad <- setdiff(measures, available)
  if (length(bad))
    .psn_stop(sprintf("%s(): unknown `measures` %s; available: %s.", fn,
                      paste(sprintf("'%s'", bad), collapse = ", "),
                      paste(available, collapse = ", ")))
  has <- vapply(measures, function(m) any(is.finite(x[[m]])), logical(1))
  if (!any(has))
    .psn_stop(sprintf("%s(): no requested measure has a finite value to plot (%s).",
                      fn, paste(measures, collapse = ", ")),
              class = c("psychnet_nothing_to_plot", "psychnet_bad_input"))
  if (!all(has))
    .psn_warn(sprintf("%s(): skipping measure(s) with no finite value: %s.", fn,
                      paste(measures[!has], collapse = ", ")),
              class = "psychnet_plot_skipped")
  measures[has]
}

# Per-element transparency (grDevices::adjustcolor() takes one alpha only).
.psn_alpha <- function(col, alpha) {
  alpha <- ifelse(is.finite(alpha), pmin(pmax(alpha, 0), 1), 0.1)
  vapply(seq_along(col), function(i)
    grDevices::adjustcolor(col[i], alpha.f = alpha[i]), character(1))
}

# Stars for a p-value, bootnet-style.
.psn_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", ""))))
}

# Point style per item for a significance flag: significant = filled circle on
# a solid segment in the emphasis colour; not significant (FALSE or NA) = open
# circle on a dashed segment in grey. NULL flag = every item emphasised.
.psn_sig_style <- function(significant, n) {
  sig <- if (is.null(significant)) rep(TRUE, n) else significant %in% TRUE
  list(col = ifelse(sig, .psn_pal$sig, .psn_pal$nonsig),
       pch = ifelse(sig, 19, 1),
       lty = ifelse(sig, 1, 2))
}

# Horizontal point-range panel: one row per item, segment = interval, dot =
# point estimate, optional reference line. With `significant`, a significant
# item is a filled point on a solid segment in the emphasis colour and a
# non-significant one an open point on a dashed segment in grey, with a key
# under the title. `point_col` / `point_pch` (one per item, input order) give
# per-item colours and shapes instead (e.g. community membership).
.psn_pointrange <- function(value, lower, upper, labels, significant = NULL,
                            main = "", xlab = "", ref = 0, sort = TRUE,
                            point_col = NULL, point_pch = NULL) {
  n <- length(value)
  ord <- if (isTRUE(sort)) order(value) else seq_len(n)
  value <- value[ord]; lower <- lower[ord]; upper <- upper[ord]
  labels <- labels[ord]
  sty <- .psn_sig_style(if (is.null(significant)) NULL else significant[ord], n)
  cols <- sty$col; pchs <- sty$pch; ltys <- sty$lty
  if (!is.null(point_col)) cols <- point_col[ord]
  if (!is.null(point_pch)) pchs <- point_pch[ord]
  xr <- range(c(lower, upper, value, ref), finite = TRUE)
  y <- seq_len(n)
  graphics::plot(NA, xlim = xr, ylim = c(0.5, n + 0.5), yaxt = "n",
                 xlab = xlab, ylab = "", main = main, bty = "n")
  if (length(ref)) graphics::abline(v = ref, col = .psn_pal$ref, lty = 2)
  graphics::segments(lower, y, upper, y, col = cols, lwd = 2, lty = ltys)
  graphics::points(value, y, pch = pchs, col = cols, lwd = 1.5,
                   cex = min(1, 0.9 * 30 / max(n, 1)) + 0.4)
  if (!is.null(significant))
    graphics::mtext("filled = significant, open = not significant",
                    side = 3, line = 0.2, cex = 0.7)
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
# pair differs significantly (diverging blue if the row item is larger, red if
# smaller) and left faint otherwise. Every off-diagonal cell also prints the
# sign of the difference as "+" / "-" (followed by significance stars), so the
# direction is readable without colour.
.psn_box_matrix <- function(items, observed, df, main = "", show_p = TRUE) {
  # Order axes by observed value (high to low), as documented; df is keyed by
  # item name so the index map below follows the reordering automatically.
  o <- order(observed, decreasing = TRUE)
  items <- items[o]; observed <- observed[o]
  k <- length(items)
  tol <- sqrt(.Machine$double.eps)
  idx <- stats::setNames(seq_len(k), items)
  D <- matrix(NA_real_, k, k); P <- matrix(NA_real_, k, k)
  S <- matrix(FALSE, k, k)
  ra <- idx[as.character(df$item1)]; rb <- idx[as.character(df$item2)]
  ok <- !is.na(ra) & !is.na(rb)
  ab <- cbind(ra[ok], rb[ok]); ba <- ab[, c(2L, 1L), drop = FALSE]
  D[ab] <- df$obs_diff[ok]; D[ba] <- -df$obs_diff[ok]
  P[ab] <- df$p_value[ok];  P[ba] <- df$p_value[ok]
  sg <- df$significant[ok] %in% TRUE
  S[ab] <- sg; S[ba] <- sg
  fin_d <- abs(D)[is.finite(D)]
  dmax <- if (length(fin_d)) max(fin_d) else 1
  if (dmax < tol) dmax <- 1
  fin_o <- abs(observed)[is.finite(observed)]
  omax <- if (length(fin_o)) max(fin_o) else 1
  if (omax < tol) omax <- 1

  graphics::plot(NA, xlim = c(0.5, k + 0.5), ylim = c(0.5, k + 0.5),
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = main,
                 bty = "n", asp = 1)
  cx <- if (k > 25) 0.5 else if (k > 12) 0.65 else 0.8
  graphics::axis(1, at = seq_len(k), labels = items, las = 2, tick = FALSE,
                 cex.axis = cx)
  graphics::axis(2, at = rev(seq_len(k)), labels = items, las = 1,
                 tick = FALSE, cex.axis = cx)

  # Diagonal: the observed value, shaded by magnitude.
  dg <- seq_len(k); ydg <- k - dg + 1
  shade <- .psn_alpha(rep(.psn_pal$fill, k), 0.15 + 0.7 * abs(observed) / omax)
  graphics::rect(dg - 0.5, ydg - 0.5, dg + 0.5, ydg + 0.5,
                 col = shade, border = "white")
  graphics::text(dg, ydg, formatC(observed, format = "f", digits = 2),
                 cex = cx * 0.9, col = "#1A1A1A")

  # Off-diagonal cells (row ri, column ci) with an observed difference.
  cell <- expand.grid(ri = seq_len(k), ci = seq_len(k))
  cell <- cell[cell$ri != cell$ci, , drop = FALSE]
  rc <- cbind(cell$ri, cell$ci)
  d <- D[rc]
  keep <- !is.na(d)
  if (!any(keep)) return(invisible(NULL))
  rc <- rc[keep, , drop = FALSE]; d <- d[keep]
  sig <- S[rc]
  yy <- k - rc[, 1L] + 1; xx <- rc[, 2L]
  positive <- d > tol; negative <- d < -tol
  base_col <- ifelse(positive, .psn_pal$pos,
              ifelse(negative, .psn_pal$neg, .psn_pal$nonsig))
  alpha <- ifelse(sig, 0.25 + 0.65 * abs(d) / dmax, 0.08)
  graphics::rect(xx - 0.5, yy - 0.5, xx + 0.5, yy + 0.5,
                 col = .psn_alpha(base_col, alpha),
                 border = "white")
  glyph <- ifelse(positive, "+", ifelse(negative, "-", "0"))
  stars <- if (isTRUE(show_p)) ifelse(sig, .psn_stars(P[rc]), "") else ""
  txt_col <- ifelse(sig & alpha >= 0.5, "white",
             ifelse(sig, "#1A1A1A", "#555555"))
  graphics::text(xx, yy, paste0(glyph, stars), cex = cx,
                 col = txt_col, font = ifelse(sig, 2, 1))
  invisible(NULL)
}

# ---- centrality ------------------------------------------------------------

# Per-measure transform for the line plot, computed on the finite values (a
# custom measure may be undefined for some nodes; those stay NA = a gap).
.psn_rescale <- function(v, scale = c("raw", "z", "relative")) {
  scale <- .psn_match_arg(scale)
  fin <- v[is.finite(v)]
  tol <- sqrt(.Machine$double.eps)
  switch(scale,
    raw = v,
    z = {
      s <- if (length(fin) > 1L) stats::sd(fin) else NA_real_
      if (is.na(s) || s < tol) v - mean(fin) else (v - mean(fin)) / s
    },
    relative = {
      rng <- range(fin); d <- diff(rng)
      if (d < tol) v * 0 else (v - rng[1]) / d
    })
}

#' @describeIn net_centralities One panel per measure (see `type`, `scale`).
#'   Raises `psychnet_bad_input` for a `measures` entry that is not a measure
#'   column, skips a column with no finite value with a `psychnet_plot_skipped`
#'   warning, and errors (`psychnet_nothing_to_plot`) when none has one.
#' @export
plot.psychnet_centrality <- function(x, type = c("bar", "line"),
                                     scale = c("raw", "z", "relative"),
                                     measures = NULL, ...) {
  type <- .psn_match_arg(type)
  scale <- .psn_match_arg(scale)
  avail <- setdiff(names(x), "node")
  avail <- avail[vapply(avail, function(m) is.numeric(x[[m]]), logical(1))]
  ms <- .psn_plot_measures(x, measures %||% avail, avail,
                           "plot.psychnet_centrality")
  nodes <- x$node
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  if (type == "line") {
    # One faceted panel per measure (qgraph centralityPlot layout): nodes share
    # a single vertical order across panels, each panel keeps its OWN x-axis, so
    # the measures never overlap and no cross-measure rescaling is forced.
    nn <- length(nodes)
    val <- matrix(vapply(ms, function(m) .psn_rescale(x[[m]], scale),
                         numeric(nn)),
                  nrow = nn, dimnames = list(nodes, ms))
    # Stable node order shared by every panel: overall importance (mean
    # z-score over the measures a node has; a missing value is a gap, not 0).
    zs <- matrix(vapply(ms, function(m) .psn_rescale(x[[m]], "z"), numeric(nn)),
                 nrow = nn)
    ord <- order(rowMeans(zs, na.rm = TRUE), nodes)
    val <- val[ord, , drop = FALSE]; lab <- nodes[ord]
    np <- length(lab); y <- seq_len(np)
    key <- .psn_group_key(length(ms))
    xlab <- switch(scale, raw = "centrality",
                   z = "standardised centrality (z)",
                   relative = "relative centrality (0-1)")
    lm <- .psn_left_margin(lab)
    graphics::par(mfrow = c(1, length(ms)), oma = c(3.5, lm, 1.5, 1))
    invisible(lapply(seq_along(ms), function(j) {
      graphics::par(mar = c(0.6, 0.4, 2.2, 1.2))
      vv <- val[, j]
      graphics::plot(NA, xlim = grDevices::extendrange(vv[is.finite(vv)], f = 0.08),
                     ylim = c(0.5, np + 0.5), yaxt = "n", bty = "n",
                     xlab = "", ylab = "", main = ms[j],
                     col.main = key$col[j], font.main = 2, cex.main = 0.95)
      graphics::segments(graphics::par("usr")[1], y, graphics::par("usr")[2], y,
                         col = "#EDEDED", lwd = 1)
      if (scale == "z") graphics::abline(v = 0, col = .psn_pal$ref, lty = 2)
      graphics::lines(vv, y, col = key$col[j], lwd = 2, lty = key$lty[j])
      graphics::points(vv, y, col = key$col[j], pch = key$pch[j],
                       cex = if (np > 20) 0.8 else 1.2)
      if (j == 1L)
        graphics::axis(2, at = y, labels = lab, las = 1, tick = FALSE,
                       cex.axis = if (np > 20) 0.65 else 0.9, xpd = NA)
    }))
    graphics::mtext(xlab, side = 1, outer = TRUE, line = 1.6, cex = 0.9)
    return(invisible(x))
  }

  graphics::par(mfrow = c(1, length(ms)),
                mar = c(4, .psn_left_margin(nodes), 3, 1), oma = c(0, 0, 0, 0))
  invisible(lapply(ms, function(m)
    .psn_pointrange(x[[m]], x[[m]], x[[m]], nodes, main = m, xlab = "",
                    ref = numeric(0), sort = TRUE)))
  invisible(x)
}

# ---- bootstrap -------------------------------------------------------------

# A bootstrap centrality `measure` must be one the object was run with.
.psn_boot_measures <- function(x, measure) {
  ms <- measure %||% x$measures
  if (!is.character(ms) || !length(ms) || anyNA(ms))
    .psn_stop("`measure` must be a character vector of bootstrapped measures.")
  bad <- setdiff(ms, x$measures)
  if (length(bad))
    .psn_stop(sprintf("`measure` %s was not bootstrapped; available: %s.",
                      paste(sprintf("'%s'", bad), collapse = ", "),
                      paste(x$measures, collapse = ", ")))
  unique(ms)
}

#' @describeIn net_boot Edge, centrality, difference or predictability
#'   intervals (see `type`). Significant edges are filled points on solid
#'   segments, non-significant ones open points on dashed segments.
#' @export
plot.psychnet_bootstrap <- function(x, type = c("edges", "centrality",
                                                 "edge_diff", "centrality_diff",
                                                 "predictability"),
                                    measure = NULL, ...) {
  type <- .psn_match_arg(type)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

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
    ms <- .psn_boot_measures(x, measure)
    graphics::par(mfrow = c(1, length(ms)),
                  mar = c(4, .psn_left_margin(x$node_labels), 3, 1))
    invisible(lapply(ms, function(m)
      .psn_pointrange(x$centrality[[m]],
                      x$centrality[[paste0(m, "_lower")]],
                      x$centrality[[paste0(m, "_upper")]],
                      x$node_labels, main = m, xlab = "", ref = numeric(0))))
    return(invisible(x))
  }

  if (type == "predictability") {
    if (is.null(x$predictability))
      .psn_stop("No predictability stored. Re-run net_boot(predictability = TRUE).",
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
  m <- .psn_boot_measures(x, measure)[1]
  df <- difference_test(x, type = m)
  graphics::par(mar = c(7, 7, 3, 1))
  .psn_box_matrix(x$node_labels, x$centrality[[m]], df,
                  main = sprintf("%s differences", m))
  invisible(x)
}

#' @describeIn difference_test The bootnet "significance box" matrix: items on
#'   both axes ordered by observed value, the diagonal showing that value, and
#'   an off-diagonal cell filled when the pair differs (diverging blue when the
#'   row item is larger, red when smaller; every cell also prints the sign as
#'   `+` / `-`, followed by significance stars). `style = "forest"` draws one
#'   row per pair instead: the difference with its interval, a reference line
#'   at zero, significant pairs as filled points on solid segments and the
#'   rest as open points on dashed segments.
#' @param x A `psychnet_difference` data frame from `difference_test()`.
#' @param style For `plot()`: `"box"` (default) or `"forest"`.
#' @param ... Unused.
#' @export
plot.psychnet_difference <- function(x, style = c("box", "forest"), ...) {
  style <- .psn_match_arg(style)
  obs <- attr(x, "observed")
  ty  <- attr(x, "diff_type")
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
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

#' @describeIn net_stability Case-dropping curves with the CS-coefficients;
#'   each measure has its own Okabe-Ito colour, line type and point shape.
#' @export
plot.psychnet_stability <- function(x, ...) {
  tab <- x$table
  ms <- unique(tab$measure)
  key <- .psn_group_key(length(ms))
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(4.5, 4.5, 3, 1))
  graphics::plot(NA, xlim = range(tab$drop_prop), ylim = c(0, 1),
                 xlab = "proportion of cases dropped",
                 ylab = "mean correlation with full sample",
                 main = "Centrality stability", bty = "n")
  graphics::abline(h = x$threshold, col = .psn_pal$ref, lty = 2)
  invisible(lapply(seq_along(ms), function(i) {
    d <- tab[tab$measure == ms[i], , drop = FALSE]
    d <- d[order(d$drop_prop), , drop = FALSE]
    if (!is.null(d$sd_cor) && any(is.finite(d$sd_cor))) {
      lo <- pmax(0, d$mean_cor - d$sd_cor); hi <- pmin(1, d$mean_cor + d$sd_cor)
      graphics::polygon(c(d$drop_prop, rev(d$drop_prop)), c(lo, rev(hi)),
                        col = grDevices::adjustcolor(key$col[i], alpha.f = 0.18),
                        border = NA)
    }
    graphics::lines(d$drop_prop, d$mean_cor, col = key$col[i], lwd = 2,
                    lty = key$lty[i])
    graphics::points(d$drop_prop, d$mean_cor, col = key$col[i],
                     pch = key$pch[i], cex = 0.9)
  }))
  leg <- sprintf("%s (CS = %.2f)", ms, x$cs[ms])
  graphics::legend("bottomleft", legend = leg, col = key$col, lwd = 2,
                   lty = key$lty, pch = key$pch, bty = "n", cex = 0.8)
  invisible(x)
}

# ---- network comparison test ----------------------------------------------

#' @describeIn net_compare Permutation null of S or M, or the per-edge
#'   differences (see `type`); significant edges are labelled `*` as well as
#'   coloured.
#' @export
plot.psychnet_nct <- function(x, type = c("strength", "structure", "edges"),
                              alpha = 0.05, ...) {
  type <- .psn_match_arg(type)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  if (type %in% c("strength", "structure")) {
    stat <- if (type == "strength") x$strength else x$structure
    lab  <- if (type == "strength") {
              if (isFALSE(x$params$abs)) "global expected influence difference (S)"
              else "global strength difference (S)"
            }
            else "maximum edge difference (M)"
    graphics::par(mar = c(4.5, 4.5, 3, 1))
    graphics::hist(stat$perm, breaks = 30, col = .psn_pal$nonsig,
                   border = "white", xlab = lab,
                   xlim = range(c(stat$perm, stat$observed)),
                   main = sprintf("NCT permutation null (p = %.3f)", stat$p_value))
    graphics::abline(v = stat$observed, col = .psn_pal$mark, lwd = 3)
    graphics::legend("topright", legend = "observed", col = .psn_pal$mark,
                     lwd = 3, bty = "n", cex = 0.85)
    return(invisible(x))
  }

  ed <- as.data.frame(x, alpha = alpha)
  # significance is marked by a direct "*" label as well as by colour
  labels <- paste0(ed$from, "--", ed$to, ifelse(ed$significant, " *", ""))
  cols <- ifelse(ed$significant, .psn_pal$sig, .psn_pal$nonsig)
  graphics::par(mar = c(4, .psn_left_margin(labels), 3, 1))
  ord <- order(ed$difference, ed$from, ed$to)
  graphics::barplot(ed$difference[ord], names.arg = labels[ord], horiz = TRUE,
                    las = 1, col = cols[ord], border = NA,
                    xlab = "|edge difference|",
                    main = sprintf("Per-edge differences (%d with p < %g, marked *)",
                                   sum(ed$significant), alpha),
                    cex.names = if (length(labels) > 20) 0.6 else 0.8)
  invisible(x)
}

# ---- network ---------------------------------------------------------------

#' @describeIn psychnet-class Draws the network with `cograph::splot()` (a
#'   Suggested package; an error says so when it is not installed).
#' @export
plot.psychnet <- function(x, ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    .psn_stop("Network plotting needs the 'cograph' package. ",
         "Install it, or use net_centralities()/net_boot() with their plot ",
         "methods for diagnostics.", call. = FALSE)
  }
  invisible(cograph::splot(x, ...))
}
