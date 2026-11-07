# Within-person centering and scaling: a data-preparation verb whose output
# feeds any estimator or framework verb (psychnet(), net_compare(), net_boot()).

#' Center and scale variables within persons
#'
#' Removes each person's own mean from every variable and, by default, divides
#' by that person's own standard deviation, so that every network estimated
#' from the result describes *within-person* associations: how a person's
#' deviations from their own typical level go together. This is the
#' within-person standardization of `mlVAR::mlVAR(scaleWithin = TRUE)`
#' (centering per person, then division by the person's SD), offered as a
#' standalone step so it can precede any analysis. The same scaling is built
#' into [psychnet()] as `within =`, which also lets [net_compare()] pair two
#' occasions by person.
#'
#' Means and standard deviations use each person's available values (the SD
#' has denominator n - 1, as [stats::sd()]). A value is set to `NA` when its
#' person has fewer than two observed values of that variable, or (with
#' `scale = TRUE`) when the person's values of it are constant: there is no
#' within-person variation to describe. Both cases raise a
#' `psychnet_constant_within` warning with the count.
#'
#' **Two occasions per person.** With exactly two observations, each person's
#' centered deviations at the two occasions are mirror images (`x1 - m = -(x2 -
#' m)`), so the networks of the two occasions are identical by construction
#' and a comparison between them cannot detect anything; after scaling, every
#' value is \eqn{\pm 1/\sqrt{2}}. When no person has more than two
#' observations a `psychnet_mirror_deviations` warning says so. Within-person
#' comparisons of occasions need three or more occasions per person.
#'
#' @param data A data frame with one row per observation (person x occasion).
#' @param id Name of the column identifying persons. Missing identifiers are an
#'   error.
#' @param vars Names of the columns to center/scale. Default `NULL`: every
#'   numeric column except `id`. Name them explicitly when the data carry other
#'   numeric columns (an occasion index, a covariate) that must stay as they
#'   are.
#' @param center Subtract each person's mean. Default `TRUE`.
#' @param scale Divide by each person's standard deviation. Default `TRUE`.
#' @return A data frame with the same rows and columns as `data`, the `vars`
#'   columns replaced by their within-person centered (and scaled) values and
#'   every other column unchanged.
#'
#'   Errors of class `psychnet_bad_input` for a missing `id` column, missing
#'   identifiers, or `vars` that are not numeric columns of `data`.
#' @references Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, D. (2018).
#'   The Gaussian graphical model in cross-sectional and time-series data.
#'   *Multivariate Behavioral Research*, 53(4), 453-480.
#'   \doi{10.1080/00273171.2018.1454823}
#' @examples
#' v <- c("efficacy", "planning", "effort", "enjoyment", "anxiety")
#' sw <- scale_within(srl_panel, id = "student", vars = v)
#' head(sw)
#'
#' # The same scaling inside the front door: one within-person network per
#' # week; compare two weeks with net_compare(nets, c(1, 4)), which pairs the
#' # students automatically (see ?psychnet).
#' psychnet(srl_panel, vars = efficacy:anxiety, group = "week",
#'          within = "student")
#' @export
scale_within <- function(data, id, vars = NULL, center = TRUE, scale = TRUE) {
  if (!is.data.frame(data))
    .psn_stop("`data` must be a data frame with one row per observation.")
  if (!(is.character(id) && length(id) == 1L && id %in% names(data)))
    .psn_stop("`id` must name one column of `data`.")
  .psn_check_flag(center, "center")
  .psn_check_flag(scale, "scale")
  ids <- data[[id]]
  if (anyNA(ids))
    .psn_stop(sprintf("`%s` has %d missing identifier(s); every row needs a person.",
                      id, sum(is.na(ids))))
  if (is.null(vars)) {
    vars <- setdiff(names(data)[vapply(data, is.numeric, logical(1))], id)
  }
  bad <- setdiff(vars, names(data))
  if (length(bad))
    .psn_stop("`vars` not in `data`: ", paste(bad, collapse = ", "), ".")
  non_numeric <- vars[!vapply(data[vars], is.numeric, logical(1))]
  if (length(non_numeric))
    .psn_stop("`vars` must be numeric columns; not numeric: ",
              paste(non_numeric, collapse = ", "), ".")
  if (id %in% vars) .psn_stop("`vars` must not include the `id` column.")
  if (!length(vars)) .psn_stop("No numeric columns to scale.")

  f <- factor(ids)
  per_person <- tabulate(as.integer(f))
  if (max(per_person) <= 2L)
    .psn_warn(paste0(
      "No person has more than two observations. With two occasions each ",
      "person's within-person deviations are mirror images, so the occasions' ",
      "networks are identical by construction (and scaled values are +/-0.707); ",
      "within-person comparisons of occasions need three or more."),
      class = "psychnet_mirror_deviations")

  scaled <- lapply(vars, function(v) {
    x <- as.numeric(data[[v]])
    n_obs <- stats::ave(!is.na(x), f, FUN = sum)
    m <- stats::ave(x, f, FUN = function(z) mean(z, na.rm = TRUE))
    s <- stats::ave(x, f, FUN = function(z) {
      zo <- z[!is.na(z)]
      if (length(zo) < 2L) NA_real_ else stats::sd(zo)
    })
    out <- if (center) x - m else x
    # no within-person variation: fewer than two values, or (scaling) a
    # constant person -- a relative tolerance, so units do not matter
    undefined <- n_obs < 2L
    if (scale) {
      undefined <- undefined | (!is.na(s) & s <= 1e-12 * pmax(abs(m), 1e-300))
      out <- out / s
    }
    undefined <- undefined & !is.na(x)
    out[undefined] <- NA_real_
    list(values = out, n_undefined = sum(undefined))
  })
  data[vars] <- lapply(scaled, `[[`, "values")
  n_undefined <- sum(vapply(scaled, `[[`, integer(1), "n_undefined"))
  if (n_undefined > 0L)
    .psn_warn(sprintf(paste0(
      "%d value(s) set to NA: their person has fewer than two observed values ",
      "of that variable%s, so there is no within-person variation."),
      n_undefined, if (scale) " or a constant one" else ""),
      class = "psychnet_constant_within")
  data
}
