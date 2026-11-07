# Shared conditions and argument checks.
#
# Every input-contract error in psychnets is a classed condition, so callers
# and tests can catch it precisely: class c("psychnet_bad_input",
# "psychnet_error") unless a more specific class is given. The message is the
# contract a user reads; the class is the contract code relies on.

# Raise a classed error. Arguments are pasted together like stop(); `call.` is
# accepted (and ignored -- the call is never shown) so stop(..., call. = FALSE)
# call sites convert one-for-one.
#' @noRd
.psn_stop <- function(..., call. = FALSE, class = "psychnet_bad_input") {
  msg <- paste0(...)
  stop(errorCondition(msg, class = unique(c(class, "psychnet_error")),
                      call = NULL))
}

# Raise a classed warning (class plus "psychnet_warning").
#' @noRd
.psn_warn <- function(..., class) {
  msg <- paste0(...)
  warning(warningCondition(msg, class = unique(c(class, "psychnet_warning")),
                           call = NULL))
}

# A single finite number, optionally bounded. `lower`/`upper` are inclusive
# unless `open = TRUE`. `null_ok` allows NULL (returned unchanged).
#' @noRd
.psn_check_number <- function(x, name, lower = -Inf, upper = Inf,
                              open = FALSE, whole = FALSE, null_ok = FALSE) {
  if (is.null(x) && null_ok) return(invisible(x))
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (ok) {
    ok <- if (open) x > lower && x < upper else x >= lower && x <= upper
  }
  if (ok && whole) ok <- abs(x - round(x)) < sqrt(.Machine$double.eps)
  if (!ok) {
    rng <- if (open) sprintf("in (%s, %s)", format(lower), format(upper))
           else sprintf("in [%s, %s]", format(lower), format(upper))
    .psn_stop(sprintf("`%s` must be a single finite %s %s.", name,
                      if (whole) "whole number" else "number", rng))
  }
  invisible(x)
}

# A single TRUE/FALSE.
#' @noRd
.psn_check_flag <- function(x, name) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x)))
    .psn_stop(sprintf("`%s` must be TRUE or FALSE.", name))
  invisible(x)
}

# Node labels: NULL, or one unique, non-missing label per variable. `p` is the
# number of variables the labels must name.
#' @noRd
.psn_check_labels <- function(labels, p) {
  if (is.null(labels)) return(invisible(NULL))
  if (length(labels) != p) {
    .psn_stop(sprintf(paste0(
      "`labels` has %d entries but the network has %d nodes (a non-numeric or ",
      "constant column may have been dropped from the data)."),
      length(labels), p))
  }
  if (anyNA(labels) || anyDuplicated(labels))
    .psn_stop("`labels` must be unique and non-missing.")
  invisible(as.character(labels))
}

# The sample size that goes with a correlation matrix: a single whole number
# of at least `min`.
#' @noRd
.psn_check_n <- function(n, min = 1, name = "n") {
  .psn_check_number(n, name, lower = min, whole = TRUE)
}

# match.arg() with a classed error. Called like match.arg(): with `choices`
# omitted, the choices are read from the calling function's formal default, as
# base match.arg() does. Partial matching and `several.ok` behave as in base.
#' @noRd
.psn_match_arg <- function(arg, choices, several.ok = FALSE) {
  arg_name <- deparse(substitute(arg))
  if (missing(choices)) {
    caller <- sys.function(sys.parent())
    choices <- eval(formals(caller)[[as.character(substitute(arg))]],
                    envir = parent.frame())
  }
  tryCatch(match.arg(arg, choices, several.ok = several.ok),
           error = function(e) {
             .psn_stop(sprintf("`%s` %s", arg_name,
                               sub("^'arg' ", "", conditionMessage(e))))
           })
}
