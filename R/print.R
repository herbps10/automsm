#' Print method for objects of class "automsm"
#'
#' @param x An object of class \code{automsm}
#' @param ... Additional arguments (not currently used)
#' @return \code{x}, invisibly
#' @importFrom cli cli_text
#' @export
print.automsm <- function(x, ...) {
  label <- automsm_label(x$estimand)

  cli::cli_text("{.strong automsm}: {label} ({.field n} = {x$n})")

  invisible(x)
}

#' Human-readable label for an estimand
#' @noRd

automsm_label <- function(estimand) {
  switch(
    estimand,
    cate = "conditional average treatment effect",
    dose_response = "dose-response",
    longitudinal_dose_response = "longitudinal dose-response",
    estimand
  )
}
