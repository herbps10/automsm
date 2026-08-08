#' Print method for objects of class "targeted_msm"
#'
#' @param x An object of class \code{targeted_msm}
#' @param ... Additional arguments (not currently used)
#' @return \code{x}, invisibly
#' @importFrom cli cli_text
#' @export
print.targeted_msm <- function(x, ...) {
  label <- targeted_msm_label(x$estimand)

  cli::cli_text("{.strong Targeted MSM}: {label} ({.field n} = {x$n})")

  invisible(x)
}

#' Human-readable label for an estimand
#' @noRd

targeted_msm_label <- function(estimand) {
  switch(
    estimand,
    cate = "conditional average treatment effect",
    dose_response = "dose-response",
    longitudinal_dose_response = "longitudinal dose-response",
    estimand
  )
}
