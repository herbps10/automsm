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
    treatment_effect_modification = "treatment effect modification",
    categorical_dose_response = "categorical dose-response function",
    longitudinal_treatment = "longitudinal treatment",
    estimand
  )
}
