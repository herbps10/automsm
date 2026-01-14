#' Summary method for class "targeted_msm"
#' @param object object of class "targeted_msm"
#' @param ... additional arguments (not currently used)
#' @export
summary.targeted_msm <- function(object, ...) {
  cat("Marginal Structural Model: ")
  if(object$estimand == "treatment_effect_modification") {
    summary_treatment_effect_modification(object)
  }
  else if(object$estimand == "categorical_dose_response") {
    summary_categorical_dose_response(object)
  }
}

#' @importFrom stringr str_pad
#' @noRd
summary_treatment_effect_modification <- function(x) {
  cat("Treatment Effect Modification\n")

  f <- function(x) stringr::str_pad(round(x, 2), 5, side = "left", pad = " ")

  cat("One-step estimator\n")
  cat("    Est      SE    2.5%   97.5%\n")
  for(index in seq_along(x$onestep$est)) {
    cat(" ", f(x$onestep$est[index]), " ", f(x$onestep$se[index]), " ", f(x$onestep$lower[index]), " ", f(x$onestep$upper[index]), "\n")
  }

  if(!is.null(x$tmle)) {
    cat("TMLE estimator\n")
    cat("    Est      SE    2.5%   97.5%\n")
    for(index in seq_along(x$tmle$est)) {
      cat(" ", f(x$tmle$est[index]), " ", f(x$tmle$se[index]), " ", f(x$tmle$lower[index]), " ", f(x$tmle$upper[index]), "\n")
    }
  }
}

#' @importFrom stringr str_pad
#' @noRd
summary_categorical_dose_response <- function(x) {
  cat("Categorical Dose-Response Function\n")

  f <- function(x) stringr::str_pad(round(x, 2), 5, side = "left", pad = " ")

  cat("One-step estimator\n")
  cat("    Est      SE    2.5%   97.5%\n")
  for(index in seq_along(x$onestep$est)) {
    cat(" ", f(x$onestep$est[index]), " ", f(x$onestep$se[index]), " ", f(x$onestep$lower[index]), " ", f(x$onestep$upper[index]), "\n")
  }

  if(!is.null(x$tmle)) {
    cat("TMLE estimator\n")
    cat("    Est      SE    2.5%   97.5%\n")
    for(index in seq_along(x$tmle$est)) {
      cat(" ", f(x$tmle$est[index]), " ", f(x$tmle$se[index]), " ", f(x$tmle$lower[index]), " ", f(x$tmle$upper[index]), "\n")
    }
  }
}
