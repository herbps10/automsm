#' Print method for objects of class "targeted_msm"
#' @param x object of class "targeted_msm"
#' @param ... additional arguments (not currently used)
#' @export
print.targeted_msm <- function(x, ...) {
  cat("MSM for ")
  if(x$estimand == "treatment_effect_modification") {
    cat("treatment effect modification.")
  }
  else if(x$estimand == "categorical_dose_response") {
    cat("categorical dose-response function.")
  }
}
