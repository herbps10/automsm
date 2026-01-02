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
