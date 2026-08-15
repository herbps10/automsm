#' Shared setup for nuisance estimation
#' @importFrom origami make_folds
#' @importFrom SuperLearner SuperLearner.CV.control predict.SuperLearner SuperLearner
#' @importFrom stats gaussian binomial
#' @noRd
nuisance_setup <- function(n, control, outcome_type) {
  cv <- control$folds
  if(is.null(cv)) {
    cv <- origami::make_folds(n, origami::folds_vfold, V = control$outer_folds)
    if (control$outer_folds == 1) {
      cv[[1]]$training_set <- cv[[1]]$validation_set
    }
  }

  list(
    cv = cv,
    cv_control = SuperLearner::SuperLearner.CV.control(V = control$inner_folds),
    outcome_family = if(outcome_type == "binomial") {
      stats::binomial()
    }
    else {
      stats::gaussian()
    }
  )
}

#' @noRd
#' @importFrom withr with_seed
with_nuisance_seed <- function(control, expr) {
  if(is.null(control$seed)) force(expr)
  else withr::with_seed(control$seed, expr)
}

#' Fit a SuperLearner
#' @noRd
sl_fit <- function(Y, X, SL.library, family, cv_control) {
  SuperLearner::SuperLearner(
    Y = Y, X = X,
    SL.library = SL.library,
    family = family,
    cvControl = cv_control,
    env = environment(SuperLearner::SuperLearner)
  )
}

#' Predict from a SuperLearner fit, optionally bounding the predictions
#' @param bounds Optional length-2 numeric passed to \code{bound()}; \code{NULL} for none
#' @importFrom SuperLearner predict.SuperLearner
#' @noRd
sl_predict <- function(fit, newdata, bounds = NULL, epsilon = 1e-5) {
  p <- as.numeric(
    SuperLearner::predict.SuperLearner(fit, newdata = newdata, onlySL = TRUE)$pred
  )
  if(!is.null(bounds)) p <- bound(p, bounds[1], bounds[2], epsilon)
  p
}

#' Fit a SuperLearner and predict on one or more new datasets
#'
#' @param newdata A named list of data frames to predict on.
#' @noRd
sl_fit_predict <- function(Y, X, newdata, SL.library, family, cv_control, bounds = NULL, epsilon = 1e-5) {
  fit <- sl_fit(Y, X, SL.library, family, cv_control)
  list(
    fit = fit,
    pred = lapply(newdata, function(nd) sl_predict(fit, nd, bounds, epsilon))
  )


}
