#' Predict from non-parametric marginal structural model
#' @param object object of class "targeted_msm"
#' @param newdata data frame of new data
#' @param estimator which estimator to use for predictions: "onestep", "bayes", or "tmle" (the default).
#' @param type type of prediction; "point" for point estimates, "joint" for draws from joint distribution.
#' @param ... additional arguments (not currently used)
#' @return For predictions of type "point", a vector of predictions is returned.
#' For predictions of type "joint", a list is returned, where each element is a vector of
#' predictions for the corresponding row of newdata.
#' @exportS3Method
predict.targeted_msm <- function(object, newdata, estimator = "tmle", type = "point", ...) {
  if(object$estimand == "treatment_effect_modification") {
    predict_treatment_effect_modification(object, ...)
  }
  else if(object$estimand == "categorical_dose_response") {
    predict_categorical_dose_response(object, ...)
  }
}

#' @importFrom stats model.matrix var
#' @importFrom mvtnorm rmvnorm
#' @noRd
predict_treatment_effect_modification <- function(fit, newdata, estimator = "tmle", type = "point") {
  design_matrix <- stats::model.matrix(fit$formula, newdata)

  n <- fit$n

  draws <- 1e3
  if(estimator %in% c("tmle", "onestep")) {
    if(estimator == "tmle") {
      est <- fit$tmle$est
      joint_draws <- mvtnorm::rmvnorm(draws, mean = fit$tmle$est, sigma = stats::var(fit$tmle$eif) / n)
    }
    else if(estimator == "onestep") {
      est <- fit$onestep$est
      joint_draws <- mvtnorm::rmvnorm(draws, mean = fit$onestep$est, sigma = stats::var(fit$onestep$eif) / n)
    }

    if(type == "point") {
      (design_matrix %*% t(t(est)))[,1]
    }
    else {
      apply(design_matrix, 1, function(x) (joint_draws %*% t(t(x)))[, 1], simplify = FALSE)
    }
  }
  else if(estimator == "bayes") {
    post <- matrix(fit$tmle$samples, nrow = 4 * dim(fit$tmle$samples)[2], ncol = dim(fit$tmle$samples)[3])
    apply(design_matrix, 1, function(x) (post %*% t(t(x)))[, 1], simplify = FALSE)
  }
}

#' @importFrom stats model.matrix var
#' @importFrom mvtnorm rmvnorm
#' @noRd
predict_categorical_dose_response <- function(fit, newdata, estimator = "tmle", type = "point") {
  design_matrix <- stats::model.matrix(fit$formula, newdata)

  n <- fit$n

  draws <- 1e3
  if(estimator %in% c("tmle", "onestep")) {
    if(estimator == "tmle") {
      est <- fit$tmle$est
      joint_draws <- mvtnorm::rmvnorm(draws, mean = fit$tmle$est, sigma = stats::var(fit$tmle$eif) / n)
    }
    else if(estimator == "onestep") {
      est <- fit$onestep$est
      joint_draws <- mvtnorm::rmvnorm(draws, mean = fit$onestep$est, sigma = stats::var(fit$onestep$eif) / n)
    }

    if(type == "point") {
      (design_matrix %*% t(t(est)))[,1]
    }
    else {
      apply(design_matrix, 1, function(x) (joint_draws %*% t(t(x)))[, 1], simplify = FALSE)
    }
  }
  else if(estimator == "bayes") {
    post <- matrix(fit$tmle$samples, nrow = 4 * dim(fit$tmle$samples)[2], ncol = dim(fit$tmle$samples)[3])
    apply(design_matrix, 1, function(x) (post %*% t(t(x)))[, 1], simplify = FALSE)
  }
}
