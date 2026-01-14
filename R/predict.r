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
    predict_treatment_effect_modification(object, newdata, estimator, type, ...)
  }
  else if(object$estimand == "categorical_dose_response") {
    predict_categorical_dose_response(object, newdata, estimator, type, ...)
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
      as.numeric(fit$working_model(torch::torch_tensor(est), torch::torch_tensor(design_matrix)))
    }
    else {
      asplit(apply(joint_draws, 1, function(x) as.matrix(fit$working_model(torch::torch_tensor(x), torch::torch_tensor(design_matrix))), simplify = TRUE), 1)
    }
  }
  else if(estimator == "bayes") {
    post <- matrix(fit$tmle$samples, nrow = 4 * dim(fit$tmle$samples)[2], ncol = dim(fit$tmle$samples)[3])
    apply(design_matrix, 1, function(x) (post %*% t(t(x)))[, 1], simplify = FALSE)
  }
}

#' @importFrom stats model.matrix var
#' @importFrom mvtnorm rmvnorm
#' @importFrom torch torch_tensor
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
      as.numeric(fit$working_model(torch::torch_tensor(est), torch::torch_tensor(design_matrix)))
    }
    else {
      asplit(apply(joint_draws, 1, function(x) as.matrix(fit$working_model(torch::torch_tensor(x), torch::torch_tensor(design_matrix))), simplify = TRUE), 1)
    }
  }
  else if(estimator == "bayes") {
    post <- matrix(fit$tmle$samples, nrow = 4 * dim(fit$tmle$samples)[2], ncol = dim(fit$tmle$samples)[3])
    apply(design_matrix, 1, function(x) (post %*% t(t(x)))[, 1], simplify = FALSE)
  }
}
