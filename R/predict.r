#' @export
predict.targeted_msm <- function(object, ...) {
  if(object$estimand == "treatment_effect_modification") {
    predict_treatment_effect_modification(object, ...)
  }
}

predict_treatment_effect_modification <- function(fit, newdata, estimator = "tmle", type = "point") {
  design_matrix <- model.matrix(fit$f, newdata)

  n <- fit$n

  draws <- 1e3
  if(estimator %in% c("tmle", "onestep")) {
    if(estimator == "tmle") {
      est <- fit$tmle$est
      joint_draws <- mvtnorm::rmvnorm(draws, mean = fit$tmle$est, sigma = var(fit$tmle$eif) / n)
    }
    else if(estimator == "onestep") {
      est <- fit$onestep$est
      joint_draws <- mvtnorm::rmvnorm(draws, mean = fit$onestep$est, sigma = var(fit$onestep$eif) / n)
    }

    if(type == "point") {
      (design_matrix %*% t(t(est)))[,1]
    }
    else {
      apply(design_matrix, 1, \(x) (joint_draws %*% t(t(x)))[, 1], simplify = FALSE)
    }
  }
  else if(estimator == "bayes") {
    post <- matrix(fit$tmle$samples, nrow = 4 * dim(fit$tmle$samples)[2], ncol = dim(fit$tmle$samples)[3])
    apply(design_matrix, 1, \(x) (post %*% t(t(x)))[, 1], simplify = FALSE)
  }
}
