#' Predict from a "targeted_msm" object
#'
#' Evaluates the fitted working model at supplied design points.
#'
#' @param object An object of class \code{"targeted_msm"}.
#' @param newdata A data frame of new data.
#' @param estimator Which estimator to use for predictions:
#'   \code{"onestep"}, \code{"bayes"},
#'   or \code{"tmle"} (default). Note that the \code{"bayes"} estimator
#'   is not available for the \code{"longitudinal_treatment"} estimand.
#' @param type Type of prediction; \code{"point"} for point estimates,
#'   \code{"joint"} for draws from joint distribution.
#' @param ... Additional arguments (not currently used).
#' @return For predictions of type \code{"point"}, a vector of predictions is
#'   returned. For predictions of type \code{"joint"}, a list is returned, where
#'   each element is a vector of predictions for the corresponding row of
#'   \code{newdata}
#' @exportS3Method
predict.targeted_msm <- function(
  object,
  newdata,
  estimator = "tmle",
  type = "point",
  ...
) {
  supported <- c(
    "treatment_effect_modification",
    "categorical_dose_response",
    "longitudinal_treatment"
  )

  if(!(object$estimand %in% supported)) {
    stop(
      "predict() is not supported for estimand '", object$estimand, "'.",
      call. = FALSE
    )
  }

  if(object$estimand == "longitudinal_treatment" && estimator == "bayes") {
    stop("The 'bayes' estimator is not available for the ",
         "'longitudinal_treatment' estimand.",
         call. = FALSE)
  }

  predict_working_model(object, newdata, estimator, type)
}

#' Shared prediction helper for targeted_msm objects
#'
#' @importFrom stats model.matrix var
#' @importFrom mvtnorm rmvnorm
#' @importFrom torch torch_tensor
#' @noRd
predict_working_model <- function(
  fit,
  newdata,
  estimator = "onestep",
  type = "point"
) {
  design_matrix <- stats::model.matrix(fit$formula, newdata)
  n <- fit$n
  draws <- 1e3

  if(estimator %in% c("tmle", "onestep")) {
    est <- fit[[estimator]]$est
    eif <- fit[[estimator]]$eif

    if(is.null(est) || is.null(eif)) {
      stop("Estimator '", estimator, "' is not available in this object.", call. = FALSE)
    }

    if(type == "point") {
      as.numeric(fit$working_model(
        torch::torch_tensor(est),
        torch::torch_tensor(design_matrix)
      ))
    }
    else {
      joint_draws <- mvtnorm::rmvnorm(
        draws,
        mean = est,
        sigma = stats::var(eif) / n
      )

      asplit(
        apply(joint_draws, 1, function(x) {
          as.matrix(fit$working_model(torch::torch_tensor(x), torch::torch_tensor(design_matrix)))
        }, simplify = TRUE),
        1
      )
    }
  }
  else if(estimator == "bayes") {
    post <- matrix(
      fit$tmle$samples,
      nrow = 4 * dim(fit$tmle$samples)[2],
      ncol = dim(fit$tmle$samples)[3]
    )
    apply(
      design_matrix,
      1,
      function(x) (post %*% t(t(x)))[, 1],
      simplify = FALSE
    )
  }
  else {
    stop(
      "Unknown estimator '", estimator, "'. ",
      "Choose one of 'tmle', 'onestep', or 'bayes'.",
      call. = FALSE
    )
  }
}
