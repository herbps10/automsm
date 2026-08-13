#' Assemble final result object for single-time-point NP-MSM estimators
#' @param estimand Character string identifying the estimand (e.g., "cate" or "dose_response")
#' @param p Number of working model coefficients
#' @param n Sample size
#' @param formula Working model formula
#' @param working_model Working model
#' @param loss loss function
#' @param learners_trt Propensity score SuperLearner libraries
#' @param learners_outcome Outcome regression SuperLearner libraries
#' @param nuisance Estimated nuisance parameters
#' @param plugin List containing plugin estimate
#' @param onestep List containing one-step estimate
#' @param tmle List containing TMLE estimate (or \code{NULL}, if TMLE was not run)
#' @param bayes_tmle List containing TMLE estimate (or \code{NULL}, if Bayesian TMLE was not run)
#' @param ... Additional arguments
#'
#' @return list of class \code{"automsm"}
#' @noRd
assemble_result <- function(estimand, p, n, formula, working_model, loss, terms, learners_trt, learners_outcome, nuisance, plugin, onestep, tmle, bayes_tmle, ...) {
  res <- list(
    estimand = estimand,
    p = p,
    n = n,
    formula = formula,
    working_model = working_model,
    loss = loss,
    terms = terms,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    nuisance = nuisance,
    plugin = plugin,
    onestep = onestep,
    tmle = tmle,
    bayes_tmle = bayes_tmle,
    ...
  )

  class(res) <- "automsm"

  res
}
