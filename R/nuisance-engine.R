#' Seed an expression with a derived offset
#' @noRd
with_seed_offset <- function(seed, offset, expr) {
  if(is.null(seed)) force(expr) else withr::with_seed(seed + offset, expr)
}

#' @param initial function() list(pi = n x tau matrix, mu = k x n x (tau + 1) array, ...)
#' @param refit_node function(t, target_train, fold_index) list(train = |T_v x K, valid = |V_v| x K)
#'   or NULL if the engine cannot refit
#' @noRd
new_nuisance_engine <- function(initial, refit_node = NULL, folds, label = "custom") {
  checkmate::assert_function(initial)
  checkmate::assert_function(refit_node, null.ok = TRUE)

  structure(
    list(
      initial = initial,
      refit_node = refit_node,
      folds = folds,
      supports_refit = !is.null(refit_node),
      label = label
    ),
    class = "automsm_nuisance_engine"
  )
}

#' @noRd
frozen_engine <- function(nuisance_estimates, folds = NULL) {
  new_nuisance_engine(initial = function() nuisance_estimates, refit_node = NULL, folds = folds, label = "frozen")
}

#' @noRd
longitudinal_sl_engine <- function(data, Ls, As, Y, regimes, control, outcome_type) {
  tau <- length(As)
  bounds <- if(identical(outcome_type, "binomial")) c(0, 1) else NULL
  learners_refit <- control$learners_refit %||% control$learners_outcome

  # Construct folds exactly once
  setup <- with_seed_offset(control$seed, 0L, nuisance_setup(nrow(data), control, outcome_type))

  cache <- new.env(parent = emptyenv())

  initial <- function() {
    if(!is.null(cache$initial)) return(cache$initial)
    cache$initial <- with_seed_offset(control$seed, 1L, list(
      pi = estimate_longitudinal_dose_response_propensity_scores(
        data, Ls, As, setup$cv, control$learners_trt, setup$cv_control, control$epsilon
      ),
      mu = estimate_ice_chain(
        data, Ls, As, Y, regimes, setup$cv, control$learners_outcome, setup$outcome_family, setup$cv_control, bounds, control$epsilon
      )
    ))
    cache$initial
  }

  refit_node <- function(t, target_train, fold_index) {
    family <- if(t == tau) setup$outcome_family else stats::gaussian()
    with_seed_offset(control$seed, 1000L + 100L * t + fold_index, fit_ice_node(t, target_train, data, Ls, As, regimes,
                 fold = setup$cv[[fold_index]],
                 learners = learners_refit, family = family, bounds = bounds,
                 epsilon = control$epsilon, cv_control = setup$cv_control))
  }

  new_nuisance_engine(
    folds = setup$cv,
    label = "SuperLearner",
    initial = initial,
    refit_node = refit_node
  )
}
