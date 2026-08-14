#' Construct an estimand specification
#'
#' @noRd
new_msm_spec <- function(estimand,
                         init_state, steps, psi_from_state,
                         make_clever, fluctuation_objective, apply_update, delta,
                         tol, sweep_criterion,
                         tmle_loss,
                         nan_guard = TRUE,
                         supports_bayes = FALSE,
                         bayes_loglik = NULL, dpsi_depsilon = NULL,
                         extra_result = function(problem, state) list()) {
  spec <- list(
    estimand = estimand,
    init_state = init_state, steps = steps, psi_from_state = psi_from_state,
    make_clever = make_clever, fluctuation_objective = fluctuation_objective, apply_update = apply_update, delta = delta,
    tol = tol, sweep_criterion = sweep_criterion,
    tmle_loss = tmle_loss,
    nan_guard = nan_guard,
    supports_bayes = supports_bayes,
    bayes_loglik = bayes_loglik, dpsi_depsilon = dpsi_depsilon,
    extra_result = extra_result
  )

  validate_spec(spec)
  structure(spec, class = "msm_spec")
}

validate_spec <- function(spec) {
  fns <- c("init_state", "steps", "psi_from_state", "make_clever", "fluctuation_objective", "apply_update", "delta", "sweep_criterion", "extra_result")
  for(f in fns) checkmate::assert_function(spec[[f]], .var.name = f)
  checkmate::assert_number(spec$tol, lower = 0, finite = TRUE)
  checkmate::assert_flag(spec$supports_bayes)
  if(spec$supports_bayes) {
    checkmate::assert_function(spec$bayes_loglik)
    checkmate::assert_function(spec$dpsi_depsilon)
  }
  TRUE
}
