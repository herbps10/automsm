#' Construct an estimand specification
#'
#' @noRd
new_msm_spec <- function(
  estimand,
  init_state, steps, psi_from_state,
  make_clever, mu_loss, apply_update, delta,
  tol,
  tmle_loss,
  nan_guard = TRUE,
  supports_bayes = FALSE,
  bayes_loglik = NULL, dpsi_depsilon = NULL,
  nuisance_contract = NULL,
  bayes_clever_scale = function(problem, state) NULL,
  sweep_criterion = default_sweep_criterion
) {
  spec <- list(
    estimand = estimand,
    init_state = init_state, steps = steps, psi_from_state = psi_from_state,
    make_clever = make_clever, mu_loss = mu_loss, apply_update = apply_update, delta = delta,
    tol = tol, sweep_criterion = sweep_criterion,
    tmle_loss = tmle_loss,
    nan_guard = nan_guard,
    supports_bayes = supports_bayes,
    bayes_loglik = bayes_loglik, dpsi_depsilon = dpsi_depsilon,
    bayes_clever_scale = bayes_clever_scale,
    nuisance_contract = nuisance_contract
  )

  validate_spec(spec)
  structure(spec, class = "msm_spec")
}

#' @noRd
default_sweep_criterion <- function(eps_list) {
  max(vapply(eps_list, function(e) max(abs(as.numeric(e))), numeric(1)))
}

validate_spec <- function(spec) {
  fns <- c("init_state", "steps", "psi_from_state", "make_clever",
           "mu_loss", "apply_update", "delta", "sweep_criterion",
           "bayes_clever_scale", "nuisance_contract")

  for(f in fns) checkmate::assert_function(spec[[f]], .var.name = f)

  checkmate::assert_number(spec$tol, lower = 0, finite = TRUE)
  checkmate::assert_flag(spec$supports_bayes)
  if(spec$supports_bayes) {
    checkmate::assert_function(spec$bayes_loglik)
    checkmate::assert_function(spec$dpsi_depsilon)
  }
  TRUE
}
