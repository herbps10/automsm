fit_msm <- function(problem, spec,
                    tmle = tmle_control(),
                    bayes = bayes_control(),
                    onestep = onestep_control()) {
  if(bayes$enabled && !tmle$enabled) {
    stop("`bayes = TRUE` requires `tmle = TRUE`.", call. = FALSE)
  }

  if(bayes$enabled && !spec$supports_bayes) {
    stop("Generalized Bayesian TMLE is not supported for estimand '", spec$estimand, "'.", call. = FALSE)
  }

  state0 <- spec$init_state(problem)
  psi0 <- spec$psi_from_state(problem, state0, detach = TRUE)

  base <- estimate_plugin_and_onestep(
    problem$Lm_fn, psi0, problem$design_matrix, problem$Q0,
    spec$delta(problem, state0), problem$p,
    joint_draws = onestep$joint_draws, seed = onestep$seed
  )

  state0$beta <- torch::torch_tensor(base$plugin$est, requires_grad = TRUE)

  tmle_res <- NULL
  fit <- NULL
  if(tmle$enabled) {
    fit <- run_tmle(problem, spec, tmle, state0)
    tmle_res <- finalize_tmle(problem, spec, fit)
  }

  bayes_res <- if(bayes$enabled) run_bayes_tmle(problem, spec, fit, bayes) else NULL

  list(base = base, fit = fit, tmle = tmle_res, bayes = bayes_res, state0 = state0)
}
