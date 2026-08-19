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

  validate_nuisance(problem, spec, bayes_enabled = bayes$enabled)

  checkmate::assert_number(bayes$prior(rep(0, problem$p)), finite = TRUE, .var.name = "bayes_control(prior)")

  problem$clamp <- if(tmle$linear) NULL else tmle$clamp %||% problem$nuisance_control$epsilon

  state0 <- spec$init_state(problem)
  psi0 <- spec$psi_from_state(problem, state0)

  base <- estimate_plugin_and_onestep(
    problem$Lm_fn, psi0, problem$design_matrix, problem$Q0,
    spec$delta(problem, state0), problem$p,
    joint_draws = onestep$joint_draws, seed = onestep$seed,
    batched_beta = problem$batched_beta
  )

  state0$beta <- torch::torch_tensor(base$plugin$est, requires_grad = TRUE)

  tmle_res <- NULL
  fit <- NULL
  if(tmle$enabled) {
    fit <- spec$driver(problem, spec, tmle, state0)
    tmle_res <- finalize_tmle(problem, spec, fit, tmle)
  }

  bayes_res <- if(bayes$enabled) {
    run_bayes_tmle(problem, spec, fit, bayes, eif = tmle_res$eif)
  } else NULL

  list(base = base, fit = fit, tmle = tmle_res, bayes = bayes_res, state0 = state0)
}

#' Shared plug-in and one-step estimation
#'
#' @noRd
estimate_plugin_and_onestep <- function(
  Lm_fn, psi, design_matrix, Q, Delta, p, joint_draws, seed, batched_beta
) {
  n <- dim(design_matrix)[1]

  plugin <- B(Lm_fn, psi, design_matrix, Q, p)
  beta <- torch::torch_tensor(plugin, requires_grad = TRUE)

  onestep_eif <- eif(Lm_fn, psi, beta, design_matrix, Q, Delta, p, batched_beta)
  onestep_est <- beta + colMeans(onestep_eif)
  onestep_se <- apply(onestep_eif, 2, sd) / sqrt(n)
  onestep_lower <- onestep_est + stats::qnorm(0.025) * onestep_se
  onestep_upper <- onestep_est + stats::qnorm(0.975) * onestep_se


  sample_joint <- function() {
    mvtnorm::rmvnorm(
      joint_draws,
      mean = as.numeric(onestep_est),
      sigma = var(as.matrix(onestep_eif)) / n
    )
  }

  if(joint_draws > 0) {
    onestep_joint <- if(is.null(seed)) sample_joint() else withr::with_seed(seed, sample_joint())
  }
  else {
    onestep_joint <- NULL
  }

  list(
    plugin = list(est = as.numeric(plugin)),
    onestep = list(
      est = as.numeric(onestep_est),
      se = as.numeric(onestep_se),
      lower = as.numeric(onestep_lower),
      upper = as.numeric(onestep_upper),
      eif = onestep_eif,
      psi = as.matrix(psi$detach()$cpu()),
      joint_draws = onestep_joint
    )
  )
}
