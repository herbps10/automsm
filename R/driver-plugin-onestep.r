#' Shared plug-in and one-step estimation
#'
#' @noRd
estimate_plugin_and_onestep <- function(
  Lm_fn, psi, design_matrix, Q, Delta, p, joint_draws, seed
) {
  n <- dim(design_matrix)[1]

  plugin <- B(Lm_fn, psi, design_matrix, Q, p)
  beta <- torch::torch_tensor(plugin, requires_grad = TRUE)

  onestep_est <- onestep(
    Lm_fn,
    psi,
    beta,
    design_matrix,
    Q,
    Delta,
    p
  )

  if(!is.null(seed)) set.seed(seed)
  onestep_joint <- mvtnorm::rmvnorm(
    joint_draws,
    mean = as.numeric(onestep_est$est),
    sigma = var(as.matrix(onestep_est$eif)) / n
  )

  list(
    plugin = list(est = as.numeric(plugin)),
    onestep = list(
      est = as.numeric(onestep_est$est),
      se = as.numeric(onestep_est$se),
      lower = as.numeric(onestep_est$lower),
      upper = as.numeric(onestep_est$upper),
      eif = onestep_est$eif,
      joint_draws = onestep_joint
    )
  )
}
