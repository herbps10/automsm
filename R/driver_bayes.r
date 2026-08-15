#' Jacobian of eps -> beta under the Bayes fluctuation model
#'
#' J(eps) = sum_k dB_dpsi\[k\]' %*% (dpsi_k/deps) + dB_dQ' %*% (dQ/deps)
#'
#' @noRd
bayes_jacobian <- function(problem, psi, Q, beta, dpsi_list, epsilon, K_Q) {
  p <- problem$p
  J_psi <- dB_dpsi(problem$Lm_fn, psi, beta, problem$design_matrix, Q, p) # (K, n, p)
  jac <- torch::torch_zeros(c(p, p))
  for(k in seq_len(problem$K)) {
    jac <- jac + J_psi[k, , ]$transpose(1, 2)$matmul(dpsi_list[[k]])
  }

  Qc <- Q$detach()$clone()$requires_grad_(TRUE)
  J_Q <- dB_dQ(problem$Lm_fn, psi, beta, problem$design_matrix, Qc, p)
  jac <- J_Q$transpose(1, 2)$matmul(dQ_fluctuation_depsilon(epsilon, K_Q, Q))
}

#' Generalized Bayesian log-density of the fluctuation parameter
#'
#' The generalized posterior is
#'   log pi(eps | O) = log f(O | eps) + log pi_beta(B(eps)) + log|det J(eps)|
#'
#' Returns the list expected by AdaptMCMC: the first element is the
#' log-density, remaining elements are recorded in `extra.values`
#' @noRd
bayes_log_density <- function(epsilon, problem, spec, state, clever, K_Q, condvar, prior) {
  step <- spec$steps(problem)[[1]]

  state <- spec$apply_update(problem, step, state, epsilon, clever, K_Q)

  psi <- spec$psi_from_state(problem, state)
  beta <- B(problem$Lm_fn, psi, problem$design_matrix, state$Q, problem$p)

  loglik <- spec$bayes_loglik(epsilon, problem, state, clever, state$Q, condvar)
  dpsi <- spec$dpsi_depsilon(problem, state, clever, epsilon)
  jac <- bayes_jacobian(problem, psi, state$Q, beta, dpsi, epsilon, K_Q)

  target <- loglik + prior(as.numeric(beta)) + as.numeric(jac$det()$abs()$log())

  list(log.density = as.numeric(target), beta = as.numeric(beta))
}

#' @noRd
run_bayes_tmle <- function(problem, spec, fit, control) {
  if(!spec$supports_bayes) {
    stop("Generalized Bayesian TMLE is not supported for estimand '", spec$estimand, "'.", call. = FALSE)
  }
  if(!fit$converged) {
    warning("TMLE did not converge; skipping the generalized Bayesian step.")
    return(list(samples = NULL, acc_rate = NA_real_))
  }

  p <- problem$p
  final <- fit$final
  condvar <- if(is.null(problem$nuisance$condvar)) NULL else as_float_tensor(problem$nuisance$condvar)

  log_dens <- function(epsilon) {
    bayes_log_density(torch::torch_tensor(epsilon), problem, spec, fit$state, final$clever, final$K_Q, condvar, control$prior)
  }

  samples <- array(dim = c(control$chains, control$draws, p))
  acc <- 0
  for(chain in seq_len(control$chains)) {
    mcmc <- adaptMCMC::MCMC(
      log_dens, n = control$draws, init = as.numeric(final$epsilon),
      adapt = TRUE, acc.rate = control$acc_rate, scale = rep(control$scale, p)
    )
    samples[chain, , ] <- matrix(unlist(mcmc$extra.values), ncol = p, nrow = control$draws, byrow = TRUE)

    acc <- acc + mcmc$acceptance.rate / control$chains
  }
  list(samples = samples, acc_rate = acc)
}
